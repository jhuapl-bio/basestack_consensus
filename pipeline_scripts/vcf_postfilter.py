#!/usr/bin/env python

import os
import vcf
import argparse
import pandas as pd
import numpy as np
import pysam
import warnings

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def collect_depths(bamfile):
    """
    Collect read depth of coverage per reference position in a BAM file
    Modified from: https://github.com/artic-network/fieldbioinformatics/ artic_make_depth_mask.py
    """
    
    # check the BAM file exists
    if not os.path.exists(bamfile):
        raise Exception("bamfile doesn't exist (%s)" % bamfile)

    # open the BAM file
    bamFile = pysam.AlignmentFile(bamfile, 'rb')

    # get the reference name from the bamfile
    refName = bamFile.get_reference_name(0)

    # create a depth vector to hold the depths at each reference position
    depths = [0] * bamFile.get_reference_length(refName)

    # generate the pileup
    for pileupcolumn in bamFile.pileup(refName, max_depth=10000, truncate=False, min_base_quality=0):

        # process the pileup column
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_refskip:
                depths[pileupcolumn.pos] += 1 #includes deletions

    # close file and return depth vector
    bamFile.close()
    return depths


def collect_position_pileup(bamfile,position):
    """
    Get full pileup at a single position
    Used to calculate depth and minor allele frequency at that position
    """
    
    # check the BAM file exists
    if not os.path.exists(bamfile):
        raise Exception("bamfile doesn't exist (%s)" % bamfile)

    # open the BAM file
    bamFile = pysam.AlignmentFile(bamfile, 'rb')

    # get the reference name from the bamfile
    refName = bamFile.get_reference_name(0)
    
    # get pileup and depth at this position
    # pileup uses 0-based indexing, so need position-1 for actual genomic position
    pileup = []
    pileup.append(0)
    for pileupcolumn in bamFile.pileup(refName, start=position-1, stop=position, max_depth=10000, min_base_quality=0):
        for pileupread in pileupcolumn.pileups:
            if pileupcolumn.pos == position and not pileupread.is_refskip:
                pileup[0] += 1
                if not pileupread.is_del:
                    pileup.append(pileupread.alignment.query_sequence[pileupread.query_position-1])
    
    bamFile.close()
    return pileup # first position is depth including deletions


def calculate_depth_threshold(ntc_bamfile,call_depth_factor):
    """
    Determine the read depth threshold to use for calling non-ambiguous bases
    based on a negative control (NTC) included on the same sequencing run
    """
    
    # get negative control depth
    cov = collect_depths(ntc_bamfile)
    median_depth = np.median(cov)
    return(call_depth_factor*median_depth)


def mask_consensus_sites(consensus,bamfile,depth_threshold,outdir,prefix):
    """
    Mask sites in the consensus genome based on the depth threshold
    calculated from the negative control on the same run
    
    This function could be replaced by adding the --depth parameter when calling artic_make_depth_mask

    """
    
    # get depth across genome
    cov = collect_depths(bamfile)
    
    # load current consensus sequence
    cons = list(SeqIO.parse(open(consensus),"fasta"))[0]
    cons = list(cons.seq.upper())
    
    assert len(cons)==len(cov)
    
    # change basecalls to N if coverage is below threshold
    # save list of newly-masked bases
    ambig=[]
    newmask=[]
    for pos,base in enumerate(cons):
        if base=='N':
            ambig.append(pos)
            continue
        else:
            rd = cov[pos]
            if rd < depth_threshold:
                cons[pos] = 'N'
                newmask.append(pos)
    
    # output newly-masked bases to file
    filename=os.path.join(outdir,prefix+'.new_masked_sites.txt')
    with open(filename, 'w') as file:
        file.write('ambiguous bases in %s.consensus.fasta\n' % prefix)
        file.write(",".join([str(x) for x in ambig]))
        file.write('\n\nadditional masked positions with depth requirement of %d:\n' % depth_threshold)
        file.write(",".join([str(x) for x in newmask]))
    
    # save new consensus sequence and output to file
    seq = ''.join(cons)
    new_record = SeqRecord(Seq(seq),id=prefix,description="")
    filepath = os.path.join(outdir,prefix+'.mask.fasta')
    SeqIO.write(new_record,filepath,"fasta")
    
    # return path to new consensus genome
    return(filepath)
    

def snp_in_nextstrain(pos,alt,vcf_nextstrain,ns_snp_threshold):
    """
    Function that returns true if a SNP has been seen in published sequences
    Requires the SNP to be found in a specific number of published sequences
    to avoid confounding with SNPs that may be a result of sequencing errors in other samples
    
    Currently the occurence is based on all alternative alleles in published sequences
    Soon to be updated using the occurence of a particular allele
    """
    
    # read in the nextstrain vcf as a dataframe
    # note: the header is hard-coded and will need to be updated if the header is altered
    ns_snps = pd.read_csv(vcf_nextstrain,sep='\t',skiprows=3)
    ns_snps = ns_snps[['POS','ALT','TOTAL_SAMPLES','OCCURENCES']]
    
    # if the position has not been variable before, return false
    if pos not in ns_snps.POS.values:
        return(False)
    
    # if the position has been variable before
    # check if the specific allele has been found
    else:
        tmp = ns_snps[ns_snps.POS==pos]
        alleles = tmp.ALT.values[0].split(',')
        
        if alt not in alleles:
            return(False)
        else:
            # if the alternate allele has been found before
            # check if it has been found enough times
            #counts = [x.split(',') if ',' in str(x) else x for x in tmp.OCCURENCES.values]
            if tmp.OCCURENCES.values[0] >= ns_snp_threshold:
                return(True)
            else:
                return(False)


def get_allele_counts(pileup,depth):
    
    # get the allele counts to output
    A_count = pileup.count('A')
    T_count = pileup.count('T')
    C_count = pileup.count('C')
    G_count = pileup.count('G')
    O_count = depth - (A_count+T_count+C_count+G_count)
    allele_string = 'A:%d:T:%d:C:%d:G:%d:O:%d' % (A_count,T_count,C_count,G_count,O_count)
    
    return(allele_string)


def refine_variant_calls(vcffile,bamfile,ntc_bamfile,consensus,coverage_flag,depth_threshold,maf_flag,snp_depth_factor,vcf_nextstrain,ns_snp_threshold,outdir,prefix):
    """
    Annotate VCF with information about confidence in variant calls:
        - note if a variant is present in the consensus genome (y/n)
        - flag a variant if the read depth at this variant position is close to the depth threshold
        - flag a variant if the cumulative frequency of all minor alleles is greater than maf_flag percent
    Additionally, impose more stringent read depth requirements on any variants also present in the negative control
    """
    
    # read vcf
    vcf_sample = vcf.Reader(filename=vcffile)
    
    # load current consensus sequence
    cons = list(SeqIO.parse(open(consensus),"fasta"))[0]
    cons = list(cons.seq.upper())
    
    # save list of all variant flag data
    var_data = []
    
    # get the coverage flag threshold
    frac = float(coverage_flag/100)
    lowend = depth_threshold - (depth_threshold*frac)
    highend = depth_threshold + (depth_threshold*frac)
    
    for record in vcf_sample:
        
        # ignore indels
        # print a warning if indels are found in the input vcf
        if len(record.REF) != len(record.ALT[0]):
            warnings.warn('Indel found at position %d\n Please note that indels are not carried over to the final VCF' % record.POS, Warning)
            continue
        
        # deal with multiple consecutive snps
        # for example at position 28881
        for i in range(len(record.REF)):
            
            # initialize an empty dictionary that we will fill in
            pos_data = {}
            
            # get the position and alternate allele for this snp
            pos = record.POS + i # this will just be the position if there is only one snp
            alt = str(record.ALT[0])[i]
            
            # get read depth and pileup at this read position
            pileup = collect_position_pileup(bamfile, pos)
            depth = pileup[0]
            pileup = pileup[1:]

            # check if coverage is close to depth threshold
            if lowend<depth<highend:
                pos_data['depth_flag'] = 'depth within %s%% of threshold' % (coverage_flag)
                    
            # check if this position is called in the consensus genome
            if cons[pos]=='N':
                pos_data['consensus_var'] = False
            else:
                pos_data['consensus_var'] = True
                    
                # apply different filters depending on whether alt allele is found in negative control
                ntc_pileup = collect_position_pileup(ntc_bamfile, pos)
                if alt in ntc_pileup:
                        
                    # require coverage at this sample to be some multiple of the negative control
                    if depth >= (snp_depth_factor * ntc_pileup[0]):
                        pos_data['consensus_var'] = True
                        # flag this position if the cumulative minor allele frequency is high
                        # here the minor allele is everything except the called ALT base
                        # this can include deletions
                        maf = float( (depth - pileup.count(alt)) / depth ) * 100
                        if maf >= maf_flag:
                            pos_data['maf_flag'] = 'MAF>%.2f' % float(maf_flag/100)
                        
                    # if this position does not pass the stringent coverage filter
                    else:
                        # change this base to 'N' in the consensus genome
                        cons[pos] = 'N'
                        pos_data['consensus_var'] = False
                        pos_data['ntc_flag'] = 'allele in NTC'
                    
                # if the allele is not the in negative control
                else:
                    # flag this position if the cumulative minor allele frequency is high
                    # here the minor allele is everything except the called ALT base
                    # this can include deletions
                    maf = float( (depth - pileup.count(alt)) / depth ) * 100
                    if maf >= maf_flag:
                        pos_data['maf_flag'] = 'MAF>%.2f' % float(maf_flag/100)
            
            # add a flag if the snp has not previously been seen before
            if not snp_in_nextstrain(pos, alt, vcf_nextstrain, ns_snp_threshold):
                pos_data['new_flag'] = 'not in nextstrain'
            
            # after checking for all flags
            # add a few values and then append this dictionary to the data list
            pos_data['pos'] = pos
            pos_data['alt'] = alt
            pos_data['ref'] = record.REF[i]
            pos_data['depth'] = depth
            pos_data['alleles'] = get_allele_counts(pileup,depth)
            var_data.append(pos_data)
    
    # after looping through all positions
    
    # output a data frame with positions and flags
    df = pd.DataFrame(var_data)
    
    # add any missing columns for consistency
    for colname in ['pos','ref','alt','consensus_var','depth','depth_thresh','alleles','depth_flag','maf_flag','ntc_flag','new_flag']:
        if colname not in df:
            df[colname] = np.nan
    
    df['flags'] = df[['depth_flag','maf_flag','ntc_flag','new_flag']].apply(lambda x: '; '.join(x.dropna()), axis=1)
    df['depth_thresh'] = depth_threshold
    df = df[['pos','ref','alt','consensus_var','depth','depth_thresh','alleles','flags']]
    
    filepath = os.path.join(outdir,prefix+'.variant_data.txt')
    df.to_csv(filepath,sep='\t',index=False)


def make_final_fasta(consensus,prefix,unambig_thresh,outdir):
    
    # load current consensus sequence
    cons = list(SeqIO.parse(open(consensus),"fasta"))[0]
    cons = list(cons.seq.upper())
    
    assert len(cons)==29903
    
    # count number of ambiguous bases in sequence
    ambig = cons.count('N')
    unambig = 29903-ambig
    
    # save new consensus sequence and output to file
    seq = ''.join(cons)
    new_record = SeqRecord(Seq(seq),id=prefix,description="")
    
    # save new output file depending on unambig threshold
    if unambig > unambig_thresh:
        filepath = os.path.join(outdir,prefix+'.complete.fasta')
    else:
        filepath = os.path.join(outdir,prefix+'.partial.fasta')
        
    SeqIO.write(new_record,filepath,"fasta")


def parse_arguments():
   parser = argparse.ArgumentParser()
   parser.add_argument('--vcffile', type=str, help='path to vcf file of sample')
   parser.add_argument('--bamfile', type=str, help='path to bam file of sample')
   parser.add_argument('--consensus', type=str, help='path to fasta file of sample consensus sequence')
   parser.add_argument('--ntc-bamfile', type=str, help='path to bam file of negative control')
   parser.add_argument('--outdir', '-o', type=str, help='directory name to write output files to')
   parser.add_argument('--vcf-nextstrain', type=str, help='path to vcf containing all nextstrain snps')
   parser.add_argument('--prefix', type=str, default='sample', help='prefix of all saved output files')
   parser.add_argument('--coverage-flag', type=int, default=20, help='flag variants with depth within this percentage of threshold')
   parser.add_argument('--maf-flag', type=int, default=10, help='flag variants with minor allele frequency with at least this value')
   parser.add_argument('--call-depth-factor', type=int, default=2, help='factor by which depth must exceed median NTC depth to call a base')
   parser.add_argument('--snp-depth-factor', type=int, default=5, help='factor by which depth must exceed NTC depth to call a variant seen in the NTC at that position')
   parser.add_argument('--unambig-threshold', type=int, default=25000, help='number of unambiguous bases required in final genome')
   parser.add_argument('--ns-snp-threshold', type=int, default=5, help='number of published samples with a particular snp needed to count it as previously seen')
   
   args = parser.parse_args()
   return(args)

if __name__ == "__main__":
    
    args = parse_arguments()
    
    depth_threshold = max(20,calculate_depth_threshold(args.ntc_bamfile, args.call_depth_factor))
    mask_cons = mask_consensus_sites(args.consensus, args.bamfile, depth_threshold, args.outdir, args.prefix)
    refine_variant_calls(args.vcffile, args.bamfile, args.ntc_bamfile, mask_cons, args.coverage_flag, depth_threshold, args.maf_flag, args.snp_depth_factor, args.vcf_nextstrain, args.ns_snp_threshold, args.outdir,args.prefix)
    make_final_fasta(mask_cons, args.prefix, args.unambig_threshold, args.outdir)