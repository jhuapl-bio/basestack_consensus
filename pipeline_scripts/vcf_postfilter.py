#!/usr/bin/env python

import os
import sys
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
    
    # mask sites at the beginning and end where amplicons do not cover
    # change basecalls to N if coverage is below threshold
    # save list of newly-masked bases
    ambig=[]
    newmask=[]
    for pos,base in enumerate(cons):
        
        # mask beginning and end
        # remember zero indexing
        if (0<=pos<=53) | (29836<=pos<=29902):
            cons[pos] = 'N'
            continue
        
        # save bases that were already 'N'
        # mask based on coverage
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
    

def depth_near_threshold(depth,pileup,depth_threshold,coverage_flag):
    """
    Function that returns a depth flag string if the read depth at a position is
    within a pre-specified percentage of the depth threshold
    """
    
    # get the coverage flag threshold
    frac = float(coverage_flag/100)
    lowend = depth_threshold - (depth_threshold*frac)
    highend = depth_threshold + (depth_threshold*frac)
    
    # check if coverage is close to depth threshold
    if lowend<depth<highend:
        return('depth within %s%% of threshold' % (coverage_flag))
    else:
        return(np.nan)


def high_minor_allele_freq(depth,pileup,alt,maf_flag):
    """
    Function that returns a MAF flag string if the cumulative minor allele frequency
    at a position is higher than a pre-specified value
    """
    
    # flag this position if the cumulative minor allele frequency is high
    # here the minor allele is everything except the called ALT base
    # this can include deletions
    maf = float( (depth - pileup.count(alt)) / depth ) * 100
    if maf >= maf_flag:
        return('MAF>%.2f' % float(maf_flag/100))
    else:
        return(np.nan)


def allele_in_ntc(pos,alt,depth,ntc_bamfile,snp_depth_factor):
    """ 
    Function that returns a flag string if the alternate allele is present in the negative control
    and the coverage in the sample is not more than snp_depth_factor * coverage in negative control
    """
    
    # get the pileup at this position in the negative control
    ntc_pileup = collect_position_pileup(ntc_bamfile, pos)
    
    if alt in ntc_pileup:
        # require coverage at this sample to be some multiple of the negative control
        if depth <= (snp_depth_factor * ntc_pileup[0]):
            return('allele in NTC')
    
    # if alt not in negative control or depth is high enough
    return(np.nan)


def snp_in_nextstrain(pos,ref,alt,vcf_nextstrain,ns_snp_threshold):
    """
    Function that returns a flag string if a SNP has not been seen in published sequences
    Requires the SNP to be found in a specific number of published sequences
    to avoid confounding with SNPs that may be a result of sequencing errors in other samples
    """
    
    # read in the nextstrain vcf as a dataframe
    # note: the header is hard-coded and will need to be updated if the header is altered
    ns_snps = pd.read_csv(vcf_nextstrain,sep='\t',skiprows=3)
    ns_snps = ns_snps[['POS','REF','ALT','TOTAL_SAMPLES','OCCURENCES']]
    
    # if the position has not been variable before, return false
    if pos not in ns_snps.POS.values:
        return('not in nextstrain')
    
    # if the position has been variable before
    # check if the specific allele has been found
    else:
        tmp = ns_snps[ns_snps.POS==pos]
        assert ref == tmp.REF.values[0]
        alleles = tmp.ALT.values[0].split(',')
        
        if alt not in alleles:
            return('not in nextstrain')
        else:
            # if the alternate allele has been found before
            # check if it has been found enough times
            idx = alleles.index(alt)
            counts = [x.split(',') if ',' in str(x) else x for x in tmp.OCCURENCES.values][0]
            if int(counts[idx]) >= ns_snp_threshold:
                return(np.nan)
            else:
                return('not in nextstrain')


def variant_caller_mismatch(supp_vec):
    """
    Function that returns a flag string if a variant has not been detected by all callers
    Currently assumes callers are: nanopolish, medaka, samtools (in that order)
    """
    
    # return different codes for different mismatch strings
    if supp_vec == '111':
        return(np.nan)
    elif supp_vec == '100':
        return('mismatch(n)')
    elif supp_vec == '010':
        return('mismatch(m)')
    elif supp_vec == '001':
        return('mismatch(s)')
    elif supp_vec == '110':
        return('mismatch(n+m)')
    elif supp_vec == '101':
        return('mismatch(n+s)')
    elif supp_vec == '011':
        return('mismatch(m+s)')
    else:
        sys.exit('%s is not a valid support vector' % supp_vec)


def ambig_in_key_position(pos,vcf_nextstrain,cons):
    """ 
    Function that returns a flag string if a position is at an important site
    but is an ambiguous base ('N') in the consensus genome
    """
    
    # read in the nextstrain vcf as a dataframe
    # note: the header is hard-coded and will need to be updated if the header is altered
    ns_snps = pd.read_csv(vcf_nextstrain,sep='\t',skiprows=3)
    ns_snps = ns_snps[['POS','CONF_FLAG']]
    
    key_snps = ns_snps[ns_snps['CONF_FLAG']=='YES']
    key_snps = list(key_snps.POS.values)
    
    # no flag needed if this position is not one of the important ones
    if pos not in key_snps:
        return(np.nan)
    
    # if it is an important position
    else:
        if cons[pos-1]=='N':
            return('ambig in key position')
        else:
            return(np.nan)


def add_key_ambiguous_positions(variants,cons,vcf_nextstrain):
    """ 
    Function that returns a dataframe of positions not called as variants in a sample
    but that are ambiguous at key positions
    """
    
    # read in the nextstrain vcf as a dataframe
    # note: the header is hard-coded and will need to be updated if the header is altered
    ns_snps = pd.read_csv(vcf_nextstrain,sep='\t',skiprows=3)
    ns_snps = ns_snps[['POS','REF','CONF_FLAG']]
    
    key_snps = ns_snps[ns_snps['CONF_FLAG']=='YES']
    key_snps = list(key_snps.POS.values)
    
    # create a dataframe in which to store any new values
    df = pd.DataFrame()
    
    # loop through important snps
    for pos in key_snps:
        if (pos not in variants) and (cons[pos-1]=='N'):
            data={}
            data['pos']=pos
            data['ref']=ns_snps[ns_snps.POS==pos].REF.values[0]
            data['alt']='N'
            data['unambig']=False
            data['key_flag']='ambig in key position'
            data = pd.DataFrame([data], columns=data.keys())
            df = pd.concat([df,data],ignore_index=True)
    
    return(df)


def get_allele_counts(pileup,depth):
    
    # get the allele counts to output
    A_count = pileup.count('A')
    T_count = pileup.count('T')
    C_count = pileup.count('C')
    G_count = pileup.count('G')
    O_count = depth - (A_count+T_count+C_count+G_count)
    allele_string = 'A:%d:T:%d:C:%d:G:%d:O:%d' % (A_count,T_count,C_count,G_count,O_count)
    
    return(allele_string)


def make_final_fasta(cons,prefix,unambig_thresh,outdir):
    
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


def main():
    
    args = parse_arguments()
    
    depth_threshold = max(20,calculate_depth_threshold(args.ntc_bamfile, args.call_depth_factor))
    mask_cons = mask_consensus_sites(args.consensus, args.bamfile, depth_threshold, args.outdir, args.prefix)
    
    # read vcf as text file
    vcf_sample = pd.read_csv(args.vcffile,sep='\t',skiprows=1)
    vcf_sample = vcf_sample[['POS','REF','ALT','INFO']]
    vcf_sample.columns = ['pos','ref','alt','info']
    
    # load current consensus sequence
    cons = list(SeqIO.parse(open(mask_cons),"fasta"))[0]
    cons = list(cons.seq.upper())
    
    # set up the dataframe to store results
    df = pd.DataFrame(
        columns=['pos','ref','alt','in_consensus','unambig','depth','depth_thresh','alleles','depth_flag','maf_flag','ntc_flag','new_flag','vc_flag','key_flag'])
        
    # loop through all positions
    for pos in vcf_sample.pos:
        
        tmp = vcf_sample[vcf_sample.pos==pos]
        
        # ignore indels
        # print a warning if indels are found in the input vcf
        if len(tmp.ref.values[0]) != len(tmp.alt.values[0]):
            warnings.warn('Indel found at position %d\n indels are not carried over to the final VCF!' % pos)
            continue
        
        # start a dictionary to store data for this sample
        data = tmp.to_dict('records')[0]
        
        # store information from info column and remove it from dictionary
        info = dict(item.split("=") for item in data['info'].split(";"))
        del data['info']
        
        # get read depth and pileup at this read position
        pileup = collect_position_pileup(args.bamfile, pos)
        depth = pileup[0]
        pileup = pileup[1:]
            
        # ignore this position if the depth is too low
        if depth < depth_threshold:
            continue
        
        # add basic data to this record
        data['depth'] = depth
        data['depth_thresh'] = depth_threshold
        data['alleles'] = get_allele_counts(pileup,depth)
        
        # add flags to this record
        data['depth_flag'] = depth_near_threshold(depth,pileup,depth_threshold,args.coverage_flag)
        data['maf_flag'] = high_minor_allele_freq(depth, pileup, data['alt'], args.maf_flag)
        data['new_flag'] = snp_in_nextstrain(pos, data['ref'], data['alt'], args.vcf_nextstrain, args.ns_snp_threshold)
        data['vc_flag'] = variant_caller_mismatch(info['SUPP_VEC'])
        data['ntc_flag'] = allele_in_ntc(pos, data['alt'], depth, args.ntc_bamfile, args.snp_depth_factor)
        
        # modify consensus genome based on ntc flag
        # remember consensus is zero-indexed but we are dealing with 1-indexed positions
        if not pd.isna(data['ntc_flag']):
            cons[pos-1]='N'
            
        # add a flag if a key position is ambiguous in the consensus
        # must happen after any masking that occurs due to NTC flags
        data['key_flag'] = ambig_in_key_position(pos, args.vcf_nextstrain, cons)
        
        # mark which positions are unambiguous in the consensus
        data['unambig'] = [False if cons[pos-1]=='N' else True][0]
        
        # mark which variants are actually in the consensus genome
        if (data['unambig']==True) and (info['SUPP_VEC'] in ['111','110','101','100']):
            assert cons[pos-1]==data['alt']
            data['in_consensus'] = True
        else:
            data['in_consensus'] = False
        
        # after adding all the flags
        # add this record to the final dataframe
        data = pd.DataFrame([data], columns=data.keys())
        df = pd.concat([df,data],ignore_index=True,sort=False)
    
    # after looping through all positions
    # add positions at key snps if not in variant list
    df = pd.concat([df,add_key_ambiguous_positions(list(df.pos.values), cons, args.vcf_nextstrain)],ignore_index=True,sort=False)
    filepath = os.path.join(args.outdir,args.prefix+'.variant_data.txt')
    df.to_csv(filepath,sep='\t',index=False)
    
    # then make the final consensus
    make_final_fasta(cons,args.prefix,args.unambig_threshold,args.outdir)


if __name__ == "__main__":
    main()