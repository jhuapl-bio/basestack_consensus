#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np
import pysam

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
    for pileupcolumn in bamFile.pileup(refName, max_depth=10000, truncate=False, min_base_quality=13):

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
    
    # fix zero versus one based indexing
    position = position-1
    
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
    for pileupcolumn in bamFile.pileup(refName, start=position, stop=position+1, max_depth=10000, min_base_quality=13):
        for pileupread in pileupcolumn.pileups:
            if pileupcolumn.pos == position and not pileupread.is_refskip:
                pileup[0] += 1
                #if not pileupread.is_del:
                if not (pileupread.query_position) is None:
                    base=pileupread.alignment.query_sequence[pileupread.query_position]
                    pileup.append(base)
    
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


def minor_allele_freq(depth,pileup,alt,maf_flag,hold_flag):
    """
    Function that returns a MAF flag string if the cumulative minor allele frequency
    at a position is higher than a pre-specified value and indicates if the position
    is a candidate within host variant, or a potentially worrisome mixed position
    """
    
    # determine the allele frequency of the listed alt allele
    AF = float(pileup.count(alt)/depth)
    
    # convert the flag thresholds to decimals
    maf = maf_flag/100.0
    hold = hold_flag/100.0
    
    # the case that there are no flags
    if AF<maf or AF>(1-maf):
        return(AF,np.nan,np.nan)
    
    # if there are flags, distinguish between the isnv and mixed scenarios
    if maf<AF<hold or (1-hold)<AF<(1-maf):
        return(AF,'%0.2f<maf<%0.2f' % (maf,hold),np.nan)
    elif hold<AF<(1-hold):
        return(AF,np.nan,'mixed position')


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
    
    If illumina data is available, the support vector will be 6 bits instead of 3
    """
    
    # here we are only interested in mismatches between nanopore variant calles
    supp_vec = supp_vec[:3]
    
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



def strand_bias_detected(strandAF,strand_threshold):
    """ 
    Function that returns a flag string if a variant is called unequally on the forward and reverse strands
    strandAF order is: positive alts, total positive reads, negative alts, total negative reads
    """
    
    # parse the strandAF string
    strandAF = [int(x) for x in strandAF.split(',')]
    
    # get frequency for each strand
    posAF = float(strandAF[0]/strandAF[1])*100
    negAF = float(strandAF[2]/strandAF[3])*100
    
    # compare frequencies to threshold
    if (posAF<strand_threshold) and (negAF<strand_threshold):
        return(np.nan) # no bias if both are low frequency
    elif posAF<strand_threshold:
        return('strand bias: low +AF')
    elif negAF<strand_threshold:
        return('strand bias: low -AF')
    else:
        return(np.nan) # no bias if both are high frequency


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


def get_allele_counts(pileup,depth):
    
    # get the allele counts to output
    A_count = pileup.count('A')
    T_count = pileup.count('T')
    C_count = pileup.count('C')
    G_count = pileup.count('G')
    O_count = depth - (A_count+T_count+C_count+G_count)
    allele_string = 'A:%d:T:%d:C:%d:G:%d:O:%d' % (A_count,T_count,C_count,G_count,O_count)
    
    return(allele_string)


def in_homopolymer_region(pos):
    """ 
    Function that reports if the position is in a known homopolymer region
    Currently uses a hard-coded list of positions, but can be expanded to take in output of other studies
    """
    
    # current homopolymer list
    hp = [241,3037,11083,12119,29700]
    
    if pos in hp:
        return(True)
    else:
        return(False)
    

def ont_illumina_mismatch(illumina,supp_vec):
    """ 
    Function that reports if there is a variant called by nanopore callers but not by
    one or more illumina variant callers
    """
    
    if not illumina:
        return(np.nan)
    else:
        assert len(supp_vec)==6
        
        # we are only interested in any potential illumina mismatches
        # we assume this position was called by at least one nanopore variant caller
        supp_vec = supp_vec[3:]
        
        if supp_vec == '111':
            return(np.nan)
        elif supp_vec == '000':
            return('mismatch()')
        elif supp_vec == '100':
            return('mismatch(f)')
        elif supp_vec == '010':
            return('mismatch(i)')
        elif supp_vec == '001':
            return('mismatch(s)')
        elif supp_vec == '110':
            return('mismatch(f+i)')
        elif supp_vec == '101':
            return('mismatch(f+s)')
        elif supp_vec == '011':
            return('mismatch(i+s)')
        else:
            sys.exit('%s is not a valid support vector' % supp_vec)
    


def add_key_ambiguous_positions(chrom,variants,cons,depth_threshold,vcf_nextstrain,bamfile):
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
            pileup=collect_position_pileup(bamfile,pos)
            depth=pileup[0]
            pileup = pileup[1:]
            data={}
            data['chrom']=chrom
            data['pos']=pos
            data['ref']=ns_snps[ns_snps.POS==pos].REF.values[0]
            data['alt']='N'
            data['consensus_base']='N'
            data['in_consensus']=False
            data['unambig']=False
            data['ont_depth']=depth
            data['ont_depth_thresh']=depth_threshold
            data['ont_alleles']=get_allele_counts(pileup,depth)
            data['key_flag']='ambig in key position'
            data['status']='Yes*'
            
            data = pd.DataFrame([data], columns=data.keys())
            df = pd.concat([df,data],ignore_index=True)
    
    return(df)


def status_by_flags(data):
    """
    Determine the status of each variant based on flags and other conditions
    """
    
    # flag keywords are 'depth','MAF','new','NTC','mismatch','key','SB'
    # these are unique words to each type of flag
    
    # situations that automatically lead to maybe
    maybe_flags = ['depth_flag','ntc_flag']
    for flag in maybe_flags:
        if not pd.isna(data[flag]):
            return('Maybe')
    
    # situations that don't depend on mismatches
    if not pd.isna(data['mixed_flag']) and pd.isna(data['illumina_vc_flag']) and data['illumina_AF']>0.85 and data['in_consensus'] and data['illumina'] and data['illumina_depth']:
        return('Yes')
        
    # specific situations with mismatches
    if not pd.isna(data['ont_vc_flag']):
        if data['ont_vc_flag']=='mismatch(s)' and data['illumina_vc_flag']=='mismatch()' and data['in_consensus']==False and data['illumina'] and data['illumina_depth']:
            return('Yes')
        elif data['ont_AF']<0.2 and data['in_consensus']==False and data['ont_vc_flag']=='mismatch(s)':
            return('Yes')
        elif data['ont_AF']<0.25 and data['in_consensus']==False and not pd.isna(data['sb_flag']):
            return('Yes')
        elif 0.2<data['ont_AF']<0.3 and data['in_consensus']==False and data['ont_vc_flag']=='mismatch(s)':
            return('Yes*')
        
        else:
            return('Maybe')
    
    # if we have gotten to this point
    # there is no mismatch within the ont data
    if data['ont_AF']<0.2 and data['in_consensus']==False:
        return('Yes')
    elif 0.2<data['ont_AF']<0.3 and data['in_consensus']==False:
        return('Yes*')
    
    # account for the unlikely case that there is no mismatch but there is strand bias
    if not pd.isna(data['sb_flag']):
        return('Maybe')
    
    # mixed or illumina mismatch flags not yet accounted for should be maybe
    if not pd.isna(data['mixed_flag']):
        return('Maybe')
    if not pd.isna(data['illumina_vc_flag']):
        return('Maybe')
    
    # other specific situations
    if data['ont_AF']<0.85 and data['in_consensus']==True:
        return('Yes*')
    if not pd.isna(data['key_flag']):
        return('Yes*')
    if not pd.isna(data['new_flag']):
        if data['illumina'] and data['illumina_depth']:
            return('Yes')
        else:
            return('Yes*')
    
    # if we make it this far there are no flags
    return('Yes')


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
   parser.add_argument('--vcf-nextstrain', type=str, help='path to vcf containing all nextstrain snps')
   
   parser.add_argument('--outdir', '-o', type=str, help='directory name to write output files to')
   parser.add_argument('--prefix', type=str, default='sample', help='prefix of all saved output files')
   
   parser.add_argument('--coverage-flag', type=int, default=20, help='flag variants with depth within this percentage of threshold')
   parser.add_argument('--maf-flag', type=int, default=15, help='flag variants with minor allele frequency with at least this value')
   parser.add_argument('--hold-flag', type=int, default=30, help='flag variants for additional validation if minor allele frequency is at least this value')
   parser.add_argument('--call-depth-factor', type=int, default=2, help='factor by which depth must exceed median NTC depth to call a base')
   parser.add_argument('--snp-depth-factor', type=int, default=5, help='factor by which depth must exceed NTC depth to call a variant seen in the NTC at that position')
   parser.add_argument('--unambig-threshold', type=int, default=25000, help='number of unambiguous bases required in final genome')
   parser.add_argument('--ns-snp-threshold', type=int, default=5, help='number of published samples with a particular snp needed to count it as previously seen')
   parser.add_argument('--strand-threshold', type=int, default=20, help='minimum minor allele frequency on each strand required for unbiased call')
   
   args = parser.parse_args()
   return(args)


def main():
    
    args = parse_arguments()
    
    depth_threshold = max(20,calculate_depth_threshold(args.ntc_bamfile, args.call_depth_factor))
    mask_cons = mask_consensus_sites(args.consensus, args.bamfile, depth_threshold, args.outdir, args.prefix)
    
    # read vcf as text file
    vcf_sample = pd.read_csv(args.vcffile,sep='\t',skiprows=2)
    vcf_sample = vcf_sample[['#CHROM','POS','REF','ALT','QUAL','INFO']]
    vcf_sample.columns = ['chrom','pos','ref','alt','nanopolish_qual','info']
    
    # get the path to illumina bam file
    illumina_bam=open(args.vcffile).readlines()[1]
    illumina_bam=illumina_bam.split('=')[1].strip() # remove whitespace characters
    if illumina_bam=='None':
        illumina=False
    else:
        assert os.path.exists(illumina_bam)
        illumina=True
    
    # load current consensus sequence
    cons = list(SeqIO.parse(open(mask_cons),"fasta"))[0]
    cons = list(cons.seq.upper())
    
    # set up the dataframe to store results
    df = pd.DataFrame(
        columns=['chrom','pos','ref','alt','consensus_base','status','homopolymer','in_consensus','unambig',
                 'ont_depth','illumina_depth','ont_depth_thresh','illumina_depth_thresh','ont_AF','illumina_AF','ont_alleles','illumina_alleles','strand_counts',
                 'medaka_qual','nanopolish_qual','illumina','min_illumina_depth','maf_flag',
                 'mixed_flag','depth_flag','ntc_flag','new_flag','ont_vc_flag','illumina_vc_flag','sb_flag','key_flag'])
        
    # loop through all positions
    for pos in vcf_sample.pos:
        
        tmp = vcf_sample[vcf_sample.pos==pos]
        
        # ignore indels
        # print a warning if indels are found in the input vcf
        if len(tmp.ref.values[0]) != len(tmp.alt.values[0]):
            print('Warning: indel found at position %d - indels are ignored!' % pos)
            continue
        
        # start a dictionary to store data for this sample
        data = tmp.to_dict('records')[0]
        
        # store information from info column and remove it from dictionary
        info = dict(item.split("=") for item in data['info'].split(";"))
        del data['info']
        
        # ignore positions with no nanopore variant calls
        if illumina:
            if info['SUPP_VEC'][:3]=='000':
                continue
        
        # get read depth and pileup at this read position
        pileup = collect_position_pileup(args.bamfile, pos)
        depth = pileup[0]
        pileup = pileup[1:]
        
        if illumina:
            illumina_pileup = collect_position_pileup(illumina_bam, pos)
            illumina_depth = illumina_pileup[0]
            illumina_pileup = illumina_pileup[1:]
            
        # ignore this position if the depth is too low
        if depth < depth_threshold:
            continue
        
        # add basic data to this record
        data['ont_depth'] = depth
        data['illumina_depth'] = [illumina_depth if illumina else np.nan][0]
        data['ont_depth_thresh'] = depth_threshold
        data['illumina_depth_thresh'] = [20 if illumina else np.nan][0]
        data['ont_alleles'] = get_allele_counts(pileup,depth)
        data['strand_counts'] = [info['STRANDAF'] if 'STRANDAF' in info.keys() else np.nan][0]
        data['medaka_qual'] = [info['pred_q'] if 'pred_q' in info.keys() else np.nan][0]
        data['illumina'] = illumina
        data['illumina_alleles'] = [get_allele_counts(illumina_pileup, illumina_depth) if illumina else np.nan][0]
        data['min_illumina_depth'] = [illumina_depth>20 if illumina else np.nan][0]
        data['illumina_AF'] = [float(illumina_pileup.count(data['alt'])/illumina_depth) if illumina else np.nan][0]
        
        # add flags to this record
        data['depth_flag'] = depth_near_threshold(depth,pileup,depth_threshold,args.coverage_flag)
        data['new_flag'] = snp_in_nextstrain(pos, data['ref'], data['alt'], args.vcf_nextstrain, args.ns_snp_threshold)
        data['ont_vc_flag'] = variant_caller_mismatch(info['SUPP_VEC'])
        data['sb_flag'] = [strand_bias_detected(info['STRANDAF'], args.strand_threshold) if 'STRANDAF' in info.keys() else np.nan][0]
        data['ntc_flag'] = allele_in_ntc(pos, data['alt'], depth, args.ntc_bamfile, args.snp_depth_factor)
        data['homopolymer'] = in_homopolymer_region(pos)
        data['ont_AF'],data['maf_flag'],data['mixed_flag'] = minor_allele_freq(depth, pileup, data['alt'], args.maf_flag, args.hold_flag)
        data['illumina_vc_flag'] = ont_illumina_mismatch(illumina, info['SUPP_VEC'])
        
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
        if (data['unambig']==True) and (info['SUPP_VEC'][:3] in ['111','110','101','100']):
            assert cons[pos-1]==data['alt']
            data['in_consensus'] = True
        else:
            data['in_consensus'] = False
        
        # report the base called in the consensus genome
        data['consensus_base'] = cons[pos-1]
        
        # double check we are reporting the consensus base correctly
        if data['in_consensus']==True and data['unambig']==True:
            assert data['consensus_base']==data['alt']
        elif data['in_consensus']==False and data['unambig']==True:
            assert data['consensus_base']==data['ref']
        elif data['unambig']==False:
            assert data['consensus_base']=='N'
        
        # after adding all the flags
        # get the status of this position
        data['status'] = status_by_flags(data)
        # add this record to the final dataframe
        data = pd.DataFrame([data], columns=data.keys())
        df = pd.concat([df,data],ignore_index=True,sort=False)
    
    # after looping through all positions
    # add positions at key snps if not in variant list
    df = pd.concat([df,add_key_ambiguous_positions(vcf_sample.chrom.values[0],list(df.pos.values), cons, depth_threshold, args.vcf_nextstrain, args.bamfile)],ignore_index=True,sort=False)
    df = df.replace(np.nan,'.')
    filepath = os.path.join(args.outdir,args.prefix+'.variant_data.txt')
    df.to_csv(filepath,sep='\t',index=False)
    
    # then make the final consensus
    make_final_fasta(cons,args.prefix,args.unambig_threshold,args.outdir)


if __name__ == "__main__":
    main()