#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

from samtools_funcs import collect_position_pileup

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
    if maf<=AF<hold or (1-hold)<AF<=(1-maf):
        return(AF,'%0.2f<maf<%0.2f' % (maf,hold),np.nan)
    elif hold<=AF<=(1-hold):
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
    if strandAF[1]>0:
        posAF = float(strandAF[0]/strandAF[1])*100
    else:
        posAF=0.0
    
    if strandAF[3]>0:
        negAF = float(strandAF[2]/strandAF[3])*100
    else:
        negAF=0.0
    
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
    

def illumina_mismatch(supp_vec):
    """ 
    Function that reports if there are mismatches in the illumina variant callers
    """
    
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


def get_illumina_support(illumina_AF,supp_vec,maf_flag):
    """
    Function that returns Yes/Maybe/Mixed/No indicating illumina support for a particular variant
    """

    # get the mismatch status for the illumina variant callers
    mismatch = illumina_mismatch(supp_vec)
    
    # get the support threshold we are interested in
    thresh = 1 - (maf_flag/100.0)
    
    # illumina support for a variant
    # use allele frequency to determine if this support is high or mixed
    if pd.isna(mismatch):
        if illumina_AF>thresh:
            return('yes')
        else:
            return('mixed')
    
    # no illumina support for a variat
    if mismatch=='mismatch()':
        return('no')
    
    # any other mismatch indicates there is unclear support
    else:
        return('maybe')