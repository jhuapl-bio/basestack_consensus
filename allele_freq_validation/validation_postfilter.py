#!/usr/bin/env python

import os
import sys
import argparse
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from variant_status import status_by_case
import variant_flags as fl


def calculate_depth_threshold(ntc_depthfile,call_depth_factor):
    """
    Determine the read depth threshold to use for calling non-ambiguous bases
    based on a negative control (NTC) included on the same sequencing run
    """
    
    # get negative control depth
    cov = pd.read_csv(ntc_depthfile,sep='\t',header=None,names=['chrom','pos','depth'])
    median_depth = cov.depth.median()
    return(call_depth_factor*median_depth)


def mask_consensus_sites(consensus,depthfile,depth_threshold,outdir,prefix):
    """
    Mask sites in the consensus genome based on the depth threshold
    calculated from the negative control on the same run
    
    This function could be replaced by adding the --depth parameter when calling artic_make_depth_mask

    """
    
    # get depth across genome
    cov = pd.read_csv(depthfile,sep='\t',header=None,names=['chrom','pos','depth'])
    cov = pd.Series(cov.depth.values,index=cov.pos).to_dict()
    
    # load current consensus sequence
    cons = list(SeqIO.parse(open(consensus),"fasta"))[0]
    cons = list(cons.seq.upper())
    
    assert len(cons)==29903
    
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
            ambig.append(pos+1)
            continue
        else:
            rd = [cov[pos+1] if (pos+1) in cov.keys() else 0][0]
            if rd < depth_threshold:
                cons[pos] = 'N'
                newmask.append(pos+1)
    
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


def parse_allele_counts(info,alt,method):
    
    # use a different info field depending on the method
    if method=='ont':
        pos_alleles = info['POSITIVE_STRAND_FREQUENCIES']
        neg_alleles = info['NEGATIVE_STRAND_FREQUENCIES']
    elif method=='illumina':
        pos_alleles = info['ILLUMINA_POSITIVE_STRAND_FREQUENCIES']
        neg_alleles = info['ILLUMINA_NEGATIVE_STRAND_FREQUENCIES']
    else:
        sys.exit('not a valid sequencing method')
    
    pos_alleles = [int(x) for x in pos_alleles.split(',')]
    neg_alleles = [int(x) for x in neg_alleles.split(',')]
    
    # calculate the total depth
    depth = sum(pos_alleles) + sum(neg_alleles)
    
    # get the allele string across both strands
    alleles = [pos_alleles[i]+neg_alleles[i] for i in range(len(pos_alleles))]
    allele_string = 'A:%d:C:%d:G:%d:T:%d:N:%d:O:%d' % (alleles[0],alleles[1],alleles[2],alleles[3],alleles[4],alleles[5])
    
    # get the alternate allele frequency
    idx = ['A','C','G','T','N','O'].index(alt)
    alt_allele_freq = [float(alleles[idx]/depth) if depth>0 else 0.0][0]
    
    return(depth,alt_allele_freq,allele_string)


def add_key_ambiguous_positions(chrom,variants,cons,depth_threshold,vcf_nextstrain,mpileup):
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
    
    # load mpileup file to get depth at this position
    names=['chrom', 'pos', 'ref', 'depth', 'pileup','qual']
    cov = pd.read_csv(mpileup,sep='\t',names=names)
    cov = pd.Series(cov.depth.values,index=cov.pos).to_dict()
    
    # create a dataframe in which to store any new values
    df = pd.DataFrame()
    
    # loop through important snps
    for pos in key_snps:
        if (pos not in variants) and (cons[pos-1]=='N'):
            data={}
            data['chrom']=chrom
            data['pos']=pos
            data['ref']=ns_snps[ns_snps.POS==pos].REF.values[0]
            data['alt']='N'
            data['consensus_base']='N'
            data['in_consensus']=False
            data['unambig']=False
            data['ont_depth']=[cov[pos] if pos in cov.keys() else 0][0]
            data['ont_depth_thresh']=depth_threshold
            data['key_flag']='ambig in key position'
            data['case']=19
            data['status']='Yes*'
            data['description']=data['key_flag']
            
            data = pd.DataFrame([data], columns=data.keys())
            df = pd.concat([df,data],ignore_index=True)
    
    return(df)


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
   parser.add_argument('--mpileup', type=str, help='path to mpileup of bam file of sample')
   parser.add_argument('--depthfile', type=str, help='path to depth of bam file of sample')
   parser.add_argument('--consensus', type=str, help='path to fasta file of sample consensus sequence')
   parser.add_argument('--ntc-bamfile', type=str, help='path to bam file of negative control')
   parser.add_argument('--ntc-depthfile', type=str, help='path to depth file of negative control')
   parser.add_argument('--vcf-nextstrain', type=str, help='path to vcf containing all nextstrain snps')
   parser.add_argument('--case-defs', type=str, help='path to csv containing case numbers and definitions')
   
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
    
    depth_threshold = max(20,calculate_depth_threshold(args.ntc_depthfile, args.call_depth_factor))
    mask_cons = mask_consensus_sites(args.consensus, args.depthfile, depth_threshold, args.outdir, args.prefix)
    
    # read vcf as text file
    vcf_sample = pd.read_csv(args.vcffile,sep='\t',skiprows=2)
    vcf_sample = vcf_sample[['#CHROM','POS','REF','ALT','QUAL','INFO']]
    vcf_sample.columns = ['chrom','pos','ref','alt','nanopolish_qual','info']
    
    # load current consensus sequence
    cons = list(SeqIO.parse(open(mask_cons),"fasta"))[0]
    cons = list(cons.seq.upper())
    
    # get the depth used for masking
    cov = pd.read_csv(args.depthfile,sep='\t',header=None,names=['chrom','pos','depth'])
    cov = pd.Series(cov.depth.values,index=cov.pos).to_dict()
    
    # set up the dataframe to store results
    df = pd.DataFrame(
        columns=['chrom','pos','ref','alt','consensus_base','case','description','status','homopolymer','in_consensus','unambig',
                 'ont_depth','illumina_depth','ont_depth_thresh','illumina_depth_thresh','ont_AF','illumina_AF','ont_alleles','illumina_alleles',
                 'ont_strand_counts','medaka_qual','nanopolish_qual','illumina_support',
                 'depth_flag','ntc_flag','vc_flag','mixed_flag','maf_flag','sb_flag','key_flag','new_flag'])
        
    # loop through all positions
    for pos in vcf_sample.pos:
        
        tmp = vcf_sample[vcf_sample.pos==pos]
        
        # ignore indels
        # print a warning if indels are found in the input vcf
        if len(tmp.ref.values[0]) != len(tmp.alt.values[0]):
            #print('Warning: indel found at position %d - indels are ignored!' % pos)
            continue
        
        # start a dictionary to store data for this sample
        data = tmp.to_dict('records')[0]
        
        # store information from info column and remove it from dictionary
        info = dict(item.split("=") for item in data['info'].split(";"))
        del data['info']
        
        # determine if there is illumina data for this sample
        if len(info['SUPP_VEC'])==6:
            illumina=True
        else:
            illumina=False
        
        # get ont read depth and pileup at this read position
        depth,alt_allele_freq,allele_string = parse_allele_counts(info, data['alt'], 'ont')
        
        # ignore this position if the depth is too low
        
        
        if cov[pos] < depth_threshold:
            continue
        
        # get illumina read depth and pileup if applicable
        if illumina:
            illumina_depth,illumina_alt_allele_freq,illumina_allele_string = parse_allele_counts(info, data['alt'], 'illumina')
        
        # add basic data to this record
        data['ont_depth'] = depth
        data['illumina_depth'] = [illumina_depth if illumina else np.nan][0]
        data['ont_depth_thresh'] = depth_threshold
        data['illumina_depth_thresh'] = [20 if illumina else np.nan][0]
        data['ont_alleles'] = allele_string
        data['ont_AF'] = alt_allele_freq
        data['medaka_qual'] = [info['pred_q'] if 'pred_q' in info.keys() else np.nan][0]
        data['illumina_alleles'] = [illumina_allele_string if illumina else np.nan][0]
        data['illumina_AF'] = [illumina_alt_allele_freq if illumina else np.nan][0]
        
        # add flags to this record
        data['depth_flag'] = fl.depth_near_threshold(depth,depth_threshold,args.coverage_flag)
        data['new_flag'] = fl.snp_in_nextstrain(pos, data['ref'], data['alt'], args.vcf_nextstrain, args.ns_snp_threshold)
        data['vc_flag'] = fl.variant_caller_mismatch(info['SUPP_VEC'])
        data['sb_flag'],data['ont_strand_counts'] = fl.strand_bias_detected(info,data['alt'],args.strand_threshold)
        data['ntc_flag'] = fl.allele_in_ntc(pos, data['alt'], depth, args.ntc_bamfile, args.snp_depth_factor)
        data['homopolymer'] = fl.in_homopolymer_region(pos)
        data['maf_flag'],data['mixed_flag'] = fl.minor_allele_freq(depth, alt_allele_freq, args.maf_flag, args.hold_flag)
        data['illumina_support'] = [fl.get_illumina_support(illumina_alt_allele_freq,info['SUPP_VEC'],args.maf_flag) if illumina and illumina_depth>=20 else np.nan][0]
        
        # modify consensus genome based on ntc flag
        # remember consensus is zero-indexed but we are dealing with 1-indexed positions
        if not pd.isna(data['ntc_flag']):
            cons[pos-1]='N'
            
        # add a flag if a key position is ambiguous in the consensus
        # must happen after any masking that occurs due to NTC flags
        data['key_flag'] = fl.ambig_in_key_position(pos, args.vcf_nextstrain, cons)
        
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
        data['case'],data['description'],data['status'] = status_by_case(data, args.case_defs)
        # add this record to the final dataframe
        data = pd.DataFrame([data], columns=data.keys())
        df = pd.concat([df,data],ignore_index=True,sort=False)
    
    # after looping through all positions
    # add positions at key snps if not in variant list
    df = pd.concat([df,add_key_ambiguous_positions(vcf_sample.chrom.values[0],list(df.pos.values), cons, depth_threshold, args.vcf_nextstrain, args.mpileup)],ignore_index=True,sort=False)
    df = df.replace(np.nan,'.')
    filepath = os.path.join(args.outdir,args.prefix+'.variant_data.txt')
    df.to_csv(filepath,sep='\t',index=False)
    
    # then make the final consensus
    make_final_fasta(cons,args.prefix,args.unambig_threshold,args.outdir)


if __name__ == "__main__":
    main()
