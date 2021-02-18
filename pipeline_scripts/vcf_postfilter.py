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
from Bio import AlignIO

from variant_status import status_by_case
import variant_flags as fl
import masking_funcs as msk

def parse_indel_data(data,info,masked_align,var_idx):
    """
    Function to create final dataframe row for indels
    Involves checking for variant caller matches and presence in consensus
    But no reporting of frequencies is performed
    """
    
    ### check if variant in consensus and add indel flag
    
    # assume not in consensus unless proven otherwise
    data['in_consensus']=False
    data['consensus_base']=data['ref']
    case='IND'
    description='indel detected'
    status='Pass'
    
    # extra check: confirm that the reference genome matches the VCF
    span = len(data['ref'])
    ins_len = abs(len(data['ref'])-len(data['alt']))
    assert ''.join([x for x in masked_align[0,var_idx:(var_idx+span)]]).upper()==data['ref']
    
    if len(data['alt'])==1: # deletion  
    
        data['indel_flag']='deletion'
        
        # determine if the sample genome has the deletion
        region_seq = ''.join(masked_align[1,(var_idx+1):(var_idx+span)])
        if region_seq==''.join(['-']*(span-1)): # deletion present
            data['in_consensus']=True
            data['consensus_base']='del'
            if (ins_len)%3==0:
                case='IND-C'; description='indel detected in consensus'; status='Check'
            else:
                case='IND-FS'; description='frameshifting indel detected in consensus'; status='Alarm'
            # note that there could also be a variant at this position
    
    elif len(data['ref'])==1: # insertion
    
        data['indel_flag']='insertion'
    
        # determine if the sample genome has the insertion
        inserted_bases = ''.join(masked_align[0,(var_idx+1):(var_idx+len(data['alt']))])
        if inserted_bases==''.join(['-']*(len(data['alt'])-1)): # insertion present
            data['in_consensus']=True
            data['consensus_base']='ins'
            if (ins_len)%3==0:
                case='IND-C'; description='indel detected in consensus'; status='Check'
            else:
                case='IND-FS'; description='frameshifting indel detected in consensus'; status='Alarm'
    
    else:
        sys.exit('position %s contains an undefined indel' % data['pos'])
    
    
    ### add other flags
    data['vc_flag'] = fl.variant_caller_mismatch(info['SUPP_VEC'])
    data['case']=case
    data['description']=description
    data['status']=status
    
    ### return data dict
    return(data)


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


def check_ambiguous_positions(cov_dict,masked_align,df,depth_threshold,masked_sites,key_vars,chrom,runid,samplename,barcode):
    """ 
    Function that returns a dataframe of positions that need to be curated
    These include:
        - positions > depth threshold
        - positions < depth threshold if they are in the list of key positions
    """
    
    # load in key variants
    key_snps = pd.read_csv(key_vars,sep='\t',header=None,names=['pos'])
    key_snps = list(key_snps.pos.values)
    
    # initialize data frame of all positions to be added
    ambig_df = pd.DataFrame()
    
    # get variant list
    variants = list(df.pos.values)
    
    # add positions to variant list that are covered by deletions
    newvars = []
    if len(variants)>0: # prevent issues due to lack of any variants in this sample
        for pos in variants:
            reflen = len(df[df.pos==pos]['ref'].values[0])
            if reflen>1:
                for i in range(1,reflen):
                    newvars.append(pos+i)
    variants = variants+newvars
    
    # find positions that are ambiguous and not in variant list
    ambig_pos = masked_align[2,np.where(masked_align[1,:]=='N')][0]
    ambig_pos = [int(x) for x in ambig_pos if x!='-'] # convert to int for comparison with variants
    ambig_pos_novar = [x for x in ambig_pos if x not in variants]
    
    # loop through ambiguous positions and note which are key positions
    for pos in ambig_pos_novar:
        
        # check depth at this position
        # include positions above depth threshold and those below if in a key position
        if (cov_dict[pos]>=depth_threshold and pos not in masked_sites.amp_mask.values) or (pos in key_snps):
            
            # determine case number and description
            if cov_dict[pos]<depth_threshold or (pos in masked_sites.amp_mask.values):
                assert pos in key_snps # only include low depth if at key position
                case='UNK-K'
                description='ambiguous base at key position'
            
            else:
                # we are only here if there is no variant at this site
                # so we can assume this is likely some sort of variant calling issue
                case='UNK-V'
                description='ambiguous base due to variant calling issue'
            
            # set up the data
            data={}; data['run_id']=runid; data['sample']=samplename; data['barcode']=barcode; data['chrom']=chrom
            data['pos']= pos; data['alt']='N'
            data['ref']= masked_align[0,np.where(masked_align[2,:]==str(pos))[0][0]].upper()
            data['consensus_base']='N'; data['in_consensus']=False; data['unambig']=False
            data['ont_depth']=cov_dict[pos]; data['ont_depth_thresh']=depth_threshold
            data['status']='Pass'; data['case']=case; data['description']=description
            
            data = pd.DataFrame([data], columns=data.keys())
            ambig_df = pd.concat([ambig_df,data],ignore_index=True)
            
    return(ambig_df)


def make_final_fasta(masked_align,samplename,unambig_thresh,outdir):
    
    """
    Get final fasta from alignment
    Assign as complete or partial given unambiguous base threshold
    """
    
    # create final consensus genome from alignment
    # do not include deletions as gaps
    masked_consensus = [base for base in list(masked_align[1,:]) if base != '-']
    
    # count number of ambiguous bases in sequence
    ambig = masked_consensus.count('N')
    unambig = len(masked_consensus)-ambig
    
    # save new consensus sequence and output to file
    cons = ''.join(masked_consensus).upper()
    new_record = SeqRecord(Seq(cons),id=samplename,description="")
    
    # save new output file depending on unambig threshold
    if unambig > unambig_thresh:
        filepath = os.path.join(outdir,samplename+'.complete.fasta')
    else:
        filepath = os.path.join(outdir,samplename+'.partial.fasta')
        
    SeqIO.write(new_record,filepath,"fasta")


def parse_arguments():
   parser = argparse.ArgumentParser()
   parser.add_argument('--vcffile', type=str, help='path to vcf file of sample')
   parser.add_argument('--depthfile', type=str, help='path to depth of bam file of sample')
   parser.add_argument('--aln-to-ref', type=str, help='path to alignment of consensus fasta to ref')
   parser.add_argument('--ntc-bamfile', type=str, help='path to bam file of negative control')
   parser.add_argument('--ntc-depthfile', type=str, help='path to depth file of negative control')
   parser.add_argument('--global-vars', type=str, help='path to tsv containing observed variants in global data')
   parser.add_argument('--key-vars', type=str, help='path to txt file containing a list of key positions')
   parser.add_argument('--homopolymers', type=str, help='path to txt file containing a list of homopolymer positions')
   parser.add_argument('--amplicons', type=str, help='path to file containing amplicon ranges')
   parser.add_argument('--case-defs', type=str, help='path to csv containing case numbers and definitions')
   
   parser.add_argument('--outdir', '-o', type=str, help='directory name to write output files to')
   parser.add_argument('--samplename', type=str, default='sample', help='prefix of all saved output files')
   
   parser.add_argument('--default-depth-threshold', type=int, default=20, help='default depth threshold')
   parser.add_argument('--coverage-flag', type=int, default=20, help='flag variants with depth within this percentage of threshold')
   parser.add_argument('--maf-flag', type=int, default=25, help='flag variants with minor allele frequency with at least this value')
   parser.add_argument('--call-depth-factor', type=int, default=2, help='factor by which depth must exceed 95th percentile of NTC depth to call a base')
   parser.add_argument('--snp-depth-factor', type=int, default=5, help='factor by which depth must exceed NTC depth to call a variant seen in the NTC at that position')
   parser.add_argument('--unambig-threshold', type=int, default=27000, help='number of unambiguous bases required in final genome')
   parser.add_argument('--ns-snp-threshold', type=int, default=1, help='fraction of published samples with a particular snp needed to count it as previously seen')
   parser.add_argument('--strand-threshold', type=int, default=10, help='minimum minor allele frequency on each strand required for unbiased call')
   
   args = parser.parse_args()
   return(args)


def main():
    
    args = parse_arguments()
    
    ### MASK AMPLICONS ###
    
    # load start and end positions of amplicons used for sequencing
    amp = pd.read_csv(args.amplicons,sep='\t')
    
    # use NTC to calculate depth threshold if present
    if not args.ntc_depthfile=="None":
        ntc_cov = pd.read_csv(args.ntc_depthfile,sep='\t',header=None,names=['chrom','pos','depth'])
        ntc_cov = pd.Series(ntc_cov.depth.values,index=ntc_cov.pos).to_dict() # depth across genome
        depth_threshold = max(args.default_depth_threshold,
                              msk.calculate_depth_threshold(ntc_cov,amp,args.call_depth_factor))
    else:
        depth_threshold = args.default_depth_threshold
    
    # load alignment file
    align = AlignIO.read(args.aln_to_ref,"fasta") # consensus aligned to reference
    align_array = np.array([list(rec) for rec in align], dtype=str, order="F")
    
    # load depth file
    cov = pd.read_csv(args.depthfile,sep='\t',header=None,names=['chrom','pos','depth'])
    del cov['chrom']
    cov_dict = pd.Series(cov.depth.values,index=cov.pos).to_dict()
    
    # adjust amplicon and depth files for insertions
    align_array = msk.add_ref_positions(align_array)
    
    # mask amplicons based on depth threshold
    masked_align,masked_sites = msk.mask_consensus_sites(align_array,cov,depth_threshold,amp,args.outdir,args.samplename)
    
    ### FLAG VARIANTS ###
    
    # load variant vcf
    vcf_sample = pd.read_csv(args.vcffile,sep='\t',skiprows=2)
    vcf_sample = vcf_sample[['#CHROM','POS','REF','ALT','INFO']]
    vcf_sample.columns = ['chrom','pos','ref','alt','info']

    # set up the dataframe to store results
    df = pd.DataFrame(
        columns=['run_id','sample','barcode','chrom','pos','ref','alt','consensus_base','case','description',
                 'status','homopolymer','in_consensus','unambig',
                 'ont_depth','illumina_depth','ont_depth_thresh','illumina_depth_thresh',
                 'ont_AF','illumina_AF','ont_alleles','illumina_alleles','ont_strand_counts',
                 'medaka_qual','nanopolish_qual','illumina_support',
                 'depth_flag','ntc_flag','indel_flag','vc_flag','mixed_flag','maf_flag',
                 'sb_flag','key_flag','new_flag'])

    # loop through all positions
    for pos in vcf_sample.pos:
        
        tmp = vcf_sample[vcf_sample.pos==pos]
        
        # get the alignment array index of the position
        var_idx = np.where(align_array[2,:]==str(pos))[0][0]
        
        # start a dictionary to store data for this sample
        data = tmp.to_dict('records')[0]
        
        # add some basic data to this record
        data['sample'] = args.samplename
        data['run_id'] = (args.depthfile).split('/')[-4]
        data['barcode'] = (args.depthfile).split('/')[-1].split('.')[0].split('_')[1]
        
        # store information from info column and remove it from dictionary
        info = dict(item.split("=") for item in data['info'].split(";"))
        del data['info']

        # determine if there is illumina data for this sample
        if len(info['SUPP_VEC'])==6:
            illumina=True
        else:
            illumina=False

        # ignore positions with no nanopore variant calls
        if illumina:
            if info['SUPP_VEC'][:3]=='000':
                continue
        
        ### SPECIAL TREAMENT FOR INSERTIONS AND DELETIONS ###
        
        # allow indels but do not worry about allele frequencies
        # check for variant caller mismatch and presence in consensus
        # but no other flags
        if len(data['alt']) != len(data['ref']):
            data = parse_indel_data(data,info,masked_align,var_idx)
            # ignore samtools-only indels
            # do not add them to final vcf
            if data['vc_flag']=="mismatch(s)":
                continue # without adding this row to the dataframe
            else:
                data = pd.DataFrame([data], columns=data.keys())
                df = pd.concat([df,data],ignore_index=True,sort=False)
                continue # after adding indel to dataframe

        ### ASSEMBLE DATA FOR SINGLE NUCLEOTIDE POLYMORPHISMS
        
        # get ont read depth and pileup at this read position
        depth,alt_allele_freq,allele_string = parse_allele_counts(info, data['alt'], 'ont')

        # ignore this position if the depth is too low
        if (cov_dict[pos] < depth_threshold) or (pos in masked_sites.amp_mask.values):
            continue

        # get illumina read depth and pileup if applicable
        if illumina:
            illumina_depth,illumina_alt_allele_freq,illumina_allele_string = parse_allele_counts(info, data['alt'], 'illumina')
        
        # add basic data to this record
        
        data['sample'] = (args.depthfile).split('/')[-1].split('.')[0].split('_')[0]
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
        data['new_flag'] = fl.snp_in_nextstrain(pos, data['ref'], data['alt'], args.global_vars, args.ns_snp_threshold)
        data['vc_flag'] = fl.variant_caller_mismatch(info['SUPP_VEC'])
        data['sb_flag'],data['ont_strand_counts'] = fl.strand_bias_detected(info,data['alt'],args.strand_threshold)
        
        data['ntc_flag'] = fl.allele_in_ntc(pos, data['alt'], depth, args.ntc_bamfile, args.snp_depth_factor)
        data['homopolymer'] = fl.in_homopolymer_region(pos,args.homopolymers)
        data['maf_flag'],data['mixed_flag'] = fl.minor_allele_freq(depth, alt_allele_freq, args.maf_flag)
        data['illumina_support'] = [fl.get_illumina_support(illumina_alt_allele_freq,info['SUPP_VEC'],args.maf_flag) if illumina and illumina_depth>=20 else np.nan][0]
        
        # update quality value if nanopolish did not call a variant
        # the QUAL field represents the medaka quality score if nanopolish didn't call a variant
        # not an issue with samtools only calls because then there is no quality score
        if (data['vc_flag']=='mismatch(m)') or (data['vc_flag']=='mismatch(m+s)'):
            data['nanopolish_qual']='.'

        # modify consensus genome based on ntc flag
        if data['ntc_flag']=='allele in NTC':
            masked_align[1,var_idx]='N'
            
        # add a flag if a key position is ambiguous in the consensus
        # must happen after any masking that occurs due to NTC flags
        data['key_flag'] = fl.ambig_in_key_position(pos,args.key_vars,masked_align,var_idx)
        
        # mark which positions are unambiguous in the consensus
        data['unambig'] = [False if (masked_align[1,var_idx]=='N') | (masked_align[1,var_idx]=='n') else True][0]
        
        # mark which variants are actually in the consensus genome
        IUPAC=['R','Y','S','W','K','M','B','D','H','V']
        if (data['unambig']==True) and masked_align[1,var_idx].upper()==data['alt']:
            data['in_consensus'] = True
        elif (data['unambig']==True) and np.isin(masked_align[1,var_idx].upper(),IUPAC):
            data['in_consensus'] = 'IUPAC'
            data['in_consensus'] = 'IUPAC'
        else:
            data['in_consensus'] = False
        
        # report the base called in the consensus genome
        data['consensus_base'] = masked_align[1,var_idx].upper()
        
        # double check we are reporting the consensus base correctly
        if data['in_consensus']==True and data['unambig']==True:
            assert data['consensus_base']==data['alt']
        elif data['in_consensus']=='IUPAC' and data['unambig']==True:
            assert data['consensus_base'] in IUPAC
        elif data['in_consensus']==False and data['unambig']==True:
            assert data['consensus_base']==data['ref']
        elif data['unambig']==False:
            assert data['consensus_base']=='N'
        
        # after adding all the flags
        # get the status of this position
        data['case'],data['description'],data['status'] = status_by_case(data, args.case_defs, args.maf_flag)
        
        # add this record to the final dataframe
        data = pd.DataFrame([data], columns=data.keys())
        df = pd.concat([df,data],ignore_index=True,sort=False)
    
    
    ### ADD POSITIONS AT KEY SNPS IF NOT IN VARIANT LIST ###
    
    runid = (args.depthfile).split('/')[-4]
    barcode = (args.depthfile).split('/')[-1].split('.')[0].split('_')[1]
    ambig_df = check_ambiguous_positions(cov_dict,masked_align,df,depth_threshold,masked_sites,args.key_vars,align[0].id,runid,args.samplename,barcode)
    df = pd.concat([df,ambig_df],ignore_index=True,sort=False)
    
    # replace nan values for ease of text parsing
    df = df.replace(np.nan,'.')
    
    # save output file
    filepath = os.path.join(args.outdir,args.samplename+'.variant_data.txt')
    df.to_csv(filepath,sep='\t',index=False)
    
    # then make the final consensus
    make_final_fasta(masked_align,args.samplename,args.unambig_threshold,args.outdir)


if __name__ == "__main__":
    main()
