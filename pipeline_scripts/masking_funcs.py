#!/usr/bin/env python

import numpy as np
import pandas as pd
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def add_ref_positions(align_array):
    """
    Function to add reference-based positions to alignment
    All computation will use reference based positions
    And indices for these positions will be searched in alignment
    
    Deletions do not affect the reference and therefore are not explicitly addressed
    """
    
    # add row of positions to alignment
    # we will use this to keep track of the reference-based position at all times
    pos = list(range(1,align_array.shape[1]+1))
    
    # detect positions with insertions
    ins_sites = list(np.where(align_array[0,:]=='-')[0])
    ins_sites = [x+1 for x in ins_sites] # convert to 1-indexing
    
    ins_counter = 0
    for ins in ins_sites:
                
        ins_adj = ins-ins_counter
        
        # fix reference position
        pos[pos.index(ins_adj):pos.index(ins_adj)]='-'
        
        ins_counter=ins_counter+1
    
    # add reference positions to alignment
    pos = pos[0:(len(pos)-ins_counter)]
    align_array = np.concatenate((align_array,np.reshape(pos,(1,np.array(pos).size))),axis=0)
    
    return(align_array)
    

def get_amp_sites(i,amp):
    """ 
    Function to return unique positions for each amplicon
    Wrapping in a function to more easily deal with first and last amplicons
    """
    
    # get number of amplicons
    last = amp.amplicon.values[-1]
    
    # deal with special first amplicon
    if i==1:
        amp_sites = range(amp['unique_start'][0],int(amp[amp.amplicon==i+1]['primer_f_stop'])+1)
    
    # deal with special last amplicon
    elif i==last:
        amp_sites = range(int(amp[amp.amplicon==i-1]['primer_r_start']),amp['primer_r_start'][last-1])
    
    # deal with all other amplicons
    else:
        prev = amp[amp.amplicon==i-1]
        subs = amp[amp.amplicon==i+1]
        amp_sites = range(int(prev['primer_r_start']),int(subs['primer_f_stop'])+1)
    
    # remember values are 1-indexed
    return(amp_sites)
    

def calculate_depth_threshold(ntc_cov,amp,call_depth_factor):
    """
    Determine the read depth threshold to use for calling non-ambiguous bases
    based on a negative control (NTC) included on the same sequencing run
    """
    
    # initialize list of average amplicon depths
    avg_amp_depths = []
    
    # loop through amplicons
    for i in amp.amplicon:
        
        # get list of depths at positions within this amplicon
        amp_sites = get_amp_sites(i,amp)
        avg_amp_depths.append(np.mean([ntc_cov[pos] for pos in amp_sites]))
    
    # return 95% quantile of this average amplicon list
    # multiplied by the call depth factor
    qvalue = np.quantile(avg_amp_depths,0.95,interpolation='lower')
    
    # round this number up to better reflect the amplicon depth mode
    qvalue = int(np.ceil(qvalue))
    return(call_depth_factor*qvalue)


def mask_failed_amplicons(align_array,cov,amp,depth_threshold):
    """ 
    Function to mask sites in the consensus genome that are in failed amplicons
    """
    
    # save list of failed amplicons and masked sites
    failed_amplicons = []
    ampmask = []
    
    # convert coverage to dictionary for faster indexing
    cov_dict = pd.Series(cov.depth.values,index=cov.pos).to_dict()
    
    # loop through amplicons
    for i in amp.amplicon:
                
        # get list of depths at positions within this amplicon
        # amplicon sites and depth positions are both reference-based
        amp_sites = get_amp_sites(i,amp)
        depths = [cov_dict[pos] for pos in amp_sites]
        
        # calculate metric to be used to assess amplicon failure
        # if any position in amplicon is less than depth threshold
        if sum(d < depth_threshold for d in depths) > 0:
            
            # if consecutive amplicons failed, mask region between them
            if (i-1) in failed_amplicons:
                amp_sites_prev = get_amp_sites(i-1,amp)
                amp_sites = range(min(amp_sites_prev),max(amp_sites)+1)
            
            # now mask all the sites in amp_sites
            # this is where we need to find the indices of the ref positions in the array
            left_mask = np.where(align_array[2,:]==str(min(amp_sites)))[0][0]
            right_mask = np.where(align_array[2,:]==str(max(amp_sites)))[0][0]
            mask_sites = range(left_mask,right_mask+1)
            align_array[1,mask_sites] = 'N'
            
            # add amplicon to list of masked amplicons
            failed_amplicons.append(i)
            ampmask=ampmask+list(amp_sites)
        
    # remove any duplicates from ampmask list
    # caused by remasking of overlapping regions
    ampmask = sorted(list(set(ampmask)))
    return(align_array,failed_amplicons,ampmask)


def mask_consensus_sites(align_array,cov,depth_threshold,amp,outdir,samplename):
    """
    Mask sites in the consensus genome based on the depth threshold
    calculated from the negative control on the same run
    
    This function could be replaced by adding the --depth parameter when calling artic_make_depth_mask
    Though artic_make_depth_threshold runs on non-primertrimmed bam files
    """
    cov_dict = pd.Series(cov.depth.values,index=cov.pos).to_dict()
    
    # all stored positions should be in reference space
    ambig=[] # store ambig positions coming out of ivar pipeline
    depthmask=[] # store all masked positions due to depth mask
    newmask=[] # store newly masked positions due to depth mask
        
    # save bases that are 'N' out of assembly
    # account for possible case differences
    ambig_idx = list(np.where((align_array[1,:]=='N') | (align_array[1,:]=='n'))[0])
    ambig = ambig + [int(x) for x in align_array[2,ambig_idx] if x!='-']
    ambig_below_threshold = [x for x in ambig if cov_dict[x]<depth_threshold]
    depthmask = depthmask + ambig_below_threshold
    
    # mask sites at beginning and end where amplicons do not cover
    # remember zero indexing of array
    first_primer_stop = amp[amp.amplicon==1]['primer_f_stop'].values[0]
    last_primer_start = (amp[amp.amplicon==amp.amplicon.values[-1]]['primer_r_start'].values[0])
    first_primer_stop_idx = np.where(align_array[2,:]==str(first_primer_stop))[0][0]
    last_primer_start_idx = np.where(align_array[2,:]==str(last_primer_start))[0][0]
    align_array[1,0:(first_primer_stop_idx+1)]='N'
    align_array[1,last_primer_start_idx:]='N'
    depthmask = depthmask + [int(x) for x in align_array[2,0:(first_primer_stop_idx+1)]]
    depthmask = depthmask + [int(x) for x in align_array[2,last_primer_start_idx:] if x!='-']
    
    # mask bases if coverage is below depth threshold
    # keep track of which were already 'N'
    current_mask = np.where((align_array[1,:]=='N') | (align_array[1,:]=='n'))[0]
    current_mask_pos = [int(x) for x in align_array[2,current_mask] if x!='-']
    
    pos_below_threshold = cov[cov.depth<depth_threshold].pos.values
    pos_below_threshold = [str(x) for x in pos_below_threshold]
    idx_below_threshold = np.nonzero(np.in1d(align_array[2,:],pos_below_threshold))[0]
    align_array[1,idx_below_threshold]='N'
    
    # convert string back into int for storage
    pos_below_threshold = [int(x) for x in pos_below_threshold]
    newmask = newmask + sorted(list(set(pos_below_threshold)-set(current_mask_pos)))
    depthmask = depthmask + newmask # capture Ns outside of amplicons
    
    # mask all bases within failed amplicons
    masked_alignment,failed_amplicons,ampmask = mask_failed_amplicons(align_array,cov,amp,depth_threshold)
    
    # get list of all depth/amplicon-masked sites
    depthmask = sorted(list(set(depthmask)))
    allmask = sorted(list(set(depthmask + ampmask)))
    # get list of all ambiguous sites not masked for coverage
    nomask = [x for x in ambig if x not in allmask]
    
    # output newly-masked bases to file
    d = dict(assembly_ambig=ambig,new_depth_mask=newmask,all_depth_mask=depthmask,amp_mask=ampmask,
             failed_amps=failed_amplicons,all_masked=allmask,non_masked_ambig=nomask)
    if any(a != [] for a in d.values()):
        masked_sites = pd.DataFrame(dict([ (k,pd.Series(v)) if v!=[] else (k,np.nan) for k,v in d.items() ]))
    else:
        masked_sites = pd.DataFrame()
    masked_sites['depth_thresh']=depth_threshold
    
    filename=os.path.join(outdir,samplename+'.new_masked_sites.txt')
    masked_sites.to_csv(filename,sep='\t',index=False)
    
    # save new consensus sequence and output to file
    masked_consensus = [base for base in list(masked_alignment[1,:]) if base != '-']
    masked_consensus = ''.join(masked_consensus).upper()

    new_record = SeqRecord(Seq(masked_consensus),id=samplename,description="")
    filepath = os.path.join(outdir,samplename+'.mask.fasta')
    SeqIO.write(new_record,filepath,"fasta")
    
    # return masked alignment and dataframe of masked sites
    return(masked_alignment,masked_sites)