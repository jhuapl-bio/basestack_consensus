#!/usr/bin/env python

import os
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO


def is_complete(rundir,prefix):
    """
    Determine if a genome has passed the filters required to be marked as complete
    Return True/False and the number of unambiguous bases
    """
    
    genome_path = os.path.join(rundir,prefix+'.complete.fasta')
    if os.path.exists(genome_path):
        cons = list(SeqIO.parse(open(genome_path),"fasta"))[0]
        cons = list(cons.seq.upper())
        assert len(cons) == 29903
        ambig = cons.count('N')
        return(True,29903-ambig)
    else:
        partial_path = os.path.join(rundir,prefix+'.partial.fasta')
        assert os.path.exists(partial_path)
        cons = list(SeqIO.parse(open(partial_path),"fasta"))[0]
        cons = list(cons.seq.upper())
        assert len(cons) == 29903
        ambig = cons.count('N')
        return(False,29903-ambig)
    

def status_by_flags(flagstr,complete):
    """
    Determine the genome status for a sample from the combined flag string
    and whether or not it passes the complete genome threshold
    """
    
    # not passing the genome threshold is an automatic no
    if not complete:
        return('No')
    else:
        # flag keywords are 'depth','MAF','nextstrain,'NTC'
        # these are unique words to each type of flag
        # flags other than MAF only put a genome in the maybe category
        if ('nextstrain' in flagstr) or ('depth' in flagstr) or ('NTC' in flagstr):
            return('Maybe')
        elif ('MAF' in flagstr):
            return('Yes*')
        else:
            return('Yes')


def generate_postfilter_summary(rundir):
    """
    Generate a summary table of postfilter results
    """
    
    # intialize lists
    samplenames = []
    alpha = []
    flags = []
    coverage = []
    status = []
    snps = []
    
    # loop through the variant data files in the postfilter run directory
    for entry in os.scandir(rundir):
        if entry.path.endswith('variant_data.txt'):
            vardata = entry.path
            prefix = entry.path.split('.')[0].split('/')[-1]
            samplename = prefix.split('_')[0].split('-',1)[1]            
            samplenames.append(samplename)
            
            # get the flags for this sample
            var = pd.read_csv(vardata,sep='\t')
            flag = []
            snp = []
            
            for pos in var.pos:
                tmp = var[var.pos==pos]
                
                snp.append(''.join([tmp.ref.values[0],str(pos),tmp.alt.values[0]]))
                
                if not pd.isna(tmp.flags.values[0]):
                    flag.append(str(pos)+':'+str(tmp.flags.values[0]))
                else:
                    flag.append(np.nan)
            
            # join the flags
            flagstr = ', '.join([x for x in flag if not pd.isna(x)])
            flags.append(flagstr)
            
            # join the snps
            snps.append(', '.join(snp))
            
            # get the percent coverage for this sample
            complete,unambig = is_complete(rundir,prefix)
            coverage.append('{:.2%}'.format(unambig/29903.0))
            
            # get the flag status for this genome
            status.append(status_by_flags(flagstr,complete))
            
            # include sample in alpha group if it passes completeness threshold
            if complete:
                alpha.append('Yes')
            else:
                alpha.append('No')  
            
    
    # make the dataframe
    df = pd.DataFrame({'Sample':samplenames,'Coverage':coverage,'Variants':snps,'Flags':flags,'Alpha':alpha,'Status':status})
    df.to_csv(os.path.join(rundir,'postfilt_summary.txt'),sep='\t',index=False)
    
    
def parse_arguments():
   parser = argparse.ArgumentParser()
   parser.add_argument('--rundir', type=str, help='path postfilter results for a particular run')
   args = parser.parse_args()
   return(args)


if __name__ == "__main__":
    
    args = parse_arguments()
    generate_postfilter_summary(args.rundir)
