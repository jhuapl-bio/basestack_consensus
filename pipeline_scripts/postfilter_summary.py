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
    

def status_by_flags(flagstr,consensus_var):
    """
    Determine the genome status for a sample from the combined flag string
    and whether or not it passes the complete genome threshold
    """
    
    # flag keywords are 'depth','MAF','nextstrain,'NTC','mistmatch'
    # these are unique words to each type of flag
    # flags other than MAF only in called bases put a genome in the maybe category
    # mismatch flags in uncalled bases put a genome in the yes* category
    
    if consensus_var == True:
        if ('nextstrain' in flagstr) or ('depth' in flagstr) or ('NTC' in flagstr) or ('mismatch' in flagstr):
            return('Maybe')
        elif ('MAF' in flagstr):
            return('Yes*')
        else:
            return('Yes')
    else:
        if 'mismatch' in flagstr:
            return('Yes*')
        else:
            return('Yes')


def generate_postfilter_summary(rundir):
    """
    Generate a summary table of postfilter results
    """
    
    alldata = pd.DataFrame()
    
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
            samplename = prefix.split('_')[0]           
            samplenames.append(samplename)
            
            # get the flags for this sample
            var = pd.read_csv(vardata,sep='\t')
            printflag = []
            snp = []
            stat = []
            
            for pos in var.pos:
                tmp = var[var.pos==pos]
                
                # get all the flags (for status determination)
                # and get only flags to print
                allflag = tmp[['depth_flag','maf_flag','ntc_flag','new_flag','vc_flag']].apply(lambda x: ', '.join(str(pos)+':'+x.dropna()), axis=1).values[0]
                stat.append(status_by_flags(allflag, tmp.consensus_var.values[0]))
                
                if tmp.consensus_var.values[0] == True:
                    if not allflag=='':
                        printflag.append(allflag)
                    if pd.isna(tmp.vc_flag.values[0]):
                        snp.append(''.join([tmp.ref.values[0],str(pos),tmp.alt.values[0]]))
                    else:
                        if not 'samtools only' in tmp.vc_flag.values[0]:
                            snp.append(''.join([tmp.ref.values[0],str(pos),tmp.alt.values[0]]))
                            # NEED TO UPDATE WITH NANOPOLISH
                            
                else:
                    if not pd.isna(tmp.vc_flag.values[0]):
                        printflag.append(str(pos)+':'+tmp.vc_flag.values[0])
                    if not pd.isna(tmp.depth_flag.values[0]):
                        printflag.append(str(pos)+':'+tmp.depth_flag.values[0])
                    snp.append(''.join([tmp.ref.values[0],str(pos),'N']))
            
            # get the percent coverage for this sample
            complete,unambig = is_complete(rundir,prefix)
            coverage.append('{:.2%}'.format(unambig/29903.0))
            
            # get the status for this genome
            if not complete:
                status.append('No')
            elif 'Maybe' in stat:
                status.append('Maybe')
            elif 'Yes*' in stat:
                status.append('Yes*')
            else:
                status.append('Yes')
            
            # include sample in alpha group if it passes completeness threshold
            if complete:
                alpha.append('Yes')
            else:
                alpha.append('No')
            
            # hide the flags and snps for discarded genomes
            if complete:
                # join the flags
                flagstr = ', '.join(printflag)
                flags.append(flagstr)
            
                # join the snps
                snps.append(', '.join(snp))
            else:
                flags.append('see full output')
                snps.append('see full output')
            
            
            
            # join this dataframe to all the others
            var['sample']=samplename
            alldata = pd.concat([alldata,var],ignore_index=True)
            
    
    # make the dataframe
    df = pd.DataFrame({'Sample':samplenames,'Coverage':coverage,'Variants':snps,'Flags':flags,'Alpha':alpha,'Status':status})
    df.to_csv(os.path.join(rundir,'postfilt_summary.txt'),sep='\t',index=False)
    
    # output the large table
    alldata.to_csv(os.path.join(rundir,'postfilt_all.txt'),sep='\t',index=False)
    
def parse_arguments():
   parser = argparse.ArgumentParser()
   parser.add_argument('--rundir', type=str, help='path postfilter results for a particular run')
   args = parser.parse_args()
   return(args)


if __name__ == "__main__":
    
    args = parse_arguments()
    generate_postfilter_summary(args.rundir)
