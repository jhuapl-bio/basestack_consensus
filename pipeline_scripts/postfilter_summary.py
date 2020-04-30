#!/usr/bin/env python

import os
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
pd.options.mode.chained_assignment = None


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
    

def status_by_flags(flagstr,in_consensus,allele_freq):
    """
    Determine the genome status for a sample from the combined flag string
    and whether or not it passes the complete genome threshold
    """
    
    # flag keywords are 'depth','MAF','new','NTC','mismatch','key','SB'
    # these are unique words to each type of flag
    
    # situations that automatically lead to maybe
    if ('depth' in flagstr) or ('NTC' in flagstr):
        return('Maybe')
        
    # situations that can lead to maybe or yes*
    if 'mismatch' in flagstr:
        if ('mismatch(s)' in flagstr) and (allele_freq<0.2) and (in_consensus==False):
            return('Yes*')
        elif ('mismatch(s)' in flagstr) and ('SB' in flagstr) and (in_consensus==False):
            return('Yes*')
        else:
            return('Maybe')
    
    # accound for the unlikely case where there is no mismatch but there is strand bias
    if 'SB' in flagstr:
        return('Maybe')
    
    # situations that lead to yes*
    # if we have gotten to this point, then the only possible flags are 'MAF','key','new'
    if ('MAF' in flagstr) or ('key' in flagstr) or ('new' in flagstr):
        return('Yes*')
    
    # if we make it this far there are no flags
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
    num_flagged = []
    
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
            
            # replace '.' empty values with np.nan
            var = var.replace('.',np.nan)
            
            for pos in var.pos:
                tmp = var[var.pos==pos]
                
                # get all the flags (for status determination)
                # and get only flags to print
                
                # replace flags with abbreviated versions
                for colname in ['depth_flag','maf_flag','ntc_flag','new_flag','sb_flag','key_flag']:
                    if not pd.isna(tmp[colname].values[0]):
                        tmp.loc[tmp[colname].str.contains('depth'), colname] = 'depth'
                        tmp.loc[tmp[colname].str.contains('MAF'), colname] = 'MAF'
                        tmp.loc[tmp[colname].str.contains('NTC'), colname] = 'NTC'
                        tmp.loc[tmp[colname].str.contains('nextstrain'), colname] = 'new'
                        tmp.loc[tmp[colname].str.contains('key'), colname] = 'key'
                        tmp.loc[tmp[colname].str.contains('strand bias'), colname] = 'SB'
                        
                
                allflag = tmp[['depth_flag','maf_flag','ntc_flag','new_flag','vc_flag','sb_flag','key_flag']].apply(lambda x: ', '.join(str(pos)+':'+x.dropna()), axis=1).values[0]
                stat.append(status_by_flags(allflag, tmp.in_consensus.values[0], tmp.allele_freq.values[0]))
                
                if tmp.unambig.values[0] == True:
                    if not allflag=='':
                        printflag.append(allflag)
                    if pd.isna(tmp.vc_flag.values[0]):
                        snp.append(''.join([tmp.ref.values[0],str(pos),tmp.alt.values[0]]))
                    else:
                        # values that mean this position is not actually a snp in the consensus
                        no_snp = ['mismatch(m)','mismatch(s)','mismatch(m+s)']
                        if any(mismatch_snp in tmp.vc_flag.values[0] for mismatch_snp in no_snp)==False:
                            snp.append(''.join([tmp.ref.values[0],str(pos),tmp.alt.values[0]]))
                            
                else:
                    if not pd.isna(tmp.vc_flag.values[0]):
                        printflag.append(str(pos)+':'+tmp.vc_flag.values[0])
                    if not pd.isna(tmp.depth_flag.values[0]):
                        printflag.append(str(pos)+':'+tmp.depth_flag.values[0])
                    if not pd.isna(tmp.key_flag.values[0]):
                        printflag.append(str(pos)+':'+tmp.key_flag.values[0])
                    snp.append(''.join([tmp.ref.values[0],str(pos),'N']))
            
            # get the percent coverage for this sample
            complete,nt = is_complete(rundir,prefix)
            coverage.append('{:.2%}'.format(nt/29903.0))
            
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
            
            # make the flags look nice
            if not printflag:
                flags.append('')
                num_flagged.append(0)
                snps.append(', '.join(snp))
            else:
                # join the flags
                flagstr = ', '.join(printflag)
                
                # format the flags
                flaglist=[item.split(":") for item in flagstr.split(", ")]
                flagdict={}
                
                for key,val in flaglist:
                    flagdict.setdefault(key, []).append(val)
                
                flagstr='; '.join("{!s}={!r}".format(key,val) for (key,val) in flagdict.items())
                flagstr=flagstr.replace('\'','')
                
                num_flagged.append(len(flagdict))
                
                if not complete:
                    flags.append('see full output')
                    snps.append('see full output')
                else:
                    # add the flags
                    flags.append(flagstr)
                    snps.append(', '.join(snp))
            
            # join this dataframe to all the others
            var['sample']=samplename
            alldata = pd.concat([alldata,var],ignore_index=True)
    
    # make the dataframe
    df = pd.DataFrame({'Sample':samplenames,'Coverage':coverage,'Variants':snps,'Flags':flags,'Flagged Positions':num_flagged,'Alpha':alpha,'Status':status})
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
