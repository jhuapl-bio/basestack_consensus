#!/usr/bin/env python

import pandas as pd
import sys
    
def case_by_flags(data,maf_flag):
    """
    Assign case numbers to specific variant situations
    """
    
    # convert the maf_flag into a fraction
    maf_flag = maf_flag/100.0
    
    # convert no NTC warning into NaN
    if data['ntc_flag']=='NTC=None':
        data['ntc_flag']=np.nan

    # situations that automatically lead to an alarm
    
    ## CASE 1
    if not pd.isna(data['depth_flag']):
        return('COV')
    
    ## CASE 2
    if not pd.isna(data['ntc_flag']):
        return('NTC')
    
    ## CASE 3
    if data['illumina_support']=='maybe':
        return('DIS-IL')
    
    ## CASE 4
    if not pd.isna(data['illumina_support']):
        if data['ont_AF']<maf_flag and data['illumina_AF']>(1-maf_flag):
            return('DIS')
        elif data['ont_AF']>(1-maf_flag) and data['illumina_AF']<maf_flag:
            return('DIS')
        
    ## CASE 5
    if data['ont_AF']<maf_flag and data['in_consensus']:
        return('LFC')
    
    ## CASE 6
    if data['ont_AF']>(1-maf_flag) and data['in_consensus']==False:
        return ('HFC')
    
    
    # situtions with mixed frequency variants
    
    ## CASE 7
    if data['illumina_support']=='mixed':
        return('MIX-IL')
    
    # at this point, illumina is yes/no/none
    
    if not pd.isna(data['mixed_flag']):
        
        ## CASE 8
        if data['in_consensus']==False and (not pd.isna(data['sb_flag'])) and (not data['illumina_support']=='yes'):
            return('SB')
        
        ## CASE 9
        elif data['illumina_support']=='yes' and data['homopolymer'] and data['in_consensus']:
            return('HP')
        
        ## CASE 10
        elif not pd.isna(data['illumina_support']):
            return('MIX-CS')
        
        ## CASE 11
        else:
            other_allele_freq = float(int(data['ont_alleles'].split(':')[11])/data['ont_depth'])
            if (other_allele_freq-0.02 <= (1-data['ont_AF']) <= other_allele_freq+0.02) and data['homopolymer'] and data['in_consensus']:
                return('MIX-HP')
        
            ## CASE 12
            else:
                return('MIX-MS')
    
    # at this point we know the frequency is either <maf_flag or >(1-maf_flag)
    # and that the in consensus status matches the high/low status
    
    # specific situations for variants not seen before
    if not pd.isna(data['new_flag']):
        
        ## CASE 13
        if data['illumina_support']=='yes' and data['in_consensus']:
            return('NEW-IL')
        
        ## CASE 14
        elif data['ont_AF']>0.9 and data['in_consensus']:
            return('NEW-HF')
        
        ## CASE 15
        elif data['in_consensus']:
            return('NEW-NI')
    
    # if there are no other worrisome flags
    # ignore low frequency variants and accept high frequency variants
    # remember we already know the consensus status matches the high/low status
    # we also know that illumina is not mixed or maybe, and that if there is illumina data, illumina status matches high/low status
    
    ## CASE 16
    if data['ont_AF']<maf_flag and data['unambig']:
        return('LFV')
    
    ## CASE 17
    if data['ont_AF']>(1-maf_flag) and data['unambig']:
        return('HFV')
    
    # at the end, ensure we catch all remaining consensus N variants
    
    ## CASE 18
    if data['unambig']==False:
        return('UNK-C')
    
    # all cases should be covered by this point
    # we should never get here
    sys.exit('you found a scenario not covered; please modify case_by_flags')


def status_by_case(variant_data,case_definitions,maf_flag):
    
    # load a text file with case definitions
    defs = pd.read_csv(case_definitions)
    
    # determine the case for this particular variant
    case = case_by_flags(variant_data,maf_flag)
    
    # get description and status for this case from definitions table
    current_case = defs[defs.case==case]
    description = current_case.description.values[0]
    status = current_case.status.values[0]
    
    return(case,description,status)