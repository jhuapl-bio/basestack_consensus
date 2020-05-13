#!/usr/bin/env python

import pandas as pd
    
def case_by_flags(data):
    """
    Assign case numbers to specific variant situations
    """
    
    # situations that automatically lead to maybe
    
    ## CASE 1
    if not pd.isna(data['depth_flag']):
        return(1)
    
    ## CASE 2
    if not pd.isna(data['ntc_flag']):
        return(2)
    
    ## CASE 3
    if data['illumina_support']=='mixed':
        return(3)
    
    ## CASE 4
    if data['illumina_support']=='maybe':
        return(4)
    
    ## CASE 5
    if pd.isna(data['illumina_support']):
        if data['ont_AF']<0.2 and data['illumina_AF']>0.8:
            return(5)
    
    
    # situations that do no depend on mismatches
    
    ## CASE 6
    if data['illumina_support']=='yes' and data['in_consensus']:
        if (not pd.isna(data['mixed_flag'])) or (not pd.isna(data['maf_flag'])):
            return(6)
        
    
    # situations with mismatches
    if not pd.isna(data['vc_flag']):
    
        ## CASE 7
        if data['vc_flag']=='mismatch(s)' and data['in_consensus']==False and data['ont_AF']<0.2:
            return(7)
        
        ## CASE 8
        elif (not pd.isna(data['sb_flag'])) and data['in_consensus']==False and data['ont_AF']<0.25:
            return(8)
        
        ## CASE 9
        elif data['vc_flag']=='mismatch(s)' and data['in_consensus']==False and data['ont_AF']<0.3 and data['illumina_support']=='no':
            return(9)
        
        ## CASE 10
        elif data['vc_flag']=='mismatch(s)' and data['in_consensus']==False and data['ont_AF']<0.3 and pd.isna(data['illumina_support']):
            return(10)
        
        ## CASE 11
        else:
            return(11)
    
    
    # situations with mixed flags (0.3<ont_AF<0.7)
    if not pd.isna(data['mixed_flag']):
        
        ## CASE 12
        if data['homopolymer'] and data['in_consensus'] and pd.isna(data['illumina_support']):
            return(12)
        
        ## CASE 13
        else:
            return(13)
        
    
    # situations with maf flags (0.15<ont_AF<0.3 or 0.7<ont_AF<0.85)
    if not pd.isna(data['maf_flag']):
        
        ## CASE 12 - REPEAT FROM ABOVE
        if data['homopolymer'] and data['in_consensus'] and pd.isna(data['illumina_support']):
            return(12)
        
        ## CASE 14
        elif data['ont_AF']<0.2 and data['in_consensus']==False:
            return(14)
        
        ## CASE 15
        elif data['ont_AF']<0.3 and data['in_consensus']:
            return(15)
        
        ## CASE 16
        elif 0.2<data['ont_AF']<0.3 and data['in_consensus']==False:
            return(16)
        
        ## CASE 17
        else:
            return(17)
    
    
    # other specific situations
    
    ## CASE 18
    if not pd.isna(data['sb_flag']):
        return(18)
    
    ## CASE 19
    if not pd.isna(data['key_flag']):
        return(19)
    
    ## CASE 20
    if data['illumina_support']=='yes' and not pd.isna(data['new_flag']):
        return(20)
    
    ## CASE 21
    if pd.isna(data['illumina_support']) and not pd.isna(data['new_flag']):
        return(21)

    ## CASE 23
    if data['illumina_support']=="illumina_only":
        return(23)
    
    
    # if we make it this far, there are no flags
    
    ## CASE 22
    flags = ['depth_flag','ntc_flag','vc_flag','mixed_flag','maf_flag','sb_flag','key_flag','new_flag']
    assert all([pd.isna(data[flag]) for flag in flags])
    return(22)


def status_by_case(variant_data,case_definitions):
    
    # load a text file with case definitions
    defs = pd.read_csv(case_definitions)
    
    # determine the case for this particular variant
    case = case_by_flags(variant_data)
    
    # get description and status for this case from definitions table
    current_case = defs[defs.case==case]
    description = current_case.description.values[0]
    status = current_case.status.values[0]
    
    return(case,description,status)
