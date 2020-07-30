#!/usr/bin/env python
"""
Parse sequencing metadata to tsv in GISAID format

"""

import pandas as pd
import numpy as np 
import sys, json , os
import argparse
import random , time
from datetime import datetime , date , timedelta
from pandas.io.json import json_normalize

# INPUT tsv file should have following columns

# Order of metadata excel sheet currently 
INPUT_TSV_FIELDS =  ['id' , 'name' , 'strain' , 'sex' ,'age' , 'state' , 'state_home', 'date_fuzzed' , 'replicate', 'plate', 'well' , 'barcode', 'ng/uL', 'ct_value', 'run_date', 'run_id' , 'depth' ]

GISAID_FIELDS = [ "Virus name",
                "Type",
                "Passage details/history",
                "Collection date",
                "Location",
                "Additional location information",
                "Host",
                "Additional host information",
                "Gender",
                "Patient age",
                "Patient status",
                "Specimen source",
                "Outbreak",
                "Last vaccinated",
                "Treatment",
                "Sequencing technology",
                "Assembly method",
                "Coverage",
                "Originating lab",
                "Address",
                "Sample ID given by the sample provider",
                "Submitting lab",
                "Address.1",
                "Sample ID given by the submitting laboratory",
                "Authors",
                "Submitter",
                "Submission date",
                "Run date",
                "Run folder" ]

#states_file = "/home/idies/workspace/covid19/ncov_reference/US_states.csv"

# Common Metadata Details
# Modify these variables to the project requirements
region  = "North America"
country = "USA"
seq_tech = "Nanopore MinION"
assembly_method = "ARTIC v1.0.0_20cab1a"
lab_name = "Johns Hopkins Hospital Department of Pathology"
lab_address = "600 N. Wolfe St; Baltimore, MD 21287"
authors = "Peter M. Thielen, Thomas Mehoke, Shirlee Wohl, Srividya Ramakrishnan, Melanie Kirsche, Amanda Ernlund, Oluwaseun Falade-Nwulia, Timothy Gilpatrick, Paul Morris, Norah Sadowski, Nídiá Trovao, Victoria Gniazdowski, Michael Schatz, Stuart C. Ray, Winston Timp, Heba Mostafa"
submitter = "JHU"
sex = {"M" : "Male" , "F" : "Female" , "U" : "Unknown"}

DROP_COLS = [ 'plate' , 'well' , 'id','replicate','barcode', 'ng/uL' , 'ct_value']


def get_fasta_lengths(fasta):
    fa_lens = {} 
    for seq_record in SeqIO.parse(fasta, "fasta"):
        if seq_record.id not in fa_lens:
           fa_lens[seq_record.id] =  len(seq_record)
        else: 
           print("INFO : Duplicate record found for " + seq_record.id )
    return fa_lens

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Parse augur json files")
    parser.add_argument("-in",dest="tsv_file",type=str, required=True, help="Path to the tsv file from sequencing runs")
    parser.add_argument("-c",dest="cov_file",type=str, required=False, help="Path to the tsv file from sequencing runs including coverage from runs")
    parser.add_argument("-y",dest="config_file",type=str, default="metadata.yaml",required=False, help="Path to the config file with common metadata columns")
    parser.add_argument("-s",dest="states_file",type=str,required=True, help="Path to the states file with mapping to 2 letter abbreviations to full state information")
    parser.add_argument("-o",dest="out", type=str, default="gisaid_metadata.tsv", help="Output file path")
    args = parser.parse_args()
    tsv_file = args.tsv_file
    out_file = args.out
    out_dir = os.path.dirname(out_file)
    out_base = os.path.basename(out_file).split(".tsv")[0] 
    out_file2 = os.path.join(out_dir,out_base + "_nexstrain.tsv")

    states_file = args.states_file

    if args.cov_file is not None:
       cov_file = args.cov_file
    else:
       cov_file = None
       coverage = None
    
    #DROP_COLS = [ 'plate' , 'well' , 'id']
    CUR_DATE =  datetime.strptime('2020-05-02', '%Y-%m-%d')
    
    data_tsv = pd.read_table(tsv_file)
    missing_fields = [ x for x in INPUT_TSV_FIELDS if x not in  data_tsv.columns.tolist() ]

    if missing_fields:
       print("ERROR : Input tsv file has missing files " + ", ".join(missing_fields))
       sys.exit()
    states_csv = pd.read_csv(states_file ,sep = ",")
    states = dict(zip(states_csv.Code,states_csv.State))
    
    # Strip white spaces from columns
    data_tsv = data_tsv.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    
    # Remove rows with strain empty
    data_tsv['strain'] = data_tsv['strain'].str.replace(" ","")
    
    data_tsv.drop(DROP_COLS, axis = 1,inplace=True) 
    
    data_tsv=data_tsv.rename(columns = { "strain" : "Virus name" ,
                                         "name" :  "Sample ID given by the sample provider",
                                         #"state" : "Location" , 
                                         "date_fuzzed" : "Collection date" , 
                                         "age" : "Patient age" ,
                                         "sex" : "Gender" ,
                                         "run_date" : "Run date" , 
                                         "run_id" : "Run folder",
                                         "depth" :  "Coverage" } ) 
    
    data_tsv["Gender"].replace(sex, inplace=True)
    data_tsv['state'].replace(states, inplace = True)
    data_tsv['state_home'].replace(states, inplace = True)
    
    data_tsv["Location"] = region + " / " + country + " / " + data_tsv["state"]
    data_tsv["Type"] = "betacoronavirus"
    data_tsv["Passage details/history"]  = "Original"
    data_tsv["Additional location information"] =  region + " / " + country + " / " + data_tsv["state_home"]
    data_tsv["Additional host information"] = ""
    data_tsv["Host"] = "Human"
    data_tsv["Patient status"] = "unknown"
    data_tsv["Specimen source"] = "Nasopharyngeal swab"
    data_tsv["Outbreak"] = ""
    data_tsv["Last vaccinated"] = ""
    data_tsv["Treatment"] = ""
    data_tsv["Sequencing technology"] = seq_tech 
    data_tsv["Assembly method"] = assembly_method
    if coverage is not None:
       data_tsv["Coverage"] = data_tsv["Virus name"].map(coverage).fillna(None) 
    data_tsv["Originating lab"] = lab_name
    data_tsv["Address"] = lab_address
    
    data_tsv[['prefix','country' , 'sample_name' , 'year']] = data_tsv["Virus name"].str.split('/',expand=True,n=4)
    data_tsv["Submitting lab"] = lab_name
    data_tsv["Address.1"] = lab_address 
    data_tsv["Sample ID given by the submitting laboratory"] = data_tsv['Sample ID given by the sample provider']
    data_tsv["Authors" ] = authors
    data_tsv["Submitter"] = submitter
    #data_tsv["Submission date"] = pd.to_datetime(CUR_DATE, format='%Y-%m-%d', errors='coerce')
    data_tsv["Submission date"] = ""
    
    # Drop newly added columns
    data_tsv.drop(['prefix','country' , 'sample_name' , 'year','state' , 'state_home'], axis = 1,inplace=True)
    
    if (len(GISAID_FIELDS) == len(data_tsv.columns.tolist())):
            data_tsv = data_tsv[GISAID_FIELDS]
            next_tsv = data_tsv[data_tsv.columns[:-4]]
            
    else:
           #print(data_tsv.columns.tolist()) 
           #print(str(len(data_tsv.columns.tolist())) + " not equal to " + str(len(GISAID_FIELDS)))
           assert len(GISAID_FIELDS) != len(data_tsv.columns.tolist()), 'Problem processing metadata fields from GISAID ; generated only :' + str(len(data_tsv.columns))
    
    not_there = [ x for x in GISAID_FIELDS if x not in  next_tsv.columns.tolist() ]
    print("INFO : Processed Samples : " +  str(data_tsv.shape[0]) + " ; Metadata fields : " + str(len(data_tsv.columns.tolist())) )    
    data_tsv.to_csv(out_file, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
    next_tsv.to_csv(out_file2, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
