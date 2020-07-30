#!/usr/bin/env python
"""
Parse excel and prepare data for submission to GISAID

"""

## Pipeline Steps :
# Step 1: Parse beta filter excel file , identify as is submissions list , modified submissions list
# Step 2: Create positions text file for samples in modified list; update the positions in the run/artic*/5-postfilter directory
# Step 3: Create a dated tab seperated report file with submission sample names , file path , STATUS (ASIS / UPDATED )
# Step 4: Make fasta file of all samples from all runs from that date
# Step 5: (Optional) perform validation of fasta and metadata
# Step 6: Rename them to be GISAID compliant , add that file to jhu_sequences/submissions/dated_sample_list.fasta
# Step 7: Check against submitted samples in submission folder ; mark them for resubmission
# Step 8: Subset metadata from master list , make dated metadata.tsv
# Step 9: Make nextstrain submission files 

import pandas as pd
import numpy as np 
import json , os
import argparse
import random , time
from datetime import datetime , date , timedelta
from pandas.io.json import json_normalize
import sys , glob ,re , tempfile ,shutil
import random , time
from datetime import datetime , date , timedelta
from pyfaidx import Fasta
from datetime import date
import warnings
from pytz import timezone
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import logging 
import time
import operator

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'


class GISAIDSubmission:
    """ 
    Process Beta filter cleared fasta for GISAID Submission
    """
    
    def __init__(self,excel_link,google_api_key,sub_manifest,SEQ_PATH,FASTA_PATH,MASTER_META,SUBMITTED_SEQ,SUBMISSIONS_PATH):
        """
        :param excel_link: Link to Primary excel sheet for beta filtering process
        :param google_api_key: Path to the google api key json
        :param sub_manifest: Submission Manifest file for a sequencing run
        :param SEQ_PATH: path to sequencing runs
        :param FASTA_PATH: path to fasta file within sequencing runs
        :param MASTER_META: path to master metadata excel sheet
        :param SUBMITTED_SEQ: path to fasta file of already submitted sequences
        :param SUBMISSIONS_PATH: path to submissions
        """
        # Path to the files needed for submission
        self.EXCEL_LINK = excel_link
        self.GOOGLE_API_KEY = google_api_key
        self.id_mapping_file = sub_manifest
        self.SEQ_PATH = SEQ_PATH
        self.FASTA_PATH = FASTA_PATH
        self.MASTER_META = MASTER_META
        self.SUBMITTED_SEQ = SUBMITTED_SEQ
        self.SUBMISSIONS_PATH = SUBMISSIONS_PATH
        self.FILE_SUFFIX = "_update_pos.txt"
        self.COMPLETE_FASTA = "*.complete.fasta"
        self.UPDATE_FASTA = "*.complete.updated.fasta"
        self.DATE_FMT ="%Y%m%d"
        self.CUR_DATE = datetime.now(timezone('US/Eastern'))
        self.CUR_DATE = self.CUR_DATE.strftime(self.DATE_FMT)
        self.PREFIX="hCoV-19/USA/" # GISAID Compliant
        self.SUFFIX = "/2020"      # GISAID Compliant
        
        # List below is just for reference of what is in the excel sheet
        self.primary_beta_fields = ["sample",
                       "position",
                       "run", 
                       "run_id" , 
                       "examiner", 
                       "proposed outcome", 
                       "change needed?", 
                       "FROM" ,
                       "TO", 
                       "comments",
                       "review complete",
                       "sample ready for modification+submission (all positions reviewed)", 
                       "marked on primary spreadsheet" ,
                       "submitted" ,
                       "submission date"] 


        self.DROP_COLS = [ 'comments' , 'examiner' , 'marked on primary spreadsheet', 'proposed outcome',  "review complete" , "run" ]
        self.RESUB_DROP_COLS = [ 'comments' , 'marked on primary spreadsheet', 'proposed outcome',  "review complete" , "run" ]
        self.SUB_DROP_COLS = [ "FROM" , "TO" ,  "run_id" , "sample ready for modification+submission (all positions reviewed)" , "change needed?"]
        self.CHANGE_DROP_COLS = ["run_id" , "sample ready for modification+submission (all positions reviewed)" , "change needed?"]
        self.SUBMIT_FILTER = "sample ready for modification+submission (all positions reviewed)"
        self.UPDATE_FILTER = "change needed?"
        self.SUB_STATUS = "submitted"
        self.submissions_dir = os.path.join(self.SUBMISSIONS_PATH, self.CUR_DATE) 
        if not os.path.exists(self.submissions_dir):
           os.makedirs(self.submissions_dir)
        self.report_file = os.path.join(self.submissions_dir, self.CUR_DATE + "_sample_submissions_report.txt")
        self.resub_report_file = os.path.join(self.submissions_dir, self.CUR_DATE + "_resubmissions_report.txt")

    def download_google_sheet(self):
        """
        Function to download a google spreadsheet
        """
        
        scope = ['https://spreadsheets.google.com/feeds',
                 'https://www.googleapis.com/auth/drive']

        credentials = ServiceAccountCredentials.from_json_keyfile_name(
        self.GOOGLE_API_KEY, scope)

        gc = gspread.authorize(credentials)

        sh = gc.open_by_url(self.EXCEL_LINK)
        #sheet_list = []
        #for spreadsheet in gc.openall():
        #     sheet_list.append(spreadsheet.title)
 
        #if "Curation" not in sheet_list or "Resubmission" not in sheet_list:
        #    self.log("ERROR: Spreadsheet  \"Curation\" or \"Resubmission\" not present !")
        #    sys.exit() 
        try:
           all_wks = sh.worksheets()
           c_wks = sh.worksheet("Curation")
           c_data = c_wks.get_all_values()
           c_headers = c_data.pop(3)
           self.EXCEL_FILE = os.path.join( self.submissions_dir, "specific_variant_validation_MASTER.xlsx")
           self.log("INFO  : Saving input excel sheet : {0}".format(self.EXCEL_FILE))
           export_file = c_wks.export(format='xlsx')
        
           #f = open(self.EXCEL_FILE , 'wb')
           #f.write(export_file)
           #f.close()
           #export_file = worksheet.export(format='xlsx')

           #c_wks = sh.worksheet("Curation")
           r_wks = sh.worksheet("Resubmission")
           r_data = r_wks.get_all_values()
           ## Third line is the header
           r_headers = r_data.pop(0)
           self.cdata_tsv = pd.DataFrame(c_data, columns=c_headers)
           self.cdata_tsv = self.cdata_tsv.iloc[3:]
           self.rdata_tsv = pd.DataFrame(r_data, columns=r_headers)
        except:
           raise

    def log(self,message):
        """ Log messages to standard output. """
        print(time.ctime() + ' --- ' + message, flush=True)

    def map_submission_names(self):
        """ 
        Extract submission names from submission_manifest.txt
        """
        sub_file = self.id_mapping_file
        sname_dict= {}
        sheader_dict = {}
        with open(sub_file,'r') as fs:
             for line in fs: 
                 run_id, sample_name, sample_header , renamed_header, coveragex = line.strip().split("\t")  
                 s_key = sample_name
                 h_key = sample_header
                 if  s_key not in sname_dict:
                     sname_dict[s_key] = renamed_header
                 if  h_key not in sheader_dict:
                     sheader_dict[h_key] = renamed_header
        return ( sname_dict , sheader_dict )
 
    def get_sample_list(self,df,sample):
        """
        get list of files from the sample columns in pandas dataframe
        """
        return df[sample].to_list()

    def get_fasta_header(self,seq):
        """
         get header from fasta file
        """
        with open(seq,'r') as fd:
           header = []
           for line in fd:
               if line[0]=='>':
                   header.append(line.strip().split(">")[1].split(" ")[0])
        return header 

    def get_submit_file_list(self,df,run_path):
        """
        get file list of complete fasta from SEQ_PATH and FASTA_PATH
    
        """
        dir_paths = df[run_path].to_list()
        c_list = []
        for i in dir_paths:
            cfile_list = glob.glob(i + self.COMPLETE_FASTA)
            if not cfile_list:
               self.log("ERROR : Fasta file not found : " + (i + self.COMPLETE_FASTA) )
               sys.exit()
            cfile = glob.glob(i + self.COMPLETE_FASTA)[0]
            c_list.append(cfile)
        return c_list

    def get_change_file_list(df,run_path):
        """
        get file list of complete fasta from SEQ_PATH and FASTA_PATH
    
        """
        dir_paths = df[run_path].to_list()
        u_list = []
        for i in dir_paths:
            ufile = glob.glob(i + self.UPDATE_FASTA)[0]
            u_list.append(ufile)
        return u_list

    def update_fasta(self,fa_file,pos_file):
        """
        Update fasta file based on positions in the text
    
        """
        out_dir =  os.path.dirname(fa_file)
        fbase = os.path.basename(fa_file)
        fa_base = fbase.split(".fasta")[0]
        out_file =  os.path.join(out_dir,fa_base + ".updated.fasta")
        shutil.copy(fa_file,out_file)
        # Get header from fasta
        fa = Fasta(fa_file)
        f = open(fa_file).readlines()
        header = f.pop(0).strip().split(">")[1].split(" ")[0]
        with open(pos_file) as mut_table:
           next(mut_table)
           # make sure you're editing a copy of the original file
           with Fasta(out_file, mutable=True) as fasta:
               for line in mut_table:
                   pos, ref , base = line.rstrip().split("\t")
                   if "#POS" in str(pos):
                      continue
                   # convert 1-based to 0-based coordinates
                   cur_base = str(fasta[header][int(pos) - 1])
                   if ( ref == cur_base ) :
                      self.log(" Updating nucleotide base : " + str(pos) + " : " +  str(fasta[header][int(pos) - 1]) + " - " + base )
                      fasta[header][int(pos) - 1] = base
                   elif ( cur_base == base ) :
                      self.log(" Reference  nucleotide base is already updated : " + str(pos) + " : " +  cur_base + " - " + base )
                   else:
                      assert ref == cur_base , " ERROR : Mismatch in consensus base reported in the tsv; " + str(pos) + ":" +  cur_base + " : "  + ref + " check file : " + pos_file 
                   #assert base != fasta[header][int(pos) - 1] , " ERROR : Consensus not updated at position  " + str(pos) + ":" +  cur_base + " - " +  base 

    def update_marked_files(self,change_dict):
        """
        Make text file of postions to be changed; call update fasta on each file
        """
        result = {}
        for d in change_dict:
            if result.get(d['sample']):
               result[d['sample']]['#POS'].append(d['position'])
               result[d['sample']]['FROM'].append(d['FROM'])
               result[d['sample']]['TO'].append(d['TO'])
            else:
               result[d['sample']] = {}
               result[d['sample']]['#POS'] = [d['position']]
               result[d['sample']]['FROM'] = [d['FROM']]
               result[d['sample']]['TO'] = [d['TO']] 
               result[d['sample']]['OUTFILE'] = d['OUTFILE']
        pos_files = []    
        fa_files = []
        up_files = []
        up_samples = []
        for key in result:
            if result[key].get("OUTFILE"):
               pos_files.append(result[key]["OUTFILE"])
               fa_dir = os.path.dirname(result[key]["OUTFILE"])
               fa_file = glob.glob(fa_dir + "/" + key + self.COMPLETE_FASTA)[0]
               up_file = fa_dir + "/" + key + self.UPDATE_FASTA
               fa_files.append(fa_file)
               up_files.append(up_file)
               up_samples.append(key)
               with open(result[key]["OUTFILE"], 'w') as f:
                    f.write("#POS\tFROM\tTO\n")
                    for i in range(len(result[key]['#POS'])):
                        f.write("{0}\t{1}\t{2}\n".format(result[key]["#POS"][i],result[key]["FROM"][i],result[key]["TO"][i]))
    
        ## Call update fasta for all complete fasta and pos_files
        for f,p in zip(fa_files,pos_files):
            self.log("INFO : Updating complete fasta :" + f + " at positions : " + p)
            self.update_fasta(f,p)
        return (up_samples,up_files)
     

    def update_resubmission_files(self,change_dict):
        """
        Make text file of postions to be changed; call update fasta on each file
        
        """
        result = {}
        for d in change_dict:
            if result.get(d['sample']):
               result[d['sample']]['position'].append(d['position'])
               result[d['sample']]['FROM'].append(d['FROM'])
               result[d['sample']]['TO'].append(d['TO'])
            else:
               result[d['sample']] = {}
               result[d['sample']]['#POS'] = [d['position']]
               result[d['sample']]['FROM'] = [d['FROM']]
               result[d['sample']]['TO'] = [d['TO']]
               result[d['sample']]['OUTFILE'] = d['OUTFILE']
        pos_files = []
        fa_files = []
        re_files = []
        re_samples = []
        for key in result:
            if result[key].get("OUTFILE"):
               pos_files.append(result[key]["OUTFILE"])
               fa_dir = os.path.dirname(result[key]["OUTFILE"])
               fa_file = glob.glob(fa_dir + "/" + key + self.COMPLETE_FASTA)[0]
               up_file = fa_dir + "/" + key + self.UPDATE_FASTA
               if glob.glob(up_file):
                  up_file =  glob.glob(up_file)[0]
                  base = os.path.basename(up_file).split(".fasta")[0]
                  temp_file = fa_dir + "/" + self.CUR_DATE +  "_" + base + ".fasta"
                  shutil.copy(up_file, temp_file) 
                  fa_file = glob.glob(up_file)[0]
               fa_files.append(fa_file)
               re_files.append(up_file)
               re_samples.append(key)
               # Check if position file exists and add only postitional updates that are not present
               if os.path.exists(result[key]["OUTFILE"]):
                  posfile = open(result[key]["OUTFILE"], 'r')
                  pos_list = posfile.readlines()
                  posfile.close()
                  append_write = 'a+' # append if already exists
               else:
                  append_write = 'w+' # make a new file if not
               #print("append_write :" + append_write )
               with open(result[key]["OUTFILE"], append_write) as f:
                    header = "#POS\tFROM\tTO"
                    if f.tell() == 0:
                       print("File : " + result[key]["OUTFILE"] + "does not exist ; creating" )
                       f.write(header+"\n")
                    print([ x.strip("\n\r") for x in f] )                    
                    for i in range(len(result[key]['#POS'])):
                            new_pos = "{0}\t{1}\t{2}".format(result[key]["#POS"][i],result[key]["FROM"][i],result[key]["TO"][i])
                            if not any(new_pos == x.strip() for x in f):  
                               f.write("{0}\t{1}\t{2}\n".format(result[key]["#POS"][i],result[key]["FROM"][i],result[key]["TO"][i]))
   
        ## Call update fasta for all complete fasta and pos_files
        for f,p in zip(fa_files,pos_files):
            self.log("INFO : Updating complete fasta :" + f + " at positions : " + p)
            self.update_fasta(f,p)
        return (re_samples,re_files)

    def write_report(self):
        """
        Function to report to get the list of samples submitted by date
        """
        pass
    
    def sed_inplace(self,filename, pattern, repl):
        """
        Perform the pure-Python equivalent of in-place `sed` substitution: e.g.,
        `sed -i -e 's/'${pattern}'/'${repl}' "${filename}"`.
        """
        # For efficiency, precompile the passed regular expression.
        pattern_compiled = re.compile(pattern)

        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp_file:
             with open(filename) as src_file:
                 for line in src_file:
                      tmp_file.write(pattern_compiled.sub(repl, line))

        # Overwrite the original file with the munged temporary file in a
        # manner preserving file attributes (e.g., permissions).
        shutil.copystat(filename, tmp_file.name)
        shutil.move(tmp_file.name, filename)
    
    def get_submission_names(self,fa_list):
        """
        Get the GISAID submission names for all the files seleceted for submission
        
        """
        meta_tsv = pd.read_csv(self.MASTER_META,sep="\t")
        sname = meta_tsv["Virus name"].str.split('/').str[2].to_list()
        sname_prefix = [ i.split("-")[0] for i in sname ]
        sname_real = [ i.split("-")[1] for i in sname ]
        sname_idx = [sname_real.index(e) for i, e in enumerate(fa_list) ]
        sprefix = [ sname_prefix[i] for i in sname_idx ]
        snames = [ sname_real[i] for i in sname_idx ]
        renamed_headers = [ "{0}{1}-{2}{3}".format(self.PREFIX,i,j,self.SUFFIX) for i , j in zip(sprefix,snames) ]
        return renamed_headers 
      
    def concat_fasta_files(self,fa_list, out):
        """
        Function to concat fasta files from a list
        """
        with open(out, 'w') as outfile:
             for fname in fa_list:
                 ffile = glob.glob(fname)[0]
                 with open(ffile) as infile:
                      for line in infile:
                           outfile.write(line)
        return out
 
    def filter_metadata(self,fa_headers,out):
        """
        Function to filter metadata from a list
        """
        meta_tsv = pd.read_table(self.MASTER_META)
        meta_tsv = meta_tsv[ meta_tsv["Virus name"].isin(fa_headers) ] 
        #if drop_cols  > 0:
        #    drop(CHANGE_DROP_COLS, axis = 1,inplace=True)
        meta_tsv.to_csv(out, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
        return out

    def rename_submission_fasta(self,sub_fa_file,rename_dict):
        """
        Get the GISAID submission names for all the files selected for submission
        
        """
        fa_file = Fasta(sub_fa_file, key_function=lambda key: rename_dict[key])
        renamed_fa_file = sub_fa_file + ".renamed"
        with open(renamed_fa_file, 'w') as renamed:
             for entry in fa_file:
                 print(">" + entry.name, file=renamed)
                 for line in entry:
                      print(line, file=renamed)
        shutil.copy(renamed_fa_file, sub_fa_file)
        if os.path.isfile(renamed_fa_file):
            os.remove(renamed_fa_file)
        fa_file = Fasta(sub_fa_file)
        return sub_fa_file
    
    def make_nextstrain_beta_files(self,all_fa_list,all_meta_list,out_fa,out_meta):
        """
        Combine the latest fasta to all fasta and make a current all fasta and metadata
        """
        self.concat_fasta_files(all_fa_list,out_fa)
        #TO DO : If files already exist in fasta update them to the current fasta
        #TO DO : if not just add the fasta files to the submitted list
        #shutil.copy(out_fa,SUBMITTED_SEQ) # Update the submission fasta file
        self.filter_metadata(all_meta_list,out_meta)

    def parse_beta_filter_resubmission(self):
        """
        Function to parse through the resubmission tab in the excel sheet and make changes to the files for resubmission
        """
        try:
           data_tsv =  self.rdata_tsv
           data_tsv  = data_tsv[data_tsv["sample"] != "#example"]
           #print(data_tsv.head())
           data_tsv.drop(self.RESUB_DROP_COLS, axis = 1,inplace=True)
           #print(data_tsv.columns.to_list())
           data_tsv = data_tsv.applymap(lambda x: x.strip() if isinstance(x, str) else x)
           data_tsv = data_tsv.loc[(data_tsv[self.SUBMIT_FILTER]  == "Yes") & (data_tsv[self.SUB_STATUS] != "Yes" )]
           
           if data_tsv.empty:
              print("ERROR: Check if Samples are marked submit ready in column ( with Yes ) : " +  self.SUBMIT_FILTER + " and they are not already submitted !")
              sys.exit()

           change_samples_tsv = data_tsv.loc[data_tsv[self.UPDATE_FILTER] == "Yes"]
           
           change_samples_tsv['RUN_PATH'] = self.SEQ_PATH + "/" +  change_samples_tsv["run_id"] + "/" + self.FASTA_PATH + "/" + change_samples_tsv["sample"]
           change_samples_tsv['OUTFILE'] = change_samples_tsv['RUN_PATH'] + self.FILE_SUFFIX
           change_samples_tsv.drop(self.CHANGE_DROP_COLS, axis = 1,inplace=True)
           change_dict = change_samples_tsv.to_dict('records')
           re_samples,re_files = self.update_resubmission_files(change_dict)
           out_resubmit_seq = os.path.join(self.submissions_dir, self.CUR_DATE + "_resubmissions.fasta")
           out_resubmit_meta = os.path.join(self.submissions_dir, self.CUR_DATE + "_resubmissions_metadata.tsv")

           self.concat_fasta_files(re_files,out_resubmit_seq)
           samples_dict , headers_dict = self.map_submission_names() 
           renamed_re_samples = [ samples_dict[x] for x in re_samples ]
           #renamed_final_samples = [ samples_dict[x] for x in final_files ]
           #print(renamed_re_samples)
           self.filter_metadata(renamed_re_samples,out_resubmit_meta)
            
           ###Change this part of code for a more stable solution for renaming for GISAID
           resub_fa_ori_headers = self.get_fasta_header(out_resubmit_seq)
           #print([ i.split("_")[0].replace("MD","").replace("-","") for i in sub_fa_ori_headers ])
           #rename_dict = dict(zip(sub_fa_ori_headers,renamed_final_samples_filt))
           rename_dict = { key : headers_dict[key] if key in headers_dict else key for key in resub_fa_ori_headers }
           self.log('INFO : Renaming headers in file {}'.format(out_resubmit_seq))
           self.log('\n'.join('Renamed {} : {}'.format(key, value) for key, value in rename_dict.items()))
           self.rename_submission_fasta(out_resubmit_seq,rename_dict)
           status = ["RESUBMISSION" ] * len(re_samples) 
           resub_report_tsv = pd.DataFrame({"#SAMPLE" : re_samples ,"FILE_PATH" : re_files , "STATUS" : status })
           resub_report_tsv.to_csv(self.resub_report_file, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
           # Update submitted seq with the resubmitted fasta files and metadata

        except:
           raise
        
    def parse_beta_filter_curation(self):
        """
        parse excel file and make submission list and update list
        """
        try:
            data_tsv =  self.cdata_tsv
            data_tsv  = data_tsv[data_tsv["sample"] != "#example"]
            #print(data_tsv.head())
            data_tsv.drop(self.DROP_COLS, axis = 1,inplace=True) 
            #print(data_tsv.columns.to_list())
            data_tsv = data_tsv.applymap(lambda x: x.strip() if isinstance(x, str) else x)
        
            # Check if submit filter is yes and already submitted
            data_tsv = data_tsv.loc[(data_tsv[self.SUBMIT_FILTER]  == "Yes") & (data_tsv[self.SUB_STATUS] != "Yes" )]
            if data_tsv.empty:
              self.log("ERROR: Check if Samples are marked submit ready in column ( with Yes ) : " +  self.SUBMIT_FILTER + " and they are not already submitted !")
              sys.exit()

            # Add samples with changes needed ; Identify positions ; Update them
            change_samples_tsv = data_tsv.loc[data_tsv[self.UPDATE_FILTER] == "Yes"] 
            change_samples_tsv['RUN_PATH'] = self.SEQ_PATH + "/" +  change_samples_tsv["run_id"] + "/" + self.FASTA_PATH + "/" + change_samples_tsv["sample"] 
            change_samples_tsv['OUTFILE'] = change_samples_tsv['RUN_PATH'] + self.FILE_SUFFIX
            change_samples_tsv.drop(self.CHANGE_DROP_COLS, axis = 1,inplace=True)
            #print(change_samples_tsv.head)
            change_dict = change_samples_tsv.to_dict('records')
            #print(change_dict) 
            u_samples,u_files = self.update_marked_files(change_dict)
        
            # Submit them if change needed is No ; get list ; get files
            submit_samples_tsv = data_tsv.loc[(data_tsv[self.UPDATE_FILTER] == "No" ) & (~data_tsv["sample"].isin(u_samples))]  
            submit_samples_tsv = submit_samples_tsv.drop_duplicates(subset='sample', keep="first")
            submit_samples_tsv['RUN_PATH'] = self.SEQ_PATH + "/" +  submit_samples_tsv["run_id"] + "/" + self.FASTA_PATH + "/" + submit_samples_tsv["sample"] 
            s_samples = self.get_sample_list(submit_samples_tsv,'sample')
            s_files = self.get_submit_file_list(submit_samples_tsv,'RUN_PATH')
            self.log("INFO : Samples identified for submission AS-IS : {0} ".format(", ".join(s_samples)))
            submit_samples_tsv.drop(self.SUB_DROP_COLS, axis = 1,inplace=True)
            #print(submit_samples_tsv.head)
        
            ### Write dated report to sequencing_runs 
            self.log("INFO : Writing report {0} ".format(self.report_file))
            final_samples = s_samples + u_samples
            final_files = s_files + u_files
         
            # Get submission names for final files
            # if in submitted , throw WARNING, Already submitted sample : Sample name , skipping
            # make submission fasta
            # Rename headers 
            samples_dict , headers_dict = self.map_submission_names() 
            renamed_final_samples = [ samples_dict[x] for x in final_samples ]
            #renamed_final_samples = [ samples_dict[x] for x in final_files ]
            #print(renamed_final_samples) 
            
            # Check already submitted and selected fasta only for new submissions ; take resubmission seperately
            existing_seq = self.get_fasta_header(self.SUBMITTED_SEQ)
            resub_indices = [i for i, e in enumerate(renamed_final_samples) if e in existing_seq]
            # Add files to the submission path for submission
            outseq = os.path.join(self.submissions_dir, self.CUR_DATE + "_jhu_sequences.fasta")
            outmeta = os.path.join(self.submissions_dir, self.CUR_DATE + "_jhu_sequences_metadata.tsv")

            renamed_final_samples_filt =  [e for i,e in enumerate(renamed_final_samples) if i not in resub_indices ]
            final_files_filt =  [e for i,e in enumerate(final_files) if i not in resub_indices ]
            final_samples_filt =  [e for i,e in enumerate(final_samples) if i not in resub_indices ]
            #print(final_samples_filt)
            self.concat_fasta_files(final_files_filt,outseq)
            #print(renamed_final_samples_filt)
            self.filter_metadata(renamed_final_samples_filt,outmeta)
            
            ###Change this part of code for a more stable solution for renaming for GISAID
            #####
            sub_fa_ori_headers = self.get_fasta_header(outseq)
            #print([ i.split("_")[0].replace("MD","").replace("-","") for i in sub_fa_ori_headers ])
            #rename_dict = dict(zip(sub_fa_ori_headers,renamed_final_samples_filt))
            rename_dict = { key : headers_dict[key] if key in headers_dict else key for key in sub_fa_ori_headers }
            self.log('INFO : Renaming headers in file {}'.format(outseq))
            self.log('\n'.join('Renamed {} : {}'.format(key, value) for key, value in rename_dict.items()))
            self.rename_submission_fasta(outseq,rename_dict)
            #sed_inplace(outseq,"MDHP-00", "MD-HP")
            #sed_inplace(outseq,">MD-HP", ">" + PREFIX + "MD-HP")
            #sed_inplace(outseq,"_[A-Z,a-z,0-9]*$", SUFFIX)
            #####
        
            status = ["AS-IS" ] * len(s_samples) +  ["UPDATED" ]* len(u_samples) 
            # If sample already submitted ; mark it for resubmission
            if len(resub_indices) > 0 :
                for index in resub_indices:
                   status[index] = "OMITTED (Already submitted file, if change needed move to resubmission tab) "
            report_tsv = pd.DataFrame({"#SAMPLE" : final_samples ,"FILE_PATH" : final_files , "STATUS" : status })
            #report_tsv = pd.DataFrame({"#SAMPLE" : final_samples_filt ,"FILE_PATH" : final_files_filt , "STATUS" : status })
            report_tsv.to_csv(self.report_file, sep='\t', encoding='utf-8',  index=False, line_terminator='\r\n')
        
            #Make nextstrain beta files
            self.log("INFO : Making files for beta nextstrain")
            all_outmeta = os.path.join(self.submissions_dir , self.CUR_DATE +"_all_submitted_jhu_sequences_metadata.tsv")
            all_outfa = os.path.join(self.submissions_dir, self.CUR_DATE + "_all_submitted_jhu_sequences.fasta")
            all_meta_list = existing_seq + renamed_final_samples_filt
            all_seq_head = [self.SUBMITTED_SEQ ,outseq] 
            self.make_nextstrain_beta_files(all_seq_head,all_meta_list,all_outfa,all_outmeta)
        except:
            raise
                  
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to Parse beta filter files from curation from a google spreadsheet")
    #parser.add_argument("-in",dest="excel_file",type=str, required=True, help="Path to the tsv file from sequencing runs")
    parser.add_argument("-in",dest="excel_link",type=str, required=True, help="Path to the google spreadsheet link file from manual curation")
    parser.add_argument("--api_key",dest="api_key",type=str, required=True, help="Path to the json file of google api key")
    parser.add_argument("--submission_manifest",dest="sub_manifest",type=str, required=True, help="Path to the tsv file with run_id,sample_name,sample_header,rename_header,coveragex")
    parser.add_argument("--seq-path",dest="SEQ_PATH",type=str, required=True, help="Path to the sequencing runs folder")
    parser.add_argument("--fasta-path",dest="FASTA_PATH",type=str, required=True, help="Path to the sequencing runs folder")
    parser.add_argument("--submission-path",dest="SUBMISSIONS_PATH",type=str, required=True, help="Path to submissions directory")
    parser.add_argument("--submitted-fasta",dest="SUBMITTED_SEQ",type=str, required=True, help="Path to a single fasta of already submitted fasta files")
    parser.add_argument("--master-meta",dest="MASTER_META",type=str, required=True, help="File containing the master metadata for GISAID submission")
    parser.add_argument("--resubmission",dest="resubmission_flag",default=False,action='store_true',required=False, help="Set this flag for doing resubmissions")


    args = parser.parse_args()
    excel_link = args.excel_link
    GOOGLE_API_KEY=args.api_key
    gisaid_sub = GISAIDSubmission(excel_link,GOOGLE_API_KEY,args.sub_manifest,args.SEQ_PATH,args.FASTA_PATH,args.MASTER_META,args.SUBMITTED_SEQ,args.SUBMISSIONS_PATH)
    try:
        gisaid_sub.download_google_sheet()
       
        if not args.resubmission_flag:
           gisaid_sub.parse_beta_filter_curation()
        else:
           #gisaid_sub.parse_beta_filter_curation()
           gisaid_sub.parse_beta_filter_resubmission()
    except:
        raise
