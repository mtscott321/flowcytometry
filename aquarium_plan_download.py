aq_username = 
aq_password = 
aq_url = "http://52.27.43.242/" #Leave as is for UW BIOFAB server. 
plan_id = 34045 
job_id = 101113 
op_id = 201958 #Find this from Aquarium. Open your plan, open the Job and find the list of operations on the left hand side. 
host_folder = '/Users/Orlando/Documents/Data/SynAgFiles' #Folder that will host a new directory with the data files and figures for this analaysis
new_folder_name = 'ODL_SynAg_190116' #Name for new directory, which will be a subdirectory of 'host folder'

import pandas as pd
import os
from FlowCytometryTools import *
from pylab import *
from FlowCytometryTools import FCPlate
import pprint
import csv

import pydent
from pydent import AqSession
from tqdm import tqdm
prod = AqSession(aq_username, aq_password, aq_url ) # Production Database

#Enter a plan ID, get a list of operations.
plan = prod.Plan.find(plan_id)
job = prod.Job.find(job_id)
# for file in job.uploads:
#     file.async_download(outdir=dir_path,overwrite=True)
cwd = os.getcwd()
dir_path= cwd + '/' + str(job.id)
os.mkdir(dir_path)

# uploads = job.uploads
job_uploads=prod.Upload.where({'job': job.id})
# prod.Upload._download_files(job_uploads, dir_path, overwrite)
for u in job_uploads:
    u.download(outdir=dir_path, filename = u.name, overwrite=False)
    
