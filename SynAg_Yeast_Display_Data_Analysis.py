#!/usr/bin/env python
# coding: utf-8

'''
This script is designed to take in a plan ID and job ID of an Aquarium protcol for SynAg surface display, 
download the data, and return graphs and tables of the mean and median fluorescence. 

It is important to make sure that the yeast are displaying correctly because in the next step of SynAg, we are using yeast
mating to compare the strength of protein/protein interactions. If the proteins are not being displayed on the surface 
of yeast, then the data from these matings won't be useful. 
'''


# In[29]:


#0. Specificications
aq_username = ""
aq_password = ""
aq_url = "http://52.27.43.242/" #Leave as is for UW BIOFAB server. 
plan_id = 34064
job_id = 101199
host_folder = '' #Folder that will host a new directory with the data files and figures for this analaysis
new_folder_name = '' #Name for new directory, which will be a subdirectory of 'host folder'


# In[30]:


#1. Import the required modules
import pandas as pd
import os
from FlowCytometryTools import *
from pylab import *
from FlowCytometryTools import FCPlate
import pprint

import pydent
from pydent import AqSession
from tqdm import tqdm
prod = AqSession(aq_username, aq_password, aq_url ) # Production Databaseå


# In[31]:


#This block is finding the operations from your plan, and then getting the well numbers for each strain.
#This script assumes 2 to 4 replicates of each strain on the plate. This is very hardwired. I should rewrite 
# it so that its easier to define the set-up used.


#Enter a plan ID, get a list of operations.
plan = prod.Plan.find(plan_id)
#designate plan for later
job = prod.Job.find(job_id)
ops = plan.operations
yd_ops = filter(lambda x: x.operation_type_id == 563 and x.status == 'done', ops)

#2. Define the metadata table
plate_metadata = []

for op in yd_ops:
    op_metadata = []
    wells = filter(lambda x: x.key == 'sample_plate_location', op.data_associations)
    for well in wells:
        op_metadata.append(well.object['sample_plate_location'])
    overnite = (filter(lambda x: x.name == 'Overnight', op.field_values))
    for on in overnite:
        op_metadata.append(on.sid)
    role =  (filter(lambda x: x.name == 'Sample_role', op.field_values))
    for r in role:
        op_metadata.append(r.value)
    plate_metadata.append(op_metadata)

for dat in plate_metadata:
    if len(dat[0][0]) < 3:
        dat[0][0] = dat[0][0][0] + '0' + dat[0][0][1]
    if len(dat[0][1]) < 3:
        dat[0][1] = dat[0][1][0] + '0' + dat[0][1][1]
    if len(dat[0]) > 2:
        if len(dat[0][2]) < 3:
            dat[0][2] = dat[0][2][0] + '0' + dat[0][2][1]
        if len(dat[0][3]) < 3:
            dat[0][3] = dat[0][3][0] + '0' + dat[0][3][1]

well_info = {}

for dat in plate_metadata:
    well_info[dat[0][0]] = [dat[1],dat[2]]
    well_info[dat[0][1]] = [dat[1],dat[2]]
    if len(dat[0]) > 2:
        well_info[dat[0][2]] = [dat[1],dat[2]]
        well_info[dat[0][3]] = [dat[1],dat[2]]
    
pprint.pprint(well_info)


# In[33]:


#1. Define the data directory and alter the filenames. Define the sample plate
#Leave as is. 
cwd = os.getcwd()
dir_path= cwd + '/' + new_folder_name
os.mkdir(dir_path)
# uploads = job.upload
job_uploads=prod.Upload.where({'job': job.id})
for u in job_uploads:
    u.download(outdir=dir_path, filename = u.name, overwrite=True)

dir_path= host_folder + '/' + new_folder_name
for file in os.listdir(dir_path):
    if '.fcs' in file and 'Well_' not in file:
        file_path = dir_path + '/' + file
        rename_file_path = dir_path + '/' + 'Well_' + file
        os.rename(file_path, rename_file_path)
        print('{} was changed to Well_{}'.format(file, file))
        
plate = FCPlate.from_dir(ID='Test Plate', path=dir_path, parser='name')
plate = plate.dropna()
print(plate)


# In[34]:


#Add info into the metadata of the wells. 
for well in plate:
    for info in well_info:
        if well == info:
            plate[well].meta["Sample ID"] = well_info[info][0]
            plate[well].meta["Role"] = well_info[info][1]
            print(well_info[info][0])
            print(well)


# In[35]:


#3. Plot scatterplots (FSC vs. SSC) for all datafiles to check if there are obviously erroneous readings to be removed
tplate = plate.transform('hlog', channels=['FSC-A', 'SSC-A'])

figure(figsize=(20,10));
title('FSC vs SSC - {}'.format(plate.ID))
tplate.plot(['FSC-A','SSC-A'], bins=100, wspace=0.2, hspace=0.2, alpha=0.9);

#Add sample names


# In[36]:


# Produce boxplots of FL1-A for all samples on the same scale 
tplate = plate.transform('hlog', channels=['FL1-A'])
figure(figsize=(20,10));
tplate.plot(['FL1-A'], xlim=(0,10000),bins = 100, color = 'green')


# In[38]:


#Calculate median FL1-A for all samples and average the triplicates. Produce a table with these values including the 
# ratio of median FL1-A to the negative control.

# sample_data = {'Sample ID': [], 'Role': [], 'FL1-A mean':[], 'FL1-A sd': []}

sample_ids = []

for well in plate:
    sample_ids.append(plate[well].meta["Sample ID"])

uniq_ids = list(set(sample_ids))

array_of_medians = []

for i in uniq_ids:
    medians = []
    for well in plate:
        if plate[well].meta["Sample ID"] == i:
            medians.append(plate[well]['FL1-A'].median())
    array_of_medians.append(medians)
    
mean_of_medians = []
stdev_of_medians = []
for m in array_of_medians:
    mean_of_medians.append(np.mean(m).round(2))
    stdev_of_medians.append(np.std(m).round(1))
    
print(mean_of_medians)
print(stdev_of_medians)

d = {'Mean': mean_of_medians,'StDev': stdev_of_medians}


medians = pd.DataFrame(data = array_of_medians, index= uniq_ids, columns=['Median FL1-A1','Median FL1-A2'])
means = pd.DataFrame(d, index= uniq_ids)


print(medians)
print(means)

figure()
means.plot.bar(y='Mean',yerr='StDev', legend = False, title = "Mean Surface Display of A and Alpha Strains of SynAg Yeast", figsize = (5,5)).set(xlabel="Sample ID", ylabel="Mean FL1-A from n = 2, bars are StDev");
savefig('19June_MeanYeastDisplay_figure.png')





