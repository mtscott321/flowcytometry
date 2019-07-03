"""
03 July 2019
Maddy Scott
This program gets data from Aquarium and fits it will Hill equations
"""
#%%
import numpy as np
from lmfit import minimize, Parameters
import scipy.integrate as spi
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from sympy import *
import xlrd
import pandas as pd
import os
from FlowCytometryTools import *
from FlowCytometryTools import FCPlate, FCMeasurement
from pylab import *
import xlwt
import pydent
from pydent import AqSession
from tqdm import tqdm

#%%
"""
This box is intended to be most of the values that are hardcoded
There are other things that, by virtue of naming conventions, etc
have to be hardcoded. They will be clearly marked.
"""
aq_username = ""
aq_password = ""
aq_url = "http://52.27.43.242/" #Leave as is for UW BIOFAB server. 
plan_id = 34178
job_id = 102045
host_folder = r'C:\Users\mtsco\OneDrive\Documents\2019 Summer Klavins Lab' 
new_folder_name = '03July_Gbox_MYC3_FL'
dir_path = host_folder + '/' + new_folder_name
sheet_name = 'Gbox_MYC3_FL'
xl_name = '03July Gbox_MYC3_FL data 14:40.xls'
excel_path = host_folder + '\\' + xl_name
gboxes = [0, 1, 3, 6, 8]

#%%
"""
Defining all the functions/methods that will be used in this script
"""

#parses the names of the yeast strains to get the number of Gboxes in it
def gbox_parse(name):
    number = name.split("GBox")[1][0]
    num = int(number)
    return num

"""
Function that defines the difference (error) between the current
estimate and the actual values. Equation to be minimized.
"""
def residual(params, x, y, stderrs):
    def f(x):
        #returning f(x), based on the current parameters
        #x is an array of all the x values
        return (params['ymax'] * (x ** params['n'])) /\
                     (params['k'] + (x ** params['n']))
    yval = f(x)
    
    errs = [ (a - b)/e for a,b,e in zip(yval, y, stderrs)]
    return errs

#plotting used after completing hill equation fitting
def plotting_vals_and_fit(xdata, ydata, result, direc):
    xvals = np.linspace(0,max(xdata),420)
    plt.plot(xvals,(result.params['ymax'] * (xvals ** result.params['n'])) /\
                         (result.params['k'] + (xvals ** result.params['n'])))
    plt.scatter(xdata,ydata)
    
    plt.title("ymax = %.2f $\pm$ %.2f\nk = %.2f $\pm$ %.2f\nn = %.2f $\pm$ %.2f" % (
        
        result.params['k'].value, result.params['k'].stderr,
        result.params['n'].value, result.params['n'].stderr,
        result.params['ymax'].value, result.params['ymax'].stderr
    ))
    plt.savefig(direc, dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None, metadata=None)

#%%
"""
Connecting to Aquarium
"""
prod = AqSession(aq_username, aq_password, aq_url ) # Production Databaseå
plan = prod.Plan.find(plan_id)
job = prod.Job.find(job_id)
ops = plan.operations

#%%
"""
Making the excel file
"""
wb = xlwt.Workbook()
ws = wb.add_sheet(sheet_name, cell_overwrite_ok=True)

#%%
"""
Getting the data from Aquarium and renaming the files
"""
plate_id = 0
os.makedirs(dir_path, exist_ok=True)

#downloading data
for op in plan.operations:
    if 'Flow Cytometry 96 well' in op.operation_type.name:
        plate_id = op.input('96 well plate').child_item_id
        for da in op.data_associations:
            if "SAMPLE_UPLOAD" in da.key:
                upload = prod.Upload.find(da.upload_id)
                up_fname = upload.upload_file_name
                upload.download(outdir=dir_path, filename=up_fname)

#changing the names
for file in os.listdir(dir_path):
    if '.fcs' in file and 'Well_' not in file:
        file_path = dir_path + '/' + file
        rename_file_path = dir_path + '/' + 'Well_' + file
        os.rename(file_path, rename_file_path)
        print('{} was changed to Well_{}'.format(file, file))

#%%
"""
Writing the excel file
"""

#writing the data from downloaded .fcs files
ws.write(0, 0, 'Well')
ws.write(0, 1, 'Median FL1-A')
ws.write(0, 2, 'Strain Name')
ws.write(0, 3, 'Gboxs')
ws.write(0, 4, 'beta-estradiol conc')
i = 1
plate = FCPlate.from_dir(ID='Test Plate', path=dir_path, parser='name')
plate = plate.dropna()
for well in plate:
    ws.write(i, 0, well)
    ws.write(i, 1, plate[well]['FL1-A'].median().item())
    i += 1
#%%
    
#writing the data from the plate item
plate_item = prod.Item.find(plate_id)
all_sample_items = plate_item.as_collection().matrix
dimensions = plate_item.as_collection().dimensions
for row in range (dimensions[0]):
    for column in range (dimensions[1]):
        sam_it = plate_item.as_collection().part(row, column)
        if sam_it is not None:
            conc = float(sam_it.data_associations[2].object.get('Inducer(s)').\
                         get('beta-estradiol').get('final_concentration').get('qty'))
            #number of columns * rows completed + columns completed in current row
            xl_row = dimensions[1]*(row) + column
            ws.write(xl_row, 2, sam_it.sample.name)
            ws.write(xl_row, 3, gbox_parse(sam_it.sample.name))
            ws.write(xl_row, 4, conc)
            
wb.save(excel_path)

#%%
"""
Reading the excel file just written and turning it into x and y value arrays 
for each type of GBox. Could just make arrays earlier, but this makes the code 
easier to use again with other excel files
"""
page = xlrd.open_workbook(excel_path).sheet_by_index(0)
results = []
for i in range(len(gboxes)):
    box = gboxes[i]
    
    #define the x and y data for this set
    xdata = []
    ydata = []
    #TODO: ERRORS
    errors = []

    #for all of the cells in this column, check if it's the correct gbox
    for j in range(1, 79):
        if page.cell(j, 3).value == box:
            #append all the x and y data associated with a certain g box
            xdata.append(float(page.cell(j, 4).value))
            ydata.append(float(page.cell(j, 1).value))
            #TEMPORARY
            errors.append(1.0)
        
    #convert to np array for later functions
    xdata = np.array(xdata)
    ydata = np.array(ydata)
    ymax = np.amax(ydata)
    
    #setting the parameters for this fit
    p = Parameters()
    p.add('k', min = 0, value=500)
    p.add('n', min = 0, value=5)
    p.add('ymax', min = 0, value=ymax)
    
    #get the result for this
    result = minimize(residual,p,args=(xdata,ydata, errors))
    results.append(result)
    plot_path = dir_path + '/' + str(box) + " Gboxes with MYC3_FL"
    plotting_vals_and_fit(xdata, ydata, result, plot_path)
