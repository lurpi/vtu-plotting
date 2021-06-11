# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 00:18:20 2021

@author: AGLUR
"""


from xml.dom import minidom

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from vtk.numpy_interface import dataset_adapter as dsa

import pandas as pd
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 #to export editable text to illustrator

#%%
# Function to change string in file (necessary to modify .pvd file fromParaview/OGS)
def ogs6_prj_timestep_update(filename, t_ini, t_fin):
        txt=""
        with open(filename, 'r') as filein:
             for line in filein:
                if "<t_initial>" in line:
                    line="                    <t_initial>"+str(t_ini)+"</t_initial> \n"
                elif "<t_end>" in line:
                    line="                    <t_end>"+str(t_fin)+"</t_end> \n"
                txt=txt+line
             return txt  


# Function to change string in file (necessary to modify .pvd file fromParaview/OGS)
def inplace_change(filename, old_string, new_string):
        s=open(filename).read()
        if old_string in s:
                print('Changing "{old_string}" to "{new_string}"'.format(**locals()))
                s=s.replace(old_string, new_string)
                f=open(filename, 'w')
                f.write(s)
                f.flush()
                f.close()
        else:
                print('No occurances of "{old_string}" found.'.format(**locals()))
    
def read_vtufiles(pvd_filename):          
        inplace_change(filename,'/>','></DataSet>')
        xmldoc = minidom.parse(filename)
        itemlist = xmldoc.getElementsByTagName('DataSet')
#         print(len(itemlist))
        #print(itemlist[0].attributes['file'].value)
        filenames=[]
        vtu_times=[]
        for s in itemlist:
                filenames.append(s.attributes['file'].value)
                vtu_times.append(s.attributes['timestep'].value)
        inplace_change(filename,'></DataSet>','/>')    
        return filenames,vtu_times

def out_arr( outfile,list):
 filename=outfile
 f = open(filename, 'w')
 f.write("\n".join(list))
 f.close()
 return;
def out_str( outfile,list):
 filename=outfile
 f = open(filename, 'w')
 f.write(list)
 f.close()
 return;

def prepare_arrows(x,y,eig_vals,eig_vecs):
    soa = np.array([[x, y, 
                 x+ eig_vals[0] * eig_vecs[0][0], 
                 y+ eig_vals[0] * eig_vecs[0][1]]])
 
    soa1 = np.array([[x, y, 
                  x+ eig_vals[1] * eig_vecs[1][0], 
                  y+ eig_vals[1] * eig_vecs[1][1]]])
    if eig_vals[0]>eig_vals[1]:
        return eig_vals[0],eig_vals[1], soa, soa1
    else:
        return eig_vals[1],eig_vals[0], soa1, soa

def prepare_arrows_decomp(x,y,S1,S2,d1x,d1y,d2x,d2y):
    soa = np.array([[x, y, 
                 x+ S1 * d1x, 
                 y+ S1 * d1y]])
 
    soa1 = np.array([[x, y, 
                  x+ S2 * d2x, 
                  y+ S2 * d2y]])
    if S1>S2:
        return S1,S2, soa, soa1
    else:
        return S2,S1, soa1, soa
        
def project_reg_grid(x_data,y_data,p_data, **kwargs):
    # define grid
    
    if len(kwargs)>0: 
        x_tmp=kwargs['x']
        y_tmp=kwargs['y']
        xg=list(set(x_tmp))
        yg=list(set(y_tmp))
    else:
        xg = np.linspace(799,3201,1251)
        yg = np.linspace(-00,-1000,1251)
    x_proj, y_proj = np.meshgrid(xg, yg)
    # grid the data
    #zg = griddata((x_data, y_data), p_data, (x_proj.ravel(), y_proj.ravel()), method='nearest')    
    zg = griddata((x_data, y_data), p_data, (xg[None,:], yg[:,None]), method='linear')    


    # reshape into meshgrid size
    p_proj = zg.reshape(x_proj.shape)
    return xg, yg,zg
#%%

#1800-3200 every 1m

#all nodes
T=[]
sig=[]
folder="res_vtk"
filename=folder+'//BGRa_with-sources.pvd'

 
vtu_list,vtu_t=(read_vtufiles(filename))
for vtu_file in vtu_list:
    pass
    #print(vtu_file)

length=len(vtu_t)
temperature= [[] for oid in range(length)]
sig= [[] for oid in range(length)]
x= [[] for oid in range(length)]
y= [[] for oid in range(length)]
ratio= [[] for oid in range(length)]
time_color= [[] for oid in range(length)]
#time= [[] for oid in range(length)]
princ_stress= [[] for oid in range(length)]
S1= [[] for oid in range(length)]
S2= [[] for oid in range(length)]
dir1x= [[] for oid in range(length)]
dir1y= [[] for oid in range(length)]
dir2x= [[] for oid in range(length)]
dir2y= [[] for oid in range(length)]
princ_dir= [[] for oid in range(length)]
steps=[]
time=[]
max_sig=0  
max_T=0      

chosen_times=range(len(vtu_list))#w[11,12,22]
for ii in chosen_times:
    x_out=[]
    y_out=[]
    x_out=[]
    vtu_file=vtu_list[ii]
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(folder+'//'+vtu_file)
    reader.Update()
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)

    sig_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("sigma")) 
    T_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("temperature"))
    for jj in range(len(sig_vtu)):
           x[ii].append(nodes_nummpy_array[:,0][jj])
           y[ii].append(nodes_nummpy_array[:,1][jj])
           #sig[ii].append(np.array([[-sig_vtu[jj][0],-sig_vtu[jj][2]],[-sig_vtu[jj][2],-sig_vtu[jj][1]]])) # pick stress component 0 x 1 y 2 z 3 xy
           sig[ii].append(np.array([[-sig_vtu[jj][0],-sig_vtu[jj][3]],[-sig_vtu[jj][3],-sig_vtu[jj][1]]])) # pick stress component 0 x 1 y 2 z 3 xy
           #sig is list of lsit of arrays: arrays are processed with eigen to extract eigenvalues (principal stresses) and eigenvectors (princ stress directions)
           temperature[ii].append(T_vtu[jj])
    time.append(float(vtu_t[ii])/86400/365.25)
    steps.append(float(vtu_t[ii])/float(vtu_t[-1]))
    if max(temperature[ii])>max_T:
        max_T=max(temperature[ii])
    ## interpolate observation points

    
#%%


min_str = np.empty(len(x[ii]))
max_str = np.empty(len(x[ii]))



#%%    
chosen_times=range(len(vtu_list))#w[11,12,22]
#chosen_times=[12]
#for ii in range(len(vtu_list)):
for ii in chosen_times:
    fig = plt.figure(figsize=[7.58, 4.10])
    ax2 = fig.add_axes([ 0.12, 0.11, 0.71, 0.05])
    ax1 = fig.add_axes([0.12, 0.25, 0.71, 0.70])
    ### distribute values on coarser meshgrid
    xi,yi,temp=project_reg_grid(x[ii], y[ii], temperature[ii])
    # xi,yi,princ_stress2=project_reg_grid(x[ii], y[ii], S2[ii])
    # xi,yi,dir1xi=project_reg_grid(x[ii], y[ii], dir1x[ii])
    # xi,yi,dir2xi=project_reg_grid(x[ii], y[ii], dir2x[ii])
    # xi,yi,dir1yi=project_reg_grid(x[ii], y[ii], dir1y[ii])
    # xi,yi,dir2yi=project_reg_grid(x[ii], y[ii], dir2y[ii])
    # ### prepare arrows
    color_inf=np.min(temperature[0])
    color_sup=max_T 
    cmap=mpl.cm.RdBu_r
    ax1.contourf(xi,yi,temp,365,vmin=color_inf, vmax=color_sup ,cmap=cmap)#,extend='both',cmap=plt.cm.RdBu)#plt.legend()
    plt.savefig("Quiver_princ-stresses-"+str(float(vtu_t[ii])/86400/365.25)+"y.png", dpi=300)
    #ax1.tick_params(axis='x')
    #ax1.tick_params(axis='y')
    ax1.grid()

    norm = mpl.colors.Normalize(vmin=color_inf, vmax=color_sup)
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                norm=norm,
                                orientation='horizontal') #extend='both',
    cb1.set_label(r'Temperature (K)')
    plt.savefig("Contourf-Temp_"+str(float(vtu_t[ii])/86400/365.25)+"y.png", dpi=300)
    plt.show()
    plt.close('all')    
    ## putting principal stresses back into vtu files
#%%
