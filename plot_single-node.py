# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 21:19:04 2021

@author: AGLUR
"""



from xml.dom import minidom

import matplotlib.pyplot as plt

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from vtk.numpy_interface import dataset_adapter as dsa

import pandas as pd
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42 #to export editable text to illustrator


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


#all nodes
nodes_list=[4, 2]# 4, 2 at (2300,-525), (2300,-475) nodes smaller interval...
T=[]
sig=[]
displ=[]
folder="res_vtk"
pvd_file="kombilager.pvd"
filename=folder+'//'+pvd_file

 
vtu_list,vtu_t=(read_vtufiles(filename))

length=len(nodes_list)
time_node= [[] for oid in range(length)]
temperature= [[] for oid in range(length)]
sig= [[] for oid in range(length)]
displ= [[] for oid in range(length)]
sigxx= [[] for oid in range(length)]
sigyy= [[] for oid in range(length)]
sigxy= [[] for oid in range(length)]
displx= [[] for oid in range(length)]
disply= [[] for oid in range(length)]
x= [[] for oid in range(length)]
y= [[] for oid in range(length)]
ratio= [[] for oid in range(length)]
time_color= [[] for oid in range(length)]
#time= [[] for oid in range(length)]
colorst= [[] for oid in range(length)]
strength= [[] for oid in range(length)]
steps=[]
time=[]
max_sig=0  
max_T=0      


for jj in range(len(nodes_list)):
    vtu_file=vtu_list[0]
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(folder+'//'+vtu_file)
    reader.Update()
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    for ii in range(len(vtu_list)):
    
           vtu_file=vtu_list[ii]
           reader = vtk.vtkXMLUnstructuredGridReader()
           reader.SetFileName(folder+'//'+vtu_file)
           reader.Update()
           nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
           nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
           #sorted and ready to accomdoate values    
           displ_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("displacement"))
           sig_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("sigma")) 
           T_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("temperature"))
           node_id=nodes_list[jj]
           #print(node_id)
           x[jj].append(nodes_nummpy_array[:,0][node_id])
           y[jj].append(nodes_nummpy_array[:,1][node_id])
           sig[jj].append(sig_vtu[node_id]) # pick stress component 0 x 1 y 2 z 3 xy
           sigxx[jj].append(-sig_vtu[node_id][0]/1.e6) # pick stress component 0 x 1 y 2 z 3 xy
           sigyy[jj].append(-sig_vtu[node_id][1]/1.e6) # pick stress component 0 x 1 y 2 z 3 xy
           sigxy[jj].append(-sig_vtu[node_id][3]/1.e6) # pick stress component 0 x 1 y 2 z 3 xy
           displ[jj].append(displ_vtu[node_id])
           displx[jj].append(displ_vtu[node_id][0])  # pick stress component 0 x 1 y
           disply[jj].append(displ_vtu[node_id][1])  # pick stress component 0 x 1 y           
           temperature[jj].append(T_vtu[node_id])
           time_node[jj].append(float(vtu_t[ii])/86400/365.25)
           if T_vtu[node_id]>max_T:
               max_T=T_vtu[node_id]
               max_time=float(vtu_t[ii])/86400/365.25
    time.append(float(vtu_t[jj])/86400/365.25)
    steps.append(float(vtu_t[jj])/float(vtu_t[-1]))

#%%

fig = plt.figure(figsize=[7.8, 3.0])
ax = fig.add_subplot(111)
for i in range(len())
ax.semilogx(time_node[0],temperature[0], ls='--',lw=3.0, label=str(x[0][0])+", "+str(y[0][0]))
ax.semilogx(time_node[1],temperature[1], ls='-',lw=3.0, label=str(x[1][0])+", "+str(y[1][0]))

ax.set_xlabel('Zeit [j]', fontsize=16)
ax.set_ylabel('Temperature [K]', fontsize=16)
norm = mpl.colors.Normalize(vmin=0, vmax=np.max(time))

plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("temp_nodes_50m.png", dpi=300)
