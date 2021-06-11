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


nodesY_385=[2048, 2049, 2050, 2051, 4, 2052, 2053, 2054, 2055]# nodes smaller interval...
T=[]
sig=[]
displ=[]
folder="res_vtk"
pvd_file="kombilager.pvd"
filename=folder+'//'+pvd_file
sig=[]
displ=[]
folder="res_vtk"
filename=folder+'//BGRa_with-sources.pvd'

 
vtu_list,vtu_t=(read_vtufiles(filename))
for vtu_file in vtu_list:
    pass
    #print(vtu_file)

length=len(vtu_t)
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


for ii in range(len(vtu_list)):
    x_out=[]
    y_out=[]
    x_out=[]
    vtu_file=vtu_list[ii]
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(folder+'//'+vtu_file)
    reader.Update()
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    #sort nodes based on their x-coord, needed otherwise messy plot
    if ii==0:
      x_reord=[]
      for jj in range(len(nodesY_385)):
        node_id=nodesY_385[jj]
        x_reord.append(nodes_nummpy_array[:,0][node_id])
    b=np.argsort(x_reord)
    nod=np.array(nodesY_385)
    nodesY_385_reordered=nod[b]
    
    #sorted and ready to accomdoate values    
    displ_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("displacement"))
    sig_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("sigma")) 
    T_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("temperature"))
    for jj in range(len(nodesY_385)):
           node_id=nodesY_385_reordered[jj]
           #print(node_id)
           x[ii].append(nodes_nummpy_array[:,0][node_id])
           y[ii].append(nodes_nummpy_array[:,1][node_id])
           sig[ii].append(sig_vtu[node_id]) # pick stress component 0 x 1 y 2 z 3 xy
           sigxx[ii].append(-sig_vtu[node_id][0]/1.e6) # pick stress component 0 x 1 y 2 z 3 xy
           sigyy[ii].append(-sig_vtu[node_id][1]/1.e6) # pick stress component 0 x 1 y 2 z 3 xy
           sigxy[ii].append(-sig_vtu[node_id][3]/1.e6) # pick stress component 0 x 1 y 2 z 3 xy
           displ[ii].append(displ_vtu[node_id]) # pick stress component 0 x 1 y 2 z 3 xy
           displx[ii].append(displ_vtu[node_id][0])
           disply[ii].append(displ_vtu[node_id][1])           
           temperature[ii].append(T_vtu[node_id])
    if max(temperature[ii])>max_T:
           max_x=x[ii]
           max_sig=sig[ii]
           max_temperature=temperature[ii]
           max_T=max(temperature[ii])
           max_time=float(vtu_t[ii])/86400/365.25
    # pp.append(pp_out)
    # temperature.append(T_out)
    # x.append(x_out)
    # y.append(y_out)
    time.append(float(vtu_t[ii])/86400/365.25)
    steps.append(float(vtu_t[ii])/float(vtu_t[-1]))

#%%
 
fig = plt.figure(figsize=[7.8, 3.0])
ax = fig.add_subplot(111)
labels = range(1,len(ratio)+1)
label_max="Max temp %.2f°K after %.1f years" %(max_T, max_time)
chosen_times=[0,11,15,36]
#for ii in chosen_times:
#   ax.plot(x[ii],temperature[ii], ls='-',lw=1.0, label=str(time[ii]))
for i in chosen_times:
   ax.plot(x[i],temperature[i], ls='--',lw=3.0, label=str(time[i])[:-2]+" years")
####
"""
to plot all times and maximum temperature, 
comment out the for cycle, then uncomment
the following lines
"""   
#ax.plot(max_x, max_temperature, lw =1.5, color='orange', label = label_max
# colormap = plt.get_cmap('coolwarm_r')
# steps = [float(i)/max(steps) for i in steps] # to normalize colors    

# #colormap = plt.get_cmap('gist_ncar_r') 
# colorst = [colormap(i) for i in steps]   
# colorst2 = [colormap(i) for i in steps]
# colorst.append((0,0,0,1))
# colorst2.append((0,0,0,1))
# for Temperature1,x1,time_color2,colorst4 in zip(temperature,x,time_color,colorst2):
   # ax.plot(x1,Temperature1,color=colorst4,linestyle='--',lw=2.5)#,label="shear stress"+time_color2 )
####"""
#ax.axvspan(2350,2775, alpha=0.1, color='green') #caverns
#ax.axvspan(1360,2250, alpha=0.1, color='purple') #tunnels
#ax.set_xlim(1000,3000)
#ax.set_ylim(300,320)
ax.set_xlabel('x [m]', fontsize=16)
ax.set_ylabel('Temperature [K]', fontsize=16)
norm = mpl.colors.Normalize(vmin=0, vmax=np.max(time))

plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("temp_horiz_385-fine.png", dpi=300)

#%%
plt.close("all")

fig = plt.figure(figsize=[7.8, 3.0])
ax = fig.add_subplot(111)
labels = range(1,len(ratio)+1)
ax.plot(x[0],sigxy[0], ls='--',lw=3.0, label=str(time[0]))
ax.plot(x[14],sigxy[14], ls='-',lw=3.0, label=str((time[14]//150)*150)[:-2]+" years")
ax.plot(x[22],sigxy[22], ls='-',lw=3.0, label=str((time[22]//950)*1000)[:-2]+" years")
ax.plot(x[38],sigxy[38], ls='-',lw=3.0, label=str((time[38]//9000)*9000)[:-2]+" years")
ax.plot(x[49],sigxy[49], ls='-',lw=3.0, label=str((time[49]//100000)*100000)[:-2]+" years")

#ax.axvspan(2350,2775, alpha=0.1, color='green') #caverns
#ax.axvspan(1360,2250, alpha=0.1, color='purple') #tunnels
#ax.set_xlim(1000,3000)
#ax.set_ylim(300,320)
ax.set_xlabel('x [m]', fontsize=16)
ax.set_ylabel('Stress [MPa]', fontsize=16)

plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("sigxy_horiz_385-fine_time.png", dpi=300)
plt.close("all")
