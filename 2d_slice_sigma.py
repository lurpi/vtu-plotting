# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 23:29:36 2016

@author: Luca
"""

from xml.dom import minidom

import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
import vtk
from vtk.util.numpy_support import vtk_to_numpy

from vtk.numpy_interface import dataset_adapter as dsa
from matplotlib.patches import Rectangle
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
    
def read_vtufiles(filename):          
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
        xg = np.linspace(-50,100,51)
        yg = np.linspace(-50,100,51)
    x_proj, y_proj = np.meshgrid(xg, yg)
    # grid the data
    zg = griddata((x_data, y_data), p_data, (x_proj.ravel(), y_proj.ravel()), method='nearest')    


    # reshape into meshgrid size
    p_proj = zg.reshape(x_proj.shape)
    return x_proj.ravel(), y_proj.ravel(), p_proj.ravel()
#%%
T=[]
sig=[]
displ=[]
folder=".//" 
pvd_filename='D2023_Step2_a_Excavation.pvd'  #folder+'//TM-BGR_no_press-l16.pvd'

folder_plot=".//"

vtu_list,vtu_t=(read_vtufiles(pvd_filename))
for vtu_file in vtu_list:
    pass 
    #print(vtu_file)

length=len(vtu_t)
temperature= [[] for oid in range(length)]
sig= [[] for oid in range(length)]
x= [[] for oid in range(length)]
y= [[] for oid in range(length)]
z= [[] for oid in range(length)]
ratio= [[] for oid in range(length)]
time_color= [[] for oid in range(length)]
#time= [[] for oid in range(length)]
princ_stress= [[] for oid in range(length)]
S1= [[] for oid in range(length)]
S2= [[] for oid in range(length)]
S3= [[] for oid in range(length)]
dir1x= [[] for oid in range(length)]
dir1y= [[] for oid in range(length)]
dir1z= [[] for oid in range(length)]
dir2x= [[] for oid in range(length)]
dir2y= [[] for oid in range(length)]
dir2z= [[] for oid in range(length)]
dir3x= [[] for oid in range(length)]
dir3y= [[] for oid in range(length)]
dir3z= [[] for oid in range(length)]
princ_dir= [[] for oid in range(length)]
steps=[]
time=[]
#%%
y_section=30.3
scale_par=10.5e7


max_sig=0  
max_T=0      
#%%
chosen_times=[1,20,51]
#chosen_times=range(len(vtu_list))#w[11,12,22]
for ii in chosen_times:
    x_out=[]
    y_out=[]
    z_out=[]
    vtu_file=vtu_list[ii]
    #print(vtu_file)
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(folder+'//'+vtu_file)
    reader.Update()
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)

    #sorted and ready to accomdoate values    
    sig_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("sigma")) 
    T_vtu= vtk_to_numpy(reader.GetOutput().GetPointData().GetArray("temperature"))
    for jj in range(len(sig_vtu)):
       #print(nodes_nummpy_array[:,1][jj],(np.abs(nodes_nummpy_array[:,1][jj] - y_section)))
       if np.abs(nodes_nummpy_array[:,1][jj] - y_section)<0.2:
        # filter out the excavation
        x_r2=(nodes_nummpy_array[:,0][jj]-25)**2
        z_r2=(nodes_nummpy_array[:,2][jj]-25)**2
        rad=np.sqrt(x_r2+z_r2)
        if rad > 1.55:
           x[ii].append(nodes_nummpy_array[:,0][jj])
           y[ii].append(nodes_nummpy_array[:,1][jj])
           z[ii].append(nodes_nummpy_array[:,2][jj])
           #sig[ii].append(np.array([[-sig_vtu[jj][0],-sig_vtu[jj][2]],[-sig_vtu[jj][2],-sig_vtu[jj][1]]])) # pick stress component 0 x 1 y 2 z 3 xy
           sig[ii].append(np.array([[-sig_vtu[jj][0],-sig_vtu[jj][3],-sig_vtu[jj][5]],
                                    [-sig_vtu[jj][3],-sig_vtu[jj][1],-sig_vtu[jj][4]],
                                    [-sig_vtu[jj][5],-sig_vtu[jj][4],-sig_vtu[jj][2]]])) # pick stress component 0 x 1 y 2 z 3 xy 4 yz 5 xz
           #sig is list of lsit of arrays: arrays are processed with eigen to extract eigenvalues (principal stresses) and eigenvectors (princ stress directions)
           temperature[ii].append(T_vtu[jj])
               ## compute principal stresses with linalg from scipy 
           eigvals, eigvecs = la.eig(sig[ii][-1])
           princ_stress[ii].append(eigvals.real)
           princ_dir[ii].append(eigvecs)
           
           
           if (eigvals.real[0] > eigvals.real[1]): 
               if (eigvals.real[0] > eigvals.real[2]):
                   i = 0
                   if (eigvals.real[1] > eigvals.real[2]):
                       j=1
                       k=2
                   else:
                       j=2
                       k=1
               else:
                    i=2
                    j=0
                    k=1
           else:
               if eigvals.real[1]> eigvals.real[2]:
                   if eigvals.real[0] < eigvals.real[2]:    
                       i=1
                       j=2
                       k=0
                   else:
                       i=1
                       j=0
                       k=2
               else:
                    i=2
                    j=1
                    k=0
           S1[ii].append(eigvals.real[i])
           S2[ii].append(eigvals.real[j])
           S3[ii].append(eigvals.real[k])
           dir1x[ii].append(eigvecs.real[0][i])
           dir2x[ii].append(eigvecs.real[1][i])
           dir3x[ii].append(eigvecs.real[2][i])
           dir1y[ii].append(eigvecs.real[0][j])
           dir2y[ii].append(eigvecs.real[1][j])
           dir3y[ii].append(eigvecs.real[2][j])
           dir1z[ii].append(eigvecs.real[0][k])
           dir2z[ii].append(eigvecs.real[1][k])
           dir3z[ii].append(eigvecs.real[2][k])
           princ_stress[ii].append(eigvals.real[i],eigvals.real[j],eigvals.real[k])
           
    time.append(float(vtu_t[ii])/86400)
    steps.append(float(vtu_t[ii])/float(vtu_t[-1]))
    ## interpolate observation points
plt.scatter(x[1],z[1])
print("filtered "+str(float(len(S1[ii])))+" nodes of a total of "+str(float(len(nodes_nummpy_array))))
#%%





    
### plot the arrow plot
### prepare arrows

min_str = np.empty(len(x[ii]))
max_str = np.empty(len(x[ii]))



#%%    
#chosen_times=range(len(vtu_list))#w[11,12,22]
#for ii in range(len(vtu_list)):
print("starting the plot of the princ.stress 1&3")    
for ii in chosen_times:
    fig = plt.figure(figsize=[9, 9])
    ax = fig.add_subplot(111)
    ### distribute values on coarser meshgrid
    xi,yi,princ_stress1=project_reg_grid(x[ii], z[ii], S1[ii])
    xi,yi,princ_stress2=project_reg_grid(x[ii], z[ii], S2[ii])
    xi,yi,princ_stress3=project_reg_grid(x[ii], z[ii], S3[ii])
    xi,yi,dir1xi=project_reg_grid(x[ii], z[ii], dir1x[ii])
    xi,yi,dir1yi=project_reg_grid(x[ii], z[ii], dir1y[ii])
    xi,yi,dir1zi=project_reg_grid(x[ii], z[ii], dir1z[ii])
    xi,yi,dir2xi=project_reg_grid(x[ii], z[ii], dir2x[ii])
    xi,yi,dir2yi=project_reg_grid(x[ii], z[ii], dir2y[ii])
    xi,yi,dir2zi=project_reg_grid(x[ii], z[ii], dir2z[ii])
    xi,yi,dir3xi=project_reg_grid(x[ii], z[ii], dir3x[ii])
    xi,yi,dir3yi=project_reg_grid(x[ii], z[ii], dir3y[ii])
    xi,yi,dir3zi=project_reg_grid(x[ii], z[ii], dir3z[ii])
    ### prepare arrows, 2d plot 
    #for ji in range(len(xi)):
    cycle=len(x[ii])#len(xi)     #
    x_ver=[]
    y_ver=[]
    for ji in range(cycle):
       jj=ji # to reduce number of points plotted, arbitrary sampling
       if jj>cycle-1:
           pass
       else:
            
        #print(jj)
        Stress1, Stress2, max_str, min_str = prepare_arrows_decomp(x[ii][jj],z[ii][jj],S2[ii][jj],S3[ii][jj],dir2x[ii][jj],dir2z[ii][jj],dir3x[ii][jj],dir3z[ii][jj])
        #Stress1, Stress2, max_str, min_str = prepare_arrows_decomp(xi[jj],yi[jj],princ_stress1[jj],princ_stress3[jj],dir1xi[jj],dir1zi[jj],dir3xi[jj],dir3zi[jj])
        X, Y, U, V = zip(*max_str)
        U2 = [i*-1 for i in U]
        V2 = [i*-1 for i in V]
        X1, Y1, U1, V1 = zip(*min_str)
        
        U3 = [i*-1 for i in U1]
        V3 = [i*-1 for i in V1]
        #arbitrary scaling factor, needed for visualization (may change with different plots)
        headl=3.75 #arbitrary size factor, needed for visualization (may change with different plots)
        Q1=ax.quiver(X, Y, U, V, pivot='tip', headlength=headl,color="red", width=21e-4,scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        ax.quiver(X, Y, U2, V2, pivot='tip',headlength=headl,color="red",width=21e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        Q2=ax.quiver(X1, Y1, U1, V1, pivot='tip',headlength=headl,color="blue", width=21E-4,scale=scale_par)#, anglines='xy', scale_units='xy', scale=1)
        ax.quiver(X1, Y1, U3, V3, pivot='tip',headlength=headl,color="blue" , width=21.e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        x_ver.append(X)
        y_ver.append(Y)
        #Q1=ax.quiver(X, Y, U, V, pivot='tip', headlength=headl,color="white", width=21e-4,scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        #ax.quiver(X, Y, U2, V2, pivot='tip',headlength=headl,color="white",width=21e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        #Q2=ax.quiver(X1, Y1, U1, V1, pivot='tip',headlength=headl,color="white", width=21E-4,scale=scale_par)#, anglines='xy', scale_units='xy', scale=1)
        #ax.quiver(X1, Y1, U3, V3, pivot='tip',headlength=headl,color="white" , width=21.e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
    ax.set_xlim(15,35)
    ax.set_ylim(15,35)
    print("prepared plot "+str(ii)+" of "+str(chosen_times))
    """
    #ax.axvspan(2350,2775, alpha=0.1, color='green') #caverns
    #ax.axvspan(1360,2250, alpha=0.1, color='purple') #tunnels    

    ann = ax.annotate('HLW',
                  xy=(1360, -525), xycoords='data',
                  xytext=(0, 00), textcoords='offset points',
                  size=100,
                  bbox=dict(boxstyle="round",
                            fc=(1.0, 0.7, 0.7),
                            ec=(1., .5, .5)))#,
                  #arrowprops=dict(arrowstyle="wedge,tail_width=1.", fc=(1.0, 0.7, 0.7), ec=(1., .5, .5),                                  patchA=None,                                  patchB=el,                                  relpos=(0.2, 0.8),                                  connectionstyle="arc3,rad=-0.1"))
    """
    #plt.grid()
    # =============================================================================
    #     ax.add_patch(Rectangle((1360, -512.5), 890, -25,
    #              edgecolor = 'red',
    #              facecolor = 'purple',
    #              fill=True,
    #              alpha=0.75,
    #              lw=3,
    #              ls='--'))
    #     ax.add_patch(Rectangle((2350, -512.5), 425, -25,
    #              edgecolor = 'blue',
    #              facecolor = 'green',
    #              fill=True,
    #              alpha=0.75,
    #              lw=3,ls='--'))
    #     ax.add_patch(Rectangle((1120, -350.0), 750, 75,
    #              edgecolor = 'black',
    #              facecolor = 'white',
    #              fill=True,
    #              lw=1,zorder=1))#ax.legend()
    # =============================================================================
    #ax.scatter(2300, -395, s=220, c='green',alpha=1 ,linewidths=2.4, edgecolors='blue', marker='X', zorder=10 )
    #ax.scatter(2300, -475, s=220, c='yellow',alpha=1 ,linewidths=2.4, edgecolors='orange', marker='X', zorder=10 )
    #ax.scatter(2300, -475, s=220, c='red', alpha=0.7 ,linewidths=2.4, edgecolors='red',marker='X')
    ax.quiverkey(Q1, X=0.20, Y=0.90, U=-1e7, label="    Eigenspannung $\sigma_1$", labelpos='E', zorder=3)#, angles='xy', scale_units='xy', scale=1)
    ax.quiverkey(Q1, X=0.20, Y=0.90, U=1e7, label="", labelpos='E', zorder=3)#, angles='xy', scale_units='xy', scale=1)
    ax.quiverkey(Q2, X=0.20, Y=0.85, U=-1e7, label="    Eigenspannung $\sigma_3$", labelpos='E', zorder=3)        
    ax.quiverkey(Q2, X=0.20, Y=0.85, U=1e7, label="", labelpos='E', zorder=3)   
    ax.set_aspect('equal', adjustable='box')
    #ax.axhspan(-375,-675, alpha=0.1, color='purple') #salzstein     
    #ax.set_title("Zeit: "+str(time[ii])[0:6]+" Jahre")
    
    ax.set_xlabel('x (m)', fontsize=16)
    ax.set_ylabel('y (m)', fontsize=16)
    plt.tight_layout()
    plt.savefig(folder_plot+"Stress_1-3_slice_y"+str(float(y_section))+"m_"+str(float(vtu_t[ii])//8640/10)+"d.png", dpi=300)
    #plt.show()
    
    plt.close('all')    
    ## putting principal stresses back into vtu files
#%%
print("starting the plot of the princ.stress 1&2")
for ii in chosen_times:
    fig = plt.figure(figsize=[9, 9])
    ax = fig.add_subplot(111)
    ### distribute values on coarser meshgrid
    xi,yi,princ_stress1=project_reg_grid(x[ii], z[ii], S1[ii])
    xi,yi,princ_stress2=project_reg_grid(x[ii], z[ii], S2[ii])
    xi,yi,princ_stress3=project_reg_grid(x[ii], z[ii], S3[ii])
    xi,yi,dir1xi=project_reg_grid(x[ii], z[ii], dir1x[ii])
    xi,yi,dir1yi=project_reg_grid(x[ii], z[ii], dir1y[ii])
    xi,yi,dir1zi=project_reg_grid(x[ii], z[ii], dir1z[ii])
    xi,yi,dir2xi=project_reg_grid(x[ii], z[ii], dir2x[ii])
    xi,yi,dir2yi=project_reg_grid(x[ii], z[ii], dir2y[ii])
    xi,yi,dir2zi=project_reg_grid(x[ii], z[ii], dir2z[ii])
    xi,yi,dir3xi=project_reg_grid(x[ii], z[ii], dir3x[ii])
    xi,yi,dir3yi=project_reg_grid(x[ii], z[ii], dir3y[ii])
    xi,yi,dir3zi=project_reg_grid(x[ii], z[ii], dir3z[ii])
    ### prepare arrows, 2d plot 
    #for ji in range(len(xi)):
    cycle=len(x[ii]) #len(xi)    
    x_ver=[]
    y_ver=[]
    for ji in range(cycle):
       jj=ji*3 # to reduce number of points plotted, arbitrary sampling
       if jj>cycle-1:
           pass
       else:
            
        #print(jj)
        Stress1, Stress2, max_str, min_str = prepare_arrows_decomp(x[ii][jj],z[ii][jj],S1[ii][jj],S2[ii][jj],dir1x[ii][jj],dir1z[ii][jj],dir2x[ii][jj],dir2z[ii][jj])
        X, Y, U, V = zip(*max_str)
        U2 = [i*-1 for i in U]
        V2 = [i*-1 for i in V]
        X1, Y1, U1, V1 = zip(*min_str)
        
        U3 = [i*-1 for i in U1]
        V3 = [i*-1 for i in V1]
        headl=3.75 #arbitrary size factor, needed for visualization (may change with different plots)
        Q1=ax.quiver(X, Y, U, V, pivot='tip', headlength=headl,color="red", width=21e-4,scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        ax.quiver(X, Y, U2, V2, pivot='tip',headlength=headl,color="red",width=21e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        Q2=ax.quiver(X1, Y1, U1, V1, pivot='tip',headlength=headl,color="green", width=21E-4,scale=scale_par)#, anglines='xy', scale_units='xy', scale=1)
        ax.quiver(X1, Y1, U3, V3, pivot='tip',headlength=headl,color="green" , width=21.e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        x_ver.append(X)
        y_ver.append(Y)
    print("prepared plot "+str(ii)+" of "+str(chosen_times))
    ax.set_xlim(15,35)
    ax.set_ylim(15,35)
    ax.quiverkey(Q1, X=0.20, Y=0.90, U=-2e7, label="    Eigenspannung $\sigma_1$", labelpos='E', zorder=3)#, angles='xy', scale_units='xy', scale=1)
    ax.quiverkey(Q1, X=0.20, Y=0.90, U=2e7, label="", labelpos='E', zorder=3)#, angles='xy', scale_units='xy', scale=1)
    ax.quiverkey(Q2, X=0.20, Y=0.85, U=-2e7, label="    Eigenspannung $\sigma_2$", labelpos='E', zorder=3)        
    ax.quiverkey(Q2, X=0.20, Y=0.85, U=2e7, label="", labelpos='E', zorder=3)   
    ax.set_xlabel('x (m)', fontsize=16)
    ax.set_ylabel('y (m)', fontsize=16)
    ax.set_aspect('equal', adjustable='box')
    plt.tight_layout()
    plt.savefig(folder_plot+"Stress_1-2_slice_y"+str(float(y_section))+"m_"+str(float(vtu_t[ii])//864/100)+"d.png", dpi=300)
    plt.close('all')    
print("Done")    
