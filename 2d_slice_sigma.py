# -*- coding: utf-8 -*-
"""
Created on Mar 09 09:29:36 2023

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

def prepare_arrows_decomp_3comp(x,y,S1,S2,S3,d1x,d1y,d2x,d2y,d3x,d3y):
    soa1 = np.array([[x, y, 
                 x+ S1 * d1x, 
                 y+ S1 * d1y]])
 
    soa2 = np.array([[x, y, 
                  x+ S2 * d2x, 
                  y+ S2 * d2y]])
    soa3 = np.array([[x, y, 
                  x+ S3 * d3x, 
                  y+ S3 * d3y]])
    
    return S1,S2,S3, soa1, soa2,soa3
        
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
scale_par=12.e7


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
        # time filter to exclude points after excavation
        if ii>y_section//1:
            rad=np.sqrt(x_r2+z_r2)
        else:    
            rad=1.e16
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
           #princ_stress[ii].append(eigvals.real)
           #princ_dir[ii].append(eigvecs)
           
           """
           principal stresses are not ordered,
           part of the routine to simply reoredr them
           """
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
           dir1x[ii].append(eigvecs.real[i][0])
           dir1y[ii].append(eigvecs.real[i][1])
           dir1z[ii].append(eigvecs.real[i][2])
           dir2x[ii].append(eigvecs.real[j][0])
           dir2y[ii].append(eigvecs.real[j][1])
           dir2z[ii].append(eigvecs.real[j][2])
           dir3x[ii].append(eigvecs.real[k][0])
           dir3y[ii].append(eigvecs.real[k][1])
           dir3z[ii].append(eigvecs.real[k][2])
           princ_stress[ii].append([eigvals.real[i],eigvals.real[j],eigvals.real[k]])
           
    time.append(float(vtu_t[ii])/86400)
    steps.append(float(vtu_t[ii])/float(vtu_t[-1]))


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
    cycle=len(x[ii]) 
    fig = plt.figure(figsize=[9, 9])
    ax = fig.add_subplot(111)
    """
    ### distribute values on coarser meshgrid if needed is done here
    """
    
    """
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
    cycle=len(x[ii])#len(xi)    
    """
    ### prepare arrows, 2d plot 
    for ji in range(cycle):
       jj=ji # to reduce number of points plotted, arbitrary sampling by multiplying this index
       if jj>cycle-1:
           pass
       else:
        Stress1, Stress2, Stress3, max_str, int_str, min_str = prepare_arrows_decomp_3comp(x[ii][jj],z[ii][jj],S1[ii][jj],S2[ii][jj],S3[ii][jj],dir1x[ii][jj],dir1z[ii][jj],dir2x[ii][jj],dir2z[ii][jj],dir3x[ii][jj],dir3z[ii][jj])
        
        # prepare arrows of maximum principal stress
        X1, Y1, U1, V1 = zip(*max_str)
        U1b = [i*-1 for i in U1]
        V1b = [i*-1 for i in V1]
        # prepare arrows of intermediate princ. stress
        X2, Y2, U2, V2 = zip(*int_str)
        U2b = [i*-1 for i in U2]
        V2b = [i*-1 for i in V2]
        # prepare arrows of minimum princ. stress
        X3, Y3, U3, V3 = zip(*min_str)
        U3b = [i*-1 for i in U3]
        V3b = [i*-1 for i in V3]
        
        headl=3.75 #arbitrary size factor, needed for visualization (may change with different plots)
        
        """
        plotting, each line plot a sets of "arrows"
        """
        Q1=ax.quiver(X1, Y1, U1, V1, pivot='tip', headlength=headl,color="red", width=21e-4,scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        ax.quiver(X1, Y1, U1b, V1b, pivot='tip',headlength=headl,color="red",width=21e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        Q2=ax.quiver(X2, Y2, U2, V2, pivot='tip',headlength=headl,color="green", width=21.E-4,scale=scale_par)#, anglines='xy', scale_units='xy', scale=1)
        ax.quiver(X2, Y2, U2b, V2b, pivot='tip',headlength=headl,color="green" , width=21.e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
        Q3=ax.quiver(X3, Y3, U3, V3, pivot='tip',headlength=headl,color="blue", width=21.E-4,scale=scale_par)#, anglines='xy', scale_units='xy', scale=1)
        ax.quiver(X3, Y3, U3b, V3b, pivot='tip',headlength=headl,color="blue" , width=21.e-4, scale=scale_par)#, angles='xy', scale_units='xy', scale=1)
    ax.set_xlim(15,35)
    ax.set_ylim(15,35)
    print("prepared plot "+str(ii)+" of "+str(chosen_times))
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
    
    ax.quiverkey(Q1, X=0.20, Y=0.90, U=-.5e7, label="    Eigenspannung $\sigma_1$", labelpos='E', zorder=3)#, angles='xy', scale_units='xy', scale=1)
    ax.quiverkey(Q1, X=0.20, Y=0.90, U=.5e7, label="", labelpos='E', zorder=3)#, angles='xy', scale_units='xy', scale=1)
    ax.quiverkey(Q2, X=0.20, Y=0.85, U=-.5e7, label="    Eigenspannung $\sigma_2$", labelpos='E', zorder=3)        
    ax.quiverkey(Q2, X=0.20, Y=0.85, U=.5e7, label="", labelpos='E', zorder=3)   
    ax.quiverkey(Q3, X=0.20, Y=0.8, U=-.5e7, label="    Eigenspannung $\sigma_3$", labelpos='E', zorder=3)        
    ax.quiverkey(Q3, X=0.20, Y=0.8, U=.5e7, label="", labelpos='E', zorder=3)   
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('x (m)', fontsize=16)
    ax.set_ylabel('z (m)', fontsize=16)
    plt.tight_layout()
    plt.savefig(folder_plot+"Stress_1-2-3_slice_y"+str(float(y_section))+"m_"+str(float(vtu_t[ii])//8640/10)+"d.png", dpi=300)
    
    #plt.show()
    
    plt.close('all')    
    ## putting principal stresses back into vtu files
