# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 14:21:32 2022

@author: Mahmoud Saad

Application of the Bayes Filter for Junction Temperature estimation
"""
import os,math
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import solve_ivp
from matplotlib.patches import Patch
from scipy.interpolate import interp1d,RectBivariateSpline,RegularGridInterpolator
from scipy.io import savemat
from mpl_toolkits.mplot3d import Axes3D
#%% Import Data

plt.close('all')

path = os.path.dirname(__file__)

data_sheet_path = path +'/data/'
file_name = 'Rds_Tj_measurement.npz'
file_name_measurement='VI_curve_measurement.npz'

data=np.load(file_name,allow_pickle=True)
data_measurement=np.load(file_name_measurement,allow_pickle=True)


#Extract Model Data into numpy array
dR_ds=np.empty([3,len(data['Tj_177_degC'])])
Tj=np.empty([3,len(data['Tj_177_degC'])])
Idt= np.array(data['Idt'])
Vds=np.array([data['Vds_30_degC'],data['Vds_125_degC'],
              data['Vds_177_degC']])
Tc=[30,125,177]


dR_ds[0]=data['Rds_on_30_degC']
dR_ds[1]=data['Rds_on_125_degC']
dR_ds[2]=data['Rds_on_177_degC']

Tj[0]=data['Tj_30_degC']
Tj[1]=data['Tj_125_degC']
Tj[2]=data['Tj_177_degC']

#Extract Measurements into numpy array

Ids_meas=np.empty((4,57),dtype=object)
Vds_meas=np.empty((4,57),dtype=object)

Ids_meas_175=np.empty((1,51),dtype=object)
Vds_meas_175=np.empty((1,51),dtype=object)

Ids_meas[0]=data_measurement['Id_25_degC_1200V']
Ids_meas[1]=data_measurement['Id_85_degC_1200V']
Ids_meas[2]=data_measurement['Id_125_degC_1200V']
Ids_meas[3]=data_measurement['Id_150_degC_1200V']
Ids_meas_175[0]=data_measurement['Id_175_degC_1200V']

Vds_meas[0]=data_measurement['Vds_25_degC_1200V']
Vds_meas[1]=data_measurement['Vds_85_degC_1200V']
Vds_meas[2]=data_measurement['Vds_125_degC_1200V']
Vds_meas[3]=data_measurement['Vds_150_degC_1200V']
Vds_meas_175[0]=data_measurement['Vds_175_degC_1200V']

Tj_meas=np.array([25,85,125,150])
Tj_anfang=np.linspace(0,180,num=15)
Tj_meas_2=np.concatenate([Tj_meas,Tj_anfang.T])


#%% Bivariate Spline Interpolation of Current Voltage

f_inter=[interp1d(Idt.T,Vds[0]) , interp1d(Idt.T,Vds[1]),
         interp1d(Idt.T,Vds[2])]


Id_new=np.linspace(0,1028,num=4500)
Vds_new=np.array([f_inter[0](Id_new),f_inter[1](Id_new),
                 f_inter[2](Id_new)])

Vds_bi=RectBivariateSpline(Id_new,np.array([Tc]),Vds_new.T,
                           bbox=[min(Id_new),max(Id_new), 0,200],kx=1, ky=1)


Tj_inter=np.linspace(0,180,num=15)
Vds_inter=Vds_bi(Id_new,Tj_inter)

dR_ds_inter=np.empty((15,4500))
dR_ds_try=np.empty((3,4500))

for i in range(0,len(dR_ds_inter)):
    dR_ds_inter[i]=np.gradient(Vds_inter.T[i])/np.gradient(Id_new.T)
    # dR_ds_inter[i]=np.gradient(Vds_inter.T[i])/np.gradient(Id_new.T)
    if i<3: 
        dR_ds_try[i]=np.gradient(Vds[i])/np.gradient(Idt)

color_list=["red","black"]

Original_Patch= Patch(color=color_list[0], label='Original Signal at Tj=30,125,175 deg')
Interpolated_Patch= Patch(color=color_list[1], label='Interpolated Signals')
patches=[Original_Patch,Interpolated_Patch]

dR_ds_2=np.concatenate([dR_ds,dR_ds_inter])
plt.figure()
plt.plot(Vds_inter,Id_new,color='black')
plt.plot(Vds.T,Idt,'red')
plt.xlabel('Vds[V]')
plt.ylabel('Id[A]')
plt.grid(True)
plt.legend(handles=patches,loc='best')

plt.figure()
plt.plot(Id_new,1000*dR_ds_inter.T,color='black')
plt.plot(Idt,1000*dR_ds.T,color='red')
plt.xlabel('Id[A]')
plt.ylabel('Rds_on[mohm]')
plt.legend(handles=patches,loc='best')
plt.grid(True)




#%% Bivariate Spline Interpolation of Junction Temperature
Tj_inter=[interp1d(Idt.T,Tj[0]) , interp1d(Idt.T,Tj[1]),
         interp1d(Idt.T,Tj[2])]


Tj_new=np.array([f_inter[0](Idt),f_inter[1](Idt),
                 f_inter[2](Idt)])

Tj_bi=RectBivariateSpline(Idt,np.array([Tc]),Tj.T,
                           bbox=[min(Idt),max(Idt), 0,200],kx=1, ky=1)


Tj_anfang=np.linspace(0,180,num=15)
Tj_inter=Tj_bi(Idt,Tj_anfang)

Tj_2=np.concatenate([Tj,Tj_inter.T])

plt.figure()
plt.plot(Idt,Tj_inter,color='black')
plt.plot(Idt,Tj.T,color='red')
plt.grid(True)


#%% Fitting of the Reference Model
def fit_Rds_Tj(Idt,dR_ds,Tj,current,l):
    
     Tj_k=np.array([lin_intrplte(Idt,Tj[i,:],current) for i in range(0,len(Tj))])
     
     Ron_k=np.array([lin_intrplte(Idt,dR_ds[i,:],current) for i in range(0,len(dR_ds))])
     
     parameters=determine_parameters(Tj_k,Ron_k)
     Tj_cont=np.linspace(0,220,num=l)
     Rds_on_fit=parameters[0]+parameters[1]*Tj_cont+parameters[2]*(np.power(Tj_cont,2))
     
     rmse=determine_rmse(Ron_k,Tj_k,parameters)*1000

   
     return (Tj_k,Ron_k,Tj_cont,Rds_on_fit,rmse)
     
def determine_parameters(Tj_k,Ron_k):  
    X=np.array([np.ones_like(Tj_k),Tj_k,np.power(Tj_k,2)],dtype='float').T
    # # F= X.T@X
    # # F=np.array(F,dtype='float')
    F=np.linalg.pinv(X)
    print(F)
    # F_inv=np.linalg.inv(F)
    # E=np.matmul(F_inv,X.T)
    parameters=F@Ron_k
    
    return parameters
def determine_rmse(Ron_k,Tj_k,parameters):
    Ron_pred=parameters[0]+parameters[1]*Tj_k+parameters[2]*(np.power(Tj_k,2))
    diff=(Ron_k-Ron_pred)/Ron_k
    normalized_square_diff=diff**2/len(Ron_k)
    MSE=sum(normalized_square_diff)
    rmse=float(math.sqrt(MSE))
    
    return rmse


def lin_intrplte(x,y,x_ref): 
     i=find_nearest(x,x_ref)
     
     if x[i]<x_ref:
         m=(y[i+1]-y[i])/(x[i+1]-x[i])
         y_ref=m*(x_ref-x[i])+y[i]
         return y_ref
     if x[i]>x_ref: 
         m=(y[i]-y[i-1])/(x[i]-x[i-1])
         y_ref=m*(x_ref-x[i-1])+y[i-1]
         return y_ref
     else: 
          y_ref=y[i]
          return y_ref
      
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def plot_for_current_list(current_list,l):
    plt.figure()
    rmse_list=[]
    for current in current_list: 
        Tj_k,Ron_k,Tj_cont,Rds_on_fit,rmse= fit_Rds_Tj(Idt,dR_ds_2,Tj_2,current,l)
        rmse_list.append(rmse)
        colors=[]
        for position in range(0,len(Tj_k)): 
            if position <3: 
                colors.append('black')
            else:
                colors.append('red')
        #plt.scatter(Tj_k,1000*Ron_k,color=colors)
        plt.plot(Tj_cont,1000*Rds_on_fit)
        
    plt.title('Rds(Tj)-Curve 750V 6xST-dies Vg=18V Model')
    plt.xlabel('Tj [°C]')
    plt.ylabel('Rds(on) [mOhm]')
    plt.grid(True)
    plt.legend(['$I_d$={0:d}A, $R_{{MSE}}$ ={1:2.3f}%'.format(current,rmse_list[current_list.index(current)]) for current in current_list] ,loc='upper left')




def Bivariate_function(Ids_fit,Rds_meas,Tj): 
    
    f_inter_rd=[interp1d(Ids_fit[i].T,Rds_meas[i]) for i in range(0,len(Rds_meas))]
    
    Ids_fit=np.array(Ids_fit,dtype='float')


    dR_ds_new=np.array([f_inter_rd[i](Ids_fit[i]) for i in range(0,len(f_inter_rd))])

    dR_ds_bi=RectBivariateSpline(Ids_fit[0],np.array([Tj]),dR_ds_new.T,
                           bbox=[min(Ids_fit[0]),max(Ids_fit[0]), 0,200],kx=1, ky=1)
    Idt=np.linspace(min(Ids_fit[0]),max(Ids_fit[0]),num=4500)
    Tj_inter=np.linspace(0,180,num=15)
    dR_ds_inter=dR_ds_bi(Idt,Tj_inter)
    
    
    dR_ds_2=np.concatenate([Rds_meas,dR_ds_inter.T])
    
    return dR_ds_2,Idt
    
    
    
    
  

#%% Model vs Measurement
def plot_R_on_model(x_state,Tj_bayes,rmse):
    str1='$I_d$={0:d}A, $R{{MSE}}$ ={1:2.3f}%'.format(current,rmse)
    plt.figure()
    plt.plot(Tj_bayes,x_state,label = str1)
    # plt.plot(Tj_bayes,x_state)
    # plt.plot(Tj_bayes,z_calibrated)
    plt.grid(True)
    plt.legend(loc='upper left')
    plt.title("Model-Rds(Tj)")
    plt.xlabel('Tj [°C]')
    plt.ylabel('Rds(on) [mOhm]')
    
    
def  seperate_R_on(Rds_on,Rds_ch,Rds_dn,Tj_bayes):
    legend=['R_on','R_channel','R_drain']
    plt.figure()
    plt.plot(Tj_bayes,Rds_on)
    plt.plot(Tj_bayes,Rds_ch)
    plt.plot(Tj_bayes,Rds_dn)
    plt.grid(True)
    plt.legend(legend,loc='upper left')
    plt.title("Drain Resistance and Channel Resistance Model")
    plt.xlabel('Tj [°C]')
    plt.ylabel('Rds(on) [mOhm]')
    
def compare_model_measurement(x_state,z_meas,z_calibrated,Tj_bayes):    
    plt.figure()
    plt.plot(Tj_bayes,(z_calibrated)/(x_state))
    plt.grid(True)
    plt.title("Measurement vs Model- Ratio between Measured and Modelled")
    plt.xlabel('Tj [°C]')
    plt.ylabel('Rds(meas)/Rds(state)')
    
    
    plt.figure()
    plt.plot(Tj_bayes,x_state-z_meas)
    plt.legend(['Rds_model-Rds_meas'],loc='upper left')
    plt.grid(True)
    plt.title("Measurement vs Model- delta measurment")
    plt.xlabel('Tj [°C]')
    plt.ylabel('Rds(on) [mOhm]')
    
#%% Program

# Plot Rds_on vs Tj for a set of given currents
l=4500
current_list=[50,200,400,800]
plot_for_current_list(current_list,l)  # Plotting Model
# Ids_fit,Rds_meas=create_measurement_data(l) #Fit Ids
#plot_meas_for_current_list(current_list,Tj_meas) # Plotting Measurement
  
    
# Return Measurement and Model for a given current
current_position=0
current=current_list[current_position]
x=fit_Rds_Tj(Idt,dR_ds_2,Tj_2,current,l)
noise_state = np.random.normal(0,0.02,4500)
noise_meas = np.random.normal(0,0.04,4500)
Tj_bayes=x[2]

Rds_on= x[3]*1000
rmse=x[4]
Rds_dn= (9.85*(10**-7)*(Tj_bayes+273.15)**2.27)
Rds_ch=Rds_on-Rds_dn

#  z_meas= z[2]*1000+noise_meas
# offset=x[3][0]-z[2][0]
# diff=x[3]-z[2]
# m=(diff[1]-diff[0])/(Tj_bayes[1]-Tj_bayes[0])
# z_calibrated=z[2]+offset+m*Tj_bayes

# Plot A set of plots quantitativly comparing the measurements and Model
plot_R_on_model(Rds_on,Tj_bayes,rmse)
seperate_R_on(Rds_on,Rds_ch,Rds_dn,Tj_bayes)

#%% Plotting 3d Plots of On-state Voltage vs Temp. vs current
fig_test=plt.figure()

ax= fig_test.add_subplot(111,projection='3d')

Tj_test=np.linspace(0,180,num=4500)

Vds_test=Vds_bi(Id_new,Tj_test)

for i in range(0,100):
    ax.plot(Id_new,Vds_test[i],Tj_test)
    ax.set_xlabel('Id[A]')
    ax.set_ylabel('Vds[V]')
    ax.set_zlabel('Tj[°C]')



#%% Saving Data

Data_dic={
    'Idt': Idt, 
    'Tj': Tj_2,
    'Rds_on':dR_ds_2,
    }
savemat('Lookup.mat',Data_dic)



      
    
