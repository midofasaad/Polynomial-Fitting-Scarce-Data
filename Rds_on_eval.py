# -*- coding: utf-8 -*-
"""
Said El-Barbari
Mahmoud Saad
EA-303
22.01.2021

"""
import os
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.integrate import solve_ivp
from matplotlib.patches import Patch


# from scipy import interpolate
# from scipy.interpolate import RegularGridInterpolator
# from scipy.interpolate import RectBivariateSpline
# from scipy.optimize import fsolve
# from mpl_toolkits.mplot3d import Axes3D
# from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
# from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
# from sklearn.preprocessing import PolynomialFeatures


#%% Import Measurement Data from npz file
plt.close('all')

path = os.path.dirname(__file__)

data_sheet_path = path +'/data/'
file_name = 'VI_curve_data.npz'
file_name_measurement='VI_curve_measurement.npz'

data=np.load(file_name,allow_pickle=True)
data_measurement=np.load(file_name_measurement,allow_pickle=True)


print("data.files: "+ str(data.files))
'''
list of the files
'Id_30degC_2500A', 'Vds_30degC_2500A',
'Id_30degC_3000A', 'Vds_30degC_3000A',
'If_30degC_3000A', 'Vsd_30degC_3000A', 't_30degC',
'Id_125degC_2500A', 'Vds_125degC_2500A', 't_125degC',
'Id_177degC_2500A', 'Vds_177degC_2500A', 'If_177degC_2500A', 'Vsd_177degC_2500A', 't_177degC'
'''



#%% Extract voltage and current data sets which start at time 0 when voltage and current start changing
#----------------------------------------------------------------------------
d_start=5000  # start of change
d_end=9500  # end before turn down

Im=np.empty((3,d_end-d_start))
Vm=np.empty((3,d_end-d_start))


t=data['t_30degC'][d_start:d_end]-data['t_30degC'][d_start:d_end].min()

Im[0]=data['Id_30degC_2500A'][d_start:d_end]
Im[1]=data['Id_125degC_2500A'][d_start:d_end]
Im[2]=data['Id_177degC_2500A'][d_start:d_end]

Vm[0]=data['Vds_30degC_2500A'][d_start:d_end]
Vm[1]=data['Vds_125degC_2500A'][d_start:d_end]
Vm[2]=data['Vds_177degC_2500A'][d_start:d_end]


#%% calculate power for varying voltage and current at temperatures 30,125 and 177
Pm=Im*Vm



#%%  Calculation of Tau based on Rth anbd Cth
nr_dies_used=6           # Number of chips MOSFET  base
nr_dies_base=10           # max number of chips per module MOSFET

dies_ratio=nr_dies_used/nr_dies_base # Chip area ratio   MOSFET

# linear scaling with the number of parallel chips n -- n up; R down; C up
# Foster model
Rth=0.881584083693822*np.array([1.832480e-02,3.819517e-02,4.085566e-03,9.939172e-04,2.168942e-2,1.785258e-03])/dies_ratio       #Skalierung R_th 6 Chips, stand 082020: Rth gesamt 0.125K/W
# Rth=0.7816712208751889*np.array([1.832480e-02,3.819517e-02,4.085566e-03,9.939172e-04,2.168942e-2,1.785258e-03])/dies_ratio    #Skalierung R_th 5 Chips, stand 082020: Rth gesamt 0.133K/W
# self.Rth=np.array([1.832480e-02,3.819517e-02,4.085566e-03,9.939172e-04,2.168942e-2,1.785258e-03])/dies_ratio
Cth=np.array([1.583077e+00,1.165177e+01,7.404443e-01,5.590657e+03,7.470273e+00,2.487444e-01 ])*dies_ratio

Tau=Rth*Cth


#%% Plot VI-Curve: Current, Voltage and average power
Tc=[30, 125,177] #Temperatures at which each voltage-current measurement took place 
fig,ax=plt.subplots(3,1,figsize=(8,9))
fig.suptitle('Measured VI-Curve')
ax[0].plot(1e6*t.T,Im.T)
ax[0].set(ylabel='Id[A]')
ax[0].legend(['Tc='+str(x)+'°C' for x in Tc])
ax[0].grid()
ax[1].plot(1e6*t.T,Vm.T)
ax[1].set(ylabel='Vds[V]')
ax[1].legend(['Tc='+str(x)+'°C' for x in Tc])
ax[1].grid()
ax[2].plot(1e6*t.T,Pm.T)
ax[2].legend(['Tc='+str(x)+'°C' for x in Tc])
ax[2].set(xlabel='t[us]', ylabel='P[W]')
ax[2].grid()

#%% Calculation of Parameters for Rds on Curves
n=len(Tc)
o=2# order of the polynomial
p=np.empty((n,o+1))
pr=np.empty((n,4+1))

Id_x = np.linspace(0,d_end-d_start,t.size)
Vd_x = np.linspace(0,10,t.size)

Vds=np.empty_like(Vm)
Rds=np.empty_like(Vm)
Rdsi=np.empty_like(Vm)


Ais=np.array([np.ones_like(t), t, t**2]).T
Bis=Im[0]
X_tr=(Ais.T)@Ais
pis=np.linalg.inv(X_tr)@(Ais.T)@Bis
Idt=np.array([np.zeros_like(t),t,t**2]).T@pis.T

ax[0].plot(1e6*t.T,Idt.T)


for idx in range(n):
    Am=np.array([np.ones_like(Im[idx]), Im[idx], Im[idx]**2]).T
    Bm=Vm[idx]
    X_tr=(Am.T)@Am
    p[idx]=np.linalg.inv(X_tr)@Am.T@Bm
    print(p[idx])
    Vds[idx]=np.array([np.ones_like(Idt),Idt,Idt**2]).T@p[idx].T
    Rds[idx]=np.array([np.zeros_like(Idt),np.ones_like(Idt),2*Idt]).T@p[idx]

    Br=Rds[idx]
    Ar=np.array([np.ones_like(Im[idx]), Im[idx], Im[idx]**2]).T
    X_tr=(Ar.T)@Ar
    pr=np.linalg.inv(X_tr)@Am.T@Br
    Rdsi[idx]=np.array([np.zeros_like(Idt),np.ones_like(Idt),2*Idt]).T@pr.T


dR_ds=np.gradient(Vds,axis=1)/np.gradient(Idt)

# str_drw=d_start

# plt.figure()
# plt.title('Rds(Id)-Curve 750V 6xST-dies Vg=18V')
# plt.plot(Im[:,str_drw:].T,1000*dR_ds[:,str_drw:].T)
# plt.plot(Id_x.T,1000*Rds.T)

# plt.plot(Im[:,str_drw::200],1000*dR_ds[:,str_drw::200])

# plt.xlabel('I [A]')
# plt.ylabel('Rds(on) [mOhm]')
# plt.legend(['Tc='+str(x)+'°C' for x in Tc])
# plt.grid(True)

#%% Plot:  Rds(Id)-Curve 750V 6xST-dies Vg=18V
plt.figure()
plt.title('Rds(Id)-Curve 750V 6xST-dies Vg=18V')
plt.plot(Im.T,1000*dR_ds.T)
plt.plot(Idt.T,1000*Rdsi.T)

#plt.plot(Im[:,0::200],1000*dR_ds[:,0::200])

plt.xlabel('I [A]')
plt.ylabel('Rds(on) [mOhm]')
plt.legend(['Tc='+str(x)+'°C' for x in Tc])
plt.grid(True)



# md=LinearRegression(fit_intercept=True)
# md_L2=Ridge(alpha=0.5,fit_intercept=True)
# md_L1=Lasso(alpha=0.5,fit_intercept=True,max_iter=5000)
# # md_L12=ElasticNet(alpha=0.5,l1_ratio=1.0,fit_intercept=False,max_iter=5000)


# md.fit(Am,Bm)
# md_L2.fit(Am,Bm)
# md_L1.fit(Am,Bm)
# md_L12.fit(A_m,b_m)

#%% Plot:  IV-Curve 750V 6xST-dies Vg=18V
plt.figure()
plt.plot(Vm.T,Im.T)
plt.plot(Vds.T,Idt.T)
plt.title('IV-Curve 750V 6xST-dies Vg=18V')
plt.xlabel('Vds [V]')
plt.ylabel('Id [A]')
#plt.legend(['Tc='+str(x)+'°C' for x in Tc])
plt.legend(np.append(['Tc='+str(x)+'°C' for x in Tc],['Tc='+str(x)+'°C Fitted' for x in Tc]))
plt.grid(True)


#%% Plot:  Rds(on)-Vds Curve 750V 6xST-dies Vg=18V
plt.figure()
plt.plot(Id_x.T,dR_ds.T)
plt.plot(Id_x.T,Rds.T)
plt.title('Rds(on)-Curve 750V 6xST-dies Vg=18V')
plt.xlabel('Vds[A]')
plt.ylabel('Rds(on) [Ohm]')
plt.legend(['Tc='+str(x)+'°C Model' for x in Tc])
plt.grid(True)


#%% Plot:  Rds(on)-time Curve 750V 6xST-dies Vg=18V
plt.figure()
plt.title('Rds(on)-Curve 750V 6xST-dies Vg=18V')
plt.plot(1e6*t.T,dR_ds.T)
plt.xlabel('t [us]')
plt.ylabel('Rds(on) [Ohm]')
plt.legend(['Tc='+str(x)+'°C Model' for x in Tc])
plt.grid(True)




# plt.figure()
# plt.plot(IU_data[:,1],Rds_on*1000,label='30°C')
# plt.plot(Vd_x,Rds_on_x,label='30°C fitted')

# plt.title('Rds(on) vs.current 750V 6xST-dies Vg=18V')
# plt.xlabel('Id [A]')
# plt.ylabel('Rds(on) [mOhm]')
# plt.legend()
# plt.grid(True)
#%%
# hst0=np.histogram(Im30,300)
# hst1=np.histogram(Vm30,300)

# plt.figure()
# plt.plot(Im30)
# # plt.plot(IU_data[:,1])
# plt.grid()

# plt.figure()
# plt.plot(Vm30)
# plt.grid()

# plt.figure()
# plt.plot(hst1[1][1:],hst1[0])
# # plt.plot(hst0[1][1:],hst0[0])

# plt.grid()




#%% Calculate Junction Temperature


def get_Tj_asc(t,Tj, *args):
    """Function to calculate Junction temperature derivative with respect to time"""
    dTjdt =[]
    P = args[0]
    for Tj,Rth,Cth in zip(Tj,args[1],args[2]):
        dTj= -Tj/(Rth*Cth) + P/Cth
        dTjdt.append(dTj)               
    return dTjdt    


## SIm30ulation loop
# t=np.linspace(0,0.001,2000)
t_points=t.size

Tj=np.zeros_like(Vm)
Tj[:,0]=Tc



for i in range(n):
    Tj_init=np.zeros_like(Cth)
    for j in tqdm(range(1, t_points)):
        td = [t[j-1], t[j]]
        # solver='LSODA'
        solver='RK45'
        # solver='BDF'
        # solver='Radau'
        # solver='DOP853'
        # solver='RK23'
        dTj = solve_ivp(get_Tj_asc,td,Tj_init , args = (Pm[i][j],Rth,Cth,), method=solver,t_eval=td, rtol=1e-8, atol=1e-6)
        Tj_init = dTj.y[:,1]
        Tj[i][j] = Tj_init.sum() + Tc[i]


#%% Plot: Tj(t)-Curve 750V 6xST-dies Vg=18V
str_drw=0
plt.figure()
plt.title('Tj-Curve 750V 6xST-dies Vg=18V')
plt.plot(1e6*t[str_drw:].T,Tj[:,str_drw:].T)
plt.xlabel('t [us]')
plt.ylabel('Tj [°C]')
plt.legend(['Tc='+str(x)+'°C' for x in Tc])
plt.grid(True)

#%% Plot: Rds(Tj)-Curve 750V 6xST-dies Vg=18V
plt.figure()
plt.title('Rds(Tj)-Curve 750V 6xST-dies Vg=18V')
plt.plot(Tj[:,str_drw:].T,1000*dR_ds[:,str_drw:].T)
plt.plot(Tj[:,1000::800],1000*dR_ds[:,1000::800])
plt.xlabel('Tj [°C]')
plt.ylabel('Rds(on) [mOhm]')
plt.legend(['Tc='+str(x)+'°C' for x in Tc])
plt.grid(True)

#%% Plot: Tj(Id)-Curve 750V 6xST-dies Vg=18V
plt.figure()
plt.title('Tj(Id)-Curve 750V 6xST-dies Vg=18V')
plt.plot(Im[:,str_drw:].T,Tj[:,str_drw:].T)
plt.plot(Idt[str_drw:].T,Tj[:,str_drw:].T)
plt.xlabel('Id[A]')
plt.ylabel('Tj [°C]')
plt.legend(['Tc='+str(x)+'°C' for x in Tc])
plt.grid(True)

#%% Tj(Pv)-Curve 750V 6xST-dies Vg=18V
plt.figure()
plt.title('Tj(Pv)-Curve 750V 6xST-dies Vg=18V')
plt.plot(Pm[:,str_drw:].T,Tj[:,str_drw:].T)

plt.xlabel('P [W]')
plt.ylabel('Tj[°C]')
plt.legend(['Tc='+str(x)+'°C' for x in Tc])
plt.grid(True)





#%% Save Data relevant to the Bayes Filter

data_dic={ 
    'Rds_on_30_degC'  :dR_ds[0],
    'Vds_30_degC':Vds[0],
    'Tj_30_degC' :Tj[0],
    
   'Rds_on_125_degC'  :dR_ds[1],
   'Vds_125_degC':Vds[1],
    'Tj_125_degC' :Tj[1],
    
    'Rds_on_177_degC'  :dR_ds[2],
    'Vds_177_degC':Vds[2],
    'Tj_177_degC' :Tj[2],
    'Idt':Idt
   }


np.savez('Rds_Tj_measurement.npz', **data_dic)


