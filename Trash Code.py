# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 19:13:15 2022

@author: Mahmoud Saad
"""

#%% Linear Regression 

Id_mean=np.empty((5,1),dtype=object)
Vd_mean=np.empty((5,1),dtype=object)
slopes= np.empty((5,1),dtype=object)
intercept=np.empty((5,1),dtype=object)

for i in range(0,len(Id_mean)):
   
    if i==4:
        l=len([Ids_meas_175[0]])
        x=Vds_meas_175[0]
        y=Ids_meas_175[0]
    else:
        l=len(Ids_meas[i])
        x=Vds_meas[i]
        y=Ids_meas[i]
    
    Id_mean[i]=np.mean(y)
    Vd_mean[i]=np.mean(x)
    
    sum_1= 0 
    sum_2= 0
    for j in range(0,l): 
        
        sum_1=sum_1+(y[j]-Id_mean[i])*(x[j]-Vd_mean[i])
        sum_2 = sum_2 +np.power(x[j]-Vd_mean[i],2)
        
    slopes[i]=sum_1/sum_2 
    intercept[i]=Id_mean[i]-slopes[i]*Vd_mean[i]




Ids_fit=[]
Vds_fit=[]
for i in range(0,len(Id_mean)): 
    
    
    Vds_fitt=np.linspace(0,4,num=1000)
    Ids_fitt=slopes[i]*Vds_fitt+intercept[i]
    
    Vds_fit.append(Vds_fitt)
    Ids_fit.append(Ids_fitt)
    
Vds_fit=np.asarray(Vds_fit)    
Ids_fit=np.asarray(Ids_fit)

#%% Bivariate Spline Interpolation of  Current Resistance
f_inter_rd=[interp1d(Idt.T,dR_ds[0]) , interp1d(Idt.T,dR_ds[1]),
         interp1d(Idt.T,dR_ds[2])]


dR_ds_new=np.array([f_inter_rd[0](Idt),f_inter_rd[1](Idt),
                 f_inter_rd[2](Idt)])

dR_ds_bi=RectBivariateSpline(Idt,np.array([Tc]),dR_ds_new.T,
                           bbox=[min(Idt),max(Idt), 0,200],kx=1, ky=1)


Tj_inter=np.linspace(0,180,num=15)
dR_ds_inter=dR_ds_bi(Idt,Tj_inter)

dR_ds_2=np.concatenate([dR_ds,dR_ds_inter.T])

plt.figure()
plt.plot(Idt,dR_ds_inter)
plt.plot(Idt,dR_ds.T,color='black')
plt.grid(True)



#%%% Unscented Kalman Filter

def fx(x, delta_t):
    xout = np.empty_like(x)
    xout = x+A*omega*np.cos(omega*t[k]+phase_shift)*delta_t  
    return xout

def hx(x):
    return x # return position [x] 

sigmas = JulierSigmaPoints(n=2, kappa=1)
ukf = UnscentedKalmanFilter(dim_x=1, dim_z=1, dt=delta_t, hx=hx, fx=fx, points=sigmas)
ukf.P *= 10
ukf.R *= .5
ukf.Q = Q_discrete_white_noise(2, dt=1., var=0.03)

zs, xs = [], []
for i in range(50):
    z = i + randn()*.5
    ukf.predict()
    ukf.update(z)
    xs.append(ukf.x[0])
    zs.append(z)
    
plt.plot(xs);
plt.plot(zs, marker='x', ls='');






n=1000 
x_KF=0
Kalman_estimate=np.empty((size,1))
Kalman_variance=np.empty((size,1))
P= 0
Q=0
R=variance
alpha=1e-3
r=0
beta=2
gamma=(alpha)**2*(1+r)
factor=np.sqrt(1+gamma)
wm=gamma/(1+gamma)
wc=gamma/(1+gamma) + (1-(alpha)**2+beta)

k=1 # step intiallization

for y in measurements: 
    if k!=1:
        wm=1/(2*(1+gamma))
        wc=wm
    x_sig=np.array([x_KF,x_KF,x_KF]) +factor*np.array([0,np.sqrt(P),-np.sqrt(P)])
    print(x_sig)
    if k==3:
        break
    X_k= x_sig+A*omega*np.cos(omega*t[k]+phase_shift)*(delta_t)
    X_k_mean=sum(X_k)*wm
    P=sum((X_k-X_k_mean)**2)*wc+Q
    X_k=np.array([X_k_mean,X_k_mean+np.sqrt((1+gamma)*P),X_k_mean-np.sqrt((1+gamma)*P)])
    y_k= X_k
    y_mean=wm*sum(y_k)
    
    Pyy=wc*sum((y_k-y_mean)**2)+R
    Pxy=wc*sum((y_k-y_mean)*(X_k-X_k_mean))
    
    K=Pxy/Pyy
    
    x_KF=X_k_mean+K*(y-y_mean)
    P=P-Pyy*(K)**2
    
    
    Kalman_estimate[k]=x_KF
    Kalman_variance[k]=x_KF
    k+=1
    if k>=1000:
         break
    
    