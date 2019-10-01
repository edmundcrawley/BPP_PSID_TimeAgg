# -*- coding: utf-8 -*-
"""
Replication of Gauss part of BPP 2008
This file replicates the minimum distance estimation
"""
import numpy as np
from numpy.linalg import inv
from scipy.optimize import minimize
import numdifftools as nd

def implied_cov_BPP(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params):
    if ma==1:
        teta = params[0] 
    else:
        teta = 0
    if taste:
        varcsi = params[ma] 
    else:
        varcsi = 0.0
    var_perm = params[ma+taste:ma+taste+perm_shk_params] 
    var_tran = params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] 
    ins_perm = params[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params] 
    ins_tran = params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params] 
    var_c_error = params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params+meas_error_params] 

    dify  =np.zeros((T,T)) #/* Income */
    difcd =np.zeros((T,T)) #/* Consumption */
    difc  =np.zeros((T,T)) #/* Consumption */
    difcme=np.zeros((T,T)) #/* Measurement error of consumption */
    difyc =np.zeros((T,T)) #/* Cov Income Consumption */
    dif   =np.zeros((2*T,2*T))

    #/* This is the variance of Income */
    dify[0,0]=var_perm[0]+var_tran[0]+(1-teta)**2*var_tran[0]+teta**2*var_tran[0]
    dify[1,1]=var_perm[0]+var_tran[1]+(1-teta)**2*var_tran[0]+teta**2*var_tran[0]
    dify[2,2]=var_perm[0]+var_tran[2]+(1-teta)**2*var_tran[1]+teta**2*var_tran[0]
    for j in np.array(range(T-6))+3:
        dify[j,j]=var_perm[j-2]+var_tran[j]+(1-teta)**2*var_tran[j-1]+teta**2*var_tran[j-2] 
    dify[T-3,T-3]=var_perm[T-5]+var_tran[T-3]+(1-teta)**2*var_tran[T-4]+teta**2*var_tran[T-5]
    dify[T-2,T-2]=var_perm[T-5]+var_tran[T-3]+(1-teta)**2*var_tran[T-3]+teta**2*var_tran[T-4]
    dify[T-1,T-1]    =var_perm[T-5]+var_tran[T-3]+(1-teta)**2*var_tran[T-3]+teta**2*var_tran[T-3]
    
    dify[0,1]=-(1-teta)**2*var_tran[0]
    ###############This should go up to  np.array(range(T-2))+2
    ###############Followine line also increased by 1 index
    for j in np.array(range(T-3))+2:
        dify[j-1,j]=-(1-teta)*var_tran[j-1]+teta*(1-teta)*var_tran[j-2]

    dify[T-2,T-3]=-(1-teta)**2*var_tran[T-3]
    
    for j in np.array(range(T-2))+2:
        dify[j-2,j]=-teta*var_tran[j-2]
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i)+i):
            dify[j,i-1]=dify[i-1,j]
            
    #/* This is the variance of Consumption */
    difcd[0,0]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[0]+varcsi
    difcd[1,1]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[1]+varcsi
    difcd[2,2]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[2]+varcsi
    for j in [3,4,5]:
        difcd[j,j]=ins_perm[0]**2*var_perm[j-2]+ins_tran[0]**2*var_tran[j]+varcsi

    for j in np.array(range(T-9))+6:
        difcd[j,j]=ins_perm[varying_ins]**2*var_perm[j-2]+ins_tran[varying_ins]**2*var_tran[j]+varcsi

    difcd[T-3,T-3]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+varcsi
    difcd[T-2,T-2]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+varcsi
    difcd[T-1,T-1]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+varcsi

    
    missing_v=(var_c_error[0]+var_c_error[1]+var_c_error[1]+var_c_error[2]+var_c_error[3]+var_c_error[4]+var_c_error[5]+var_c_error[6]+var_c_error[7]+var_c_error[8])/9.0
    difcme[0,0]=2*var_c_error[0]
    for j in np.array(range(6))+1:
        difcme[j,j]=var_c_error[j]+var_c_error[j-1]

    difcme[7,7]=var_c_error[6]+missing_v
    for j in np.array(range(3))+8:
        difcme[j,j]=2*missing_v

    difcme[11,11]=var_c_error[7]+missing_v
    difcme[12,12]=var_c_error[8]+var_c_error[7]
    difcme[13,13]=2*var_c_error[8]
    
    
    for j in np.array(range(8)):
        difcme[j,j+1]=-var_c_error[j]
    for j in np.array(range(3))+8:
        difcme[j,j+1]=-missing_v

    for j in [11,12]:
        difcme[j,j+1]=-var_c_error[j-4]
 
    difc=difcme+difcd
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i))+i:
            difc[j,i-1]=difc[i-1,j]

    #/* This is the Covariance of Income and Consumption */
    
    difyc[0,0]=ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[0]
    difyc[1,1]=ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[1]
    difyc[2,2]=ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[2]

    for j in [3,4,5]:
        difyc[j,j]=ins_perm[0]*var_perm[j-2]+ins_tran[0]*var_tran[j]

    for j in np.array(range(T-9))+6:
        difyc[j,j]=ins_perm[varying_ins]*var_perm[j-2]+ins_tran[varying_ins]*var_tran[j]

    difyc[T-3,T-3]=ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]
    difyc[T-2,T-2]=ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]
    difyc[T-1,T-1]=ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]

    for j in [1,2,3,4,5,6]:
        difyc[j-1,j]=-(1-teta)*ins_tran[0]*var_tran[j-1];
    
    for j in np.array(range(T-8))+7:
        difyc[j-1,j]=-(1-teta)*ins_tran[varying_ins]*var_tran[j-1]

    difyc[T-2,T-1]=-(1-teta)*ins_tran[varying_ins]*var_tran[T-3]

    for j in [2,3,4,5,6,7]:
        difyc[j-2,j]=-teta*ins_tran[0]*var_tran[j-2]
    
    for j in np.array(range(T-8))+8:
        difyc[j-2,j]=-teta*ins_tran[varying_ins]*var_tran[j-2]

    #/* Final matrix */
    
    dif[0:T,0:T]            =difc
    dif[T:2*(T),0:T]        =difyc
    dif[0:T,T:2*(T)]        =np.transpose(difyc)
    dif[T:2*(T),T:2*(T)]    =dify
    
    difa1 = np.concatenate((dif[0:8,:],dif[11:2*T,:]),0)
    difa2 = np.concatenate((difa1[:,0:8],difa1[:,11:2*T]),1)
    
    vech_indicies = np.tril_indices(np.shape(difa2)[0])
    fm=difa2[vech_indicies]

    return fm

def implied_cov_TimeAgg(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params):
    if ma==1:
        teta = params[0] 
    else:
        teta = 0
    if taste ==1:
        var_taste = params[ma] 
    else:
        var_taste = 0.0
    var_perm = params[ma+taste:ma+taste+perm_shk_params] 
    var_tran = params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] 
    ins_perm = params[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+1] 
    ins_tran = params[ma+taste+perm_shk_params+tran_shk_params+1:ma+taste+perm_shk_params+tran_shk_params+2] 
    var_c_error = params[ma+taste+perm_shk_params+tran_shk_params+2:ma+taste+perm_shk_params+tran_shk_params+2+meas_error_params] 

    dify  =np.zeros((T,T)) #/* Income */
    difcd =np.zeros((T,T)) #/* Consumption */
    difc  =np.zeros((T,T)) #/* Consumption */
    difcme=np.zeros((T,T)) #/* Measurement error of consumption */
    difyc =np.zeros((T,T)) #/* Cov Income Consumption */
    dif   =np.zeros((2*T,2*T))
    
    ##########################################
    #/* This is the variance of Income */
    dify[0,0]=2.0/3.0*var_perm[0]+  var_tran[0]+  (1-teta)**2*var_tran[0]+    teta**2*var_tran[0]
    dify[1,1]=2.0/3.0*var_perm[0]+  var_tran[1]+  (1-teta)**2*var_tran[0]+    teta**2*var_tran[0]
    dify[2,2]=2.0/3.0*var_perm[0]+  var_tran[2]+  (1-teta)**2*var_tran[1]+    teta**2*var_tran[0]
    for j in np.array(range(T-6))+3:
        dify[j,j]=1.0/3.0*var_perm[j-2]+1.0/3.0*var_perm[j-3]+  var_tran[j]+    (1-teta)**2*var_tran[j-1]+  teta**2*var_tran[j-2]
    dify[T-3,T-3]=1.0/3.0*var_perm[T-5]+1.0/3.0*var_perm[T-6]+  var_tran[T-3]+  (1-teta)**2*var_tran[T-4]+  teta**2*var_tran[T-5]
    dify[T-2,T-2]=2.0/3.0*var_perm[T-5]+                        var_tran[T-3]+  (1-teta)**2*var_tran[T-3]+  teta**2*var_tran[T-4]
    dify[T-1,T-1]=2.0/3.0*var_perm[T-5]+                        var_tran[T-3]+  (1-teta)**2*var_tran[T-3]+  teta**2*var_tran[T-3]
    
    dify[0,1]=1.0/6.0*var_perm[0] -(1-teta)*var_tran[0]+    teta*(1-teta)*var_tran[0]
    dify[1,2]=1.0/6.0*var_perm[0] -(1-teta)*var_tran[1]+    teta*(1-teta)*var_tran[0]
    dify[2,3]=1.0/6.0*var_perm[0] -(1-teta)*var_tran[2]+    teta*(1-teta)*var_tran[1]
    for j in np.array(range(T-5))+3:
        dify[j-1,j]=1.0/6.0*var_perm[j-3]   -(1-teta)*var_tran[j-1]+    teta*(1-teta)*var_tran[j-2]
    dify[T-3,T-2]=1.0/6.0*var_perm[T-5]     -(1-teta)*var_tran[T-3]+    teta*(1-teta)*var_tran[T-4]
    dify[T-2,T-1]=1.0/6.0*var_perm[T-5]     -(1-teta)*var_tran[T-3]+    teta*(1-teta)*var_tran[T-3]
    
    for j in np.array(range(T-2))+2:
        dify[j-2,j]=-teta*var_tran[j-2]
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i)+i):
            dify[j,i-1]=dify[i-1,j]
            
    ##########################################
    #/* This is the variance of Consumption */
    difcd[0,0]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[0]+var_taste
    difcd[1,1]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[1]+var_taste
    difcd[2,2]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[2]+var_taste
    for j in [3,4,5]:
        difcd[j,j]=ins_perm[0]**2*var_perm[j-2]+ins_tran[0]**2*var_tran[j]+var_taste

    for j in np.array(range(T-9))+6:
        difcd[j,j]=ins_perm[varying_ins]**2*var_perm[j-2]+ins_tran[varying_ins]**2*var_tran[j]+var_taste

    difcd[T-3,T-3]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste
    difcd[T-2,T-2]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste
    difcd[T-1,T-1]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste

    missing_v=(var_c_error[0]+var_c_error[1]+var_c_error[2]+var_c_error[3]+var_c_error[4]+var_c_error[5]+var_c_error[6]+var_c_error[7]+var_c_error[8])/9.0
    difcme[0,0]=2*var_c_error[0]
    for j in np.array(range(6))+1:
        difcme[j,j]=var_c_error[j]+var_c_error[j-1]

    difcme[7,7]=var_c_error[6]+missing_v
    for j in np.array(range(3))+8:
        difcme[j,j]=2*missing_v

    difcme[11,11]=var_c_error[7]+missing_v
    difcme[12,12]=var_c_error[8]+var_c_error[7]
    difcme[13,13]=2*var_c_error[8]
    
    
    for j in np.array(range(8)):
        difcme[j,j+1]=-var_c_error[j]
    for j in np.array(range(3))+8:
        difcme[j,j+1]=-missing_v

    for j in [11,12]:
        difcme[j,j+1]=-var_c_error[j-4]
 
    difc=difcme+difcd
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i))+i:
            difc[j,i-1]=difc[i-1,j]
            
    ##########################################
    #/* This is the Covariance of Income and Consumption */
    difyc[0,0]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[0]
    difyc[1,1]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[1]
    difyc[2,2]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[2]

    for j in [3,4,5]:
        difyc[j,j]=1.0/2.0*ins_perm[0]*var_perm[j-2]+ins_tran[0]*var_tran[j]

    for j in np.array(range(T-9))+6:
        difyc[j,j]=1.0/2.0*ins_perm[varying_ins]*var_perm[j-2]+ins_tran[varying_ins]*var_tran[j]

    difyc[T-3,T-3]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]
    difyc[T-2,T-2]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]
    difyc[T-1,T-1]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]

    difyc[0,1]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-teta)*ins_tran[0]*var_tran[0]
    difyc[1,2]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-teta)*ins_tran[0]*var_tran[1]
    difyc[2,3]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-teta)*ins_tran[0]*var_tran[2]
    for j in [4,5,6]:
        difyc[j-1,j]=1.0/2.0*ins_perm[0]*var_perm[j-3] - (1-teta)*ins_tran[0]*var_tran[j-1]
    
    for j in np.array(range(T-8))+7:
        difyc[j-1,j]=1.0/2.0*ins_perm[varying_ins]*var_perm[j-3] - (1-teta)*ins_tran[varying_ins]*var_tran[j-1]

    difyc[T-2,T-1]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5] - (1-teta)*ins_tran[varying_ins]*var_tran[T-3]
    
    for j in [2,3,4,5,6,7]:
        difyc[j-2,j]=-teta*ins_tran[0]*var_tran[j-2]
    
    for j in np.array(range(T-8))+8:
        difyc[j-2,j]=-teta*ins_tran[varying_ins]*var_tran[j-2]
        
    ##########################################
            
    #/* Final matrix */
    dif[0:T,0:T]            =difc
    dif[T:2*(T),0:T]        =difyc
    dif[0:T,T:2*(T)]        =np.transpose(difyc)
    dif[T:2*(T),T:2*(T)]    =dify
    
    difa1 = np.concatenate((dif[0:8,:],dif[11:2*T,:]),0)
    difa2 = np.concatenate((difa1[:,0:8],difa1[:,11:2*T]),1)
    
    vech_indicies = np.tril_indices(np.shape(difa2)[0])
    fm=difa2[vech_indicies]

    return fm


# This assumes transitory shocks last for period tau and decay linearly in this time
def implied_cov_TimeAgg3(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params):
    if ma==1:
        tau = params[0] 
    else:
        tau = 0
    if taste ==1:
        var_taste = params[ma] 
    else:
        var_taste = 0.0
    var_perm = params[ma+taste:ma+taste+perm_shk_params] 
    var_tran = params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] 
    ins_perm = params[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+1] 
    ins_tran = params[ma+taste+perm_shk_params+tran_shk_params+1:ma+taste+perm_shk_params+tran_shk_params+2] 
    var_c_error = params[ma+taste+perm_shk_params+tran_shk_params+2:ma+taste+perm_shk_params+tran_shk_params+2+meas_error_params] 

    dify  =np.zeros((T,T)) #/* Income */
    difcd =np.zeros((T,T)) #/* Consumption */
    difc  =np.zeros((T,T)) #/* Consumption */
    difcme=np.zeros((T,T)) #/* Measurement error of consumption */
    difyc =np.zeros((T,T)) #/* Cov Income Consumption */
    dif   =np.zeros((2*T,2*T))
    
    ##########################################
    #/* This is the variance of Income */
    dify[0,0]=2.0/3.0*var_perm[0]+  (1.0-7.0/15.0*tau)*var_tran[0] + (1.0-8.0/15.0*tau)*var_tran[0]  +(1.0/5.0*tau)*var_tran[0]   
    dify[1,1]=2.0/3.0*var_perm[0]+  (1.0-7.0/15.0*tau)*var_tran[1] + (1.0-8.0/15.0*tau)*var_tran[0]  +(1.0/5.0*tau)*var_tran[0]  
    dify[2,2]=2.0/3.0*var_perm[0]+  (1.0-7.0/15.0*tau)*var_tran[2] + (1.0-8.0/15.0*tau)*var_tran[1]  +(1.0/5.0*tau)*var_tran[0]  
    for j in np.array(range(T-6))+3:
        dify[j,j]=1.0/3.0*var_perm[j-2]+1.0/3.0*var_perm[j-3]+  (1.0-7.0/15.0*tau)*var_tran[j]   + (1.0-8.0/15.0*tau)*var_tran[j-1]  +(1.0/5.0*tau)*var_tran[j-2]   
    dify[T-3,T-3]=1.0/3.0*var_perm[T-5]+1.0/3.0*var_perm[T-6]+  (1.0-7.0/15.0*tau)*var_tran[T-3] + (1.0-8.0/15.0*tau)*var_tran[T-4]  +(1.0/5.0*tau)*var_tran[T-5]  
    dify[T-2,T-2]=2.0/3.0*var_perm[T-5]+                        (1.0-7.0/15.0*tau)*var_tran[T-3] + (1.0-8.0/15.0*tau)*var_tran[T-3]  +(1.0/5.0*tau)*var_tran[T-4]  
    dify[T-1,T-1]=2.0/3.0*var_perm[T-5]+                        (1.0-7.0/15.0*tau)*var_tran[T-3] + (1.0-8.0/15.0*tau)*var_tran[T-3]  +(1.0/5.0*tau)*var_tran[T-3]  
    
    dify[0,1]=1.0/6.0*var_perm[0] +(-2.0/5.0*tau - (1-tau))*var_tran[0] - 1.0/15.0*tau*var_tran[0]
    dify[1,2]=1.0/6.0*var_perm[0] +(-2.0/5.0*tau - (1-tau))*var_tran[1] - 1.0/15.0*tau*var_tran[0]
    dify[2,3]=1.0/6.0*var_perm[0] +(-2.0/5.0*tau - (1-tau))*var_tran[2] - 1.0/15.0*tau*var_tran[1]
    for j in np.array(range(T-5))+3:
        dify[j-1,j]=1.0/6.0*var_perm[j-3]   +(-2.0/5.0*tau - (1-tau))*var_tran[j-1] - 1.0/15.0*tau*var_tran[j-1]
    dify[T-3,T-2]=1.0/6.0*var_perm[T-5]     +(-2.0/5.0*tau - (1-tau))*var_tran[T-3] - 1.0/15.0*tau*var_tran[T-4]
    dify[T-2,T-1]=1.0/6.0*var_perm[T-5]     +(-2.0/5.0*tau - (1-tau))*var_tran[T-3] - 1.0/15.0*tau*var_tran[T-3]
    
    for j in np.array(range(T-2))+2:
        dify[j-2,j]=-2.0/15.0*tau*var_tran[j-2]
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i)+i):
            dify[j,i-1]=dify[i-1,j]
            
    ##########################################
    #/* This is the variance of Consumption */
    difcd[0,0]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[0]+var_taste
    difcd[1,1]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[1]+var_taste
    difcd[2,2]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[2]+var_taste
    for j in [3,4,5]:
        difcd[j,j]=ins_perm[0]**2*var_perm[j-2]+ins_tran[0]**2*var_tran[j]+var_taste

    for j in np.array(range(T-9))+6:
        difcd[j,j]=ins_perm[varying_ins]**2*var_perm[j-2]+ins_tran[varying_ins]**2*var_tran[j]+var_taste

    difcd[T-3,T-3]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste
    difcd[T-2,T-2]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste
    difcd[T-1,T-1]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste

    missing_v=(var_c_error[0]+var_c_error[1]+var_c_error[2]+var_c_error[3]+var_c_error[4]+var_c_error[5]+var_c_error[6]+var_c_error[7]+var_c_error[8])/9.0
    difcme[0,0]=2*var_c_error[0]
    for j in np.array(range(6))+1:
        difcme[j,j]=var_c_error[j]+var_c_error[j-1]

    difcme[7,7]=var_c_error[6]+missing_v
    for j in np.array(range(3))+8:
        difcme[j,j]=2*missing_v

    difcme[11,11]=var_c_error[7]+missing_v
    difcme[12,12]=var_c_error[8]+var_c_error[7]
    difcme[13,13]=2*var_c_error[8]
    
    
    for j in np.array(range(8)):
        difcme[j,j+1]=-var_c_error[j]
    for j in np.array(range(3))+8:
        difcme[j,j+1]=-missing_v

    for j in [11,12]:
        difcme[j,j+1]=-var_c_error[j-4]
 
    difc=difcme+difcd
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i))+i:
            difc[j,i-1]=difc[i-1,j]
            
    ##########################################
    #/* This is the Covariance of Income and Consumption */
    difyc[0,0]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[0]*(1.0-1.0/3.0*tau)
    difyc[1,1]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[1]*(1.0-1.0/3.0*tau)
    difyc[2,2]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[2]*(1.0-1.0/3.0*tau)

    for j in [3,4,5]:
        difyc[j,j]=1.0/2.0*ins_perm[0]*var_perm[j-2]+ins_tran[0]*var_tran[j]*(1.0-1.0/3.0*tau)

    for j in np.array(range(T-9))+6:
        difyc[j,j]=1.0/2.0*ins_perm[varying_ins]*var_perm[j-2]+ins_tran[varying_ins]*var_tran[j]*(1.0-1.0/3.0*tau)

    difyc[T-3,T-3]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]*(1.0-1.0/3.0*tau)
    difyc[T-2,T-2]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]*(1.0-1.0/3.0*tau)
    difyc[T-1,T-1]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]*(1.0-1.0/3.0*tau)

    difyc[0,1]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-2.0/3.0*tau)*ins_tran[0]*var_tran[0]
    difyc[1,2]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-2.0/3.0*tau)*ins_tran[0]*var_tran[1]
    difyc[2,3]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-2.0/3.0*tau)*ins_tran[0]*var_tran[2]
    for j in [4,5,6]:
        difyc[j-1,j]=1.0/2.0*ins_perm[0]*var_perm[j-3] - (1-2.0/3.0*tau)*ins_tran[0]*var_tran[j-1]
    
    for j in np.array(range(T-8))+7:
        difyc[j-1,j]=1.0/2.0*ins_perm[varying_ins]*var_perm[j-3] - (1-2.0/3.0*tau)*ins_tran[varying_ins]*var_tran[j-1]

    difyc[T-2,T-1]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5] - (1-2.0/3.0*tau)*ins_tran[varying_ins]*var_tran[T-3]
    
    for j in [2,3,4,5,6,7]:
        difyc[j-2,j]=-1.0/5.0*tau*ins_tran[0]*var_tran[j-2]
    
    for j in np.array(range(T-8))+8:
        difyc[j-2,j]=-1.0/5.0*tau*ins_tran[varying_ins]*var_tran[j-2]
        
    ##########################################
            
    #/* Final matrix */
    dif[0:T,0:T]            =difc
    dif[T:2*(T),0:T]        =difyc
    dif[0:T,T:2*(T)]        =np.transpose(difyc)
    dif[T:2*(T),T:2*(T)]    =dify
    
    difa1 = np.concatenate((dif[0:8,:],dif[11:2*T,:]),0)
    difa2 = np.concatenate((difa1[:,0:8],difa1[:,11:2*T]),1)
    
    vech_indicies = np.tril_indices(np.shape(difa2)[0])
    fm=difa2[vech_indicies]

    return fm


# This assumes transitory shocks last for period tau and are uniformly distributed in this time
def implied_cov_TimeAgg2(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params):
    if ma==1:
        tau = params[0] 
    else:
        tau = 0
    if taste ==1:
        var_taste = params[ma] 
    else:
        var_taste = 0.0
    var_perm = params[ma+taste:ma+taste+perm_shk_params] 
    var_tran = params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] 
    ins_perm = params[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+1] 
    ins_tran = params[ma+taste+perm_shk_params+tran_shk_params+1:ma+taste+perm_shk_params+tran_shk_params+2] 
    var_c_error = params[ma+taste+perm_shk_params+tran_shk_params+2:ma+taste+perm_shk_params+tran_shk_params+2+meas_error_params] 

    dify  =np.zeros((T,T)) #/* Income */
    difcd =np.zeros((T,T)) #/* Consumption */
    difc  =np.zeros((T,T)) #/* Consumption */
    difcme=np.zeros((T,T)) #/* Measurement error of consumption */
    difyc =np.zeros((T,T)) #/* Cov Income Consumption */
    dif   =np.zeros((2*T,2*T))
    
    ##########################################
    #/* This is the variance of Income */
    dify[0,0]=2.0/3.0*var_perm[0]+  (1.0/3.0*tau +1.0-tau)*var_tran[0] + (1.0/3.0*tau +1.0-tau)*var_tran[0]  +(1.0/3.0*tau)*var_tran[0]   
    dify[1,1]=2.0/3.0*var_perm[0]+  (1.0/3.0*tau +1.0-tau)*var_tran[1] + (1.0/3.0*tau +1.0-tau)*var_tran[0]  +(1.0/3.0*tau)*var_tran[0]  
    dify[2,2]=2.0/3.0*var_perm[0]+  (1.0/3.0*tau +1.0-tau)*var_tran[2] + (1.0/3.0*tau +1.0-tau)*var_tran[1]  +(1.0/3.0*tau)*var_tran[0]  
    for j in np.array(range(T-6))+3:
        dify[j,j]=1.0/3.0*var_perm[j-2]+1.0/3.0*var_perm[j-3]+  (1.0/3.0*tau +1.0-tau)*var_tran[j]   + (1.0/3.0*tau +1.0-tau)*var_tran[j-1]  +(1.0/3.0*tau)*var_tran[j-2]   
    dify[T-3,T-3]=1.0/3.0*var_perm[T-5]+1.0/3.0*var_perm[T-6]+  (1.0/3.0*tau +1.0-tau)*var_tran[T-3] + (1.0/3.0*tau +1.0-tau)*var_tran[T-4]  +(1.0/3.0*tau)*var_tran[T-5]  
    dify[T-2,T-2]=2.0/3.0*var_perm[T-5]+                        (1.0/3.0*tau +1.0-tau)*var_tran[T-3] + (1.0/3.0*tau +1.0-tau)*var_tran[T-3]  +(1.0/3.0*tau)*var_tran[T-4]  
    dify[T-1,T-1]=2.0/3.0*var_perm[T-5]+                        (1.0/3.0*tau +1.0-tau)*var_tran[T-3] + (1.0/3.0*tau +1.0-tau)*var_tran[T-3]  +(1.0/3.0*tau)*var_tran[T-3]  
    
    dify[0,1]=1.0/6.0*var_perm[0] +(-1.0/6.0*tau - (1-tau))*var_tran[0] - 1.0/6.0*tau*var_tran[0]
    dify[1,2]=1.0/6.0*var_perm[0] +(-1.0/6.0*tau - (1-tau))*var_tran[1] - 1.0/6.0*tau*var_tran[0]
    dify[2,3]=1.0/6.0*var_perm[0] +(-1.0/6.0*tau - (1-tau))*var_tran[2] - 1.0/6.0*tau*var_tran[1]
    for j in np.array(range(T-5))+3:
        dify[j-1,j]=1.0/6.0*var_perm[j-3]   +(-1.0/6.0*tau - (1-tau))*var_tran[j-1] - 1.0/6.0*tau*var_tran[j-1]
    dify[T-3,T-2]=1.0/6.0*var_perm[T-5]     +(-1.0/6.0*tau - (1-tau))*var_tran[T-3] - 1.0/6.0*tau*var_tran[T-4]
    dify[T-2,T-1]=1.0/6.0*var_perm[T-5]     +(-1.0/6.0*tau - (1-tau))*var_tran[T-3] - 1.0/6.0*tau*var_tran[T-3]
    
    for j in np.array(range(T-2))+2:
        dify[j-2,j]=-1.0/6.0*tau*var_tran[j-2]
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i)+i):
            dify[j,i-1]=dify[i-1,j]
            
    ##########################################
    #/* This is the variance of Consumption */
    difcd[0,0]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[0]+var_taste
    difcd[1,1]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[1]+var_taste
    difcd[2,2]=ins_perm[0]**2*var_perm[0]+ins_tran[0]**2*var_tran[2]+var_taste
    for j in [3,4,5]:
        difcd[j,j]=ins_perm[0]**2*var_perm[j-2]+ins_tran[0]**2*var_tran[j]+var_taste

    for j in np.array(range(T-9))+6:
        difcd[j,j]=ins_perm[varying_ins]**2*var_perm[j-2]+ins_tran[varying_ins]**2*var_tran[j]+var_taste

    difcd[T-3,T-3]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste
    difcd[T-2,T-2]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste
    difcd[T-1,T-1]=ins_perm[varying_ins]**2*var_perm[T-5]+ins_tran[varying_ins]**2*var_tran[T-3]+var_taste

    missing_v=(var_c_error[0]+var_c_error[1]+var_c_error[2]+var_c_error[3]+var_c_error[4]+var_c_error[5]+var_c_error[6]+var_c_error[7]+var_c_error[8])/9.0
    difcme[0,0]=2*var_c_error[0]
    for j in np.array(range(6))+1:
        difcme[j,j]=var_c_error[j]+var_c_error[j-1]

    difcme[7,7]=var_c_error[6]+missing_v
    for j in np.array(range(3))+8:
        difcme[j,j]=2*missing_v

    difcme[11,11]=var_c_error[7]+missing_v
    difcme[12,12]=var_c_error[8]+var_c_error[7]
    difcme[13,13]=2*var_c_error[8]
    
    
    for j in np.array(range(8)):
        difcme[j,j+1]=-var_c_error[j]
    for j in np.array(range(3))+8:
        difcme[j,j+1]=-missing_v

    for j in [11,12]:
        difcme[j,j+1]=-var_c_error[j-4]
 
    difc=difcme+difcd
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i))+i:
            difc[j,i-1]=difc[i-1,j]
            
    ##########################################
    #/* This is the Covariance of Income and Consumption */
    difyc[0,0]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[0]*(1.0-0.5*tau)
    difyc[1,1]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[1]*(1.0-0.5*tau)
    difyc[2,2]=1.0/2.0*ins_perm[0]*var_perm[0]+ins_tran[0]*var_tran[2]*(1.0-0.5*tau)

    for j in [3,4,5]:
        difyc[j,j]=1.0/2.0*ins_perm[0]*var_perm[j-2]+ins_tran[0]*var_tran[j]*(1.0-0.5*tau)

    for j in np.array(range(T-9))+6:
        difyc[j,j]=1.0/2.0*ins_perm[varying_ins]*var_perm[j-2]+ins_tran[varying_ins]*var_tran[j]*(1.0-0.5*tau)

    difyc[T-3,T-3]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]*(1.0-0.5*tau)
    difyc[T-2,T-2]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]*(1.0-0.5*tau)
    difyc[T-1,T-1]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5]+ins_tran[varying_ins]*var_tran[T-3]*(1.0-0.5*tau)

    difyc[0,1]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-tau)*ins_tran[0]*var_tran[0]
    difyc[1,2]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-tau)*ins_tran[0]*var_tran[1]
    difyc[2,3]=1.0/2.0*ins_perm[0]*var_perm[0] - (1-tau)*ins_tran[0]*var_tran[2]
    for j in [4,5,6]:
        difyc[j-1,j]=1.0/2.0*ins_perm[0]*var_perm[j-3] - (1-tau)*ins_tran[0]*var_tran[j-1]
    
    for j in np.array(range(T-8))+7:
        difyc[j-1,j]=1.0/2.0*ins_perm[varying_ins]*var_perm[j-3] - (1-tau)*ins_tran[varying_ins]*var_tran[j-1]

    difyc[T-2,T-1]=1.0/2.0*ins_perm[varying_ins]*var_perm[T-5] - (1-tau)*ins_tran[varying_ins]*var_tran[T-3]
    
    for j in [2,3,4,5,6,7]:
        difyc[j-2,j]=-0.5*tau*ins_tran[0]*var_tran[j-2]
    
    for j in np.array(range(T-8))+8:
        difyc[j-2,j]=-0.5*tau*ins_tran[varying_ins]*var_tran[j-2]
        
    ##########################################
            
    #/* Final matrix */
    dif[0:T,0:T]            =difc
    dif[T:2*(T),0:T]        =difyc
    dif[0:T,T:2*(T)]        =np.transpose(difyc)
    dif[T:2*(T),T:2*(T)]    =dify
    
    difa1 = np.concatenate((dif[0:8,:],dif[11:2*T,:]),0)
    difa2 = np.concatenate((difa1[:,0:8],difa1[:,11:2*T]),1)
    
    vech_indicies = np.tril_indices(np.shape(difa2)[0])
    fm=difa2[vech_indicies]

    return fm


def Parameter_estimation(model, c_vector, omega, T, ma=1, taste=1, varying_ins=0):
    '''
    Replicates table 6 from BPP
    
    Parameters
    ----------
    model   : string
        takes values 'BPP' to replicate BPP method exactly, or 'TimeAgg' to do
        the time aggregated version
    c_vector : np.array
        Vector containing empirical moments
    omega : np.array
        Empirical covariance matrix for the empirical moments
    T : int
        Length of panel
    ma : int
        1 -> include moving average component, 0->don't
    taste : int
        1 -> include taste shocks, 0->don't
    varying_ins : int
        1 -> allow for insurance parameters to change in 1985, 0->don't
    Returns
    -------
    var_perm : np.array
        Array of permanent shock variances
    var_perm_se : np.array
        Array of standard errors for permanent shock variances
    var_tran : np.array
        Array of transitory shock variances
    var_tran_se : np.array
        Array of standard errors for transitory shock variances
    ins_perm : np.array
        Array of permanent shock insurance
    ins_perm_se : np.array
        Array of standard errors for permanent insurance
    ins_tran : np.array
        Array of transitory shock insurance
    ins_tran_se : np.array
        Array of standard errors for transitory shock insurance
    var_c_error : np.array
        Array of consumption measurement error variances
    var_c_error_se : np.array
        Array of standard errors for consumption measurement error variances
    '''
    if model=='BPP':
        implied_cov = lambda params, ma, taste, varying_ins, T, perm_shk_params, \
                            tran_shk_params, perm_ins_params,tran_ins_params,\
                            meas_error_params : implied_cov_BPP(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
    if model=='TimeAgg':
        implied_cov = lambda params, ma, taste, varying_ins, T, perm_shk_params, \
                            tran_shk_params, perm_ins_params,tran_ins_params,\
                            meas_error_params : implied_cov_TimeAgg3(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
        
    if model=='TimeAgg_twoshot':
        implied_cov = lambda params, ma, taste, varying_ins, T, perm_shk_params, \
                            tran_shk_params, perm_ins_params,tran_ins_params,\
                            meas_error_params : implied_cov_TimeAgg(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)

    if model=='TimeAgg_uniform':
        implied_cov = lambda params, ma, taste, varying_ins, T, perm_shk_params, \
                            tran_shk_params, perm_ins_params,tran_ins_params,\
                            meas_error_params : implied_cov_TimeAgg2(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)

    if model=='TimeAgg_lineardecay':
        implied_cov = lambda params, ma, taste, varying_ins, T, perm_shk_params, \
                            tran_shk_params, perm_ins_params,tran_ins_params,\
                            meas_error_params : implied_cov_TimeAgg3(params, ma, taste, varying_ins, T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)

    
    # get number of parameters of each type
    perm_shk_params = T-4  # time varying permanent shock variance
    tran_shk_params = T-2
    perm_ins_params = 1+varying_ins
    tran_ins_params = 1+varying_ins
    meas_error_params = 9   #time-varying measurement error variance
    
    num_params = ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params+meas_error_params
    
    init_params = np.zeros(num_params)
    
    if ma==1:
        init_params[0] = 0.1   #teta, ma component of income process
    if taste:
        init_params[ma] = 0.01  #variance of taste shocks
    init_params[ma+taste:ma+taste+perm_shk_params] = 0.03*np.ones(perm_shk_params)
    init_params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] = 0.03*np.ones(tran_shk_params)
    init_params[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params] = 1.0*np.ones(perm_ins_params)
    init_params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params] = 0.3*np.ones(tran_ins_params)
    init_params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params+meas_error_params] = 0.06*np.ones(meas_error_params)
    
    def objectiveFun(params, ma, taste, varying_ins, T, empirical_cov, weight_matrix):
        model_cov = implied_cov(params, ma, taste, varying_ins,T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
        distance = np.dot(np.dot((model_cov-empirical_cov), weight_matrix),(model_cov-empirical_cov))
        distance = distance #+ 10000*np.std(params[ma+taste:ma+taste+perm_shk_params]) #add in this to keep same variance for permanent shocks over the whole time period
        return distance
    
    # Define the weight matrix as Equal Weight Minimum Distance
    weight_matrix = np.diag(np.diag(omega)**(-1))
    
    ret = objectiveFun(init_params, ma, taste, varying_ins, T, c_vector, weight_matrix)
    
    #Solve with one method, reset init_params, then solve again. Seems to converge OK.
    solved_objective1 = minimize(objectiveFun, init_params, args=(ma, taste, varying_ins, T, c_vector, weight_matrix))  
    init_params2 = solved_objective1.x
    solved_objective = minimize(objectiveFun, init_params2, args=(ma, taste, varying_ins, T, c_vector, weight_matrix),method='Nelder-Mead')
     
    solved_params = solved_objective.x
    
    fun_for_jacob = lambda params: implied_cov(params, ma, taste, varying_ins,T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
    
    jacob = nd.Jacobian(fun_for_jacob)(solved_params)
    
    Sandwich1 = inv(np.dot(np.transpose(jacob),np.dot(weight_matrix,jacob)))
    Sandwich2 = np.dot(np.transpose(jacob),np.dot(weight_matrix,np.dot(omega,np.dot(weight_matrix,jacob))))
    cov_params = np.dot(Sandwich1,np.dot(Sandwich2,Sandwich1))
    standard_errors = np.diag(cov_params)**0.5
    
    #extract relevant numbers
    
    if ma==1:
        teta = solved_params[0] 
        teta_se = standard_errors[0] 
    else:
        teta = 0.0
        teta_se = 0.0
    if taste:
        varcsi = solved_params[ma] 
        varcsi_se = standard_errors[ma] 
    else:
        varcsi = 0.0
        varcsi_se = 0.0 
    var_perm = solved_params[ma+taste:ma+taste+perm_shk_params] 
    var_tran = solved_params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] 
    ins_perm = solved_params[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params] 
    ins_tran = solved_params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params] 
    var_c_error = solved_params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params+meas_error_params] 
    
    var_perm_se = standard_errors[ma+taste:ma+taste+perm_shk_params] 
    var_tran_se = standard_errors[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] 
    ins_perm_se = standard_errors[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params] 
    ins_tran_se = standard_errors[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params] 
    var_c_error_se = standard_errors[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params+meas_error_params] 
    
    return var_perm, var_perm_se, var_tran, var_tran_se, ins_perm, ins_perm_se, ins_tran, ins_tran_se, var_c_error, var_c_error_se, teta, teta_se, varcsi, varcsi_se
    
    
