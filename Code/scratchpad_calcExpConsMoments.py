# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 10:43:33 2019

@author: edmun
"""
num_periods = 20000
num_divisions = 200
total_obs = num_periods*num_divisions
max_decay = 100000
#np.random.seed(seed=3)
#########################Parameters
theta = 0.75
Omega = 0.25
ins_tran_test = 0.5
#########################

#tran_shocks = (var_tran_test*num_divisions)**0.5*np.random.normal(size=total_obs)
tran_shocks = np.random.normal(size=total_obs)/num_divisions**0.5

tran_cons = np.zeros_like(tran_shocks[1:])
d_tran_inc = np.zeros_like(tran_cons)



decay     = np.exp(-np.arange(total_obs)*theta/num_divisions)
decay_inc = np.exp(-np.arange(total_obs)*Omega/num_divisions)

observed_cons = np.zeros(num_periods)
for T1 in range(num_periods):
    observed_cons[T1] =  np.sum(theta*ins_tran_test/(1.0-np.exp(-theta)) * decay[0:min((T1+1)*num_divisions,max_decay)] * np.flip(tran_shocks[max(0,(T1+1)*num_divisions-max_decay):(T1+1)*num_divisions],0))

d_tran_inc_0_1 = np.zeros(num_periods)*np.nan
d_tran_inc_1_2 = np.zeros(num_periods)*np.nan
d_tran_inc_1_2_a = np.zeros(num_periods)*np.nan
d_tran_inc_1_2_b = np.zeros(num_periods)*np.nan
d_tran_inc_2_3 = np.zeros(num_periods)*np.nan
d_tran_inc_2_infty = np.zeros(num_periods)*np.nan
d_tran_inc_3_4 = np.zeros(num_periods)*np.nan
d_tran_inc_3_infty = np.zeros(num_periods)*np.nan
d_tran_inc_4_infty = np.zeros(num_periods)*np.nan
observed_inc = np.zeros(num_periods)*np.nan

tran_inc_0_1 = np.zeros(num_periods)*np.nan
tran_inc_1_2 = np.zeros(num_periods)*np.nan
tran_inc_2_3 = np.zeros(num_periods)*np.nan
tran_inc_2_infty = np.zeros(num_periods)*np.nan
tran_inc_3_4 = np.zeros(num_periods)*np.nan
for T1 in range(num_periods):
    d_tran_inc_0_1[T1]     = 1.0/(1-np.exp(-Omega))*np.sum((1-np.flip(decay_inc[0:num_divisions] ,0)) * tran_shocks[T1*num_divisions:(T1+1)*num_divisions])
    if T1>=1:
        d_tran_inc_1_2[T1] = 1.0/(1-np.exp(-Omega))*np.sum( ((2-np.exp(-Omega))*np.flip(decay_inc[0:num_divisions] ,0)  -1)  * tran_shocks[(T1-1)*num_divisions:(T1)*num_divisions])
        d_tran_inc_1_2_a[T1] = 1.0/(1-np.exp(-Omega))*np.sum( ((2-np.exp(-Omega))*np.flip(decay_inc[0:num_divisions] ,0)  )  * tran_shocks[(T1-1)*num_divisions:(T1)*num_divisions])
        d_tran_inc_1_2_b[T1] = 1.0/(1-np.exp(-Omega))*np.sum(   -1  * tran_shocks[(T1-1)*num_divisions:(T1)*num_divisions])
        tran_inc_1_2[T1] = np.sum( np.flip(decay_inc[0:num_divisions] ,0)   * tran_shocks[(T1-1)*num_divisions:(T1)*num_divisions])
    if T1>=2:
        d_tran_inc_2_infty[T1] = -(1-np.exp(-Omega))*np.sum( np.flip(decay_inc[0:min(max_decay,(T1-1)*num_divisions)] ,0)   * tran_shocks[max(0,(T1-1)*num_divisions-max_decay):(T1-1)*num_divisions])
        d_tran_inc_2_3[T1] = -(1-np.exp(-Omega))*np.sum( np.flip(decay_inc[0:num_divisions]                      ,0)   * tran_shocks[(T1-2)*num_divisions:(T1-1)*num_divisions])
        observed_inc[T1] = d_tran_inc_0_1[T1] + np.sum( np.flip(decay_inc[0:min(max_decay,T1*num_divisions)] ,0)   * tran_shocks[max(0,T1*num_divisions-max_decay):(T1)*num_divisions])
    if T1>=3:
        d_tran_inc_3_infty[T1] = -(1-np.exp(-Omega))*np.sum( np.flip(decay_inc[0:min(max_decay,(T1-2)*num_divisions)] ,0)   * tran_shocks[max(0,(T1-2)*num_divisions-max_decay):(T1-2)*num_divisions])*np.exp(-Omega)
        d_tran_inc_3_4[T1] = -(1-np.exp(-Omega))*np.sum( np.flip(decay_inc[0:num_divisions]                      ,0)   * tran_shocks[(T1-3)*num_divisions:(T1-2)*num_divisions])*np.exp(-Omega)
    if T1>=4:
        d_tran_inc_4_infty[T1] = -(1-np.exp(-Omega))*np.sum( np.flip(decay_inc[0:min(max_decay,(T1-3)*num_divisions)] ,0)   * tran_shocks[max(0,(T1-3)*num_divisions-max_decay):(T1-3)*num_divisions])*np.exp(-2*Omega)
tran_inc_0_1 = d_tran_inc_0_1

var_d_tran_inc_0_1 = np.nanmean(d_tran_inc_0_1**2)
var_d_tran_inc_1_2 = np.nanmean(d_tran_inc_1_2**2)
var_d_tran_inc_2_infty = np.nanmean(d_tran_inc_2_infty**2)
foo = [var_d_tran_inc_0_1,var_d_tran_inc_1_2,var_d_tran_inc_2_infty]
delta_obs_inc = observed_inc[1:]-observed_inc[:-1]
delta_obs_cons = observed_cons[1:]-observed_cons[:-1]

simulated_y_cov = np.array([np.nanmean(delta_obs_inc**2), 
np.nanmean(delta_obs_inc[1:]*delta_obs_inc[:-1]), 
np.nanmean(delta_obs_inc[2:]*delta_obs_inc[:-2]), 
np.nanmean(delta_obs_inc[3:]*delta_obs_inc[:-3]), 
np.nanmean(delta_obs_inc[4:]*delta_obs_inc[:-4]),
np.nanmean(delta_obs_inc[5:]*delta_obs_inc[:-5]),
np.nanmean(delta_obs_inc[6:]*delta_obs_inc[:-6]),
np.nanmean(delta_obs_inc[7:]*delta_obs_inc[:-7])])

e1t = np.exp(-theta)
e2t = np.exp(-2.0*theta)
e3t = np.exp(-3.0*theta)
e1o = np.exp(-Omega)
e2o = np.exp(-2.0*Omega)
e3o = np.exp(-3.0*Omega)
theoretical_y_cov = np.array([2.0/((1.0-e1o)**2) - (3-e1o)/(Omega*(1-e1o)),
1.0/(2*Omega) * (2.0-e1o) - 1.0/(1-e1o)**2 * (1-(1-e1o)/Omega),
-1.0/Omega*(1-e1o)         +1.0/(2.0*Omega)*(1.0-e2o)         - (2-e1o)*e1o     /(2.0*Omega)*(1.0-e2o) + 1.0/Omega*e1o     *(1.0-e1o)   + (1.0 -e1o)**2/(2.0*Omega)*e2o,
-1.0/Omega*(1-e1o)*e1o     +1.0/(2.0*Omega)*(1.0-e2o)*e1o     - (2-e1o)*e2o     /(2.0*Omega)*(1.0-e2o) + 1.0/Omega*e2o     *(1.0-e1o)   + (1.0 -e1o)**2/(2.0*Omega)*e3o,
-1.0/Omega*(1-e1o)*e2o     +1.0/(2.0*Omega)*(1.0-e2o)*e2o     - (2-e1o)*e3o     /(2.0*Omega)*(1.0-e2o) + 1.0/Omega*e3o     *(1.0-e1o)   + (1.0 -e1o)**2/(2.0*Omega)*e3o*e1o,
-1.0/Omega*(1-e1o)*e3o     +1.0/(2.0*Omega)*(1.0-e2o)*e3o     - (2-e1o)*e1o*e3o /(2.0*Omega)*(1.0-e2o) + 1.0/Omega*e3o*e1o *(1.0-e1o)   + (1.0 -e1o)**2/(2.0*Omega)*e3o*e2o,
-1.0/Omega*(1-e1o)*e3o*e1o +1.0/(2.0*Omega)*(1.0-e2o)*e3o*e1o - (2-e1o)*e2o*e3o /(2.0*Omega)*(1.0-e2o) + 1.0/Omega*e3o*e2o *(1.0-e1o)   + (1.0 -e1o)**2/(2.0*Omega)*e3o*e3o,
-1.0/Omega*(1-e1o)*e3o*e2o +1.0/(2.0*Omega)*(1.0-e2o)*e3o*e2o - (2-e1o)*e3o*e3o /(2.0*Omega)*(1.0-e2o) + 1.0/Omega*e3o*e3o *(1.0-e1o)   + (1.0 -e1o)**2/(2.0*Omega)*e3o*e3o*e1o])

simulated_c_cov = np.array([np.nanmean(delta_obs_cons**2), 
np.nanmean(delta_obs_cons[1:]*delta_obs_cons[:-1]), 
np.nanmean(delta_obs_cons[2:]*delta_obs_cons[:-2]), 
np.nanmean(delta_obs_cons[3:]*delta_obs_cons[:-3]), 
np.nanmean(delta_obs_cons[4:]*delta_obs_cons[:-4]),
np.nanmean(delta_obs_cons[5:]*delta_obs_cons[:-5]),
np.nanmean(delta_obs_cons[6:]*delta_obs_cons[:-6]),
np.nanmean(delta_obs_cons[7:]*delta_obs_cons[:-7])])

theoretical_c_cov = np.array([ins_tran_test**2*theta/(1-e1t),
                        -ins_tran_test**2*theta/2.0,
                        -ins_tran_test**2*theta/2.0 * e1t,
                        -ins_tran_test**2*theta/2.0 * e2t,
                        -ins_tran_test**2*theta/2.0 * e3t,
                        -ins_tran_test**2*theta/2.0 * e3t*e1t,
                        -ins_tran_test**2*theta/2.0 * e3t*e2t,
                        -ins_tran_test**2*theta/2.0 * e3t*e3t])

simulated_cy_cov = np.array([np.nanmean(delta_obs_cons[3:]*delta_obs_inc[:-3]), 
np.nanmean(delta_obs_cons[2:]*delta_obs_inc[:-2]), 
np.nanmean(delta_obs_cons[1:]*delta_obs_inc[:-1]), 
np.nanmean(delta_obs_cons[:]*delta_obs_inc[:]),
np.nanmean(delta_obs_cons[:-1]*delta_obs_inc[1:]),
np.nanmean(delta_obs_cons[:-2]*delta_obs_inc[2:]),
np.nanmean(delta_obs_cons[:-3]*delta_obs_inc[3:])])

theoretical_cy_cov = np.flip(np.array([-ins_tran_test*theta*(1-e1o)/(1-e1t)*e1o/(theta+Omega)*(1-e1t*e1o) + ins_tran_test*theta*(1-e1o)*e2o/(theta+Omega),
                               -ins_tran_test*theta*(1-e1o)/(1-e1t)    /(theta+Omega)*(1-e1t*e1o) + ins_tran_test*theta*(1-e1o)*e1o/(theta+Omega),
                               ins_tran_test*theta/((1-e1o)*(1-e1t))* ((2-e1o)/(Omega+theta)*(1-e1o*e1t) - (1-e1t)/theta) + ins_tran_test*theta*(1-e1o)/(theta+Omega),
                               ins_tran_test*theta/((1.0-e1o)*(1.0-e1t))*((1-e1t)/theta - (1-e1t*e1o)/(theta+Omega)) - ins_tran_test*theta/(1-e1o)*( (2-e1o)/(theta+Omega)*(1-e1o*e1t) - (1-e1t)/theta ) + ins_tran_test*theta*(1-e1o)*e1t/(theta+Omega),
                               -ins_tran_test*theta/(1-e1o)    *((1-e1t)/theta - (1-e1o*e1t)/(theta+Omega)  ) - ins_tran_test*theta/(1-e1o)*e1t*( (2-e1o)/(theta+Omega)*(1-e1o*e1t) - (1-e1t)/theta ) + ins_tran_test*theta*(1-e1o)*e2t    /(theta+Omega),
                               -ins_tran_test*theta/(1-e1o)*e1t*((1-e1t)/theta - (1-e1o*e1t)/(theta+Omega)  ) - ins_tran_test*theta/(1-e1o)*e2t*( (2-e1o)/(theta+Omega)*(1-e1o*e1t) - (1-e1t)/theta ) + ins_tran_test*theta*(1-e1o)*e3t    /(theta+Omega),
                               -ins_tran_test*theta/(1-e1o)*e2t*((1-e1t)/theta - (1-e1o*e1t)/(theta+Omega)  ) - ins_tran_test*theta/(1-e1o)*e3t*( (2-e1o)/(theta+Omega)*(1-e1o*e1t) - (1-e1t)/theta ) + ins_tran_test*theta*(1-e1o)*e3t*e1t/(theta+Omega)]),0)


