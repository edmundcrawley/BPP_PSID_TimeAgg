# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 10:43:33 2019

@author: edmun
"""
num_periods = 20000
num_divisions = 200
total_obs = num_periods*num_divisions
max_decay = 10000
#np.random.seed(seed=3)

theta_test = 0.9
theta_inc_test = 1.0
ins_tran_test = 0.3

#tran_shocks = (var_tran_test*num_divisions)**0.5*np.random.normal(size=total_obs)
tran_shocks = np.random.normal(size=total_obs)/num_divisions**0.5

tran_cons = np.zeros_like(tran_shocks[1:])
tran_inc = np.zeros_like(tran_cons)



decay = np.exp(-np.arange(total_obs)*theta_test/num_divisions)
decay_inc = np.exp(-np.arange(total_obs)*theta_inc_test/num_divisions)

tran_inc_0_1 = np.zeros(num_periods)*np.nan
tran_inc_1_2 = np.zeros(num_periods)*np.nan
tran_inc_1_2_a = np.zeros(num_periods)*np.nan
tran_inc_1_2_b = np.zeros(num_periods)*np.nan
tran_inc_2_3 = np.zeros(num_periods)*np.nan
tran_inc_2_infty = np.zeros(num_periods)*np.nan
tran_inc_3_4 = np.zeros(num_periods)*np.nan
tran_inc_3_infty = np.zeros(num_periods)*np.nan
tran_inc_4_infty = np.zeros(num_periods)*np.nan
observed_inc = np.zeros(num_periods)*np.nan
for T1 in range(num_periods):
    tran_inc_0_1[T1]     = 1.0/(1-np.exp(-theta_inc_test))*np.sum((1-np.flip(decay_inc[0:num_divisions] ,0)) * tran_shocks[T1*num_divisions:(T1+1)*num_divisions])
    if T1>=1:
        tran_inc_1_2[T1] = 1.0/(1-np.exp(-theta_inc_test))*np.sum( ((2-np.exp(theta_inc_test))*np.flip(decay_inc[0:num_divisions] ,0)  -1)  * tran_shocks[(T1-1)*num_divisions:(T1)*num_divisions])
        tran_inc_1_2_a[T1] = 1.0/(1-np.exp(-theta_inc_test))*np.sum( ((2-np.exp(theta_inc_test))*np.flip(decay_inc[0:num_divisions] ,0)  )  * tran_shocks[(T1-1)*num_divisions:(T1)*num_divisions])
        tran_inc_1_2_b[T1] = 1.0/(1-np.exp(-theta_inc_test))*np.sum(   -1  * tran_shocks[(T1-1)*num_divisions:(T1)*num_divisions])
    if T1>=2:
        tran_inc_2_infty[T1] = -(1-np.exp(-theta_inc_test))*np.sum( np.flip(decay_inc[0:min(max_decay,(T1-1)*num_divisions)] ,0)   * tran_shocks[max(0,(T1-1)*num_divisions-max_decay):(T1-1)*num_divisions])
        tran_inc_2_3[T1] = -(1-np.exp(-theta_inc_test))*np.sum( np.flip(decay_inc[0:num_divisions]                      ,0)   * tran_shocks[(T1-2)*num_divisions:(T1-1)*num_divisions])
        observed_inc[T1] = tran_inc_0_1[T1] + np.sum( np.flip(decay_inc[0:min(max_decay,T1*num_divisions)] ,0)   * tran_shocks[max(0,T1*num_divisions-max_decay):(T1)*num_divisions])
    if T1>=3:
        tran_inc_3_infty[T1] = -(1-np.exp(-theta_inc_test))*np.sum( np.flip(decay_inc[0:min(max_decay,(T1-2)*num_divisions)] ,0)   * tran_shocks[max(0,(T1-2)*num_divisions-max_decay):(T1-2)*num_divisions])*np.exp(-theta_inc_test)
        tran_inc_3_4[T1] = -(1-np.exp(-theta_inc_test))*np.sum( np.flip(decay_inc[0:num_divisions]                      ,0)   * tran_shocks[(T1-3)*num_divisions:(T1-2)*num_divisions])*np.exp(-theta_inc_test)
    if T1>=4:
        tran_inc_4_infty[T1] = -(1-np.exp(-theta_inc_test))*np.sum( np.flip(decay_inc[0:min(max_decay,(T1-3)*num_divisions)] ,0)   * tran_shocks[max(0,(T1-3)*num_divisions-max_decay):(T1-3)*num_divisions])*np.exp(-2*theta_inc_test)


var_tran_inc_0_1 = np.nanmean(tran_inc_0_1**2)
var_tran_inc_1_2 = np.nanmean(tran_inc_1_2**2)
var_tran_inc_2_infty = np.nanmean(tran_inc_2_infty**2)
foo = [var_tran_inc_0_1,var_tran_inc_1_2,var_tran_inc_2_infty]
delta_obs_inc = observed_inc[1:]-observed_inc[:-1]

[np.nanmean(delta_obs_inc**2), 
np.nanmean(delta_obs_inc[1:]*delta_obs_inc[:-1]), 
np.nanmean(delta_obs_inc[2:]*delta_obs_inc[:-2]), 
np.nanmean(delta_obs_inc[3:]*delta_obs_inc[:-3]), 
np.nanmean(delta_obs_inc[4:]*delta_obs_inc[:-4]),
np.nanmean(delta_obs_inc[5:]*delta_obs_inc[:-5]),
np.nanmean(delta_obs_inc[6:]*delta_obs_inc[:-6]),
np.nanmean(delta_obs_inc[7:]*delta_obs_inc[:-7])]

foo2=[np.nanmean(tran_inc_0_1[0:-1]*tran_inc_1_2[1:]),
np.nanmean(tran_inc_1_2[0:-1]*tran_inc_2_3[1:]),
np.nanmean(tran_inc_2_3[0:-1]*tran_inc_3_infty[1:])]
foo2

np.nanmean(tran_inc_0_1[0:-1]*tran_inc_1_2[1:])+np.nanmean(tran_inc_1_2[0:-1]*tran_inc_2_3[1:])+np.nanmean(tran_inc_2_3[0:-1]*tran_inc_3_infty[1:])


foo3=[np.nanmean(tran_inc_0_1[0:-2]*tran_inc_2_3[2:]),
np.nanmean(tran_inc_1_2_a[0:-2]*tran_inc_3_4[2:]),
np.nanmean(tran_inc_1_2_b[0:-2]*tran_inc_3_4[2:]),
np.nanmean(tran_inc_2_3[0:-2]*tran_inc_4_infty[2:])]
foo3
np.sum(foo3)

delta_obs_inc
tran_inc_0_1+tran_inc_1_2+tran_inc_2_infty
