#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:33:11 2020

@author: m1esc01
"""
c_vector, omega, T = create_moment_vector(Path("./InputFiles/CohA.csv"))


def Parameter_estimation(model, c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="None"):
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
    init_params[ma+taste:ma+taste+perm_shk_params] = 0.0275*np.ones(perm_shk_params)
    init_params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params] = 0.03*np.ones(tran_shk_params)
    init_params[ma+taste+perm_shk_params+tran_shk_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params] = 1.0*np.ones(perm_ins_params)
    init_params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params] = 0.3*np.ones(tran_ins_params)
    init_params[ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params:ma+taste+perm_shk_params+tran_shk_params+perm_ins_params+tran_ins_params+meas_error_params] = 0.06*np.ones(meas_error_params)
    
    def objectiveFun(params, ma, taste, varying_ins, T, empirical_cov, weight_matrix,fix):
        model_cov = implied_cov(params, ma, taste, varying_ins,T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
        distance = np.dot(np.dot((model_cov-empirical_cov), weight_matrix),(model_cov-empirical_cov))
        if fix == "perm_var":
            perm_var=params[ma+taste+perm_shk_params]
            params_new=np.ones_like(params)*params
            params_new[ma+taste:ma+taste+perm_shk_params]=perm_var
            model_cov = implied_cov(params_new, ma, taste, varying_ins,T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
            distance = np.dot(np.dot((model_cov-empirical_cov), weight_matrix),(model_cov-empirical_cov))
            distance = distance + 10000*np.sum((params[ma+taste+1:ma+taste+perm_shk_params])**2) #add in this to keep same variance for permanent shocks over the whole time period
        elif fix == "tran_var":
            distance = distance + 10000*(np.std(params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params]))
        elif fix == "all_var":
            distance = distance + 10000*(np.std(params[ma+taste:ma+taste+perm_shk_params])+10*np.std(params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params]))
        if fix == "mid_perm_var":
            params_new=np.ones_like(params)*params
            fix_perm_var=params[ma+taste+perm_shk_params+1]
            params_new[ma+taste+1:ma+taste+perm_shk_params-1]=fix_perm_var
            model_cov = implied_cov(params_new, ma, taste, varying_ins,T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
            distance = np.dot(np.dot((model_cov-empirical_cov), weight_matrix),(model_cov-empirical_cov))
            distance = distance + 10000*np.sum((params[ma+taste+2:ma+taste+perm_shk_params-1])**2)
        elif fix == "fix0275":
            distance = distance + 10000*(100*np.sum((params[ma+taste:ma+taste+perm_shk_params]-0.0275)**2)+10*np.std(params[ma+taste+perm_shk_params:ma+taste+perm_shk_params+tran_shk_params]))
        return distance
    
    def orig_objectiveFun(params, ma, taste, varying_ins, T, empirical_cov, weight_matrix,fix):
        model_cov = implied_cov(params, ma, taste, varying_ins,T, perm_shk_params, tran_shk_params, perm_ins_params,tran_ins_params, meas_error_params)
        distance = np.dot(np.dot((model_cov-empirical_cov), weight_matrix),(model_cov-empirical_cov))
        return distance
    
    # Define the weight matrix as Equal Weight Minimum Distance
    weight_matrix = np.diag(np.diag(omega)**(-1))
    
    ret = objectiveFun(init_params, ma, taste, varying_ins, T, c_vector, weight_matrix,fix)
    
    #Solve with one method, reset init_params, then solve again. Seems to converge OK.
    solved_objective1 = minimize(objectiveFun, init_params, args=(ma, taste, varying_ins, T, c_vector, weight_matrix,fix))  
    init_params2 = solved_objective1.x
    solved_objective = minimize(objectiveFun, init_params2, args=(ma, taste, varying_ins, T, c_vector, weight_matrix,fix),method='Nelder-Mead')
     
    solved_params = solved_objective.x
    print(orig_objectiveFun(solved_params, ma, taste, varying_ins, T, c_vector, weight_matrix,fix))
    
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


print('Replicate BPP')
var_perm_fix_midpermvar, var_perm_se_fix_midpermvar, var_tran_fix_midpermvar, var_tran_se_fix_midpermvar, ins_perm_fix_midpermvar, \
 ins_perm_se_fix_midpermvar, ins_tran_fix_midpermvar, ins_tran_se_fix_midpermvar, var_c_error_fix_midpermvar, \
 var_c_error_se_fix_midpermvar, teta_fix_midpermvar, teta_se_fix_midpermvar, varcsi_fix_midpermvar, varcsi_se_fix_midpermvar \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="mid_perm_var") 


print('Replicate BPP')
var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
 ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
 var_c_error_se_BPP, teta_BPP, teta_se_BPP, varcsi_BPP, varcsi_se_BPP \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="None") 
  
print('Replicate BPP')
var_perm_fix_permvar, var_perm_se_fix_permvar, var_tran_fix_permvar, var_tran_se_fix_permvar, ins_perm_fix_permvar, \
 ins_perm_se_fix_permvar, ins_tran_fix_permvar, ins_tran_se_fix_permvar, var_c_error_fix_permvar, \
 var_c_error_se_fix_permvar, teta_fix_permvar, teta_se_fix_permvar, varcsi_fix_permvar, varcsi_se_fix_permvar \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="perm_var") 
  
print('Replicate BPP')
var_perm_fix_tranvar, var_perm_se_fix_tranvar, var_tran_fix_tranvar, var_tran_se_fix_tranvar, ins_perm_fix_tranvar, \
 ins_perm_se_fix_tranvar, ins_tran_fix_tranvar, ins_tran_se_fix_tranvar, var_c_error_fix_tranvar, \
 var_c_error_se_fix_tranvar, teta_fix_tranvar, teta_se_fix_tranvar, varcsi_fix_tranvar, varcsi_se_fix_tranvar \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="tran_var") 
  
print('Replicate BPP')
var_perm_fix_allvar, var_perm_se_fix_allvar, var_tran_fix_allvar, var_tran_se_fix_allvar, ins_perm_fix_allvar, \
 ins_perm_se_fix_allvar, ins_tran_fix_allvar, ins_tran_se_fix_allvar, var_c_error_fix_allvar, \
 var_c_error_se_fix_allvar, teta_fix_allvar, teta_se_fix_allvar, varcsi_fix_allvar, varcsi_se_fix_allvar \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="all_var") 

print('Replicate BPP')
var_perm_fix0275, var_perm_se_fix0275, var_tran_fix0275, var_tran_se_fix0275, ins_perm_fix0275, \
 ins_perm_se_fix0275, ins_tran_fix0275, ins_tran_se_fix0275, var_c_error_fix0275, \
 var_c_error_se_fix0275, teta_fix0275, teta_se_fix0275, varcsi_fix0275, varcsi_se_fix0275 \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="fix0275") 


plt.plot(var_perm_fix_allvar)
plt.plot(var_perm_fix_permvar)
plt.plot(var_perm_fix_midpermvar)
plt.plot(var_perm_fix0275)

plt.plot(var_tran_fix_allvar)
plt.plot(var_tran_fix_midpermvar)

plt.plot(var_c_error_fix_allvar)
plt.plot(var_c_error_fix_midpermvar)
  
print('Replicate BPP')
var_perm_fix_midpermvar, var_perm_se_fix_midpermvar, var_tran_fix_midpermvar, var_tran_se_fix_midpermvar, ins_perm_fix_midpermvar, \
 ins_perm_se_fix_midpermvar, ins_tran_fix_midpermvar, ins_tran_se_fix_midpermvar, var_c_error_fix_midpermvar, \
 var_c_error_se_fix_midpermvar, teta_fix_midpermvar, teta_se_fix_midpermvar, varcsi_fix_midpermvar, varcsi_se_fix_midpermvar \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0,fix="mid_perm_var") 
 