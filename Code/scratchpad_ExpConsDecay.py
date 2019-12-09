"""

"""
import numpy as np
from numpy.linalg import inv
from scipy.optimize import minimize
import numdifftools as nd
import matplotlib.pyplot as plt

implied_cov = implied_cov_ExpConsDecay(init_params, 1, T)

# This assumes all parameters are constant, and consumption follows a random walk
def implied_cov_simple(params, taste, T):
    if taste ==1:
        var_taste = params[0] 
    else:
        var_taste = 0.0
    var_perm = params[taste] 
    var_tran = params[taste+1] 
    ins_perm = params[taste+2] 
    ins_tran = params[taste+3] 
    var_c_error = params[taste+4] 

    dify  =np.zeros((T,T)) #/* Income */
    difcd =np.zeros((T,T)) #/* Consumption */
    difc  =np.zeros((T,T)) #/* Consumption */
    difcme=np.zeros((T,T)) #/* Measurement error of consumption */
    difyc =np.zeros((T,T)) #/* Cov Income Consumption */
    dif   =np.zeros((2*T,2*T))
    
    ##########################################
    #/* This is the variance of Income */
    for j in np.array(range(T)):
        dify[j,j]=2.0/3.0*var_perm +  2.0*var_tran
    for j in np.array(range(T-1))+1:
        dify[j-1,j]=1.0/6.0*var_perm - var_tran
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i)+i):
            dify[j,i-1]=dify[i-1,j]
            
    ##########################################
    #/* This is the variance of Consumption */
    for j in np.array(range(T)):
        difcd[j,j]=ins_perm**2*var_perm + ins_tran**2*var_tran + var_taste
    for j in np.array(range(T)):
        difcme[j,j]=2.0*var_c_error
    for j in np.array(range(T-1)):
        difcme[j,j+1]=-var_c_error

    difc=difcme+difcd
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i))+i:
            difc[j,i-1]=difc[i-1,j]
            
    ##########################################
    #/* This is the Covariance of Income and Consumption */
    for j in np.array(range(T)):
        difyc[j,j]   = 1.0/2.0*ins_perm*var_perm + ins_tran*var_tran
    for j in np.array(range(T-1))+1:
        difyc[j-1,j] = 1.0/2.0*ins_perm*var_perm - ins_tran*var_tran
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

# This assumes all parameters are constant, and consumption follows a random walk
def implied_cov_ExpConsDecay(params, taste, T):
    if taste ==1:
        var_taste = params[0] 
    else:
        var_taste = 0.0
    var_perm = params[taste] 
    var_tran = params[taste+1] 
    ins_perm = params[taste+2] 
    ins_tran = params[taste+3] 
    var_c_error = params[taste+4] 
    theta = params[taste+5] 

    dify  =np.zeros((T,T)) #/* Income */
    difcd =np.zeros((T,T)) #/* Consumption */
    difc  =np.zeros((T,T)) #/* Consumption */
    difcme=np.zeros((T,T)) #/* Measurement error of consumption */
    difyc =np.zeros((T,T)) #/* Cov Income Consumption */
    dif   =np.zeros((2*T,2*T))
    
    ##########################################
    #/* This is the variance of Income */
    for j in np.array(range(T)):
        dify[j,j]=2.0/3.0*var_perm +  2.0*var_tran
    for j in np.array(range(T-1))+1:
        dify[j-1,j]=1.0/6.0*var_perm - var_tran
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i)+i):
            dify[j,i-1]=dify[i-1,j]
            
    ##########################################
    #/* This is the variance of Consumption */
    zeta2 = theta*ins_tran**2/(2*(1.0-np.exp(-theta))**2)   #useful quantity that comes up a lot
    for j in np.array(range(T)):
        difcd[j,j]=ins_perm**2*var_perm + zeta2*((1-np.exp(-2*theta))  +  np.exp(-2.0*theta)*(1.0-np.exp(theta))**2  )*var_tran + var_taste
    for i in np.array(range(T-1)):
        for j in np.array(range(T-1-i)):
            difcd[i, i+j+1] = zeta2* (np.exp(-(j+1)*theta)*(1-np.exp(-2*theta))*(1-np.exp(theta))  + np.exp(-theta*(j+3))*(1-np.exp(theta))**2   ) *var_tran
            difcd[i+j+1,i] = difcd[i, i+j+1]
    
    for j in np.array(range(T)):
        difcme[j,j]=2.0*var_c_error
    for j in np.array(range(T-1)):
        difcme[j,j+1]=-var_c_error

    difc=difcme+difcd
    
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i))+i:
            difc[j,i-1]=difc[i-1,j]
            
    ##########################################
    #/* This is the Covariance of Income and Consumption */
    for j in np.array(range(T)):
        difyc[j,j]   = 1.0/2.0*ins_perm*var_perm + ins_tran*(2.0-np.exp(-theta) )*var_tran
    for j in np.array(range(T-1))+1:
        difyc[j-1,j] = 1.0/2.0*ins_perm*var_perm - ins_tran*var_tran
    for i in np.array(range(T-1)):
        for j in np.array(range(T-1-i)):
            difyc[i+j+1,i] = ins_tran*( np.exp(-(j+1)*theta) *(1.0 -np.exp(theta) )*(1.0 -np.exp(-theta) )  )*var_tran
    ##########################################
            
    #/* Final matrix */
    dif[0:T,0:T]            =difc
    dif[T:2*(T),0:T]        =difyc
    dif[0:T,T:2*(T)]        =np.transpose(difyc)
    dif[T:2*(T),T:2*(T)]    =dify
    
    # remove missing consumption values
    if True:
        difa1 = np.concatenate((dif[0:8,:],dif[11:2*T,:]),0)
        difa2 = np.concatenate((difa1[:,0:8],difa1[:,11:2*T]),1)
    
    vech_indicies = np.tril_indices(np.shape(difa2)[0])
    fm=difa2[vech_indicies]

    return fm


def Parameter_estimation_simple(model, c_vector, omega, T, taste=1, theta=0):
    if model=='Simple':
        implied_cov = lambda params, taste, T : implied_cov_simple(params, taste, T)
        num_params = taste+5
    if model=='ExpConsDecay':
        implied_cov = lambda params, taste, T : implied_cov_ExpConsDecay(params, taste, T)
        num_params = taste+6
    if model=='ExpConsDecay_fix_theta':
        implied_cov = lambda params, taste, T, theta : implied_cov_ExpConsDecay(np.concatenate((params,[theta])), taste, T)
        num_params = taste+5
    
    init_params = np.zeros(num_params)
    if taste:
        init_params[0] = 0.01  #variance of taste shocks
    init_params[taste] = 0.03
    init_params[taste+1] = 0.03
    init_params[taste+2] = 0.5
    init_params[taste+3] = 0.25
    init_params[taste+4] = 0.5
    if model=='ExpConsDecay':
        init_params[taste+5] = 0.5
        
    def objectiveFun(params, taste, T, theta, empirical_cov, weight_matrix):
        if model=='ExpConsDecay_fix_theta':
            model_cov = implied_cov(params, taste,T, theta)
        else:
            model_cov = implied_cov(params, taste,T)
        distance = np.dot(np.dot((model_cov-empirical_cov), weight_matrix),(model_cov-empirical_cov))
        distance = distance #+ 10000*np.std(params[ma+taste:ma+taste+perm_shk_params]) #add in this to keep same variance for permanent shocks over the whole time period
        return distance
    
    # Define the weight matrix as Equal Weight Minimum Distance
    weight_matrix = np.diag(np.diag(omega)**(-1))
    
    ret = objectiveFun(init_params, taste, T, theta, c_vector, weight_matrix)
    
    #Solve with one method, reset init_params, then solve again. Seems to converge OK.
    solved_objective1 = minimize(objectiveFun, init_params, args=(taste, T, theta, c_vector, weight_matrix))  
    init_params2 = solved_objective1.x
    solved_objective = minimize(objectiveFun, init_params2, args=(taste, T, theta, c_vector, weight_matrix),method='Nelder-Mead')
     
    solved_params = solved_objective.x
    
    if model=='ExpConsDecay_fix_theta':
        fun_for_jacob = lambda params: implied_cov(params, taste, T,theta)
    else:
        fun_for_jacob = lambda params: implied_cov(params, taste, T)
    jacob = nd.Jacobian(fun_for_jacob)(solved_params)
    
    Sandwich1 = inv(np.dot(np.transpose(jacob),np.dot(weight_matrix,jacob)))
    Sandwich2 = np.dot(np.transpose(jacob),np.dot(weight_matrix,np.dot(omega,np.dot(weight_matrix,jacob))))
    cov_params = np.dot(Sandwich1,np.dot(Sandwich2,Sandwich1))
    standard_errors = np.diag(cov_params)**0.5
    
    #extract relevant numbers
    if taste:
        varcsi = solved_params[0] 
        varcsi_se = standard_errors[0] 
    else:
        varcsi = 0.0
        varcsi_se = 0.0 
    var_perm = solved_params[taste] 
    var_tran = solved_params[taste+1] 
    ins_perm = solved_params[taste+2] 
    ins_tran = solved_params[taste+3] 
    var_c_error = solved_params[taste+4] 
    
    var_perm_se = standard_errors[taste] 
    var_tran_se = standard_errors[taste+1] 
    ins_perm_se = standard_errors[taste+2] 
    ins_tran_se = standard_errors[taste+3] 
    var_c_error_se = standard_errors[taste+4] 
    if model=='ExpConsDecay':
        theta = solved_params[taste+5] 
        theta_se = standard_errors[taste+5] 
    else:
        theta = 0.0
        theta_se = 0.0
    return var_perm, var_perm_se, var_tran, var_tran_se, ins_perm, ins_perm_se, ins_tran, ins_tran_se, var_c_error, var_c_error_se, theta, theta_se, varcsi, varcsi_se
    
#implied_cov_ExpConsDecay(init_params, taste, T)
#implied_cov_simple(init_params, taste, T)

var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
 ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
 var_c_error_se_BPP, theta_BPP, theta_se_BPP, varcsi_BPP, varcsi_se_BPP \
  = Parameter_estimation_simple('ExpConsDecay', c_vector, omega, T, taste=1) 
#        
# Create some simulated data to test 
var_perm_test = 0.05
var_tran_test = 0.04
ins_perm_test = 0.5
ins_tran_test = 0.3
theta_test = 0.5

#var_perm_test = 0.0
#var_tran_test = 1.0
#ins_perm_test = 0.0
#ins_tran_test = 1.0
#theta_test = 0.2

#var_perm_test = 0.5
#var_tran_test = 1.0
#ins_perm_test = 1.0
#ins_tran_test = 0.5
#theta_test = 0.000001


num_periods = 4000
num_divisions = 20
total_obs = num_periods*num_divisions
np.random.seed(seed=6)
perm_shocks = (var_perm_test/num_divisions)**0.5*np.random.normal(size=total_obs)
tran_shocks = (var_tran_test*num_divisions)**0.5*np.random.normal(size=total_obs)
perm_inc  = np.cumsum(perm_shocks)
perm_cons = ins_perm_test*perm_inc
tran_cons = np.zeros_like(perm_cons)
decay = np.exp(-np.arange(total_obs)*theta_test/num_divisions)
for t in range(total_obs):
    tran_cons[t] = np.sum(theta_test*ins_tran_test/(1.0-np.exp(-theta_test)) * decay[0:t+1] * np.flip(tran_shocks[0:t+1],0))/num_divisions

#tran_cons = ins_tran_test*np.cumsum(tran_shocks)/num_divisions  # This is a random walk

total_inc  = perm_inc  + tran_shocks
total_cons = perm_cons + tran_cons
observed_inc = np.zeros(num_periods)
observed_cons = np.zeros(num_periods)
for T1 in range(num_periods):
    observed_inc[T1] = np.sum(total_inc[num_divisions*T1:num_divisions*(T1+1)])/num_divisions
    observed_cons[T1] = total_cons[(T1+1)*num_divisions-1]
#for T in range(num_periods):    
#    for i in range(np.min([(T+1)*num_divisions,50])):
#        observed_cons[T] += theta_test*ins_tran_test/(1.0-np.exp(-theta_test)) * np.exp(-i*theta_test/num_divisions) * tran_shocks[(T+1)*num_divisions-i-1]
delta_inc = np.diff(observed_inc)
delta_cons = np.diff(observed_cons)

T=3
#Now create a c_vector matrix
ignore = 20
delta_offset = np.zeros((num_periods-T-ignore,2*T))
for t in np.array(range(T)):
    delta_offset[:,t] = delta_cons[ignore+T-t-2:-1-t]
    delta_offset[:,t+T] = delta_inc[ignore+T-t-2:-1-t]
cov_matrix_test = np.cov(np.transpose(delta_offset))    
    
difa1 = np.concatenate((cov_matrix_test[0:8,:],cov_matrix_test[11:2*T,:]),0)
difa2 = np.concatenate((difa1[:,0:8],difa1[:,11:2*T]),1)

vech_indicies = np.tril_indices(np.shape(difa2)[0])
cov_vector_test=difa2[vech_indicies]
omega_test = np.eye(cov_vector_test.shape[0])
    
var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
 ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
 var_c_error_se_BPP, theta_BPP, theta_se_BPP, varcsi_BPP, varcsi_se_BPP \
  = Parameter_estimation_simple('ExpConsDecay', cov_vector_test, omega_test, T, taste=1) 
  
#var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
# ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
# var_c_error_se_BPP, theta_BPP, theta_se_BPP, varcsi_BPP, varcsi_se_BPP \
#  = Parameter_estimation_simple('Simple', cov_vector_test, omega_test, T, taste=1) 
  
init_params = np.array([0.0,var_perm_test,var_tran_test,ins_perm_test,ins_tran_test,0.0,theta_test])
implied_cov = implied_cov_ExpConsDecay(init_params, 1, T)

T=14
num_thetas=10
var_perm_fix_theta = np.zeros(num_thetas+1)
var_tran_fix_theta = np.zeros(num_thetas+1)
ins_perm_fix_theta = np.zeros(num_thetas+1)
ins_tran_fix_theta = np.zeros(num_thetas+1)
theta_array =  np.concatenate(([0],np.linspace(0.01,5,num_thetas)))
var_perm_fix_theta[0], var_perm_se_BPP, var_tran_fix_theta[0], var_tran_se_BPP, ins_perm_fix_theta[0], \
 ins_perm_se_BPP, ins_tran_fix_theta[0], ins_tran_se_BPP, var_c_error_BPP, \
 var_c_error_se_BPP, theta_BPP, theta_se_BPP, varcsi_BPP, varcsi_se_BPP \
  = Parameter_estimation_simple('Simple', c_vector, omega, T, taste=1)
for i in np.array(range(num_thetas))+1:
    theta = theta_array[i]
    var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
     ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
     var_c_error_se_BPP, theta_BPP, theta_se_BPP, varcsi_BPP, varcsi_se_BPP \
      = Parameter_estimation_simple('ExpConsDecay_fix_theta', c_vector, omega, T, theta=theta, taste=1)
    var_perm_fix_theta[i] = var_perm_BPP
    var_tran_fix_theta[i] = var_tran_BPP
    ins_perm_fix_theta[i] = ins_perm_BPP
    ins_tran_fix_theta[i] = ins_tran_BPP
    
    
plt.plot(theta_array,ins_perm_fix_theta)
plt.plot(theta_array,ins_tran_fix_theta)