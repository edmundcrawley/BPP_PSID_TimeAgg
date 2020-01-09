"""

"""
import numpy as np
from numpy.linalg import inv
from scipy.optimize import minimize
import numdifftools as nd
import matplotlib.pyplot as plt
from pathlib import Path
from create_moments import create_moment_vector


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

# This assumes all parameters are constant, and consumption decays exponentially
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

# This assumes all parameters are constant, and consumption AND income decay exponentially
def implied_cov_ExpConsIncDecay(params, taste, T):
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
    Omega = params[taste+6] 

    dify  =np.zeros((T,T)) #/* Income */
    difcd =np.zeros((T,T)) #/* Consumption */
    difc  =np.zeros((T,T)) #/* Consumption */
    difcme=np.zeros((T,T)) #/* Measurement error of consumption */
    difyc =np.zeros((T,T)) #/* Cov Income Consumption */
    dif   =np.zeros((2*T,2*T))
    
    # Useful quantities to store
    e1t = np.exp(-theta)
    e1o = np.exp(-Omega)
    
    ##########################################
    #/* This is the variance of Income */
    for j in np.array(range(T)):
        dify[j,j]=2.0/3.0*var_perm +  (2.0/(1-e1o)**2 - (3-e1o)/(Omega*(1-e1o)))*var_tran
    for j in np.array(range(T-1))+1:
        dify[j-1,j]=1.0/6.0*var_perm + ( (2-e1o)/(2*Omega) - (1.0-(1.0-e1o)/Omega )/(1-e1o)**2   )*var_tran
    for M in np.array(range(T-2))+2:
        for j in np.array(range(T-M))+M:
            dify[j-M,j] = ( -(1-e1o)*e1o**(M-2)/Omega  +(1-e1o**2)*e1o**(M-2)/(2*Omega) - (2-Omega)*e1o**(M-1) *(1-e1o**2)/(2*Omega) + e1o**(M-1)*(1-e1o)/Omega + (1-e1o)**2 *e1o**M /(2*Omega)   ) * var_tran
    #symetric
    for i in np.array(range(T-1))+1:
        for j in np.array(range(T-i)+i):
            dify[j,i-1]=dify[i-1,j]
            
    ##########################################
    #/* This is the variance of Consumption */
    for j in np.array(range(T)):
        difcd[j,j]=ins_perm**2*var_perm + ins_tran**2*theta/(1-e1t)*var_tran + var_taste
    for M in np.array(range(T-1))+1:
        for j in np.array(range(T-M))+M:
            difcd[j-M,j] = (-ins_tran**2*theta/2.0*e1t**(M-1) )*var_tran
            difcd[j,j-M] = difcd[j-M,j]
    
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
        difyc[j,j]   = 1.0/2.0*ins_perm*var_perm + ins_tran*theta*( ((1-e1t)/theta - (1-e1t*e1o)/(theta+Omega))/((1-e1o)*(1-e1t))  - ((2-e1o)*(1-e1o*e1t)/(theta+Omega) - (1-e1t)/theta )/(1-e1o) + (1-e1o)*e1t/(theta+Omega) )*var_tran
    for j in np.array(range(T-1))+1:
        difyc[j-1,j] = 1.0/2.0*ins_perm*var_perm + ins_tran*theta*( ( (2-e1o)*(1-e1o*e1t)/(theta+Omega)  -(1-e1t)/theta )/((1-e1o)*(1-e1t)) + (1-e1o)/(theta+Omega) )*var_tran
    for M in np.array(range(T-2))+2:
        for j in np.array(range(T-M))+M:
            difyc[j-M,j] = ins_tran*theta*( -(1-e1o)/(1-e1t)*e1o**(M-2)*(1-e1o*e1t)/(Omega+theta) + (1-e1o)*e1o**(M-1)/(theta+Omega)  )*var_tran
    for M in np.array(range(T-1))+1:
        for j in np.array(range(T-M))+M:
            difyc[j,j-M] = ins_tran*theta*( -e1t**(M-1)*( (1-e1t)/theta - (1-e1o*e1t)/(theta+Omega)  )/(1-e1o) - e1t**M*( (2-e1o)*(1-e1o*e1t)/(theta+Omega) -(1-e1t)/theta )/(1-e1o)  + (1-e1o)*e1t**(M+1)/(theta+Omega) )*var_tran
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
    if model=='ExpConsIncDecay':
        implied_cov = lambda params, taste, T : implied_cov_ExpConsIncDecay(params, taste, T)
        num_params = taste+7
    
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
    if model=='ExpConsIncDecay':
        init_params[taste+5] = 0.5
        init_params[taste+6] = 2.0
        
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
        Omega = 0.0 
        Omega_se = 0.0
    elif model=='ExpConsIncDecay':
        theta = solved_params[taste+5] 
        theta_se = standard_errors[taste+5] 
        Omega = solved_params[taste+6] 
        Omega_se = standard_errors[taste+6] 
    else:
        theta = 0.0
        theta_se = 0.0
        Omega = 0.0 
        Omega_se = 0.0
    return var_perm, var_perm_se, var_tran, var_tran_se, ins_perm, ins_perm_se, ins_tran, ins_tran_se, var_c_error, var_c_error_se, theta, theta_se, varcsi, varcsi_se, Omega, Omega_se
    
###########################################################################
# Start estimations
###########################################################################
# Calculate empirical moments
print('Create Moment Vector from Inputs for All')
c_vector, omega, T = create_moment_vector(Path("./InputFiles/CohA.csv"))

var_perm_ExpsDecay, var_perm_se_ExpsDecay, var_tran_ExpsDecay, var_tran_se_ExpsDecay, ins_perm_ExpsDecay, \
 ins_perm_se_ExpsDecay, ins_tran_ExpsDecay, ins_tran_se_ExpsDecay, var_c_error_ExpsDecay, \
 var_c_error_se_ExpsDecay, theta_ExpsDecay, theta_se_ExpsDecay, varcsi_ExpsDecay, varcsi_se_ExpsDecay, Omega_ExpsDecay, Omega_se_ExpsDecay \
  = Parameter_estimation_simple('ExpConsDecay', c_vector, omega, T, taste=1) 
  
var_perm_ExpsIncDecay, var_perm_se_ExpsIncDecay, var_tran_ExpsIncDecay, var_tran_se_ExpsIncDecay, ins_perm_ExpsIncDecay, \
 ins_perm_se_ExpsIncDecay, ins_tran_ExpsIncDecay, ins_tran_se_ExpsIncDecay, var_c_error_ExpsIncDecay, \
 var_c_error_se_ExpsIncDecay, theta_ExpsIncDecay, theta_se_ExpsIncDecay, varcsi_ExpsIncDecay, varcsi_se_ExpsIncDecay, Omega_ExpsIncDecay, Omega_se_ExpsIncDecay \
  = Parameter_estimation_simple('ExpConsIncDecay', c_vector, omega, T, taste=1) 
#        

T=14
num_thetas=10
var_perm_fix_theta = np.zeros(num_thetas+1)
var_tran_fix_theta = np.zeros(num_thetas+1)
ins_perm_fix_theta = np.zeros(num_thetas+1)
ins_tran_fix_theta = np.zeros(num_thetas+1)
var_perm_se_fix_theta = np.zeros(num_thetas+1)
var_tran_se_fix_theta = np.zeros(num_thetas+1)
ins_perm_se_fix_theta = np.zeros(num_thetas+1)
ins_tran_se_fix_theta = np.zeros(num_thetas+1)
theta_array =  np.concatenate(([0],np.linspace(0.01,5,num_thetas)))
var_perm_fix_theta[0], var_perm_se_fix_theta[0], var_tran_fix_theta[0], var_tran_se_fix_theta[0], ins_perm_fix_theta[0], \
 ins_perm_se_fix_theta[0], ins_tran_fix_theta[0], ins_tran_se_fix_theta[0], var_c_error_fix_theta, \
 var_c_error_se_fix_theta, theta_fix_theta, theta_se_fix_theta, varcsi_fix_theta, varcsi_se_fix_theta \
  = Parameter_estimation_simple('Simple', c_vector, omega, T, taste=1)
for i in np.array(range(num_thetas))+1:
    theta = theta_array[i]
    var_perm_fix_theta[i], var_perm_se_fix_theta[i], var_tran_fix_theta[i], var_tran_se_fix_theta[i], ins_perm_fix_theta[i], \
     ins_perm_se_fix_theta[i], ins_tran_fix_theta[i], ins_tran_se_fix_theta[i], var_c_error_fix_theta, \
     var_c_error_se_fix_theta, theta_fix_theta, theta_se_fix_theta, varcsi_fix_theta, varcsi_se_fix_theta \
      = Parameter_estimation_simple('ExpConsDecay_fix_theta', c_vector, omega, T, theta=theta, taste=1)
#    var_perm_fix_theta[i] = var_perm_fix_theta
#    var_tran_fix_theta[i] = var_tran_fix_theta
#    ins_perm_fix_theta[i] = ins_perm_fix_theta
#    ins_tran_fix_theta[i] = ins_tran_fix_theta
#    var_perm_se_fix_theta[i] = var_perm_fix_theta
#    var_tran_se_fix_theta[i] = var_tran_fix_theta
#    ins_perm_se_fix_theta[i] = ins_perm_fix_theta
#    ins_tran_se_fix_theta[i] = ins_tran_fix_theta
    
    
plt.plot(theta_array,ins_perm_fix_theta,color="blue")
plt.plot(theta_array,ins_tran_fix_theta,color="red")
plt.plot(theta_array,ins_perm_fix_theta+1.96*ins_perm_se_fix_theta,color="blue",linestyle="--")
plt.plot(theta_array,ins_perm_fix_theta-1.96*ins_perm_se_fix_theta,color="blue",linestyle="--")
plt.plot(theta_array,ins_tran_fix_theta+1.96*ins_tran_se_fix_theta,color="red",linestyle="--")
plt.plot(theta_array,ins_tran_fix_theta-1.96*ins_tran_se_fix_theta,color="red",linestyle="--")


###############################################################################
# Look just at the perm insurance moments
###############################################################################
#recreate covariance matrix from vector
c_matrix1 = np.zeros((25,25))
vech_indicies = np.tril_indices(25)
c_matrix1[vech_indicies] = cov_vector_test
c_matrix1 = np.transpose(c_matrix1)
c_matrix1[vech_indicies] = cov_vector_test
#add in missing consumption as N/As
c_matrix2 = np.concatenate((c_matrix1[0:8,:],np.nan*np.zeros((3,25)),c_matrix1[8:25,:]),0)
c_matrix = np.concatenate((c_matrix2[:,0:8],np.nan*np.zeros((28,3)),c_matrix2[:,8:25]),1)

# Now look at equation 9
T=14
phi_estimate = np.zeros(12)
nominator     = np.zeros(12)
denominator   = np.zeros(12)
for t in np.array(range(12))+1:
    nominator[t-1]   = c_matrix[t,T+t-1] + c_matrix[t,T+t] + c_matrix[t,T+t+1]
    denominator[t-1] = c_matrix[T+t,T+t-1] + c_matrix[T+t,T+t] + c_matrix[T+t,T+t+1]
phi_estimate = nominator/denominator
mean_phi_estimate = np.nanmean(nominator)/np.nanmean(denominator+0.0*nominator)
    



##############################################################################
# TESTING: Create some simulated data to test the moment calculation 
# returns unbiased estimates
##############################################################################
#var_perm_test = 0.02
#var_tran_test = 0.04
#ins_perm_test = 0.20685882126750385
#ins_tran_test = 0.1759568432842964
#theta_test = 0.9036180417684739
#Omega_test = 4.0

var_perm_test = var_perm_ExpsIncDecay
var_tran_test = var_tran_ExpsIncDecay
ins_perm_test = ins_perm_ExpsDecay
ins_tran_test = ins_tran_ExpsDecay
theta_test = theta_ExpsDecay
Omega_test = Omega_ExpsIncDecay


num_periods = 4000
num_divisions = 20
total_obs = num_periods*num_divisions
np.random.seed(seed=6)
perm_shocks = (var_perm_test/num_divisions)**0.5*np.random.normal(size=total_obs)
tran_shocks = (var_tran_test*num_divisions)**0.5*np.random.normal(size=total_obs)
perm_inc  = np.cumsum(perm_shocks)
perm_cons = ins_perm_test*perm_inc
tran_cons = np.zeros_like(perm_cons)
tran_inc = np.zeros_like(perm_cons)
decay = np.exp(-np.arange(total_obs)*theta_test/num_divisions)
decay_inc = np.exp(-np.arange(total_obs)*Omega_test/num_divisions)
for t in range(total_obs):
    tran_cons[t] = np.sum(theta_test*ins_tran_test/(1.0-np.exp(-theta_test)) * decay[0:t+1] * np.flip(tran_shocks[0:t+1],0))/num_divisions
    tran_inc[t]  = np.sum(Omega_test              /(1.0-np.exp(-Omega_test)) * decay_inc[0:t+1] * np.flip(tran_shocks[0:t+1],0))/num_divisions

#tran_cons = ins_tran_test*np.cumsum(tran_shocks)/num_divisions  # This is a random walk
#tran_inc = tran_shocks
total_inc  = perm_inc  + tran_inc
total_cons = perm_cons + tran_cons
observed_inc = np.zeros(num_periods)
observed_tran_inc = np.zeros(num_periods)
observed_tran_shocks = np.zeros(num_periods)
observed_cons = np.zeros(num_periods)
for T1 in range(num_periods):
    observed_inc[T1] = np.sum(total_inc[num_divisions*T1:num_divisions*(T1+1)])/num_divisions
    observed_cons[T1] = total_cons[(T1+1)*num_divisions-1]
    observed_tran_inc[T1] = np.sum(tran_inc[num_divisions*T1:num_divisions*(T1+1)])/num_divisions
    observed_tran_shocks[T1] = np.sum(tran_shocks[num_divisions*T1:num_divisions*(T1+1)])/num_divisions
#for T in range(num_periods):    
#    for i in range(np.min([(T+1)*num_divisions,50])):
#        observed_cons[T] += theta_test*ins_tran_test/(1.0-np.exp(-theta_test)) * np.exp(-i*theta_test/num_divisions) * tran_shocks[(T+1)*num_divisions-i-1]
delta_inc = np.diff(observed_inc)
delta_cons = np.diff(observed_cons)

T=14
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
 var_c_error_se_BPP, theta_BPP, theta_se_BPP, varcsi_BPP, varcsi_se_BPP , Omega_BPP, Omega_se_BPP\
  = Parameter_estimation_simple('ExpConsIncDecay', cov_vector_test, omega_test, T, taste=1) 
  
#var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
# ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
# var_c_error_se_BPP, theta_BPP, theta_se_BPP, varcsi_BPP, varcsi_se_BPP \
#  = Parameter_estimation_simple('Simple', cov_vector_test, omega_test, T, taste=1) 
  
init_params = np.array([0.0,var_perm_test,var_tran_test,ins_perm_test,ins_tran_test,0.0,theta_test])
implied_cov = implied_cov_ExpConsDecay(init_params, 1, T)

from min_distance_replication import Parameter_estimation
#Next replicate BPP
print('Replicate BPP')
var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
 ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
 var_c_error_se_BPP, teta_BPP, teta_se_BPP, varcsi_BPP, varcsi_se_BPP \
  = Parameter_estimation('BPP', cov_vector_test, omega_test, T, ma=1, taste=1, varying_ins=0) 
