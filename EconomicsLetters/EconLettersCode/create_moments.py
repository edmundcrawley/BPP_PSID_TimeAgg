# -*- coding: utf-8 -*-
"""
Replication of Gauss part of BPP 2008
This file constains a function that creates the moments and variance matrix of the moments
"""

import numpy as np
import pandas as pd


def create_moment_vector(empirical_input_file, only_income=False):
    '''
    Makes the empirical covariance moments and associated covariance matrix
    
    Parameters
    ----------
    empirical_input_file : file
        File containing income and consumption changes
    only_income : boolean
        Determines whether the output is just for income or for both the income and consumption process
    Returns
    -------
    c_vector : np.array
        Vector of covariances
    omega : np.array
        Covariance matrix for the vector of covariances
    T   : int
        Integer representing the length of the panel
    '''
    
    all_data =  np.genfromtxt(empirical_input_file, delimiter=',',skip_header =1)
    
    
    first_diff=1				#Tells whether moments are in levels (0) or FD (1)
    T1=0
    
    col_id     =0                      #/*colum where id is */ 
    col_year   =1                      #/*column where year is*/
    col_coh    =2                      #/*columns where dummies for cohorts are*/
    col_deld   =3                      #/*columns where dummies for the missing years in FD consumption*/
    coly_dif   =4                      #/*column where residuals in FD of income are*/
    coldy_dif  =5                      #/*column where dummies for no-missing residuals in FD are*/
    colc_dif   =6                      #/*column where residuals in FD of consumption are*/
    coldc_dif  =7                      #/*column where dummies for no-missing residuals in FD are*/
    col_missd  =8                      #/*number of missing consumption years in FD*/
    
    #**********First create moments just based on income
    
    # length of panel. Assumes all individuals appear in all years
    T  =int(np.max(all_data[:,col_year])-np.min(all_data[:,col_year])+1 )
    
    dif    =np.zeros([T,T])
    d_dif  =np.zeros([T,T])
    y      =np.shape(all_data)[0] 
    
    for k in range(int(y/T)):
        i=k*T
        dif_j = np.outer(all_data[i:i+T,coly_dif],all_data[i:i+T,coly_dif])
        d_dif_j=np.outer(all_data[i:i+T,coldy_dif],all_data[i:i+T,coldy_dif])
        dif = dif+dif_j
        d_dif = d_dif+d_dif_j
    dif = dif/d_dif
    income_cov = dif
    income_var = np.diag(dif)
    
    # create vector from the lower triangular elements of the covariance matrix (it is a symetric matrix)
    c_matrix=dif
    vech_indicies = np.tril_indices(np.shape(c_matrix)[0])
    c_vector=c_matrix[vech_indicies]
    dim = int((np.shape(c_matrix)[0]*(np.shape(c_matrix)[0]+1))/2)
    
    # create variance matrix for all these moments
    omega  =np.zeros([dim,dim])
    d_matrix=d_dif
    d_vector=d_matrix[vech_indicies]
    for k in range(int(y/T)):
            i=k*T
            dif_j  =np.outer(all_data[i:i+T,coly_dif],all_data[i:i+T,coly_dif])
            d_dif_j=np.outer(all_data[i:i+T,coldy_dif],all_data[i:i+T,coldy_dif])
            c_matrix_j=dif_j
            d_matrix_j=d_dif_j;
            c_vector_j=c_matrix_j[vech_indicies]
            d_vector_j=d_matrix_j[vech_indicies]
            omega=omega+(np.outer((c_vector_j-c_vector),(c_vector_j-c_vector))*(np.outer(d_vector_j,d_vector_j)))
    c_vector_income_only = c_vector
    omega_income_only=omega/np.outer(d_vector,d_vector)
    
    #**********Now create both income and consumption moments
    
    # stack consumption and income
    consumption_rows = np.concatenate((all_data[:,[colc_dif,coldc_dif,col_deld,col_id,col_year]],np.ones_like(all_data[:,[1]])),1)
    income_rows = np.concatenate((all_data[:,[coly_dif,coldy_dif]],np.zeros_like(all_data[:,[1]]),all_data[:,[col_id,col_year]],2*np.ones_like(all_data[:,[1]])),1)
    c_y_rows = np.concatenate((consumption_rows,income_rows))
    #remove rows for which consumption is missing for everyone
    c_y_rows = c_y_rows[c_y_rows[:,2]<0.5]
    # sort by id, cons/income, year (do backwards and make sure sorts are stable)
    c_y_rows = c_y_rows[c_y_rows[:,4].argsort()] # First sort doesn't need to be stable.
    c_y_rows = c_y_rows[c_y_rows[:,5].argsort(kind='mergesort')]
    c_y_rows = c_y_rows[c_y_rows[:,3].argsort(kind='mergesort')]
    
    initial_year = int(np.min(all_data[:,col_year]))
    final_year = int(np.max(all_data[:,col_year]))
    T = int(final_year-initial_year+1 )
    T1= int(T+T-all_data[0,col_missd])
    
    dif    =np.zeros([T1,T1])
    d_dif  =np.zeros([T1,T1])
    y      =np.shape(c_y_rows)[0]   
    
    for k in range(int(y/T1)):
        i=k*T1
        dif_j = np.outer(c_y_rows[i:i+T1,0],c_y_rows[i:i+T1,0])
        d_dif_j = np.outer(c_y_rows[i:i+T1,1],c_y_rows[i:i+T1,1])
        dif = dif+dif_j
        d_dif = d_dif+d_dif_j
    dif = dif/d_dif
    # create vector from the lower triangular elements of the covariance matrix (it is a symetric matrix)
    c_matrix=dif
    vech_indicies = np.tril_indices(np.shape(c_matrix)[0])
    c_vector=c_matrix[vech_indicies]
    dim = int((np.shape(c_matrix)[0]*(np.shape(c_matrix)[0]+1))/2)
    # create variance matrix for all these moments
    omega  =np.zeros([dim,dim])
    d_matrix=d_dif
    d_vector=d_matrix[vech_indicies]
    
    for k in range(int(y/T1)):
            i=k*T1
            dif_j  =np.outer(c_y_rows[i:i+T1,0],c_y_rows[i:i+T1,0])
            d_dif_j=np.outer(c_y_rows[i:i+T1,1],c_y_rows[i:i+T1,1])
            c_matrix_j=dif_j
            d_matrix_j=d_dif_j
            c_vector_j=c_matrix_j[vech_indicies]
            d_vector_j=d_matrix_j[vech_indicies]
            omega=omega+(np.outer((c_vector_j-c_vector),(c_vector_j-c_vector))*(np.outer(d_vector_j,d_vector_j)))
    c_vector_both = c_vector
    omega_both=omega/np.outer(d_vector,d_vector)
    
    
    if only_income:
        return c_vector_income_only, omega_income_only, T
    else:
        return c_vector_both, omega_both, T
    
    
    
    
    
    
    
    
    
