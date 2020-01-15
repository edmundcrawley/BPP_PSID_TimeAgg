"""
Does direct estimation of phi and psi using 8 and 9
"""

import numpy as np
import pandas as pd

empirical_input_file= Path("./InputFiles/CohA.csv")

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
T      =int(np.max(all_data[:,col_year])-np.min(all_data[:,col_year])+1 )
y      =np.shape(all_data)[0] 
num_hh = int(y/T)

#replace zeros with NAN when missing
for j in range(y):
    if all_data[j,coldc_dif]==0.0:
        all_data[j,colc_dif]=np.nan
for j in range(y):
    if all_data[j,coldy_dif]==0.0:
        all_data[j,coly_dif]=np.nan

# first row is delta_c(t)
# second is    delta_y(-1) + delta_y(0) + delta_y(1) 
# third row is delta_y(t)
# fourth row is delta_y(t+1)
deltas = np.zeros((y,4))*np.nan
deltas[:,0] = all_data[:,colc_dif]
deltas[:,2] = all_data[:,coly_dif]
deltas[:-1,3] = all_data[1:,coly_dif]
for i in range(num_hh):
    for t in np.array(range(T-2))+1:
        deltas[i*T+t,1] = all_data[i*T+t-1,coly_dif] + all_data[i*T+t,coly_dif] + all_data[i*T+t+1,coly_dif]
        
deltas_no_nan = deltas[~np.isnan(deltas).any(axis=1)]

cov_matrix = np.cov(np.transpose(deltas_no_nan))
phi_BPP =  cov_matrix[0,1]/cov_matrix[2,1]
psi_BPP =  cov_matrix[0,3]/cov_matrix[2,3]


