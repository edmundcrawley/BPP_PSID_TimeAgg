import statsmodels.api as sm
import pandas as pd

all_data = pd.read_csv(Path("./InputFiles/CohA.csv"), delimiter=',')

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
col_uy     =9                      #/*Residual of log income
col_uc     =10                      #/*Residual of log consumption

# length of panel. Assumes all individuals appear in all years
T  =int(np.max(all_data['year'])-np.min(all_data['year'])+1 )
num_hh      =int(np.shape(all_data)[0] /T)

inc_growth = np.zeros((num_hh,1))
con_growth = np.zeros((num_hh,1))
beta_n = np.zeros((T-1,1))
for n in range(T-1):
    for i in range(num_hh):
            inc_growth[i] = all_data.iat[i*T + n+1, col_uy] - all_data.iat[i*T, col_uy]
            con_growth[i] = all_data.iat[i*T + n+1, col_uc] - all_data.iat[i*T, col_uc]
    # Note the difference in argument order
    if n!=7 and n!=8:
        model = sm.OLS(con_growth, inc_growth,missing='drop').fit()
        beta_n[n] = model.params[0]