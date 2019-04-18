"""
Code to replicate Table 6 from BPP using both their methodology and the updated one
"""
import sys 
import os
sys.path.insert(0, os.path.abspath('../'))
from create_moments import create_moment_vector
from min_distance_replication import Parameter_estimation
import matplotlib.pyplot as plt
import numpy as np

time_agg_taste = 1
time_agg_ma = 0

###############################################################################
#First create the empirical moments for whole sample
print('Create Moment Vector from Inputs for All')
c_vector, omega, T = create_moment_vector(".\InputFiles\CohA.csv")

#Next replicate BPP
print('Replicate BPP')
var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
 ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
 var_c_error_se_BPP, teta_BPP, teta_se_BPP, varcsi_BPP, varcsi_se_BPP \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0) 
 
#Then do time aggregated version
print('Do Time Aggregated Version')
var_perm_TimeAgg, var_perm_se_TimeAgg, var_tran_TimeAgg, var_tran_se_TimeAgg, ins_perm_TimeAgg, \
 ins_perm_se_TimeAgg, ins_tran_TimeAgg, ins_tran_se_TimeAgg, var_c_error_TimeAgg, \
 var_c_error_se_TimeAgg, teta_TimeAgg, teta_se_TimeAgg, varcsi_TimeAgg, varcsi_se_TimeAgg \
  = Parameter_estimation('TimeAgg', c_vector, omega, T, ma=time_agg_ma, taste=time_agg_taste, varying_ins=0) 
  
#Replicate BPP with no transitory persistence
print('Replicate BPP with no transitory persistence')
var_perm_BPP_nopers, var_perm_se_BPP_nopers, var_tran_BPP_nopers, var_tran_se_BPP_nopers, ins_perm_BPP_nopers, \
 ins_perm_se_BPP_nopers, ins_tran_BPP_nopers, ins_tran_se_BPP_nopers, var_c_error_BPP_nopers, \
 var_c_error_se_BPP_nopers, teta_BPP_nopers, teta_se_BPP_nopers, varcsi_BPP_nopers, varcsi_se_BPP_nopers \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=0, taste=1, varying_ins=0) 
 
#Then do time aggregated version with uniformly distributed transitory persistence
print('Do Time Aggregated Version with Uniform Transitory Persistence')
var_perm_TimeAgg_unif, var_perm_se_TimeAgg_unif, var_tran_TimeAgg_unif, var_tran_se_TimeAgg_unif, ins_perm_TimeAgg_unif, \
 ins_perm_se_TimeAgg_unif, ins_tran_TimeAgg_unif, ins_tran_se_TimeAgg_unif, var_c_error_TimeAgg_unif, \
 var_c_error_se_TimeAgg_unif, teta_TimeAgg_unif, teta_se_TimeAgg_unif, varcsi_TimeAgg_unif, varcsi_se_TimeAgg_unif \
  = Parameter_estimation('TimeAgg_uniform', c_vector, omega, T, ma=1, taste=time_agg_taste, varying_ins=0) 
  
#Then do time aggregated version with linearly decaying  transitory persistence
print('Do Time Aggregated Version with Linear Decay Transitory Persistence')
var_perm_TimeAgg_lineardecay, var_perm_se_TimeAgg_lineardecay, var_tran_TimeAgg_lineardecay, var_tran_se_TimeAgg_lineardecay, ins_perm_TimeAgg_lineardecay, \
 ins_perm_se_TimeAgg_lineardecay, ins_tran_TimeAgg_lineardecay, ins_tran_se_TimeAgg_lineardecay, var_c_error_TimeAgg_lineardecay, \
 var_c_error_se_TimeAgg_lineardecay, teta_TimeAgg_lineardecay, teta_se_TimeAgg_lineardecay, varcsi_TimeAgg_lineardecay, varcsi_se_TimeAgg_lineardecay \
  = Parameter_estimation('TimeAgg_lineardecay', c_vector, omega, T, ma=1, taste=time_agg_taste, varying_ins=0) 
  
print('Do Time Aggregated Version with discrete two-shot persistence')
var_perm_TimeAgg_twoshot, var_perm_se_TimeAgg_twoshot, var_tran_TimeAgg_twoshot, var_tran_se_TimeAgg_twoshot, ins_perm_TimeAgg_twoshot, \
 ins_perm_se_TimeAgg_twoshot, ins_tran_TimeAgg_twoshot, ins_tran_se_TimeAgg_twoshot, var_c_error_TimeAgg_twoshot, \
 var_c_error_se_TimeAgg_twoshot, teta_TimeAgg_twoshot, teta_se_TimeAgg_twoshot, varcsi_TimeAgg_twoshot, varcsi_se_TimeAgg_twoshot \
  = Parameter_estimation('TimeAgg_twoshot', c_vector, omega, T, ma=1, taste=time_agg_taste, varying_ins=0) 
###############################################################################   
###############################################################################
#Empirical moments for non-college
print('Replicate BPP for Non-College')
c_vector_NC, omega_NC, T = create_moment_vector(".\InputFiles\CohA_nocollege.csv")

#Next replicate BPP
var_perm_BPP_NC, var_perm_se_BPP_NC, var_tran_BPP_NC, var_tran_se_BPP_NC, ins_perm_BPP_NC, \
 ins_perm_se_BPP_NC, ins_tran_BPP_NC, ins_tran_se_BPP_NC, var_c_error_BPP_NC, \
 var_c_error_se_BPP_NC, teta_BPP_NC, teta_se_BPP_NC, varcsi_BPP_NC, varcsi_se_BPP_NC \
  = Parameter_estimation('BPP', c_vector_NC, omega_NC, T, ma=1, taste=1, varying_ins=0) 
 
#Then do time aggregated version
print('Do Time Aggregated for Non-College')
var_perm_TimeAgg_NC, var_perm_se_TimeAgg_NC, var_tran_TimeAgg_NC, var_tran_se_TimeAgg_NC, ins_perm_TimeAgg_NC, \
 ins_perm_se_TimeAgg_NC, ins_tran_TimeAgg_NC, ins_tran_se_TimeAgg_NC, var_c_error_TimeAgg_NC, \
 var_c_error_se_TimeAgg_NC, teta_TimeAgg_NC, teta_se_TimeAgg_NC, varcsi_TimeAgg_NC, varcsi_se_TimeAgg_NC \
  = Parameter_estimation('TimeAgg', c_vector_NC, omega_NC, T, ma=time_agg_ma, taste=time_agg_taste, varying_ins=0) 
###############################################################################
#Empirical moments for college graduates
print('Replicate BPP for College Graduates')
c_vector_C, omega_C, T = create_moment_vector(".\InputFiles\CohA_college.csv")

#Next replicate BPP
var_perm_BPP_C, var_perm_se_BPP_C, var_tran_BPP_C, var_tran_se_BPP_C, ins_perm_BPP_C, \
 ins_perm_se_BPP_C, ins_tran_BPP_C, ins_tran_se_BPP_C, var_c_error_BPP_C, \
 var_c_error_se_BPP_C, teta_BPP_C, teta_se_BPP_C, varcsi_BPP_C, varcsi_se_BPP_C \
  = Parameter_estimation('BPP', c_vector_C, omega_C, T, ma=1, taste=1, varying_ins=0) 
 
#Then do time aggregated version
print('Do Time Aggregated for College Graduates')
var_perm_TimeAgg_C, var_perm_se_TimeAgg_C, var_tran_TimeAgg_C, var_tran_se_TimeAgg_C, ins_perm_TimeAgg_C, \
 ins_perm_se_TimeAgg_C, ins_tran_TimeAgg_C, ins_tran_se_TimeAgg_C, var_c_error_TimeAgg_C, \
 var_c_error_se_TimeAgg_C, teta_TimeAgg_C, teta_se_TimeAgg_C, varcsi_TimeAgg_C, varcsi_se_TimeAgg_C \
  = Parameter_estimation('TimeAgg', c_vector_C, omega_C, T, ma=time_agg_ma, taste=time_agg_taste, varying_ins=0) 
###############################################################################
 
###############################################################################
#Make Table6 Replication
print('Replicate Table 6 from BPP')
def mystr1(number):
    if not np.isnan(number):
        out = "{:.4f}".format(number)
    else:
        out = ''
    return out

output = "\\begin{table}  \n"
output += "\caption{Minimum-Distance Partial Insurance and Variance Estimates}  \n"
output += "\label{table:ReplicationTable}  \n"
output += "\\begin{center}  \n"
output += "\\newsavebox{\ReplicationTable}  \n"
output += "\\resizebox{!}{.35\paperheight}{  \n"
output += "\\begin{tabular}{cccc|cc|cc}  \n"
output += "\\toprule  \n"
output += "& &  \multicolumn{2}{c}{Whole Sample} &  \multicolumn{2}{c}{No College} &  \multicolumn{2}{c}{College}  \n"
output += "\\\\ \\hline  \n"
output += "& & BPP & Time Agg.  & BPP & Time Agg. & BPP & Time Agg. \n"
output += "\\\\ \\hline  \n"
output += " $\sigma^2_{P,T}$ & 1979-1981 & " +mystr1(var_perm_BPP[0])+    " &   "+mystr1(var_perm_TimeAgg[0])+ "& " +mystr1(var_perm_BPP_NC[0])+    " &   "+mystr1(var_perm_TimeAgg_NC[0])+ "& " +mystr1(var_perm_BPP_C[0])+    " &   "+mystr1(var_perm_TimeAgg_C[0])+ " \n"
output += "\\\\ (Variance perm. shock) &     & ("+mystr1(var_perm_se_BPP[0])+ ") & ("+mystr1(var_perm_se_TimeAgg[0])+ ") & ("+mystr1(var_perm_se_BPP_NC[0])+ ") & ("+mystr1(var_perm_se_TimeAgg_NC[0])+ ") & ("+mystr1(var_perm_se_BPP_C[0])+ ") & ("+mystr1(var_perm_se_TimeAgg_C[0])+ ") \n"
for i in np.array(range(8))+1:
    output += "\\\\  & "+'{:d}'.format(1981+i)+" & " +mystr1(var_perm_BPP[i])+    " &   "+mystr1(var_perm_TimeAgg[i])+ " & " +mystr1(var_perm_BPP_NC[i])+    " &   "+mystr1(var_perm_TimeAgg_NC[i])+ " & " +mystr1(var_perm_BPP_C[i])+    " &   "+mystr1(var_perm_TimeAgg_C[i])+ " \n"
    output += "\\\\  &                    & ("+mystr1(var_perm_se_BPP[i])+ ") & ("+mystr1(var_perm_se_TimeAgg[i])+ ")  & ("+mystr1(var_perm_se_BPP_NC[i])+ ") & ("+mystr1(var_perm_se_TimeAgg_NC[i])+ ")  & ("+mystr1(var_perm_se_BPP_C[i])+ ") & ("+mystr1(var_perm_se_TimeAgg_C[i])+ ") \n"
output += "\\\\  & 1990-92 & " +mystr1(var_perm_BPP[9])+    " &   "+mystr1(var_perm_TimeAgg[9])+ " & " +mystr1(var_perm_BPP_NC[9])+    " &   "+mystr1(var_perm_TimeAgg_NC[9])+ " & " +mystr1(var_perm_BPP_C[9])+    " &   "+mystr1(var_perm_TimeAgg_C[9])+ " \n"
output += "\\\\  &         & ("+mystr1(var_perm_se_BPP[9])+ ") & ("+mystr1(var_perm_se_TimeAgg[9])+ ") & ("+mystr1(var_perm_se_BPP_NC[9])+ ") & ("+mystr1(var_perm_se_TimeAgg_NC[9])+ ") & ("+mystr1(var_perm_se_BPP_C[9])+ ") & ("+mystr1(var_perm_se_TimeAgg_C[9])+ ") \n"

output += "\\\\ \\hline  \n"
output += " $\sigma^2_{Q,T}$ & 1979      & " +mystr1(var_tran_BPP[0])+    " &   "+mystr1(var_tran_TimeAgg[0])+ " & " +mystr1(var_tran_BPP_NC[0])+    " &   "+mystr1(var_tran_TimeAgg_NC[0])+ " & " +mystr1(var_tran_BPP_C[0])+    " &   "+mystr1(var_tran_TimeAgg_C[0])+ " \n"
output += "\\\\ (Variance trans. shock) &     & ("+mystr1(var_tran_se_BPP[0])+ ") & ("+mystr1(var_tran_se_TimeAgg[0])+ ") & ("+mystr1(var_tran_se_BPP_NC[0])+ ") & ("+mystr1(var_tran_se_TimeAgg_NC[0])+ ") & ("+mystr1(var_tran_se_BPP_C[0])+ ") & ("+mystr1(var_tran_se_TimeAgg_C[0])+ ") \n"
for i in np.array(range(10))+1:
    output += "\\\\  & "+'{:d}'.format(1979+i)+" & " +mystr1(var_tran_BPP[i])+    " &   "+mystr1(var_tran_TimeAgg[i])+ " & " +mystr1(var_tran_BPP_NC[i])+    " &   "+mystr1(var_tran_TimeAgg_NC[i])+ " & " +mystr1(var_tran_BPP_C[i])+    " &   "+mystr1(var_tran_TimeAgg_C[i])+ "\n"
    output += "\\\\  &                    & ("+mystr1(var_tran_se_BPP[i])+ ") & ("+mystr1(var_tran_se_TimeAgg[i])+ ")  & ("+mystr1(var_tran_se_BPP_NC[i])+ ") & ("+mystr1(var_tran_se_TimeAgg_NC[i])+ ")  & ("+mystr1(var_tran_se_BPP_C[i])+ ") & ("+mystr1(var_tran_se_TimeAgg_C[i])+ ") \n"
output += "\\\\  & 1990-92 & " +mystr1(var_tran_BPP[11])+    " &   "+mystr1(var_tran_TimeAgg[11])+ " & " +mystr1(var_tran_BPP_NC[11])+    " &   "+mystr1(var_tran_TimeAgg_NC[11])+ " & " +mystr1(var_tran_BPP_C[11])+    " &   "+mystr1(var_tran_TimeAgg_C[11])+ " \n"
output += "\\\\  &         & ("+mystr1(var_tran_se_BPP[11])+ ") & ("+mystr1(var_tran_se_TimeAgg[11])+ ") & ("+mystr1(var_tran_se_BPP_NC[11])+ ") & ("+mystr1(var_tran_se_TimeAgg_NC[11])+ ") & ("+mystr1(var_tran_se_BPP_C[11])+ ") & ("+mystr1(var_tran_se_TimeAgg_C[11])+ ") \n"
output += "\\\\ \\hline  \n"

output += " $\\theta$ &     & " +mystr1(teta_BPP)+    " &   "+"N/A"+ " & " +mystr1(teta_BPP_NC)+    " &   "+"N/A"+ " & " +mystr1(teta_BPP_C)+    " &   "+"N/A"+ " \n"
output += "\\\\ (Serial correl. trans. shock) &     & ("+mystr1(teta_se_BPP)+ ") &  & ("+mystr1(teta_se_BPP_NC)+ ") &  & ("+mystr1(teta_se_BPP_C)+ ") &  \n"
output += "\\\\ $\sigma^2_{\\xi}$ &     & " +mystr1(varcsi_BPP)+    " &   "+mystr1(varcsi_TimeAgg)+ " & " +mystr1(varcsi_BPP_NC)+    " &   "+mystr1(varcsi_TimeAgg_NC)+ " & " +mystr1(varcsi_BPP_C)+    " &   "+mystr1(varcsi_TimeAgg_C)+ " \n"
output += "\\\\ (Variance unobs. slope heterog.) &     & ("+mystr1(varcsi_se_BPP)+ ") & ("+mystr1(varcsi_se_TimeAgg)+ ") & ("+mystr1(varcsi_se_BPP_NC)+ ") & ("+mystr1(varcsi_se_TimeAgg_NC)+ ") & ("+mystr1(varcsi_se_BPP_C)+ ") & ("+mystr1(varcsi_se_TimeAgg_C)+ ") \n"
output += "\\\\ \\hline  \n"

output += " $\\phi$ &     & " +mystr1(ins_perm_BPP[0])+    " &   "+mystr1(ins_perm_TimeAgg[0])+ " & " +mystr1(ins_perm_BPP_NC[0])+    " &   "+mystr1(ins_perm_TimeAgg_NC[0])+ " & " +mystr1(ins_perm_BPP_C[0])+    " &   "+mystr1(ins_perm_TimeAgg_C[0])+ " \n"
output += "\\\\ (Partial insurance perm. shock) &     & ("+mystr1(ins_perm_se_BPP[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg[0])+ ") & ("+mystr1(ins_perm_se_BPP_NC[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg_NC[0])+ ") & ("+mystr1(ins_perm_se_BPP_C[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg_C[0])+ ") \n"
output += "\\\\ $\\psi$ &     & " +mystr1(ins_tran_BPP[0])+    " &   "+mystr1(ins_tran_TimeAgg[0])+ " & " +mystr1(ins_tran_BPP_NC[0])+    " &   "+mystr1(ins_tran_TimeAgg_NC[0])+ " & " +mystr1(ins_tran_BPP_C[0])+    " &   "+mystr1(ins_tran_TimeAgg_C[0])+ " \n"
output += "\\\\ (Partial insurance trans. shock) &     & ("+mystr1(ins_tran_se_BPP[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg[0])+ ") & ("+mystr1(ins_tran_se_BPP_NC[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg_NC[0])+ ") & ("+mystr1(ins_tran_se_BPP_C[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg_C[0])+ ") \n"
output += "\\\\ \\hline  \n"

output += " \end{tabular}   \n"
output += " } \n "
output += "\\usebox{\ReplicationTable}  \n"
output += "\settowidth\TableWidth{\\usebox{\ReplicationTable}} % Calculate width of table so notes will match  \n"
#output += "\medskip\medskip \\vspace{0.0cm} \parbox{\TableWidth}{\small  \n"
#output += "\\textbf{Notes}: The table reports the DWMD results of the parameters of interest. I also calculate time-varying measurement error in consumption (results not reported for brevity). Standard errors in parentheses.  \n"
#output += "}  \n"
output += "\end{center}  \n"
output += "\end{table}  \n"

with open('./Tables/RepTable6.tex','w') as f:
    f.write(output)
    f.close()
###############################################################################
print('Replicate Shock Variance Graph from BPP')
#Create some graphs of Transitory and Permanent Shocks over time
plt.figure(figsize=(12,5))
plt.subplot(1, 2, 1)
var_perm_BPP_graph = np.append(np.array(2*[var_perm_BPP[0]]), np.append(var_perm_BPP, np.array(2*[var_perm_BPP[-1]])))
var_perm_se_BPP_graph = np.append(np.array(2*[var_perm_se_BPP[0]]), np.append(var_perm_se_BPP, np.array(2*[var_perm_se_BPP[-1]])))
var_perm_TimeAgg_graph = np.append(np.array(2*[var_perm_TimeAgg[0]]), np.append(var_perm_TimeAgg, np.array(2*[var_perm_TimeAgg[-1]])))
var_perm_se_TimeAgg_graph = np.append(np.array(2*[var_perm_se_TimeAgg[0]]), np.append(var_perm_se_TimeAgg, np.array(2*[var_perm_se_TimeAgg[-1]])))
years = np.array(range(14))+1979
plt.plot(years, var_perm_TimeAgg_graph, label='Time Agg.', color='r')
plt.plot(years, var_perm_BPP_graph, label='BPP', color='c',linestyle='--')
#plt.plot(years, var_perm_BPP_graph+1.96*var_perm_se_BPP_graph, color='r',linestyle='--')
#plt.plot(years, var_perm_BPP_graph-1.96*var_perm_se_BPP_graph, color='r',linestyle='--')
#plt.plot(years, var_perm_TimeAgg_graph+1.96*var_perm_se_TimeAgg_graph, color='c',linestyle='--')
#plt.plot(years, var_perm_TimeAgg_graph-1.96*var_perm_se_TimeAgg_graph, color='c',linestyle='--')
plt.ylim(0.0,0.06)
plt.title("Variance of Permanent Shocks")
plt.legend(loc='lower left',ncol=2)

plt.subplot(1, 2, 2)
var_tran_BPP_graph =  np.append(var_tran_BPP, np.array(2*[var_tran_BPP[-1]]))
var_tran_se_BPP_graph = np.append(var_tran_se_BPP, np.array(2*[var_tran_se_BPP[-1]]))
var_tran_TimeAgg_graph = np.append(var_tran_TimeAgg, np.array(2*[var_tran_TimeAgg[-1]]))
var_tran_se_TimeAgg_graph = np.append(var_tran_se_TimeAgg, np.array(2*[var_tran_se_TimeAgg[-1]]))
years = np.array(range(14))+1979
plt.plot(years, var_tran_TimeAgg_graph, label='Time Agg.', color='r')
plt.plot(years, var_tran_BPP_graph, label='BPP', color='c',linestyle='--')
#plt.plot(years, var_tran_BPP_graph+1.96*var_tran_se_BPP_graph, color='r',linestyle='--')
#plt.plot(years, var_tran_BPP_graph-1.96*var_tran_se_BPP_graph, color='r',linestyle='--')
#plt.plot(years, var_tran_TimeAgg_graph+1.96*var_tran_se_TimeAgg_graph, color='c',linestyle='--')
#plt.plot(years, var_tran_TimeAgg_graph-1.96*var_tran_se_TimeAgg_graph, color='c',linestyle='--')
plt.ylim(0.0,0.06)
plt.title("Variance of Transitory Shocks")
#plt.legend(loc='lower left',ncol=2)

plt.tight_layout()
plt.savefig('./Figures/ShockVariances1980s.pdf')


###############################################################################
#Find parameter estimates for Tables 7 and 8
print('Calculate Parameters to replicate tables 7 and 8')
input_files = ['\CohA_Table7Col2.csv','\CohA_Table7Col3.csv','\CohA_Table8Col2.csv','\CohA_Table8Col3.csv','\CohA_Table8Col4.csv','\CohA_Table8Col5.csv','\CohA_Table8Col6.csv']
ins_perm_BPP_othertables = np.zeros((2,7))+np.nan
ins_tran_BPP_othertables = np.zeros((2,7))+np.nan
ins_perm_TimeAgg_othertables = np.zeros((2,7))+np.nan
ins_tran_TimeAgg_othertables = np.zeros((2,7))+np.nan
for i in range(7):
    c_vector, omega, T = create_moment_vector(".\InputFiles"+input_files[i])
    
    #Next replicate BPP
    a, b, c, d, ins_perm_BPP_othertables[0,i], \
      ins_perm_BPP_othertables[1,i],  ins_tran_BPP_othertables[0,i],  ins_tran_BPP_othertables[1,i], e, \
     f, g, h, ii, jj \
      = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0) 
     
    #Then do time aggregated version
    a, b, c, d, ins_perm_TimeAgg_othertables[0,i], \
     ins_perm_TimeAgg_othertables[1,i], ins_tran_TimeAgg_othertables[0,i], ins_tran_TimeAgg_othertables[1,i], e, \
     f, g, h, ii, jj \
      = Parameter_estimation('TimeAgg', c_vector, omega, T, ma=time_agg_ma, taste=time_agg_taste, varying_ins=0) 
###############################################################################   
#Replicate table 7 
print('Replicate table 7')
output = "\\begin{table}  \n"
output += "\caption{Minimum-Distance Partial Insurance and Variance Estimates}  \n"
output += "\label{table:ReplicationTable7}  \n"
output += "\\begin{center}  \n"
output += "\\newsavebox{\ReplicationTableSeven}  \n"
output += "\\resizebox{0.7\paperwidth}{!}{  \n"
output += "\\begin{tabular}{lcc|cc|cc}  \n"
output += "\\toprule  \n"
output += "Consumption: & \multicolumn{2}{c}{Nondurable} & \multicolumn{2}{c}{Nondurable} & \multicolumn{2}{c}{Nondurable} \n"
output += "\\\\ Income: & \multicolumn{2}{c}{net income} & \multicolumn{2}{c}{earnings only} & \multicolumn{2}{c}{male earnings} \n"
output += "\\\\ Sample: & \multicolumn{2}{c}{baseline} & \multicolumn{2}{c}{baseline} & \multicolumn{2}{c}{baseline} \n"
output += "\\\\ \\hline  \n"
output += "&  BPP & Time Agg.  & BPP & Time Agg. & BPP & Time Agg. \n"
output += "\\\\ \\hline  \n"

output += " $\\phi$ &      " +mystr1(ins_perm_BPP[0])+    " &   "+mystr1(ins_perm_TimeAgg[0])+ " & " +mystr1(ins_perm_BPP_othertables[0,0])+    " &   "+mystr1(ins_perm_TimeAgg_othertables[0,0])+ " & " +mystr1(ins_perm_BPP_othertables[0,1])+    " &   "+mystr1(ins_perm_TimeAgg_othertables[0,1])+ " \n"
output += "\\\\ (Partial insurance perm. shock) &      ("+mystr1(ins_perm_se_BPP[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg[0])+ ") & ("+mystr1(ins_perm_BPP_othertables[1,0])+ ") & ("+mystr1(ins_perm_TimeAgg_othertables[1,0])+ ") & ("+mystr1(ins_perm_BPP_othertables[1,1])+ ") & ("+mystr1(ins_perm_TimeAgg_othertables[1,1])+ ") \n"
output += "\\\\ $\\psi$      &" +mystr1(ins_tran_BPP[0])+    " &   "+mystr1(ins_tran_TimeAgg[0])+ " & " +mystr1(ins_tran_BPP_othertables[0,0])+    " &   "+mystr1(ins_tran_TimeAgg_othertables[0,0])+ " & " +mystr1(ins_tran_BPP_othertables[0,1])+    " &   "+mystr1(ins_tran_TimeAgg_othertables[0,1])+ " \n"
output += "\\\\ (Partial insurance trans. shock)     & ("+mystr1(ins_tran_se_BPP[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg[0])+ ") & ("+mystr1(ins_tran_BPP_othertables[1,0])+ ") & ("+mystr1(ins_tran_TimeAgg_othertables[1,0])+ ") & ("+mystr1(ins_tran_BPP_othertables[1,1])+ ") & ("+mystr1(ins_tran_TimeAgg_othertables[1,1])+ ") \n"
output += "\\\\ \\hline  \n"

output += " \end{tabular}   \n"
output += " } \n "
output += "\\usebox{\ReplicationTableSeven}  \n"
output += "\settowidth\TableWidth{\\usebox{\ReplicationTableSeven}} % Calculate width of table so notes will match  \n"
#output += "\medskip\medskip \\vspace{0.0cm} \parbox{\TableWidth}{\small  \n"
#output += "\\textbf{Notes}: The table reports the DWMD results of the parameters of interest. I also calculate time-varying measurement error in consumption (results not reported for brevity). Standard errors in parentheses.  \n"
#output += "}  \n"
output += "\end{center}  \n"
output += "\end{table}  \n"

with open('./Tables/RepTable7.tex','w') as f:
    f.write(output)
    f.close()
    
###############################################################################
#Replicate table 8 
print('Replicate table 8')
output = "\\begin{table}  \n"
output += "\caption{Minimum-Distance Partial Insurance and Variance Estimates}  \n"
output += "\label{table:ReplicationTable8}  \n"
output += "\\begin{center}  \n"
output += "\\newsavebox{\ReplicationTableEight}  \n"
output += "\\resizebox{0.7\paperwidth}{!}{  \n"
output += "\\begin{tabular}{lcc|cc|cc}  \n"
output += "\\toprule  \n"
output += "Consumption: & \multicolumn{2}{c}{Nondurable} & \multicolumn{2}{c}{Nondurable} & \multicolumn{2}{c}{Nondurable} \n"
output += "\\\\ Income: & \multicolumn{2}{c}{net income} & \multicolumn{2}{c}{excluding help} & \multicolumn{2}{c}{net income} \n"
output += "\\\\ Sample: & \multicolumn{2}{c}{baseline} & \multicolumn{2}{c}{baseline} & \multicolumn{2}{c}{low wealth} \n"
output += "\\\\ \\hline  \n"
output += "&  BPP & Time Agg.  & BPP & Time Agg. & BPP & Time Agg. \n"
output += "\\\\ \\hline  \n"
output += " $\\phi$ &      " +mystr1(ins_perm_BPP[0])+    " &   "+mystr1(ins_perm_TimeAgg[0])+ " & " +mystr1(ins_perm_BPP_othertables[0,2])+    " &   "+mystr1(ins_perm_TimeAgg_othertables[0,2])+ " & " +mystr1(ins_perm_BPP_othertables[0,3])+    " &   "+mystr1(ins_perm_TimeAgg_othertables[0,3])+ " \n"
output += "\\\\ (Partial insurance perm. shock) &      ("+mystr1(ins_perm_se_BPP[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg[0])+ ") & ("+mystr1(ins_perm_BPP_othertables[1,2])+ ") & ("+mystr1(ins_perm_TimeAgg_othertables[1,2])+ ") & ("+mystr1(ins_perm_BPP_othertables[1,3])+ ") & ("+mystr1(ins_perm_TimeAgg_othertables[1,3])+ ") \n"
output += "\\\\ $\\psi$      &" +mystr1(ins_tran_BPP[0])+    " &   "+mystr1(ins_tran_TimeAgg[0])+ " & " +mystr1(ins_tran_BPP_othertables[0,2])+    " &   "+mystr1(ins_tran_TimeAgg_othertables[0,2])+ " & " +mystr1(ins_tran_BPP_othertables[0,3])+    " &   "+mystr1(ins_tran_TimeAgg_othertables[0,3])+ " \n"
output += "\\\\ (Partial insurance trans. shock)     & ("+mystr1(ins_tran_se_BPP[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg[0])+ ") & ("+mystr1(ins_tran_BPP_othertables[1,2])+ ") & ("+mystr1(ins_tran_TimeAgg_othertables[1,2])+ ") & ("+mystr1(ins_tran_BPP_othertables[1,3])+ ") & ("+mystr1(ins_tran_TimeAgg_othertables[1,3])+ ") \n"

output += "\\\\ \\hline  \n"
output += "\\\\  \n"

output += "\\toprule  \n"
output += "Consumption: & \multicolumn{2}{c}{Nondurable} & \multicolumn{2}{c}{Total} & \multicolumn{2}{c}{Nondurable} \n"
output += "\\\\ Income: & \multicolumn{2}{c}{net income} & \multicolumn{2}{c}{net income} & \multicolumn{2}{c}{net income} \n"
output += "\\\\ Sample: & \multicolumn{2}{c}{high wealth} & \multicolumn{2}{c}{low wealth} & \multicolumn{2}{c}{baseline+SEO} \n"
output += "\\\\ \\hline  \n"
output += "&  BPP & Time Agg.  & BPP & Time Agg. & BPP & Time Agg. \n"
output += "\\\\ \\hline  \n"

output += " $\\phi$ &      " +mystr1(ins_perm_BPP_othertables[0,4])+    " &   "+mystr1(ins_perm_TimeAgg_othertables[0,4])+ " & " +mystr1(ins_perm_BPP_othertables[0,5])+    " &   "+mystr1(ins_perm_TimeAgg_othertables[0,5])+ " & " +mystr1(ins_perm_BPP_othertables[0,6])+    " &   "+mystr1(ins_perm_TimeAgg_othertables[0,6])+ " \n"
output += "\\\\ (Partial insurance perm. shock) &      ("+mystr1(ins_perm_BPP_othertables[1,4])+ ") & ("+mystr1(ins_perm_TimeAgg_othertables[1,4])+") & ("+mystr1(ins_perm_BPP_othertables[1,5])+ ") & ("+mystr1(ins_perm_TimeAgg_othertables[1,5])+ ") & ("+mystr1(ins_perm_BPP_othertables[1,6])+ ") & ("+mystr1(ins_perm_TimeAgg_othertables[1,6])+ ") \n"
output += "\\\\ $\\psi$      &"  +mystr1(ins_tran_BPP_othertables[0,4])+    " &   "+mystr1(ins_tran_TimeAgg_othertables[0,4])+ " & " +mystr1(ins_tran_BPP_othertables[0,5])+    " &   "+mystr1(ins_tran_TimeAgg_othertables[0,5])+ " & " +mystr1(ins_tran_BPP_othertables[0,6])+    " &   "+mystr1(ins_tran_TimeAgg_othertables[0,6])+ " \n"
output += "\\\\ (Partial insurance trans. shock)     & ("+mystr1(ins_tran_BPP_othertables[1,4])+ ") & ("+mystr1(ins_tran_TimeAgg_othertables[1,4])+ ") & ("+mystr1(ins_tran_BPP_othertables[1,5])+ ") & ("+mystr1(ins_tran_TimeAgg_othertables[1,5])+ ") & ("+mystr1(ins_tran_BPP_othertables[1,6])+ ") & ("+mystr1(ins_tran_TimeAgg_othertables[1,6])+ ") \n"
output += "\\\\ \\hline  \n"

output += " \end{tabular}   \n"

output += " } \n "
output += "\\usebox{\ReplicationTableEight}  \n"
output += "\settowidth\TableWidth{\\usebox{\ReplicationTableEight}} % Calculate width of table so notes will match  \n"
#output += "\medskip\medskip \\vspace{0.0cm} \parbox{\TableWidth}{\small  \n"
#output += "\\textbf{Notes}: The table reports the DWMD results of the parameters of interest. I also calculate time-varying measurement error in consumption (results not reported for brevity). Standard errors in parentheses.  \n"
#output += "}  \n"
output += "\end{center}  \n"
output += "\end{table}  \n"

with open('./Tables/RepTable8.tex','w') as f:
    f.write(output)
    f.close()
    
###############################################################################   
#Table to show effect of transitory persistence
print('Table to show effect of transitory persistence')
output = "\\begin{table}  \n"
output += "\caption{Minimum-Distance Partial Insurance and Variance Estimates}  \n"
output += "\label{table:Persistence}  \n"
output += "\\begin{center}  \n"
output += "\\newsavebox{\Persistence}  \n"
output += "\\resizebox{0.7\paperwidth}{!}{  \n"
output += "\\begin{tabular}{lcc|cccc}  \n"
output += "\\toprule  \n"
output += " & \multicolumn{2}{c}{BPP} & \multicolumn{4}{c}{Time Agg.}  \n"
output += "\\\\ \\hline  \n"
output += "\\\\ Persistence Type:  & None & MA(1) & None & Two-shot & Uniform & Linear Decay \n"
output += "\\\\ \\hline  \n"

output += " $\\phi$ &                               " +mystr1(ins_perm_BPP_nopers[0])+    " &   "+mystr1(ins_perm_BPP[0])+ " & " +mystr1(ins_perm_TimeAgg[0])+   " &   "+mystr1(ins_perm_TimeAgg_twoshot[0])+  " &   "+mystr1(ins_perm_TimeAgg_unif[0])+ " & " +mystr1(ins_perm_TimeAgg_lineardecay[0])+    "  \n"
output += "\\\\ (Partial insurance perm. shock) &      ("+mystr1(ins_perm_se_BPP_nopers[0])+ ") & ("+mystr1(ins_perm_se_BPP[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg_twoshot[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg_unif[0])+ ") & ("+mystr1(ins_perm_se_TimeAgg_lineardecay[0])+ ") \n"

output += "\\\\ $\\psi$ &                               " +mystr1(ins_tran_BPP_nopers[0])+    " &   "+mystr1(ins_tran_BPP[0])+ " & " +mystr1(ins_tran_TimeAgg[0])+   " &   "+mystr1(ins_tran_TimeAgg_twoshot[0])+  " &   "+mystr1(ins_tran_TimeAgg_unif[0])+ " & " +mystr1(ins_tran_TimeAgg_lineardecay[0])+    "  \n"
output += "\\\\ (Partial insurance tran. shock) &      ("+mystr1(ins_tran_se_BPP_nopers[0])+ ") & ("+mystr1(ins_tran_se_BPP[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg_twoshot[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg_unif[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg_lineardecay[0])+ ") \n"

output += "\\\\ $\\theta$ or $\\tau$ &            " +"N/A"+    " &   "+mystr1(teta_BPP)+ " & " +"N/A"+   " &   "+mystr1(teta_TimeAgg_twoshot)+  " &   "+mystr1(teta_TimeAgg_unif)+ " & " +mystr1(teta_TimeAgg_lineardecay)+    "  \n"
output += "\\\\ (Degree of Persistence) &      ("+mystr1(teta_se_BPP_nopers)+ ") & ("+mystr1(teta_se_BPP)+ ") & ("+mystr1(teta_se_TimeAgg)+ ") & ("+mystr1(teta_se_TimeAgg_twoshot)+ ") & ("+mystr1(teta_se_TimeAgg_unif)+ ") & ("+mystr1(teta_se_TimeAgg_lineardecay)+ ") \n"


output += "\\\\ \\hline  \n"

output += " \end{tabular}   \n"
output += " } \n "
output += "\\usebox{\Persistence}  \n"
output += "\settowidth\TableWidth{\\usebox{\Persistence}} % Calculate width of table so notes will match  \n"
#output += "\medskip\medskip \\vspace{0.0cm} \parbox{\TableWidth}{\small  \n"
#output += "\\textbf{Notes}: The table reports the DWMD results of the parameters of interest. I also calculate time-varying measurement error in consumption (results not reported for brevity). Standard errors in parentheses.  \n"
#output += "}  \n"
output += "\end{center}  \n"
output += "\end{table}  \n"

with open('./Tables/Persistence.tex','w') as f:
    f.write(output)
    f.close()