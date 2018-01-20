"""
Code to replicate Table 6 from BPP using both their methodology and the updated one
"""
import sys 
import os
sys.path.insert(0, os.path.abspath('../'))
from create_moments import create_moment_vector
from min_distance_replication import Parameter_estimation

###############################################################################
#First create the empirical moments for whole sample
c_vector, omega, T = create_moment_vector(".\InputFiles\CohA.csv")

#Next replicate BPP
var_perm_BPP, var_perm_se_BPP, var_tran_BPP, var_tran_se_BPP, ins_perm_BPP, \
 ins_perm_se_BPP, ins_tran_BPP, ins_tran_se_BPP, var_c_error_BPP, \
 var_c_error_se_BPP, teta_BPP, teta_se_BPP, varcsi_BPP, varcsi_se_BPP \
  = Parameter_estimation('BPP', c_vector, omega, T, ma=1, taste=1, varying_ins=0) 
 
#Then do time aggregated version
var_perm_TimeAgg, var_perm_se_TimeAgg, var_tran_TimeAgg, var_tran_se_TimeAgg, ins_perm_TimeAgg, \
 ins_perm_se_TimeAgg, ins_tran_TimeAgg, ins_tran_se_TimeAgg, var_c_error_TimeAgg, \
 var_c_error_se_TimeAgg, teta_TimeAgg, teta_se_TimeAgg, varcsi_TimeAgg, varcsi_se_TimeAgg \
  = Parameter_estimation('TimeAgg', c_vector, omega, T, ma=1, taste=1, varying_ins=0) 
###############################################################################
#Empirical moments for non-college
c_vector_NC, omega_NC, T = create_moment_vector(".\InputFiles\CohA_nocollege.csv")

#Next replicate BPP
var_perm_BPP_NC, var_perm_se_BPP_NC, var_tran_BPP_NC, var_tran_se_BPP_NC, ins_perm_BPP_NC, \
 ins_perm_se_BPP_NC, ins_tran_BPP_NC, ins_tran_se_BPP_NC, var_c_error_BPP_NC, \
 var_c_error_se_BPP_NC, teta_BPP_NC, teta_se_BPP_NC, varcsi_BPP_NC, varcsi_se_BPP_NC \
  = Parameter_estimation('BPP', c_vector_NC, omega_NC, T, ma=1, taste=1, varying_ins=0) 
 
#Then do time aggregated version
var_perm_TimeAgg_NC, var_perm_se_TimeAgg_NC, var_tran_TimeAgg_NC, var_tran_se_TimeAgg_NC, ins_perm_TimeAgg_NC, \
 ins_perm_se_TimeAgg_NC, ins_tran_TimeAgg_NC, ins_tran_se_TimeAgg_NC, var_c_error_TimeAgg_NC, \
 var_c_error_se_TimeAgg_NC, teta_TimeAgg_NC, teta_se_TimeAgg_NC, varcsi_TimeAgg_NC, varcsi_se_TimeAgg_NC \
  = Parameter_estimation('TimeAgg', c_vector_NC, omega_NC, T, ma=1, taste=1, varying_ins=0) 
###############################################################################
#Empirical moments for college graduates
c_vector_C, omega_C, T = create_moment_vector(".\InputFiles\CohA_college.csv")

#Next replicate BPP
var_perm_BPP_C, var_perm_se_BPP_C, var_tran_BPP_C, var_tran_se_BPP_C, ins_perm_BPP_C, \
 ins_perm_se_BPP_C, ins_tran_BPP_C, ins_tran_se_BPP_C, var_c_error_BPP_C, \
 var_c_error_se_BPP_C, teta_BPP_C, teta_se_BPP_C, varcsi_BPP_C, varcsi_se_BPP_C \
  = Parameter_estimation('BPP', c_vector_C, omega_C, T, ma=1, taste=1, varying_ins=0) 
 
#Then do time aggregated version
var_perm_TimeAgg_C, var_perm_se_TimeAgg_C, var_tran_TimeAgg_C, var_tran_se_TimeAgg_C, ins_perm_TimeAgg_C, \
 ins_perm_se_TimeAgg_C, ins_tran_TimeAgg_C, ins_tran_se_TimeAgg_C, var_c_error_TimeAgg_C, \
 var_c_error_se_TimeAgg_C, teta_TimeAgg_C, teta_se_TimeAgg_C, varcsi_TimeAgg_C, varcsi_se_TimeAgg_C \
  = Parameter_estimation('TimeAgg', c_vector_C, omega_C, T, ma=1, taste=1, varying_ins=0) 
###############################################################################
 
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
output += "\\begin{tabular}{cc|cc|cc|cc}  \n"
output += "\\toprule  \n"
output += "& &  \multicolumn{2}{c}{Whole Sample} &  \multicolumn{2}{c}{No College} &  \multicolumn{2}{c}{College}  \n"
output += "\\\\ \\hline  \n"
output += "& & BPP & Time Agg.  & BPP & Time Agg. & BPP & Time Agg. \n"
output += "\\\\ \\hline  \n"
output += "\\\\ $\sigma^2_{P,T}$ & 1979-1981 & " +mystr1(var_perm_BPP[0])+    " &   "+mystr1(var_perm_TimeAgg[0])+ "& " +mystr1(var_perm_BPP_NC[0])+    " &   "+mystr1(var_perm_TimeAgg_NC[0])+ "& " +mystr1(var_perm_BPP_C[0])+    " &   "+mystr1(var_perm_TimeAgg_C[0])+ " \n"
output += "\\\\ (Variance perm. shock) &     & ("+mystr1(var_perm_se_BPP[0])+ ") & ("+mystr1(var_perm_se_TimeAgg[0])+ ") & ("+mystr1(var_perm_se_BPP_NC[0])+ ") & ("+mystr1(var_perm_se_TimeAgg_NC[0])+ ") & ("+mystr1(var_perm_se_BPP_C[0])+ ") & ("+mystr1(var_perm_se_TimeAgg_C[0])+ ") \n"
for i in np.array(range(8))+1:
    output += "\\\\  & "+'{:d}'.format(1981+i)+" & " +mystr1(var_perm_BPP[i])+    " &   "+mystr1(var_perm_TimeAgg[i])+ " & " +mystr1(var_perm_BPP_NC[i])+    " &   "+mystr1(var_perm_TimeAgg_NC[i])+ " & " +mystr1(var_perm_BPP_C[i])+    " &   "+mystr1(var_perm_TimeAgg_C[i])+ " \n"
    output += "\\\\  &                    & ("+mystr1(var_perm_se_BPP[i])+ ") & ("+mystr1(var_perm_se_TimeAgg[i])+ ")  & ("+mystr1(var_perm_se_BPP_NC[i])+ ") & ("+mystr1(var_perm_se_TimeAgg_NC[i])+ ")  & ("+mystr1(var_perm_se_BPP_C[i])+ ") & ("+mystr1(var_perm_se_TimeAgg_C[i])+ ") \n"
output += "\\\\  & 1990-92 & " +mystr1(var_perm_BPP[9])+    " &   "+mystr1(var_perm_TimeAgg[9])+ " & " +mystr1(var_perm_BPP_NC[9])+    " &   "+mystr1(var_perm_TimeAgg_NC[9])+ " & " +mystr1(var_perm_BPP_C[9])+    " &   "+mystr1(var_perm_TimeAgg_C[9])+ " \n"
output += "\\\\  &         & ("+mystr1(var_perm_se_BPP[9])+ ") & ("+mystr1(var_perm_se_TimeAgg[9])+ ") & ("+mystr1(var_perm_se_BPP_NC[9])+ ") & ("+mystr1(var_perm_se_TimeAgg_NC[9])+ ") & ("+mystr1(var_perm_se_BPP_C[9])+ ") & ("+mystr1(var_perm_se_TimeAgg_C[9])+ ") \n"

output += "\\\\ \\hline  \n"
output += "\\\\ $\sigma^2_{Q,T}$ & 1979      & " +mystr1(var_tran_BPP[0])+    " &   "+mystr1(var_tran_TimeAgg[0])+ " & " +mystr1(var_tran_BPP_NC[0])+    " &   "+mystr1(var_tran_TimeAgg_NC[0])+ " & " +mystr1(var_tran_BPP_C[0])+    " &   "+mystr1(var_tran_TimeAgg_C[0])+ " \n"
output += "\\\\ (Variance trans. shock) &     & ("+mystr1(var_tran_se_BPP[0])+ ") & ("+mystr1(var_tran_se_TimeAgg[0])+ ") & ("+mystr1(var_tran_se_BPP_NC[0])+ ") & ("+mystr1(var_tran_se_TimeAgg_NC[0])+ ") & ("+mystr1(var_tran_se_BPP_C[0])+ ") & ("+mystr1(var_tran_se_TimeAgg_C[0])+ ") \n"
for i in np.array(range(10))+1:
    output += "\\\\  & "+'{:d}'.format(1979+i)+" & " +mystr1(var_tran_BPP[i])+    " &   "+mystr1(var_tran_TimeAgg[i])+ " & " +mystr1(var_tran_BPP_NC[i])+    " &   "+mystr1(var_tran_TimeAgg_NC[i])+ " & " +mystr1(var_tran_BPP_C[i])+    " &   "+mystr1(var_tran_TimeAgg_C[i])+ "\n"
    output += "\\\\  &                    & ("+mystr1(var_tran_se_BPP[i])+ ") & ("+mystr1(var_tran_se_TimeAgg[i])+ ")  & ("+mystr1(var_tran_se_BPP_NC[i])+ ") & ("+mystr1(var_tran_se_TimeAgg_NC[i])+ ")  & ("+mystr1(var_tran_se_BPP_C[i])+ ") & ("+mystr1(var_tran_se_TimeAgg_C[i])+ ") \n"
output += "\\\\  & 1990-92 & " +mystr1(var_tran_BPP[11])+    " &   "+mystr1(var_tran_TimeAgg[11])+ " & " +mystr1(var_tran_BPP_NC[11])+    " &   "+mystr1(var_tran_TimeAgg_NC[11])+ " & " +mystr1(var_tran_BPP_C[11])+    " &   "+mystr1(var_tran_TimeAgg_C[11])+ " \n"
output += "\\\\  &         & ("+mystr1(var_tran_se_BPP[11])+ ") & ("+mystr1(var_tran_se_TimeAgg[11])+ ") & ("+mystr1(var_tran_se_BPP_NC[11])+ ") & ("+mystr1(var_tran_se_TimeAgg_NC[11])+ ") & ("+mystr1(var_tran_se_BPP_C[11])+ ") & ("+mystr1(var_tran_se_TimeAgg_C[11])+ ") \n"
output += "\\\\ \\hline  \n"

output += "\\\\ $\\theta$ &     & " +mystr1(teta_BPP)+    " &   "+"N/A"+ " \n"
output += "\\\\ (Serial correl. trans. shock) &     & ("+mystr1(teta_se_BPP)+ ") &  \n"
output += "\\\\ $\sigma^2_{\\xi}$ &     & " +mystr1(varcsi_BPP)+    " &   "+mystr1(varcsi_TimeAgg)+ " \n"
output += "\\\\ (Variance unobs. slope heterog.) &     & ("+mystr1(varcsi_se_BPP)+ ") & ("+mystr1(varcsi_se_TimeAgg)+ ") \n"
output += "\\\\ \\hline  \n"

output += "\\\\ $\\phi$ &     & " +mystr1(ins_perm_BPP[0])+    " &   "+mystr1(ins_perm_TimeAgg[0])+ " \n"
output += "\\\\ (Partial insurance perm. shock) &     & ("+mystr1(ins_perm_se_BPP[0])+ ") & "+mystr1(ins_perm_se_TimeAgg[0])+ " \n"
output += "\\\\ $\\psi$ &     & " +mystr1(ins_tran_BPP[0])+    " &   "+mystr1(ins_tran_TimeAgg[0])+ " \n"
output += "\\\\ (Partial insurance trans. shock) &     & ("+mystr1(ins_tran_se_BPP[0])+ ") & ("+mystr1(ins_tran_se_TimeAgg[0])+ ") \n"
output += "\\\\ \\hline  \n"


output += " \end{tabular}   \n"
output += " } \n "
output += "\usebox{\ReplicationTable}  \n"
output += "\settowidth\TableWidth{\usebox{\ReplicationTable}} % Calculate width of table so notes will match  \n"
#output += "\medskip\medskip \\vspace{0.0cm} \parbox{\TableWidth}{\small  \n"
#output += "\\textbf{Notes}: The table reports the DWMD results of the parameters of interest. I also calculate time-varying measurement error in consumption (results not reported for brevity). Standard errors in parentheses.  \n"
#output += "}  \n"
output += "\end{center}  \n"
output += "\end{table}  \n"

with open('./Tables/RepTable6.tex','w') as f:
    f.write(output)
    f.close()
