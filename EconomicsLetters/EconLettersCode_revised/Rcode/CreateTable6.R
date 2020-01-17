# This file does the minimum distance optimization for BPP to create Table 6 in their paper

setwd("C:/Users/edmun/OneDrive/Documents/Research/BPP_PSID_TimeAgg/Code/Rcode")

source("./create_moments.r")
source("./min_distance_replication.r")

###############################################################################
#First create the empirical moments for whole sample
moments <- create_moments("../InputFiles/CohA.csv")

c_vector <- moments[["c_vector"]]
omega <- moments[["omega"]]
T <- moments[["T"]]


#Next replicate BPP
BPP_output = BPP_parameter_estimation(c_vector, omega, T, ma=1, taste=1, varying_ins=0) 




