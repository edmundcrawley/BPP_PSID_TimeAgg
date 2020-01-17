
**************************CREATION OF CEX FILE ********************************************************
do readmemb_AER.do
do readfmly_AER.do
do create_bigfmly_AER.do
do readmtab_AER.do
do create_bigmtab_AER.do
do cexall_new_AER.do
*******************************************************************************************************

**************************CREATION OF PSID FILE *******************************************************
do create1_AER			/*Extracts variables from the family files PSID*/
do create2_AER			/*Extracts variables from the individual file PSID*/
do adjust_AER			/*Select sample PSID*/
*******************************************************************************************************

**************************IMPUTATION, GRAPHS, SUMMARY STAtS *******************************************
do impute_AER			/*Imputes consumption into PSID from CEX demand estimation + Table 2*/
do graphs_AER			/*Graphs of the paper*/
do mean_comparison_AER		/*Table 1 of the paper*/
*******************************************************************************************************

**************************RESULTS *********************************************************************
do mindist_AER			/*Table 3,4,5,6,7*/
do mindist_extra_AER		/*Table 8*/
*******************************************************************************************************

