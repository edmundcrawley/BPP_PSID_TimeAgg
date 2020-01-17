***Table 1

cd C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files
clear
set more off
cap log close
log using C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\mean_comparison,replace t

u data3,clear
keep if seo==0
keep if age>=30 & age<=65

gen ndur=exp(lc)
gen income_at=y-ftax 
ren y income_bt
gen totfood=fout+food

su age ncomp kid white edd* regd* part* income_at income_bt totfood if year==1980
su age ncomp kid white edd* regd* part* income_at income_bt totfood if year==1983
su age ncomp kid white edd* regd* part* income_at income_bt totfood if year==1986
su age ncomp kid white edd* regd* part* income_at income_bt totfood if year==1989
su age ncomp kid white edd* regd* part* income_at income_bt totfood if year==1992


u cexall,clear

keep if age>=30 & age<=65
gen     educ=1 if educh==7| educh ==1 | educh ==2 
replace educ=2 if educh ==3
replace educ=3 if educh ==4| educh ==5 | educh ==6
gen    white=race =="1" 				/* dummy for white husband */
tab educ  ,gen(edd)
tab region,gen(regd)

gen     coh=1  if yb>=1950 & yb<=1959
replace coh=2  if yb>=1940 & yb<=1949
replace coh=3  if yb>=1930 & yb<=1939
replace coh=4  if yb>=1920 & yb<=1929
tab coh,gen(cohd)
drop if coh==.

replace weekh =0 if weekh ==.
replace hourh =0 if hourh ==.
gen yhourh=weekh*hourh
gen parth=yhourh>0

replace weekw =0 if weekw ==.
replace hourw =0 if hourw ==.
gen yhourw=weekw*hourw
gen partw=yhourw>0

gen totfood=food+fout
gen lc=ln(ndur)

su age ncomp kid white edd* regd* part* income_at income totfood if year==1980
su age ncomp kid white edd* regd* part* income_at income totfood if year==1983
su age ncomp kid white edd* regd* part* income_at income totfood if year==1986
su age ncomp kid white edd* regd* part* income_at income totfood if year==1989
su age ncomp kid white edd* regd* part* income_at income totfood if year==1992

log close
clear
