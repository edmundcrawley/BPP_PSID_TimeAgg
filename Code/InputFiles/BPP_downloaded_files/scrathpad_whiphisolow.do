
***************************************************************************
* Experiment with regressing consumption on income in PSID
***************************************************************************
u cohA_scratch, clear
sort person year
fillin person year
gen yduy=duy!=.
replace  duy=0 if duy==.
gen yduc=duc!=.        
replace  duc=0 if duc==.
egen maxy=max(year)
egen miny=min(year)
scalar nyears=maxy-miny+1
drop maxy miny
scalar nind=_N/nyears
sort person year
capture drop id
gen id=group(nind)
compress

egen nobsdif=sum(yduy),by(year)	
drop if nobsdif<50
drop nobsdif
egen miny=min(year)
egen maxy=max(year)
egen temp1=sum(yduc),by(id)
egen nmissd=max(temp1)
replace nmissd=(maxy-miny+1)-nmissd
gen ndrod=year==1987|year==1988|year==1989
drop maxy miny temp1

sort id year
gen log_food = log(food)

reg lc logy
//get 0.43

reg uc uy
//get 0.33

reg duc duy
//get 0.11

reg log_food logy
reg lc log_food

//regress changes over increasing differences
capture drop dlogy_* 
capture drop dlc_*
capture drop duy_*
capture drop duc_*
local start_year = 2 //start in 1979 since 1978 consumption is messed up
foreach n of numlist `start_year'/13 {
	by id: gen dlogy_`n' = logy[`n'+1]-logy[`start_year']
	by id: gen dlc_`n' = lc[`n'+1]-lc[`start_year']
	by id: gen duy_`n' = uy[`n'+1]-uy[`start_year']
	by id: gen duc_`n' = uc[`n'+1]-uc[`start_year']
}
matrix beta_n_lc_logy = J(14,2,.)
matrix beta_n_duc_duy = J(14,2,.)
foreach n of numlist `start_year'/7 {
	reg dlc_`n'  dlogy_`n' 
	matrix foo = e(b)
	matrix beta_n_lc_logy[`n',1] = foo[1,1]
	matrix beta_n_lc_logy[`n',2] = foo[1,2]
}
// n=8 and 9 are missing
foreach n of numlist 10/13 {
	reg dlc_`n'  dlogy_`n'
	matrix foo = e(b)
	matrix beta_n_lc_logy[`n',1] = foo[1,1]
	matrix beta_n_lc_logy[`n',2] = foo[1,2]
}
foreach n of numlist `start_year'/7 {
	reg duc_`n'  duy_`n' 
	matrix foo = e(b)
	matrix beta_n_duc_duy[`n',1] = foo[1,1]
	matrix beta_n_duc_duy[`n',2] = foo[1,2]
}
// n=8 and 9 are missing
foreach n of numlist 10/13 {
	reg duc_`n'  duy_`n' 
	matrix foo = e(b)
	matrix beta_n_duc_duy[`n',1] = foo[1,1]
	matrix beta_n_duc_duy[`n',2] = foo[1,2]
}

***************************************************************************
* Experiment with regressing consumption on income in CEX data
***************************************************************************
set mem 48m
set matsize 800
set more off 


u cexall,clear

reg ndur income_at
reg ndur income_at if ndur<60000 & income_at<60000 
reg ndur income 

reg ndurplus income_at
reg ndurplus income_at if ndur<60000 & income_at<60000 
reg ndurplus income 

gen log_ndurplus = log(ndurplus)
gen log_ndur = log(ndur)
gen log_food = log(food)
gen log_income = log(income)
gen log_income_at = log(income_at)
reg ndurplus ndur
reg log_ndurplus log_ndur

reg ndur food
reg log_ndur log_food

reg log_ndur log_income_at
reg log_food log_income_at

ivreg2 log_ndur (log_income = educ year )
ivreg2 log_ndur (lwh        = educ year )

ivreg2 lq (lx =lwhi ) 
ivreg2 lq (lx =lwhi lwwi ) 

#delimit;
ivreg2 log_income (lx lxed2 lxed3 lxkid2 lxkid3 lxkid4 lxy* =lwhi* lwwi*) age age2
                 lp lpalc lpfut lptra edd2 edd3 regd1 regd2 regd3 cohd1-cohd7 
		     kidd2-kidd4 ncomp white if complete==1;					
#delimit cr

egen decile_income = xtile(income_at), n(10)
bys decile_income: egen decile_inc = mean(income_at)
bys decile_income: egen decile_ndur = mean(ndur)
bys decile_income: egen decile_ndurplus = mean(ndurplus)
bys decile_income: egen decile_cons1 = mean(cons1)
bys decile_income: egen decile_cons2 = mean(cons2)
bys decile_income: egen decile_food = mean(food)
bys decile_income: egen decile_fout = mean(fout)
bys decile_income: egen decile_wageh = mean(wageh)
bys decile_income: egen decile_wagehw = mean(wageh+wagew)

scatter decile_ndur decile_inc if rand<0.02 & age>30 & age<60
reg decile_cons2 decile_inc if rand<0.02 & age>30 & age<60 & decile_income > 1.5

scatter decile_ndurplus decile_inc if rand<0.02 & age>30 & age<60
scatter decile_cons1 decile_inc if rand<0.02 & age>30 & age<60
scatter decile_cons2 decile_inc if rand<0.02 & age>30 & age<60
scatter decile_food decile_inc if rand<0.02 & age>30 & age<60
scatter decile_fout decile_inc if rand<0.02 & age>30 & age<60
scatter decile_wageh decile_inc if rand<0.02 & age>30 & age<60
scatter decile_wagehw decile_inc if rand<0.02 & age>30 & age<60
scatter decile_ndur decile_inc if rand<0.02 & age>30 & age<60

reg decile_wagehw decile_inc if rand<0.02 & age>30 & age<60

gen foo1 = log(decile_ndur)
gen foo2 = log(decile_inc)
scatter foo1 foo2 if rand<0.02 & age>30 & age<60
reg foo1 foo2 if rand<0.02 & age>30 & age<60 & decile_income > 1.5
