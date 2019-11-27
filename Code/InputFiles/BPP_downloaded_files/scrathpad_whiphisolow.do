
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

reg lc logy
//get 0.43

reg uc uy
//get 0.33

reg duc duy
//get 0.11

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
