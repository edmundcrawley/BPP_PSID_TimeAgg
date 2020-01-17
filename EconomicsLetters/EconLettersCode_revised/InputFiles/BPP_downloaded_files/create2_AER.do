clear
cap log close
set mem 200m
set more off

cd C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files
cap log close
log using C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\create2_AER,replace t

*This is the individual data set
*keep first family number, then sequence number, then relationship to head
*Note: there's no relationship to head in the first year of data


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\ind6817er,clear
rename *, lower
rename er* v*
#delimit;
keep  v30001 v30020 v30043 v30067 v30091 v30117 v30138 v30160 v30188 v30217 v30246 v30283 v30313 v30343 
      v30373 v30399 v30429 v30463 v30498 v30535 v30570 v30606 v30642 v30689 v30733 v30806 v33101 v33201
	v33301 v33401 		/* Interview number */
             v30021 v30044 v30068 v30092 v30118 v30139 v30161 v30189 v30218 v30247 v30284 v30314 v30344 
      v30374 v30400 v30430 v30464 v30499 v30536 v30571 v30607 v30643 v30690 v30734 v30807 v33102 v33202
	v33302 v33402	 	/* Sequence number */
      	 v30003 v30022 v30045 v30069 v30093 v30119 v30140 v30162 v30190 v30219 v30248 v30285 v30315 
	v30345 v30375 v30401 v30431 v30465 v30500 v30537 v30572 v30608 v30644 v30691 v30735 v30808 v33103
	v33203 v33303 v33403;	 /* Relationship to Head */

drop if v30001>=7001 & v30001<=9308; /*Drop latino households*/

gen seo=v30001>=5000 & v30001<=7000;	

drop if   v30003!=1  & v30022!=1  & v30045!=1  & v30069!=1  & v30093!=1 
        & v30119!=1  & v30140!=1  & v30162!=1  & v30190!=1  & v30219!=1 
        & v30248!=1  & v30285!=1  & v30315!=1  & v30345!=1  & v30375!=1 
        & v30401!=10 & v30431!=10 & v30465!=10 & v30500!=10 & v30537!=10 
        & v30572!=10 & v30608!=10 & v30644!=10 & v30691!=10 & v30735!=10
        & v30808!=10 & v33103!=10 & v33203!=10 & v33303!=10 & v33403!=10;	/*Drop those who are never heads*/

#delimit cr
ren v30001 id67
ren v30020 id68
ren v30043 id69
ren v30067 id70 
ren v30091 id71 
ren v30117 id72 
ren v30138 id73
ren v30160 id74
ren v30188 id75
ren v30217 id76
ren v30246 id77
ren v30283 id78
ren v30313 id79
ren v30343 id80
ren v30373 id81
ren v30399 id82
ren v30429 id83
ren v30463 id84
ren v30498 id85
ren v30535 id86
ren v30570 id87
ren v30606 id88
ren v30642 id89
ren v30689 id90
ren v30733 id91
ren v30806 id92
ren v33101 id93
ren v33201 id94
ren v33301 id95
ren v33401 id96

gen pid67=v30003==1
gen pid68=v30022==1  & v30021>=1 & v30021<=20
gen pid69=v30045==1  & v30044>=1 & v30044<=20
gen pid70=v30069==1  & v30068>=1 & v30068<=20
gen pid71=v30093==1  & v30092>=1 & v30092<=20
gen pid72=v30119==1  & v30118>=1 & v30118<=20 
gen pid73=v30140==1  & v30139>=1 & v30139<=20 
gen pid74=v30162==1  & v30161>=1 & v30161<=20 
gen pid75=v30190==1  & v30189>=1 & v30189<=20 
gen pid76=v30219==1  & v30218>=1 & v30218<=20 
gen pid77=v30248==1  & v30247>=1 & v30247<=20 
gen pid78=v30285==1  & v30284>=1 & v30284<=20 
gen pid79=v30315==1  & v30314>=1 & v30314<=20 
gen pid80=v30345==1  & v30344>=1 & v30344<=20 
gen pid81=v30375==1  & v30374>=1 & v30374<=20 
gen pid82=v30401==10 & v30400>=1 & v30400<=20
gen pid83=v30431==10 & v30430>=1 & v30430<=20
gen pid84=v30465==10 & v30464>=1 & v30464<=20
gen pid85=v30500==10 & v30499>=1 & v30499<=20
gen pid86=v30537==10 & v30536>=1 & v30536<=20
gen pid87=v30572==10 & v30571>=1 & v30571<=20
gen pid88=v30608==10 & v30607>=1 & v30607<=20
gen pid89=v30644==10 & v30643>=1 & v30643<=20
gen pid90=v30691==10 & v30690>=1 & v30690<=20
gen pid91=v30735==10 & v30734>=1 & v30734<=20
gen pid92=v30808==10 & v30807>=1 & v30807<=20
gen pid93=v33103==10 & v33102>=1 & v33102<=20
gen pid94=v33203==10 & v33202>=1 & v33202<=20
gen pid95=v33303==10 & v33302>=1 & v33302<=20
gen pid96=v33403==10 & v33402>=1 & v33402<=20

keep id* pid* seo

gen person=1
replace person=sum(person)
reshape long id pid,i(person seo) j(year)
replace year=1900+year

*drop people who are no heads in that year*
drop if pid==0

save       person,replace

clear


u fam,clear
sort id year
save fam, replace

u person,clear
compress
sort id year
merge id year using fam
tab _merge
drop if _merge!=3			/*_merge==1 is 1993-96, _merge==2 is latino subsample*/
drop _merge
save data_additional_vars, replace

// Use data file from original BPP, but add in additional vars
u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\dataBPP2008, clear
merge 1:1 id year using data_additional_vars
drop _merge
save data, replace

log close
clear
