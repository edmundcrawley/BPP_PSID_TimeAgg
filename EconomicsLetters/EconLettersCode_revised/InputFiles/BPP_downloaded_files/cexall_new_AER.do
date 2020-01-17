cd c:\papers\partial\revise2

clear
capture log close
set mem 200m
set more off
log using c:\papers\partial\revision\cexall_gorbacheva,replace t


u c:\papers\partial\revise2\mtabdata,clear

gen refmo=real(REF_MO)
gen year=real(REF_YR)
replace year=1900+year if year<1900
drop REF*

drop if year<1980			/*There are no data before 1980 in officially CEX released data*/

ren cons totnd
ren educ educa
gen other=durables+health+educa
keep year family intno refmo foodin foodout totnd educa health other
sort family
merge family using c:\cex\johnson						
tab _merge
drop _merge

gen tot1=totnd+educa+health	/*Non durable + exp on health and education*/
gen tot2=totnd+other		/*Total expenditure*/
drop educa health other
compress
drop if foodin==.|foodin==0
drop if totnd ==.|totnd ==0

/*Takes maximum consumption for those with duplicate ref.months (8 observations overall)*/

egen fin =max(foodin) ,by(family year refmo)	
egen fou =max(foodout),by(family year refmo)
egen nd  =max(totnd)  ,by(family year refmo)
egen co1 =max(tot1)   ,by(family year refmo)
egen co2 =max(tot2)   ,by(family year refmo)

sort family refmo
drop if family==family[_n-1] & year==year[_n-1] & refmo==refmo[_n-1]	/*drop those with duplicate ref.months*/
drop foodin foodout totnd tot1 tot2

/* Sum the consumption observations over the period a family is observed */

sort family year
egen food   =sum(fin) , by(family)
egen fout   =sum(fou) , by(family)
egen ndur   =sum(nd)  , by(family)
egen cons1  =sum(co1) , by(family)
egen cons2  =sum(co2) , by(family)

gen ndurplus=ndur+imp_hs+imp_vh			/*imp_hs and imp_vh are the variables contained in the file johnson.dta */

drop fin fou nd co1 co2
compress

/* Drop those with less than 12 months of records*/

egen n=sum(family!=.),by(family)
drop if n<12
drop n

/*Allocates to the year with the most observations*/

gen origyear=year
egen maxyear  =max(year), by(family)
egen minyear  =min(year), by(family)
gen        aux=year==maxyear
egen      maux=sum(aux), by(family)
sort family year refmo

replace year=minyear if maux<=6
replace year=maxyear if maux> 6

gen     o_year=year    if maux==12
replace o_year=maxyear if maux<=6
replace o_year=minyear if maux>6

drop if year>1992				/*There's no PSID data after 1992*/
collapse ndur ndurplus cons1 cons2 food fout o_year maux,by(family year)
sort family year
save cexfood,replace

u c:\papers\partial\revise2\fmlydata,clear
gen complete=RESPSTAT=="1"
gen region=real(REGION)
ren FAM_S ncomp
gen famtype=real(FAM_TY)
ren FINCB incbtax
ren FINCA incatax
ren PERS comp18
ren SEX_REF sex
gen educh=real(school)
gen educw=real(schoolw)
ren num_wk weekh 
ren num_wkw weekw 
ren wkhrs hourh
ren wkhrsw hourw
ren wages salaryh
ren wages salaryw

egen inc_ass=rsum(FININCX INTEARNX)

#delimit;
keep family intno region ncomp age agew race racew sex educh educw famtype incbtax incatax weekh weekw hourh hourw
     comp18 salaryh salaryw complete inc_ass;
#delimit cr

gen inc_netass=incatax-inc_ass

sort family
merge family using cexfood
tab _merge
drop if _merge!=3
drop _merge

egen t=sum(incbtax==0),by(family)
egen n=sum(family!=.),by(family)

drop if n==t					/*drop those who always report zero income*/
drop t n 

							/*assign income at the last interview, or the one before if this is zero	*/
							/*Everybody here has 5 interviews, so maxint=4 throughout				*/
							/*Here we are assign income at the 2nd interview which is carried through	*/
egen maxint=max(intno),by(family)
gen maxint1=maxint-1
sort family intno
qui by family:gen     inc=incbtax   if intno==maxint  & incbtax!=0
qui by family:replace inc=incbtax   if intno==maxint1 & inc[_n+1]==. & intno!=intno[_n+1]

qui by family:gen     inc2=incatax   if intno==maxint  & incatax!=0
qui by family:replace inc2=incatax   if intno==maxint1 & inc2[_n+1]==. & intno!=intno[_n+1]

qui by family:gen     inc3=inc_netass  if intno==maxint  & inc_netass!=0
qui by family:replace inc3=inc_netass  if intno==maxint1 & inc3[_n+1]==. & intno!=intno[_n+1]

gen salh=salaryh if intno==maxint		/*We use only the 5th interview - it's ok to have a 0 wage*/
gen salw=salaryw if intno==maxint

egen wageh=mean(salh),by(family)
egen wagew=mean(salw),by(family)

egen income    =max(inc),by(family)
egen income_at =max(inc2),by(family)
egen income_net=max(inc3),by(family)

drop maxint* inc inc2 inc3 salh salw

/*Drop those with missing values on region of residence*/
egen n=sum(family!=.),by(family)
egen r=sum(region==.),by(family)
drop if r==n

ren comp18 kid

/*drop those with some changes in family composition - we want to keep only married couples throughout to match selection in PSID*/
/*These are either singles, or lone parents, or "other non-nuclear HH types"*/
egen unstable=sum(famtype==6|famtype==7|famtype==8|famtype==9),by(family)
keep if unstable==0
drop unstable

egen maxint=max(intno),by(family)
keep if intno==maxint

#delimit;
keep region ncomp age* sex race famtype ndur ndurplus cons1 cons2 food fout income income_at income_net educh educw family year  
     weekh weekw hourh hourw kid maux o_year wageh wagew complete;
#delimit cr

sort year
merge year using c:\papers\partial\dta\natpr
tab _merge
drop if _merge!=3
drop _merge

gen     pfin =pfoodin                                      if year==o_year
replace pfin =((12-maux)/12)*pfoodin+(maux/12)*pfoodinlead if year <o_year
replace pfin =(maux/12)*pfoodin+((12-maux)/12)*pfoodinlag  if year >o_year
gen     pfaw =pfout                                        if year==o_year
replace pfaw =((12-maux)/12)*pfout+(maux/12)*pfoutlead     if year <o_year
replace pfaw =(maux/12)*pfout+((12-maux)/12)*pfoutlag      if year >o_year
gen     pf   =pfood                                        if year==o_year
replace pf   =((12-maux)/12)*pfood+(maux/12)*pfoodlead     if year <o_year
replace pf   =(maux/12)*pfood+((12-maux)/12)*pfoodlag      if year >o_year

drop maux o_year cpil* pfoodl* pfoodinl* pfoutl*
compress

gen yb=year-age
drop if yb<=1919|yb>=1960
drop if age<30|age>65

drop if income<food & food!=.

save cexall,replace
erase cexfood.dta

log close
clear
