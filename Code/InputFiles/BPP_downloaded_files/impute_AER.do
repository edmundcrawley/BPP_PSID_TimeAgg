cd C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files

/* This program runs the structural demand regression using
the CEX and then imputes the total non-durable consumption
in the PSID data */

clear
capture log close
log using C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files,replace t

set mem 48m
set matsize 800
set more off 


u cexall,clear

gen     educ=1 if educh==7| educh ==1 | educh ==2 
replace educ=2 if educh ==3
replace educ=3 if educh ==4| educh ==5 | educh ==6

gen    white=race =="1" 				/* dummy for white husband */
tab educ  ,gen(edd)
tab region,gen(regd)

gen     coh=1  if yb>=1955 & yb<=1959
replace coh=2  if yb>=1950 & yb<=1954
replace coh=3  if yb>=1945 & yb<=1949
replace coh=4  if yb>=1940 & yb<=1944
replace coh=5  if yb>=1935 & yb<=1939
replace coh=6  if yb>=1930 & yb<=1934
replace coh=7  if yb>=1925 & yb<=1929
replace coh=8  if yb>=1920 & yb<=1924
tab coh,gen(cohd)
drop if coh==.

replace pf  =pf/100

gen lq  =ln((food+fout)/pf)
gen lp  =ln(pf)
gen lx  =ln(ndur)
gen lxs =ln(ndurplus)
gen lX  =ln(cons2)
gen lXs =ln(cons2+(ndurplus-ndur))

replace weekh =0 if weekh ==.
replace hourh =0 if hourh ==.
gen yhourh=weekh*hourh
gen parth=yhourh>0

replace weekw =0 if weekw ==.
replace hourw =0 if hourw ==.
gen yhourw=weekw*hourw
gen partw=yhourw>0

gen lxed2=lx*edd2
gen lxed3=lx*edd3
gen lxsed2=lxs*edd2
gen lxsed3=lxs*edd3
gen lXed2=lX*edd2
gen lXed3=lX*edd3
gen lXsed2=lXs*edd2
gen lXsed3=lXs*edd3

gen kidcat=recode(kid,0,1,2,4)
tab kidcat,gen(kidd)

gen lxkid2=lx*kidd2
gen lxkid3=lx*kidd3
gen lxkid4=lx*kidd4
gen lxskid2=lxs*kidd2
gen lxskid3=lxs*kidd3
gen lxskid4=lxs*kidd4
gen lXkid2=lX*kidd2
gen lXkid3=lX*kidd3
gen lXkid4=lX*kidd4
gen lXskid2=lXs*kidd2
gen lXskid3=lXs*kidd3
gen lXskid4=lXs*kidd4

tab year,gen(yrd)

cap program drop doit
program def doit
        local i=1
        while `i' < 13 {
		gen lxy`i'=lx*yrd`i'
		gen lxsy`i'=lxs*yrd`i'
		gen lXy`i'=lX*yrd`i'
		gen lXsy`i'=lXs*yrd`i'
        local i=`i'+1               
        }
end
qui doit

sort year
merge year using finprice
keep if _merge==3
drop _merge

gen lpalc=ln(palc/100) 
gen lpfut=ln(pfutil/100)
gen lptra=ln(ptran/100)

gen lwh=ln(wageh/(hourh*weekh))
gen lww=ln(wagew/(hourw*weekw))

gen age2=age^2

egen lwhi=mean(lwh),by(coh educ year)
egen lwwi=mean(lww),by(coh educ year)

//egen lwhi=mean(ln(income_at)),by(coh educ year)
//egen lwwi=mean(ln(income_at)),by(coh educ year)

//gen lwhi=lwh
//gen lwwi=lww

//gen lwhi=income_at
//gen lwwi=income_at

cap program drop doit
program def doit
        local i=1
        while `i' < 13 {
		gen lwhiy`i'=lwhi*yrd`i'
		gen lwwiy`i'=lwwi*yrd`i'
        local i=`i'+1               
        }
end
qui doit

gen lwhik2 =lwhi*kidd2
gen lwhik3 =lwhi*kidd3
gen lwhik4 =lwhi*kidd4
gen lwhied2=lwhi*edd2
gen lwhied3=lwhi*edd3
gen lwwik2 =lwwi*kidd2
gen lwwik3 =lwwi*kidd3
gen lwwik4 =lwwi*kidd4
gen lwwied2=lwwi*edd2
gen lwwied3=lwwi*edd3


***DEMAND STRATEGY IV: non durables
/*Note: All regressions below are estimated only among complete income respondents*/

#delimit;
ivreg2 lq (lx lxed2 lxed3 lxkid2 lxkid3 lxkid4 lxy* =lwhi* lwwi*) age age2
                 lp lpalc lpfut lptra edd2 edd3 regd1 regd2 regd3 cohd1-cohd7 
		     kidd2-kidd4 ncomp white if complete==1;					
#delimit cr
testparm lxy*

gen b0_lx    =_b[lx]                     

gen b0_lxed2  =_b[lxed2]
gen b0_lxed3  =_b[lxed3]

gen b0_lxkid2 =_b[lxkid2]
gen b0_lxkid3 =_b[lxkid3]
gen b0_lxkid4 =_b[lxkid4]

cap program drop doit
program def doit
        local i=1
        while `i' < 13 {
		gen b0_lxy`i'=_b[lxy`i']
        local i=`i'+1               
        }
end
qui doit

gen b0_age    =_b[age]
gen b0_age2   =_b[age2]

gen b0_lp     =_b[lp]
gen b0_lpalc  =_b[lpalc]
gen b0_lpfut  =_b[lpfut]
gen b0_lptra  =_b[lptra]

gen b0_edd2   =_b[edd2]
gen b0_edd3   =_b[edd3] 
gen b0_regd1  =_b[regd1] 
gen b0_regd2  =_b[regd2] 
gen b0_regd3  =_b[regd3] 
gen b0_cohd1  =_b[cohd1] 
gen b0_cohd2  =_b[cohd2] 
gen b0_cohd3  =_b[cohd3] 
gen b0_cohd4  =_b[cohd4] 
gen b0_cohd5  =_b[cohd5] 
gen b0_cohd6  =_b[cohd6] 
gen b0_cohd7  =_b[cohd7] 
gen b0_kid2   =_b[kidd2]
gen b0_kid3   =_b[kidd3]
gen b0_kid4   =_b[kidd4]

gen b0_ncomp  =_b[ncomp] 
gen b0_dumwh  =_b[white]

gen b0_cons   =_b[_cons]

***DEMAND STRATEGY IV: total expenditure
#delimit;
ivreg2 lq (lX lXed2 lXed3 lXkid2 lXkid3 lXkid4 lXy*=lwhi* lwwi* ) age age2
                 lp lpalc lpfut lptra edd2 edd3 regd1 regd2 regd3 cohd1-cohd7
		     kidd2-kidd4 ncomp white if complete==1;
#delimit cr

gen b1_lx    =_b[lX]                     

gen b1_lxed2 =_b[lXed2]
gen b1_lxed3 =_b[lXed3]

gen b1_lxkid2=_b[lXkid2]
gen b1_lxkid3=_b[lXkid3]
gen b1_lxkid4=_b[lXkid4]

cap program drop doit
program def doit
        local i=1
        while `i' < 13 {
		gen b1_lxy`i'=_b[lXy`i']
        local i=`i'+1               
        }
end
qui doit

gen b1_age  =_b[age]
gen b1_age2 =_b[age2]

gen b1_lp    =_b[lp]
gen b1_lpalc =_b[lpalc]
gen b1_lpfut =_b[lpfut]
gen b1_lptra =_b[lptra]

gen b1_edd2  =_b[edd2]
gen b1_edd3  =_b[edd3] 
gen b1_regd1 =_b[regd1] 
gen b1_regd2 =_b[regd2] 
gen b1_regd3 =_b[regd3] 
gen b1_cohd1 =_b[cohd1] 
gen b1_cohd2 =_b[cohd2] 
gen b1_cohd3 =_b[cohd3] 
gen b1_cohd4 =_b[cohd4] 
gen b1_cohd5 =_b[cohd5] 
gen b1_cohd6 =_b[cohd6] 
gen b1_cohd7 =_b[cohd7] 
gen b1_kid2   =_b[kidd2]
gen b1_kid3   =_b[kidd3]
gen b1_kid4   =_b[kidd4]

gen b1_ncomp =_b[ncomp] 
gen b1_dumwh =_b[white]

gen b1_cons  =_b[_cons]


***DEMAND STRATEGY IV: non durables+services from durables
#delimit;
ivreg2 lq (lxs lxsed2 lxsed3 lxskid2 lxskid3 lxskid4 lxsy* =lwhi* lwwi*) age age2
                 lp lpalc lpfut lptra edd2 edd3 regd1 regd2 regd3 cohd1-cohd7 
		     kidd2-kidd4 ncomp white if complete==1;
#delimit cr

gen b2_lx    =_b[lxs]                     

gen b2_lxed2  =_b[lxsed2]
gen b2_lxed3  =_b[lxsed3]

gen b2_lxkid2 =_b[lxskid2]
gen b2_lxkid3 =_b[lxskid3]
gen b2_lxkid4 =_b[lxskid4]

cap program drop doit
program def doit
        local i=1
        while `i' < 13 {
		gen b2_lxy`i'=_b[lxsy`i']
        local i=`i'+1               
        }
end
qui doit

gen b2_age    =_b[age]
gen b2_age2   =_b[age2]

gen b2_lp     =_b[lp]
gen b2_lpalc  =_b[lpalc]
gen b2_lpfut  =_b[lpfut]
gen b2_lptra  =_b[lptra]

gen b2_edd2   =_b[edd2]
gen b2_edd3   =_b[edd3] 
gen b2_regd1  =_b[regd1] 
gen b2_regd2  =_b[regd2] 
gen b2_regd3  =_b[regd3] 
gen b2_cohd1  =_b[cohd1] 
gen b2_cohd2  =_b[cohd2] 
gen b2_cohd3  =_b[cohd3] 
gen b2_cohd4  =_b[cohd4] 
gen b2_cohd5  =_b[cohd5] 
gen b2_cohd6  =_b[cohd6] 
gen b2_cohd7  =_b[cohd7] 
gen b2_kid2   =_b[kidd2]
gen b2_kid3   =_b[kidd3]
gen b2_kid4   =_b[kidd4]

gen b2_ncomp  =_b[ncomp] 
gen b2_dumwh  =_b[white]

gen b2_cons   =_b[_cons]

***DEMAND STRATEGY IV: tot exp+services from durables
#delimit;
ivreg2 lq (lXs lXsed2 lXsed3 lXskid2 lXskid3 lXskid4 lXsy* =lwhi* lwwi*) age age2
                 lp lpalc lpfut lptra edd2 edd3 regd1 regd2 regd3 cohd1-cohd7 
		     kidd2-kidd4 ncomp white if complete==1;
#delimit cr

gen b3_lx    =_b[lXs]                     

gen b3_lxed2  =_b[lXsed2]
gen b3_lxed3  =_b[lXsed3]

gen b3_lxkid2 =_b[lXskid2]
gen b3_lxkid3 =_b[lXskid3]
gen b3_lxkid4 =_b[lXskid4]

cap program drop doit
program def doit
        local i=1
        while `i' < 13 {
		gen b3_lxy`i'=_b[lXsy`i']
        local i=`i'+1               
        }
end
qui doit

gen b3_age    =_b[age]
gen b3_age2   =_b[age2]

gen b3_lp     =_b[lp]
gen b3_lpalc  =_b[lpalc]
gen b3_lpfut  =_b[lpfut]
gen b3_lptra  =_b[lptra]

gen b3_edd2   =_b[edd2]
gen b3_edd3   =_b[edd3] 
gen b3_regd1  =_b[regd1] 
gen b3_regd2  =_b[regd2] 
gen b3_regd3  =_b[regd3] 
gen b3_cohd1  =_b[cohd1] 
gen b3_cohd2  =_b[cohd2] 
gen b3_cohd3  =_b[cohd3] 
gen b3_cohd4  =_b[cohd4] 
gen b3_cohd5  =_b[cohd5] 
gen b3_cohd6  =_b[cohd6] 
gen b3_cohd7  =_b[cohd7] 
gen b3_kid2   =_b[kidd2]
gen b3_kid3   =_b[kidd3]
gen b3_kid4   =_b[kidd4]

gen b3_ncomp  =_b[ncomp] 
gen b3_dumwh  =_b[white]

gen b3_cons   =_b[_cons]

cap program drop doit
program def doit
        local i=1
        while `i' < 13 {
		gen lqy`i'=lq*yrd`i'
        local i=`i'+1               
        }
end
qui doit

keep in 1
collapse b0* b1* b2* b3*
save b0,replace

clear

************************************NOW IMPUTE IN PSID *******************************
set more off

u data2,clear

sort id year
merge id year using tax9192			/*TAXSIM data for 1991-92*/
tab _merge
drop if _merge==2
drop _merge

replace ftax=fiitax if year>1990		/*Tax data for 91-92 obtained from TAXSIM*/
drop fiitax

*****Impute for those few observations that don't give federal tax from taxsim****
gen y2=y^2
gen y3=y^3
gen y4=y^4
gen y5=y^5
ren fsize ncomp
gen ncomp2=ncomp^2
gen asset2=asset^2
gen post86=year>1986
reg ftax y y2 y3 y4 y5 ncomp ncomp2 post86 kids self asset if year<1991
predict ftaxp
drop y2-asset2 post86
replace ftax=ftaxp if year>=1991 & ftax==.

sort yb

/* Generate the cohort groups */
gen     coh=1  if yb>=1955 & yb<=1959
replace coh=2  if yb>=1950 & yb<=1954
replace coh=3  if yb>=1945 & yb<=1949
replace coh=4  if yb>=1940 & yb<=1944
replace coh=5  if yb>=1935 & yb<=1939
replace coh=6  if yb>=1930 & yb<=1934
replace coh=7  if yb>=1925 & yb<=1929
replace coh=8  if yb>=1920 & yb<=1924
tab coh,gen(cohd)
drop if coh==.

tab region, gen(regd)
tab educ, gen(edd)
gen white=race==1
gen age2=age^2
replace hours=0 if hours==.
gen parth=hours>0 & hours!=.
gen partw=hourw>0 & hourw!=.

/* Now I merge the PSID data with the price data*/

sort year
merge year using natpr
tab _merge
drop if _merge!=3
drop _merge

gen pf   =pfood/100
gen lq  =ln((food+fout)/pf)
gen lp  =ln(pf)

gen kidcat=recode(kids,0,1,2,4)
tab kidcat,gen(kidd)	

sort year
merge year using finprice
keep if _merge==3
drop _merge

gen lpalc=ln(palc/100) 
gen lpfut=ln(pfutil/100)
gen lptra=ln(ptran/100)

merge using b0
drop _merge

tab year,gen(yrd)

cap program drop doit
program def doit
        local i=0
        while `i' < 4 {
		#delimit;
		gen lc`i'=(lq-(
		b`i'_age[1]    *age+
		b`i'_age2[1]   *age2+
		
		b`i'_lp[1]     *lp+
		b`i'_lpalc[1]  *lpalc+
		b`i'_lpfut[1]  *lpfut+
		b`i'_lptra[1]  *lptra+
		
		b`i'_cohd1[1]  *cohd1+
		b`i'_cohd2[1]  *cohd2+
		b`i'_cohd3[1]  *cohd3+
		b`i'_cohd4[1]  *cohd4+
		b`i'_cohd5[1]  *cohd5+
		b`i'_cohd6[1]  *cohd6+
		b`i'_cohd7[1]  *cohd7+
		b`i'_edd2[1]   *edd2+
		b`i'_edd3[1]   *edd3+ 
		b`i'_regd1[1]  *regd1+ 
		b`i'_regd2[1]  *regd2+ 
		b`i'_regd3[1]  *regd3+ 
		b`i'_kid2[1]   *kidd2+
		b`i'_kid3[1]   *kidd3+
		b`i'_kid4[1]   *kidd4+
		
		b`i'_ncomp[1]  *ncomp+ 
		b`i'_dumwh[1]  *white+
		
		b`i'_cons[1]  ))
		/(b`i'_lx[1]+b`i'_lxed2[1]*edd2+b`i'_lxed3[1]*edd3+b`i'_lxkid2[1]*kidd2+b`i'_lxkid3[1]*kidd3+b`i'_lxkid4[1]*kidd4
		+b`i'_lxy1[1]*yrd3+b`i'_lxy2[1]*yrd4+b`i'_lxy3[1]*yrd5+b`i'_lxy4[1]*yrd6+b`i'_lxy5[1]*yrd7+b`i'_lxy6[1]*yrd8+b`i'_lxy7[1]*yrd9
		+b`i'_lxy8[1]*yrd10+b`i'_lxy9[1]*yrd11+b`i'_lxy10[1]*yrd12+b`i'_lxy11[1]*yrd13+b`i'_lxy12[1]*yrd14);
		#delimit cr
        local i=`i'+1               
        }
end
qui doit

ren lc0 lc
ren lc1 lC
ren lc2 lcs
ren lc3 lCs

drop b0* b1* b2* b3* yrd*
drop coh cohd* kidcat-lptra

compress

gen     coh=1  if yb>=1950 & yb<=1959
replace coh=2  if yb>=1940 & yb<=1949
replace coh=3  if yb>=1930 & yb<=1939
replace coh=4  if yb>=1920 & yb<=1929

save data3,replace

log close
clear

