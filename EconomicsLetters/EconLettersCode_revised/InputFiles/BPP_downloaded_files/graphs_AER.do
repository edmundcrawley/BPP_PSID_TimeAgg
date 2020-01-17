cd C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files
clear
set mem 24m
cap log close
set more off
log using graphs_variance,replace

u data3,clear

******Does not use osbervations with topcoded income or financial income or federal taxes paid******
replace asset=. if trunca==1
replace y    =. if truncy==1
replace ftax =. if ftax==99999 & year>=1982 & year<=1985

							/*Y=LY+TR+FI (LY=lab inc, TR=Transfers, FI=Financ inc) */
gen ratio_ya=(y-asset)/y			/*Ratio of non-financial income to income, i.e. (LY+TR)/Y */
replace ratio_ya=0 if ratio_ya<0

keep if seo==0
drop if age<30|age>65

gen nas_scale=sqrt(ncomp)
replace lc=lc-ln(nas_scale)
							/*Two definitions of after-tax income*/
							/*Def #1: Total income - Taxes on total income */
							/*Def #2: Non-financial income - Taxes on non-financial income */
gen     ly_alt=ln(y-ftax)-ln(nas_scale)				
replace ly    =ln(y-asset-ratio_ya*ftax)-ln(nas_scale)

gen lf=ln(food+fout)-ln(nas_scale)
collapse (mean) m_psd=lc m_yp=ly f_p=lf (sd) v_psd=lc v_yp=ly v_yp_alt=ly_alt v_fp=lf,by(year)

replace v_yp=v_yp^2
replace v_yp_alt=v_yp_alt^2
replace v_f=v_f^2
replace v_f=. if v_f==0
replace m_psd=. if m_psd==0
replace v_psd=. if v_psd==0

replace v_psd=v_psd^2

sort year

keep m_* v_* year f_p
sort year
save mean_psid,replace

u cexall,clear

gen     coh=1  if yb>=1950 & yb<=1959
replace coh=2  if yb>=1940 & yb<=1949
replace coh=3  if yb>=1930 & yb<=1939
replace coh=4  if yb>=1920 & yb<=1929
tab coh,gen(cohd)
drop if coh==.
drop if age <30|age >65

gen nas_scale=sqrt(ncomp)
gen lx=ln(ndur)-ln(nas_scale)
gen ly  =ln(income)				/*This is not used*/

gen lf=ln(food+fout)-ln(nas_scale)
collapse (mean) m_cex=lx m_yc=ly f_c=lf (sd) v_cex=lx v_yc=ly v_fc=lf,by(year)

replace v_yc=v_yc^2
replace v_fc=v_fc^2
replace m_cex=. if m_cex==0
replace v_cex=. if v_cex==0

replace v_cex=v_cex^2

keep m_* v_* f_c year
sort year
save mean_cex,replace

u mean_cex,clear
merge year using mean_psid
tab _merge
drop _merge
save mean,replace

erase mean_cex.dta
erase mean_psid.dta

sort year
keep if year>=1980 & year<=1992

gen m_psd1=m_psd if year<=1986
gen m_psd2=m_psd if year>=1989
gen f_p1=f_p if year<=1986
gen f_p2=f_p if year>=1989
lab var m_psd1 "Mean of log(C) PSID, 80-86"
lab var m_psd2 "Mean of log(C) PSID, 89-92"
lab var f_p1 "Mean of log(f) PSID, 80-86"
lab var f_p2 "Mean of log(f) PSID, 89-92"
lab var m_cex "Mean of log(C) CEX"
lab var f_c "Mean of log(f) CEX"
#delimit;
gr7 m_psd1 m_psd2 m_cex year,c(lll[-]) sy(opi) xlabel(1980 1982 to 1992) ylabel(8.6 8.8 to 9.4) l1("log(C)") 
				saving(C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\mean_wrong,replace) title(Panel A);
gr7 f_p1 f_p2 f_c year,c(lll[-]) sy(opi) xlabel(1980 1982 to 1992) ylabel(7.5 7.7 to 8.3) l1("log(food)") 
				saving(C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\food,replace) title(Panel B);
#delimit cr

merge using b0
drop _merge
tab year,gen(yrd)
#delimit;
gen gamma=b0_lx[1]+b0_lxy1[1]*yrd1+b0_lxy2[1]*yrd2+b0_lxy3[1]*yrd3+b0_lxy4[1]*yrd4+b0_lxy5[1]*yrd5+b0_lxy6[1]*yrd6+b0_lxy7[1]*yrd7
+b0_lxy8[1]*yrd8+b0_lxy9[1]*yrd9+b0_lxy10[1]*yrd10+b0_lxy11[1]*yrd11+b0_lxy12[1]*yrd12;
#delimit cr
drop b0* yrd*

gen diff_f=(1/gamma)*(f_p-f_c)
gen diff_c=m_psd-m_cex
lab var diff_c "log(C,PSID)-log(C,CEX)"
lab var diff_f "(log(f,PSID)-log(f,CEX))/beta"
gen m_psd_c=m_psd-diff_f

gen m_psd_c1=m_psd_c if year<=1986
gen m_psd_c2=m_psd_c if year>=1989
lab var m_psd_c1 "Mean of log(C) PSID 80-86, corr."
lab var m_psd_c2 "Mean of log(C) PSID 89-92, corr."

#delimit;
gr7 diff_f diff_c year,c(ll) sy(oi) xlabel(1980 1982 to 1992) ylabel(-.05 0 to .15) l1(" ") 
				saving(C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\diff,replace) title(Panel C);
gr7 m_psd_c1 m_psd_c2 m_cex year,	c(lll[-]) s(opi) xla(1980 1982 to 1992) yla(8.6 8.8 to 9.4) l1("Corrected log(C)")
				saving(C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\mean_right,replace) title(Panel D);
gr7 using	C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\mean_wrong.gph C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\food.gph
		C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\diff.gph C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\mean_right.gph, 
		saving(C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\mean,replace);
#delimit cr

****************The graph above is in the Annex ********************

lab var v_psd "Var. of log(C) PSID"
lab var v_cex "Var. of log(C) CEX"

*****FIGURE 3 **********
#delimit;
twoway (scatter v_psd year, c(l) clp(-) yaxis(1) msymbol(oh)) (scatter v_cex year,c(l) yaxis(2) msymbol(o)) if year>=1980, 
	ylabel(.18(.02).30, axis(1) nogrid) ylabel(.11(.02).23,axis(2) nogrid) xlabel(1980(2)1992) 
	ytitle("PSID", axis(1)) ytitle("CEX",axis(2)) xtitle(Year) graphregion(fcolor(white))	
	saving(C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\var_right_sm,replace);
#delimit cr

***Income inequality and consumption inequality

gen svy=v_yp_alt if year==1980
gen svc=v_cex if year==1980
egen tempy=mean(svy)
egen tempc=mean(svc)
gen normv_cex=v_cex+tempy-tempc
drop svy-tempc

lab var v_yp_alt "Var. of log(Y) PSID"
lab var normv_cex "Var. of log(C) CEX, normalized"
gr7 normv_cex v_yp year,	c(ll[-]) s(oi) xla(1980 1982 to 1992) yla(.24(.02).38) l1(" ")
gr7 v_cex v_yp year    ,	c(ll[-]) s(oi) xla(1980 1982 to 1992) yla(.12(.02).38) l1(" ")
gr7 v_yp_alt v_cex year    ,	c(ll[-]) s(oi) xla(1980 1982 to 1992) rla(.13(.02).25) yla(.25(.02).37) rescale

smooth 4253eh,twice v_yp_alt,gen(vyp_sm)
smooth 4253eh,twice v_cex,gen(vcc_sm)

lab var vyp_sm "Var. of log(Y) PSID, smoothed"
lab var vcc_sm "Var. of log(C) CEX, smoothed"

*****FIGURE 1 *****
#delimit;
twoway (scatter vyp_sm v_yp_alt year, c(l) yaxis(1) msymbol(i o)) (scatter vcc_sm v_cex year,c(l) clp(-) yaxis(2) msymbol(i sh)) if year>=1980, 
	ylabel(.26(.02).36, axis(1) nogrid) ylabel(.145(.02).245,axis(2) nogrid) xlabel(1980(2)1992) ytitle("Var(log(Y)) PSID", axis(1)) 
	ytitle("Var(log(C)) CEX",axis(2)) xtitle(Year) graphregion(fcolor(white));
#delimit cr


***** This produces Figure 4        ********
***** It uses MD estimation results ********
u C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\cy_var_vs_y_var.dta,clear
lab var vp_cy "Using cons. and income data"
lab var vp_y "Using only income data"
#delimit;
scatter vp_cy vp_y year,c(l l) ylabel(0(.01).04, nogrid) xla(1980(5)1990) ytitle(" ") 
	graphregion(fcolor(white)) clp(solid dash);
#delimit cr
clear

*****This produces Figure 6
***** It uses MD estimation results ********
u C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\trans_var_using_e_y,clear
smooth 4253eh,twice earnings_at,gen(e)
smooth 4253eh,twice income  ,gen(y)
lab var e "Var(log(male earnings))"
lab var y "Var(log(net family income))"
#delimit;
twoway (scatter e year, c(l) yaxis(1) msymbol(i o)) (scatter y year,c(l) clp(-) yaxis(2) msymbol(i sh)) if year>=1980, 
	ylabel(.06(.01).10, axis(1) nogrid) ylabel(.03(.01).07,axis(2) nogrid) xlabel(1980(2)1992) ytitle("Male earnings", axis(1)) 
	ytitle("Net family income",axis(2)) xtitle(Year) graphregion(fcolor(white));
#delimit cr

u C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\fig6.dta,clear
lab var zt_cy2 "Baseline sample"
lab var zt_no_college2 "No college education"
lab var zt_college2 "College education"
lab var zt_cy_seo2 "Adding SEO"
lab var zt_total_cons2 "Total exp."

*****This produces Figure 7
***** It uses MD estimation results ********
#delimit;
scatter zt_cy2 year if year<1990, ylabel(0(.01).04,nogrid) c(l) graphregion(fcolor(white)) saving(temp1,replace);
scatter zt_no_college2 zt_college2 year if year<1990, ylabel(0(.01).04,nogrid) c(l l) clp(-) graphregion(fcolor(white)) 
	ytitle(" ")saving(temp2,replace);
scatter zt_cy2 zt_cy_seo2 year if year<1990, ylabel(0(.01).04,nogrid) c(l l) clp(-) graphregion(fcolor(white)) 
	ytitle(" ")saving(temp3,replace);
scatter zt_cy2 zt_total_cons2 year if year<1990, ylabel(0(.01).04,nogrid) c(l l) clp(-) graphregion(fcolor(white)) 
	ytitle(" ")saving(temp4,replace);
#delimit cr

graph combine temp1.gph temp2.gph temp3.gph temp4.gph, iscale(.70) graphregion(fcolor(white))

*****This produces Figure 5
***** It uses MD estimation results ********
u C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\fit_model
#delimit;
twoway (scatter act_c pred_c year, c(l l) clp(solid dash) yaxis(1) msymbol(i i)) if year>=1980, 
	ylabel(.1(.04).22, axis(1) nogrid) xlabel(1980(2)1992) title(Var(Dlog(c))) legend(off)
	xtitle(Year) graphregion(fcolor(white)) saving(temp1,replace);
twoway (scatter act_y pred_y year, c(l l) clp(solid dash) yaxis(1) msymbol(i i)) if year>=1980, 
	ylabel(.06(.02).12, axis(1) nogrid) xlabel(1980(2)1992) title(Var(Dlog(y))) legend(off)
	xtitle(Year) graphregion(fcolor(white)) saving(temp2,replace);
twoway (scatter act_cy pred_cy year, c(l l) clp(solid dash) yaxis(1) msymbol(i i)) if year>=1980, 
	ylabel(0(.01).03, axis(1) nogrid) xlabel(1980(2)1992) title(Cov(Dlog(c),Dlog(y))) legend(off)
	xtitle(Year) graphregion(fcolor(white)) saving(temp3,replace);
#delimit cr

graph combine temp1.gph temp2.gph temp3.gph, graphregion(fcolor(white))


*****This produces Figure 2
cd C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files

u data3,clear

keep if seo==0
drop if age<30|age>65

replace ly=ln(y)

gen nas_scale=sqrt(ncomp)
replace lc=lc-ln(nas_scale)
replace ly=ln(y-ftax)-ln(nas_scale)

gen lf=ln(food+fout)-ln(nas_scale)
collapse (mean) age m_psd=lc m_yp=ly f_p=lf (sd) v_psd=lc v_yp=ly v_fp=lf,by(year coh)

replace v_y=v_y^2
replace v_f=v_f^2
replace v_f=. if v_f==0
replace m_psd=. if m_psd==0
replace v_psd=. if v_psd==0

replace v_psd=v_psd^2

lab var year     "Year"

keep m_* v_* year f_p coh age
sort coh year
save mean_psid_coh,replace

u cexall,clear

gen     coh=1  if yb>=1950 & yb<=1959
replace coh=2  if yb>=1940 & yb<=1949
replace coh=3  if yb>=1930 & yb<=1939
replace coh=4  if yb>=1920 & yb<=1929
tab coh,gen(cohd)
drop if coh==.

drop if age <30|age >65

gen nas_scale=sqrt(ncomp)
gen lx=ln(ndur)-ln(nas_scale)
gen ly  =ln(income)

gen lf=ln(food+fout)-ln(nas_scale)
collapse (mean) age m_cex=lx m_yc=ly f_c=lf (sd) v_cex=lx v_yc=ly v_fc=lf,by(year coh)

replace v_yc=v_yc^2
replace v_fc=v_fc^2
replace m_cex=. if m_cex==0
replace v_cex=. if v_cex==0

replace v_cex=v_cex^2

keep m_* v_* year f_c coh age
sort coh year
save mean_cex_coh,replace

u mean_cex_coh,clear
merge coh year using mean_psid_coh
tab _merge
drop _merge
save mean_coh,replace

erase mean_cex_coh.dta
erase mean_psid_coh.dta

sort coh year

gen vc1=v_cex if coh==1
gen vc2=v_cex if coh==2
gen vc3=v_cex if coh==3
gen vc4=v_cex if coh==4

gen vy1=v_yp if coh==1
gen vy2=v_yp if coh==2
gen vy3=v_yp if coh==3
gen vy4=v_yp if coh==4

smooth 4253eh,twice vc1    ,gen(vc1sm)
smooth 4253eh,twice vc2    ,gen(vc2sm)
smooth 4253eh,twice vc3    ,gen(vc3sm)
smooth 4253eh,twice vc4    ,gen(vc4sm)

smooth 4253eh,twice vy1    ,gen(vy1sm)
smooth 4253eh,twice vy2    ,gen(vy2sm)
smooth 4253eh,twice vy3    ,gen(vy3sm)
smooth 4253eh,twice vy4    ,gen(vy4sm)

lab var age "Age"
lab var vc1sm "Born 1950s"
lab var vc2sm "Born 1940s"
lab var vc3sm "Born 1930s"
lab var vc4sm "Born 1920s"
lab var vy1sm "Born 1950s"
lab var vy2sm "Born 1940s"
lab var vy3sm "Born 1930s"
lab var vy4sm "Born 1920s"

sort coh year

#delimit;
gr7 vc1sm vc2sm vc3sm vc4sm age if age>=30 & age<=65, xla(30(5)65) yla(.1(.05).25) c(llll) 
	title("The variance of log consumption over the life-cycle, CEX");
gr7 vy1sm vy2sm vy3sm vy4sm age if age>=30 & age<=65, xla(30(5)65) yla(.15(.05).5) c(llll) 
	title("The variance of log income over the life-cycle, PSID");
#delimit cr

#delimit;
scatter vc1sm vc2sm vc3sm vc4sm age if age>=30 & age<=65, xla(30(5)65) yla(.1(.05).25,nogrid) c(l l l l) graphregion(fcolor(white)) s(oh s x i) ytitle("");
#delimit cr


log close
clear
