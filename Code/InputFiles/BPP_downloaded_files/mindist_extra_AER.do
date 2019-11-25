cd C:\Users\edmun\OneDrive\Documents\Research\Norway\BPP_project\BP_AER_download
clear
cap log close
set mem 36m
set matsize 800
set more off

*TABLE 6 *
/*BASELINE:	010010*/
*TABLE 8 *
/*COL. 3:	010017*/  		/*Low wealth*/
/*COL. 4:	010018*/  		/*High wealth*/ 
/*COL. 5:	011017*/		/*Low wealth + total cons*/ 

version 6.0
				scalar incdef  =0		/* 0=earnings+transfers, 1=earnings+transfer+financial income		*/
								/* 2=earnings, 3=earnings+financial income, 4=male earnings			*/
								/* 5=male+female earnings, 6=income net of interg transf			*/
				scalar aftertax=1		/* 0=income before tax, 1=income after tax 					*/
				scalar totexp  =0		/* 0=non dur cons, 1=tot exp, 2=food, 3=non dur cons+ services from dur */
								/* 4=tot exp+services from dur							*/
				scalar SEO     =0		/* 0=representative sample, 1=representative sample+SEO, 2=Just SEO 	*/
				scalar agesel  =1		/* 0=20-65, 1=30-65, 2=25-65								*/
				scalar group   =8		/* 0=whole sample, 1=ww always, 2= ww sometimes, 3=unionized, 4=not un.
								   5=low wealth, 6=high wealth, 7=low wealth using 1979 data		
								   8=high wealth using 1979 data							*/
u data3,clear

******Does not use osbervations with topcoded income or financial income or federal taxes paid******
replace asset=. if trunca==1
replace y    =. if truncy==1
replace ftax =. if ftax==99999 & year>=1982 & year<=1985

gen ratio_ya=(y-asset)/y
replace ratio_ya=0 if ratio_ya<0

if SEO==0		{keep if seo==0}      
else if SEO==1	{keep if seo==0|seo==1}
else			{keep if seo==1}

if agesel==0	{drop if age<20|age>65}      
else if agesel==1	{drop if age<30|age>65}
else if agesel==2	{drop if age<25|age>65}
ren ncomp fsize

*Add a CPI price index and transform variables in real term
*This series is drawn from http://woodrow.mpls.frb.fed.us/economy/calc/hist1913.html
*CPI Base year is chained; 1982-1984 = 100
*2000 U.S. Department Of Labor, Bureau of Labor Statistics, Washington, D.C. 20212
*Consumer Price Index All Urban Consumers - (CPI-U) U.S. city average All items 1982-84=100

gen     price=1.403 if year==1992
replace price=1.362 if year==1991
replace price=1.307 if year==1990
replace price=1.240 if year==1989
replace price=1.183 if year==1988
replace price=1.136 if year==1987
replace price=1.096 if year==1986
replace price=1.076 if year==1985
replace price=1.039 if year==1984
replace price=0.996 if year==1983
replace price=0.965 if year==1982
replace price=0.909 if year==1981
replace price=0.824 if year==1980
replace price=0.726 if year==1979
replace price=0.652 if year==1978
replace price=0.606 if year==1977
replace price=0.569 if year==1976
replace price=0.538 if year==1975
replace price=0.493 if year==1974
replace price=0.444 if year==1973
replace price=0.418 if year==1972
replace price=0.405 if year==1971
replace price=0.388 if year==1970
replace price=0.367 if year==1969
replace price=0.348 if year==1968
replace price=0.334 if year==1967

*****Interest rate series ******
gen     rate=0.0764 if year==1978
replace rate=0.1056 if year==1979
replace rate=0.1210 if year==1980
replace rate=0.1460 if year==1981
replace rate=0.1094 if year==1982
replace rate=0.0899 if year==1983
replace rate=0.0990 if year==1984
replace rate=0.0771 if year==1985
replace rate=0.0609 if year==1986
replace rate=0.0588 if year==1987
replace rate=0.0694 if year==1988
replace rate=0.0844 if year==1989
replace rate=0.0769 if year==1990
replace rate=0.0543 if year==1991
replace rate=0.0348 if year==1992

gen tot_ass=((asset/rate)/price)+(house/price)
egen miny=min(year),by(person)
gen temp=tot_ass if year==miny
egen init_ass=mean(temp),by(person)

if incdef==0		{replace     y=(y-asset)/price}
else	if incdef==1	{replace     y=(y)/price}
else	if incdef==2	{replace     y=(y-trhw-troth-asset)/price}
else	if incdef==3	{replace     y=(y-trhw-troth)/price}
else	if incdef==4	{replace     y=ly/price}
else	if incdef==5	{replace     y=(ly+wly)/price}
else				{replace     y=(y-asset-help)/price}

replace    ly=ly   /price
replace   wly=wly  /price
replace  ftax=ftax /price

replace lc  =lc-ln(price)
replace lC  =lC-ln(price)
replace lcs =lcs-ln(price)
replace lCs =lCs-ln(price)

if totexp==0    		{drop lC lcs lCs}      
else if totexp==1       {drop lc lcs lCs
                 		 ren lC lc}
else if totexp==2		{drop lc lC lcs lCs
				 gen totf=food+fout+fstmp
				 gen lc=ln(totf/pf)}
else if totexp==3		{drop lc lC lCs
				 ren lcs lc}
else if totexp==4		{drop lc lC lcs
				 ren lCs lc}

if      aftertax==0  	{gen logy=ln(y)}      
else			      {gen logy=ln((y-ratio_ya*ftax))}

tab year,gen(yrd)

gen black =race==2
gen other =race>=3
gen selfe =self==3
gen empl  =empst==1
gen unempl=empst==2
gen retir =empst==3
gen kidsout   =outkid==1
gen bigcity=smsa==1|smsa==2
tab yb,gen(ybd)
tab fsize,gen(fd)
tab kids,gen(chd)
gen extra=(tyoth)>0

#delimit;
xi	i.educ*i.year i.white*i.year i.black*i.year i.other*i.year i.kidsout*i.year i.region*i.year i.extra*i.year
	i.year*i.bigcity i.year*i.kids i.fsize*i.year i.empl*i.year i.unempl*i.year i.retir*i.year;

/*Controls present:	year dummies, year of birth dummies, education dummies, race dummies,# of kids dummies, #of family member dummies,
				empl. status dummies, dummy for income recipient other than h/w, region dummies, dummy for big city, dummy for kids not in FU*/
/*Interactions present: educ*year race dummies*year employed*year unemployed*year retired*year region*year big city*year*/	
				 	
qui reg logy	yrd* ybd* edd*  white black other fd* chd* empl  unempl retir extra regd* bigcity kidsout 
			          Iey_* Iwy_* Iby_* Ioy_*          Ieyb* Iuy_*  Irya*       Iry_* Iyb_*;
predict uy if e(sample),res;
qui reg lc	      yrd* ybd* edd*  white black other fd* chd* empl  unempl retir extra regd* bigcity kidsout 
			          Iey_* Iwy_* Iby_* Ioy_*          Ieyb* Iuy_*  Irya*       Iry_* Iyb_*;
predict uc if e(sample),res;
#delimit cr

drop I*

xtile tempa=init_ass if year==miny,nq(10)
egen assq=mean(tempa),by(person)

if group==1		{keep if wly>0}
else if group==2	{keep if wly==0}
else if group==3	{keep if nn>0}
else if group==4	{keep if nn==0}
else if group==7	{keep if assq<=2}
else if group==8	{keep if assq>2}

sort person year
qui by person:gen duc=uc-uc[_n-1]
qui by person:gen duy=uy-uy[_n-1]
keep person year duy duc
save cohA,replace

         u cohA,clear
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
         gen id=group(nind)
         compress
	   gen coh=0
         keep  id year coh duy yduy duc yduc
         order id year coh duy yduy duc yduc
         save cohA,replace

*Set a minimum of 50 observations per year

         u cohA,clear
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
	   saveold cohA,replace

u cohA,clear
sort id year
keep  id year coh ndrod duy yduy duc yduc nmissd
order id year coh ndrod duy yduy duc yduc nmissd
saveold cohA,replace
stcmd cohA.dta cohA.dat /y
erase c:\papers\partial\revise2\cohA.dta

cd c:\gausswin
!gauss c:\papers\partial\revise2\cmV_AER.run
!exit

!gauss c:\papers\partial\revise2\MD_AER.prg
!exit
clear
