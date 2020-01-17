******DOES SAMPLE SELECTION *************

cd C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files
clear
set mem 100m
set more off
cap log close
log using C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\adjust_AER,replace t

u data,clear
drop pid
sort person year

sort person year
qby person: gen dyear=year-year[_n-1]
egen todrop=sum(dyear>1 & dyear!=.),by(person)		/*Drop those with intermittent "headship" */
egen n=sum(person!=.),by(person)				/*Drop those appearing only once		*/
replace todrop=1 if n==1
drop if todrop>0
drop todrop dyear n

drop if year<1978							/*We don't use data before 78*/

*Takes into account the retrospective nature of the data
replace age=age-1

*Makes the variable race consistent over the years
*Now 1 is white, 2 black, 3 others

gen     rc=race if race<=2
replace rc=3    if race>=3 & race<=7
drop race
ren rc race

/*
Makes the variable education consistent over the years
Now 1 is 0-11 grades (includes those with no schooling, i.e. can't read/write people); 
    2 is 12 grades, i.e. high school w or w/o nonacademic training;
    3 is college dropout, college degree, or collage & advanced/professional degree;
    a missing point denotes na/dk.
*/

gen     sc=1 if educ<=3 & year<=1989                  	/*0-11  grades*/
replace sc=2 if (educ==4|educ==5) & year<=1989		      /*High school or 12 grades+nonacademic training*/
replace sc=3 if educ>=6 & educ<=8 & year<=1989            	/*College dropout, BA degree, or collage & adv./prof. degree*/
replace sc=. if educ>8 & year<=1989                   	/*Missing, NA, DK*/

replace sc=1 if educ>=0  & educ<=11 & year>1989           	/*0-11  grades*/
replace sc=2 if educ==12                        		/*High school or 12 grades+nonacademic training*/
replace sc=3 if educ>=13 & educ<=17 & year>1989           	/*College dropout, BA degree, or collage & adv./prof. degree*/
replace sc=. if educ>17 & year>1989                   	/*Missing, NA, DK*/
drop educ
ren sc educ

/*Make some manual adjustment*/

sort person year
qui by person:replace educ=educ[_n-1] if educ==.
gsort person -year
qui by person:replace educ=educ[_n-1] if educ==.
sort person year

*Take the maximum grade achieved as the relevant education level

egen maxed=max(educ),by(person)
gen educ2=educ
replace educ=maxed
drop maxed 


*Design a demographically unstable household as one where some family composition change (apart from changes in
*people other than head-wife takes place).
*Then drop households with unstable demographical pattern

/*
fchg=0 (No change), 1 (Change in members other than Head or Wife), 2 (Head same, wife left,die, or is new)
3 (wife is now head), 4 (ex female head married and huisband is now head), 5 (some sample member other than ex Head or Wife is now head)
7 (ex wife head because husband in inst., now husband back and is head), 8 (Other) 
*/

egen miny=min(year),by(person)
gen todrop1=fchg>1
replace todrop1=0 if year==miny & todrop1==1	/*starting household structure is when they first appear in sample*/
egen todrop2=sum(todrop1),by(person)
drop if todrop2!=0
drop todrop* miny


/*Drop female household heads*/
drop if sex==2

*Drop households with missing observations on race, education, or state of residence

drop if race==.

gen m=educ==.
egen mm=sum(m),by(person)
drop if mm>0
drop m mm

****Some manual adjustments****
replace age = 30 if person==5084 & age==98 
replace age = 31 if person==5084 & age==. 
replace agew =77 if person==6927 & agew==.

/*Recode age so that there is no gap or jump*/

egen lasty=max(year),by(person)
gen lastage=age if year==lasty
gen b=year-lastage
replace b=0 if b==.
egen yb=sum(b),by(person)
replace age=year-yb
drop lasty lastage yb b

/* Generate the regions                          */

egen n=sum(state==.|state==0|state==99),by(person)
drop if n>0
drop n

# delimit;
gen     region=1 if state==6  | state==18 | state==20 | state==28 | state==29 | state==31 
                              | state==37 | state==38 | state==44; 						 /* North East*/
replace region=2 if state==12 | state==13 | state==14 | state==15 | state==21 | state==22
                              | state==24 | state==26 | state==33 | state==34 | state==40 | state==48; /*Midwest*/
replace region=3 if state==1  | state==3  | state==7  | state==8  | state==9  | state==10 | state==16  
                              | state==17 | state==19 | state==23 | state==32 | state==35 
                              | state==39 | state==41 | state==42 | state==45 | state==47; 		 /*South*/
replace region=4 if state==2  | state==4  | state==5  | state==11 | state==25 | state==27 | state==30 
                              | state==36 | state==43 | state==46 | state==49 | state==50 | state==51; /*West*/                              
#delimit cr

gen n=marit==1
egen nn=sum(n),by(person)
gen b=1
egen bb=sum(b),by(person)
keep if nn==bb				/*Keep those who are continuously married*/
drop nn bb n b

/*Drop those with extreme (+ve or -ve) growth (outliers)*/

sort person year
qui by person:gen gy=(y-y[_n-1])/(y[_n-1])
gen m=gy>5 & gy!=.|gy<-0.8 & gy!=.|y<=100			/*drop those whose income grows more than 500%, falls by more than 80%, or is below $100*/
egen mm=sum(m),by(person)
lab var mm "0=no income growth outlier"
tab mm
drop if mm>0
drop m gy mm


sort person year

gen yb=year-age

drop if yb<=1919|yb>=1960
egen help=rsum(hhelp whelp ohelp)	/*Create variable help received, from husb, wife or others*/
drop hhelp whelp ohelp

save data2,replace

log close
clear
