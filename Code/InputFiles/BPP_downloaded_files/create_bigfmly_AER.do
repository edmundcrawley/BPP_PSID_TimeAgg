clear
set mem 250m
cd c:\papers\partial\revise2

cap program drop doit
program def doit
     local i=80
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\famyq`i'`j'.dta,clear
				egen assets_1=rsum(SAVACCTX CKBKACTX USBNDX   SECESTX  MONYOWDX)
				egen diff    =rsum(COMPSAVX COMPCKGX COMPBNDX COMPSECX COMPOWDX)
				gen  assets_0=assets_1-diff
				gen from=`i'`j'
				save c:\papers\partial\revise2\fq`i'`j'.dta,replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit


				u c:\papers\partial\revise2\famyq861old.dta,clear
				egen assets_1=rsum(SAVACCTX CKBKACTX USBNDX   SECESTX  MONYOWDX)
				egen diff    =rsum(COMPSAVX COMPCKGX COMPBNDX COMPSECX COMPOWDX)
				gen  assets_0=assets_1-diff
				gen from=860
				save c:\papers\partial\revise2\fq861old,replace
				u c:\papers\partial\revise2\famyq961old.dta,clear
				egen assets_1=rsum(SAVACCTX CKBKACTX USBNDX   SECESTX  MONYOWDX)
				egen diff    =rsum(COMPSAVX COMPCKGX COMPBNDX COMPSECX COMPOWDX)
				gen  assets_0=assets_1-diff
				gen from=960
				save c:\papers\partial\revise2\fq961old,replace
				u c:\papers\partial\revise2\famyq1041.dta,clear
				egen assets_1=rsum(SAVACCTX CKBKACTX USBNDX   SECESTX  MONYOWDX)
				egen diff    =rsum(COMPSAVX COMPCKGX COMPBNDX COMPSECX COMPOWDX)
				gen  assets_0=assets_1-diff
				gen from=1041
				save c:\papers\partial\revise2\fq1041.dta,replace
clear

cap program drop doit
program def doit
     local i=80
	   while `i' < 82 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\fq`i'`j'.dta,clear
				gen str1 temp="7" if EDUCA2=="00"
				replace  temp="1" if EDUCA2=="01"|EDUCA2=="02"|EDUCA2=="03"|EDUCA2=="04"|EDUCA2=="05"|EDUCA2=="06"|EDUCA2=="07"|EDUCA2=="08"
				replace  temp="2" if EDUCA2=="09"|EDUCA2=="10"|EDUCA2=="11"
				replace  temp="3" if EDUCA2=="12"
				replace  temp="4" if EDUCA2=="13"|EDUCA2=="14"|EDUCA2=="15"
				replace  temp="5" if EDUCA2=="16"
				replace  temp="6" if EDUCA2=="17"|EDUCA2=="18"
				replace EDUCA2=temp
				save,replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=82
	   while `i' < 88 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\fq`i'`j'.dta,clear
				gen str1 temp="7" if EDUCA2=="00"
				replace  temp="1" if EDUCA2=="01"|EDUCA2=="02"|EDUCA2=="03"|EDUCA2=="04"|EDUCA2=="05"|EDUCA2=="06"|EDUCA2=="07"|EDUCA2=="08"
				replace  temp="2" if EDUCA2=="09"|EDUCA2=="10"|EDUCA2=="11"
				replace  temp="3" if EDUCA2=="12"
				replace  temp="4" if EDUCA2=="21"|EDUCA2=="22"|EDUCA2=="23"
				replace  temp="5" if EDUCA2=="24"
				replace  temp="6" if EDUCA2=="31"|EDUCA2=="32"
				replace EDUCA2=temp
				save,replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

				u c:\papers\partial\revise2\fq861old.dta,clear
				gen str1 temp="7" if EDUCA2=="00"
				replace  temp="1" if EDUCA2=="01"|EDUCA2=="02"|EDUCA2=="03"|EDUCA2=="04"|EDUCA2=="05"|EDUCA2=="06"|EDUCA2=="07"|EDUCA2=="08"
				replace  temp="2" if EDUCA2=="09"|EDUCA2=="10"|EDUCA2=="11"
				replace  temp="3" if EDUCA2=="12"
				replace  temp="4" if EDUCA2=="21"|EDUCA2=="22"|EDUCA2=="23"
				replace  temp="5" if EDUCA2=="24"
				replace  temp="6" if EDUCA2=="31"|EDUCA2=="32"
				replace EDUCA2=temp
				save,replace

cap program drop doit
program def doit
     local i=96
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\fq`i'`j'.dta,clear
				gen str1 temp="7" if EDUCA2=="00"
				replace  temp="1" if EDUCA2=="10"
				replace  temp="2" if EDUCA2=="11"
				replace  temp="3" if EDUCA2=="12"
				replace  temp="4" if EDUCA2=="13"|EDUCA2=="14"
				replace  temp="5" if EDUCA2=="15"
				replace  temp="6" if EDUCA2=="16"|EDUCA2=="17"
				replace EDUCA2=temp
				gen str1 temp2="7" if EDUC_REF=="00"
				replace  temp2="1" if EDUC_REF=="10"
				replace  temp2="2" if EDUC_REF=="11"
				replace  temp2="3" if EDUC_REF=="12"
				replace  temp2="4" if EDUC_REF=="13"|EDUC_REF=="14"
				replace  temp2="5" if EDUC_REF=="15"
				replace  temp2="6" if EDUC_REF=="16"|EDUC_REF=="17"
				replace EDUC_REF=temp2
				save,replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

				u c:\papers\partial\revise2\fq1041.dta,clear
				gen str1 temp="7" if EDUCA2=="00"
				replace  temp="1" if EDUCA2=="10"
				replace  temp="2" if EDUCA2=="11"
				replace  temp="3" if EDUCA2=="12"
				replace  temp="4" if EDUCA2=="13"|EDUCA2=="14"
				replace  temp="5" if EDUCA2=="15"
				replace  temp="6" if EDUCA2=="16"|EDUCA2=="17"
				replace EDUCA2=temp
				gen str1 temp2="7" if EDUC_REF=="00"
				replace  temp2="1" if EDUC_REF=="10"
				replace  temp2="2" if EDUC_REF=="11"
				replace  temp2="3" if EDUC_REF=="12"
				replace  temp2="4" if EDUC_REF=="13"|EDUC_REF=="14"
				replace  temp2="5" if EDUC_REF=="15"
				replace  temp2="6" if EDUC_REF=="16"|EDUC_REF=="17"
				replace EDUC_REF=temp2
				save,replace


cap program drop doit
program def doit
     local i=80
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\fq`i'`j'.dta,clear
				gen     age=AGE_REF
				replace age=AGE2 		if SEX_REF=="2" & SEX2=="1"
				gen     school=EDUC_REF
				replace school=EDUCA2 	if SEX_REF=="2" & SEX2=="1"
				gen     race=REF_RACE
				replace race=RACE2 	if SEX_REF=="2" & SEX2=="1"
				gen     wkhrs=INC_HRS1
				replace wkhrs=INC_HRS2 	if SEX_REF=="2" & SEX2=="1"
				gen     num_wk=INCWEEK1
				replace num_wk=INCWEEK2 if SEX_REF=="2" & SEX2=="1"
				gen     wages=salary1
				replace wages=salary2 	if SEX_REF=="2" & SEX2=="1"
				gen     agew=AGE2 		if SEX2=="2"	
				replace agew=AGE_REF 		if SEX2=="1"	
				gen     schoolw=EDUCA2 		if SEX2=="2"	
				replace schoolw=EDUC_REF 	if SEX2=="1"	
				gen     racew=RACE2 		if SEX2=="2"	
				replace racew=REF_RACE 		if SEX2=="1"	
				gen     wkhrsw=INC_HRS2 	if SEX2=="2"	
				replace wkhrsw=INC_HRS1 	if SEX2=="1"	
				gen     num_wkw=INCWEEK2 	if SEX2=="2"	
				replace num_wkw=INCWEEK1 	if SEX2=="1"	
				gen     wagesw=salary2 		if SEX2=="2"	
				replace wagesw=salary1 		if SEX2=="1"	
				#delimit;
				keep  NEWID intno age* school* race* FAM_* FINCATAX FINCBTAX SEX_REF RESPSTAT QINT* REGION from INTEARNX FININCX
					  BLS_URBN wkhrs* PERSLT18 assets* num_wk* wages* BUILDING RENTEQVX CUTENURE;
				#delimit cr
				sort NEWID intno
				save c:\papers\partial\revise2\data`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

				u c:\papers\partial\revise2\fq861old.dta,clear
				gen     age=AGE_REF
				replace age=AGE2 		if SEX_REF=="2" & SEX2=="1"
				gen     school=EDUC_REF
				replace school=EDUCA2 	if SEX_REF=="2" & SEX2=="1"
				gen     race=REF_RACE
				replace race=RACE2 	if SEX_REF=="2" & SEX2=="1"
				gen     wkhrs=INC_HRS1
				replace wkhrs=INC_HRS2 	if SEX_REF=="2" & SEX2=="1"
				gen     num_wk=INCWEEK1
				replace num_wk=INCWEEK2 if SEX_REF=="2" & SEX2=="1"
				gen     wages=salary1
				replace wages=salary2 	if SEX_REF=="2" & SEX2=="1"
				gen     agew=AGE2 		if SEX2=="2"	
				replace agew=AGE_REF 		if SEX2=="1"	
				gen     schoolw=EDUCA2 		if SEX2=="2"	
				replace schoolw=EDUC_REF 	if SEX2=="1"	
				gen     racew=RACE2 		if SEX2=="2"	
				replace racew=REF_RACE 		if SEX2=="1"	
				gen     wkhrsw=INC_HRS2 	if SEX2=="2"	
				replace wkhrsw=INC_HRS1 	if SEX2=="1"	
				gen     num_wkw=INCWEEK2 	if SEX2=="2"	
				replace num_wkw=INCWEEK1 	if SEX2=="1"	
				gen     wagesw=salary2 		if SEX2=="2"	
				replace wagesw=salary1 		if SEX2=="1"	
				#delimit;
				keep  NEWID intno age* school* race* FAM_* FINCATAX FINCBTAX SEX_REF RESPSTAT QINT* REGION from INTEARNX FININCX
					 BLS_URBN wkhrs* PERSLT18 assets* num_wk* wages* BUILDING RENTEQVX CUTENURE;
				#delimit cr
				sort NEWID intno
				save c:\papers\partial\revise2\data861old,replace
				u c:\papers\partial\revise2\fq961old.dta,clear
				gen     age=AGE_REF
				replace age=AGE2 		if SEX_REF=="2" & SEX2=="1"
				gen     school=EDUC_REF
				replace school=EDUCA2 	if SEX_REF=="2" & SEX2=="1"
				gen     race=REF_RACE
				replace race=RACE2 	if SEX_REF=="2" & SEX2=="1"
				gen     wkhrs=INC_HRS1
				replace wkhrs=INC_HRS2 	if SEX_REF=="2" & SEX2=="1"
				gen     num_wk=INCWEEK1
				replace num_wk=INCWEEK2 if SEX_REF=="2" & SEX2=="1"
				gen     wages=salary1
				replace wages=salary2 	if SEX_REF=="2" & SEX2=="1"
				gen     agew=AGE2 		if SEX2=="2"	
				replace agew=AGE_REF 		if SEX2=="1"	
				gen     schoolw=EDUCA2 		if SEX2=="2"	
				replace schoolw=EDUC_REF 	if SEX2=="1"	
				gen     racew=RACE2 		if SEX2=="2"	
				replace racew=REF_RACE 		if SEX2=="1"	
				gen     wkhrsw=INC_HRS2 	if SEX2=="2"	
				replace wkhrsw=INC_HRS1 	if SEX2=="1"	
				gen     num_wkw=INCWEEK2 	if SEX2=="2"	
				replace num_wkw=INCWEEK1 	if SEX2=="1"	
				gen     wagesw=salary2 		if SEX2=="2"	
				replace wagesw=salary1 		if SEX2=="1"	
				#delimit;
				keep  NEWID intno age* school* race* FAM_* FINCATAX FINCBTAX SEX_REF RESPSTAT QINT* REGION from INTEARNX FININCX
					 BLS_URBN wkhrs* PERSLT18 assets* num_wk* wages* BUILDING RENTEQVX CUTENURE;
				#delimit cr
				sort NEWID intno
				save c:\papers\partial\revise2\data961old,replace
				u c:\papers\partial\revise2\fq1041.dta,clear
				gen     age=AGE_REF
				replace age=AGE2 		if SEX_REF=="2" & SEX2=="1"
				gen     school=EDUC_REF
				replace school=EDUCA2 	if SEX_REF=="2" & SEX2=="1"
				gen     race=REF_RACE
				replace race=RACE2 	if SEX_REF=="2" & SEX2=="1"
				gen     wkhrs=INC_HRS1
				replace wkhrs=INC_HRS2 	if SEX_REF=="2" & SEX2=="1"
				gen     num_wk=INCWEEK1
				replace num_wk=INCWEEK2 if SEX_REF=="2" & SEX2=="1"
				gen     wages=salary1
				replace wages=salary2 	if SEX_REF=="2" & SEX2=="1"
				gen     agew=AGE2 		if SEX2=="2"	
				replace agew=AGE_REF 		if SEX2=="1"	
				gen     schoolw=EDUCA2 		if SEX2=="2"	
				replace schoolw=EDUC_REF 	if SEX2=="1"	
				gen     racew=RACE2 		if SEX2=="2"	
				replace racew=REF_RACE 		if SEX2=="1"	
				gen     wkhrsw=INC_HRS2 	if SEX2=="2"	
				replace wkhrsw=INC_HRS1 	if SEX2=="1"	
				gen     num_wkw=INCWEEK2 	if SEX2=="2"	
				replace num_wkw=INCWEEK1 	if SEX2=="1"	
				gen     wagesw=salary2 		if SEX2=="2"	
				replace wagesw=salary1 		if SEX2=="1"	
				#delimit;
				keep  NEWID intno age* school* race* FAM_* FINCATAX FINCBTAX SEX_REF RESPSTAT QINT* REGION from INTEARNX FININCX
					 BLS_URBN wkhrs* PERSLT18 assets* num_wk* wages* BUILDING RENTEQVX CUTENURE;
				#delimit cr
				sort NEWID intno
				save c:\papers\partial\revise2\data1041,replace



cap program drop doit
program def doit
     local i=80
	   while `i' < 86 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\data`i'`j'.dta,clear
				gen family=NEWID
				sort family intno
				save,replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=86
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\data`i'`j'.dta,clear
				gen family=NEWID+70000
				sort family intno
				save,replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

				u c:\papers\partial\revise2\data961.dta,clear
				keep family intno
				sort family intno
				save c:\papers\partial\revise2\data961_f.dta,replace

				u c:\papers\partial\revise2\data861old.dta,clear
				gen family=NEWID
				sort family intno
				save,replace										/*The file 861old has no duplicates ID*/
				u c:\papers\partial\revise2\data961old.dta,clear
				gen family=NEWID+70000
				sort family intno
				merge family intno using c:\papers\partial\revise2\data961_f
				keep if _merge==1
				drop _merge
				sort family intno
				save,replace										/*The file 961old has duplicates ID - we're taking only the nonduplicates here*/ 
				u c:\papers\partial\revise2\data1041.dta,clear
				gen family=NEWID+70000
				sort family intno
				save,replace


clear
set mem 350m

u            c:\papers\partial\revise2\data801,clear
append using c:\papers\partial\revise2\data802
append using c:\papers\partial\revise2\data803
append using c:\papers\partial\revise2\data804
cap program drop doit
program def doit
     local i=81
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				append using c:\papers\partial\revise2\data`i'`j'
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit
append using c:\papers\partial\revise2\data861old
append using c:\papers\partial\revise2\data961old
append using c:\papers\partial\revise2\data1041
compress
save c:\papers\partial\revise2\fmlydata,replace


cap program drop doit
program def doit
     local i=80
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				erase c:\papers\partial\revise2\data`i'`j'.dta
				erase c:\papers\partial\revise2\fq`i'`j'.dta
				erase c:\papers\partial\revise2\famyq`i'`j'.dta
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit
				erase c:\papers\partial\revise2\data861old.dta
				erase c:\papers\partial\revise2\data961old.dta
				erase c:\papers\partial\revise2\data1041.dta
				erase c:\papers\partial\revise2\fq861old.dta
				erase c:\papers\partial\revise2\fq961old.dta
				erase c:\papers\partial\revise2\fq1041.dta
				erase c:\papers\partial\revise2\famyq861old.dta
				erase c:\papers\partial\revise2\famyq961old.dta
				erase c:\papers\partial\revise2\famyq1041.dta

				erase c:\papers\partial\revise2\data961_f.dta

clear


