clear
set mem 100m
cd c:\papers\partial\revise2

cap program drop doit
program def doit
     local i=80
	   while `i' < 82 {
		   local j=1	
		   while `j' < 5 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						JOTAXNET MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 POCC_REF 
						PROPVALX PTAXRFDX QCURRIN1 QINTRVMO QINTRVYR RACE2 REF_RACE REGION           RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX;
				ren POCC_REF OCCUP;
				ren QCURRIN1 INDUS;
				gen RENTEQVX=99999999 if CUTENURE=="1"|CUTENURE=="2";
				drop PROPVALX;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=82
	   while `i' < 84 {
		   local j=1	
		   while `j' < 5 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						JOTAXNET MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX QPREVIN1 QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX;
				ren OCCUPRE1 OCCUP;
				ren QPREVIN1 INDUS;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=84
	   while `i' < 85 {
		   local j=1	
		   while `j' < 5 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						JOTAXNET MARITAL1          MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX QPREVIN1 QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX;
				ren OCCUPRE1 OCCUP;
				ren QPREVIN1 INDUS;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=85
	   while `i' < 86 {
		   local j=1	
		   while `j' < 5 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						JOTAXNET MARITAL1          MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX QPREVIN1 QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUPRE1 OCCUP;
				ren QPREVIN1 INDUS;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

				#delimit;
				u c:\cex\cex\85\fmlyq861.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						JOTAXNET MARITAL1          MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX QPREVIN1 QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUPRE1 OCCUP;
				ren QPREVIN1 INDUS;
				#delimit cr
			save c:\papers\partial\revise2\famyq861old,replace


cap program drop doit
program def doit
     local i=86
	   while `i' < 88 {
		   local j=1	
		   while `j' < 5 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						JOTAXNET MARITAL1         MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX QPREVIN1 QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUPRE1 OCCUP;
				ren QPREVIN1 INDUS;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit


cap program drop doit
program def doit
     local i=88
	   while `i' < 93 {
		   local j=1	
		   while `j' < 5 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX QPREVIN1 QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUPRE1 OCCUP;
				ren QPREVIN1 INDUS;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=93
	   while `i' < 94 {
		   local j=1	
		   while `j' < 3 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX          QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUPRE1 OCCUP;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=93
	   while `i' < 94 {
		   local j=3	
		   while `j' < 5 {
				#delimit;
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX CUREMPL1 CUREMPL2 CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUPRE1 
							   PTAXRFDX          QINTRVMO QINTRVYR RACE2 REF_RACE REGION          RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUPRE1 OCCUP;
				gen RENTEQVX=99999999 if CUTENURE=="1"|CUTENURE=="2";
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=94
	   while `i' < 95 {
		   local j=1	
		   while `j' < 5 {
				#delimit;	
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX                   CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUCOD1 
							   PTAXRFDX          QINTRVMO QINTRVYR RACE2 REF_RACE REGION          RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUCOD1 OCCUP;
				gen RENTEQVX=99999999 if CUTENURE=="1"|CUTENURE=="2";
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=95
	   while `i' < 103 {
		   local j=1	
		   while `j' < 5 {
				#delimit;	
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX                   CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUCOD1 
							   PTAXRFDX          QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX JFDSTMPA;
				ren OCCUCOD1 OCCUP;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

cap program drop doit
program def doit
     local i=103
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				#delimit;	
				u c:\cex\cex\\`i'\fmlyq`i'`j'.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX                   CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUCOD1 
							   PTAXRFDX          QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX ;
				ren OCCUCOD1 OCCUP;
				#delimit cr
			save c:\papers\partial\revise2\famyq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

				#delimit;
				u c:\cex\cex\95\fmlyq961.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX                   CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUCOD1 
							   PTAXRFDX          QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX ;
				ren OCCUCOD1 OCCUP;
				#delimit cr
			save c:\papers\partial\revise2\famyq961old,replace

				#delimit;
				u c:\cex\cex\103\fmlyq1041.dta,clear;
				gen intno=round(((NEWID/10)-int(NEWID/10))*10);
				replace NEWID=(NEWID-intno)/10;
				keep  	NEWID intno AGE_REF AGE2 BLS_URBN CKBKACTX COMPBND COMPBNDX COMPCKGX COMPOWDX COMPSAVX BUILDING 
						COMPSECX                   CUTENURE EARNINCX EDUC_REF EDUCA2 FAM_SIZE FAM_TYPE FAMTFEDX INTEARNX FININCX 
						FEDRFNDX FEDTAXX FINCATAX FINCBTAX FSLTAXX  INC_HRS1 INC_HRS2 INCWEEK1 INCWEEK2 
						         MARITAL1 MISCTAXX MONYOWDX OTHRFNDX PERSLT18 PERSOT64 OCCUCOD1 
							   PTAXRFDX          QINTRVMO QINTRVYR RACE2 REF_RACE REGION RENTEQVX RESPSTAT  
						SAVACCTX COMPSEC SECESTX SEX_REF SEX2 SLOCTAXX SLRFUNDX TAXPROPX TOTTXPDX USBNDX
						FSSIX UNEMPLX COMPENSX WELFAREX ;
				ren OCCUCOD1 OCCUP;
				#delimit cr
			save c:\papers\partial\revise2\famyq1041,replace


*************MERGE MEMBER AND FAMILY FILES *************************************

cap program drop doit
program def doit
     local i=80
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\famyq`i'`j'.dta,clear
				sort NEWID intno
				merge NEWID intno using c:\papers\partial\revise2\mq`i'`j'
				drop if _merge!=3
				drop _merge
				save,replace
				erase c:\papers\partial\revise2\mq`i'`j'.dta
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit


				u c:\papers\partial\revise2\famyq861old.dta,clear
				sort NEWID intno
				merge NEWID intno using c:\papers\partial\revise2\mq861old
				drop if _merge!=3
				drop _merge
				save,replace
				erase c:\papers\partial\revise2\mq861old.dta
				u c:\papers\partial\revise2\famyq961old.dta,clear
				sort NEWID intno
				merge NEWID intno using c:\papers\partial\revise2\mq961old
				drop if _merge!=3
				drop _merge
				save,replace
				erase c:\papers\partial\revise2\mq961old.dta
				u c:\papers\partial\revise2\famyq1041.dta,clear
				sort NEWID intno
				merge NEWID intno using c:\papers\partial\revise2\mq1041
				drop if _merge!=3
				drop _merge
				save,replace
				erase c:\papers\partial\revise2\mq1041.dta

clear




