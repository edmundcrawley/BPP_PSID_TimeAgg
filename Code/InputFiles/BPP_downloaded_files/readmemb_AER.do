clear
set mem 250m
cd c:\papers\partial\revise2

cap program drop doit
program def doit
         local i=80
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
			u c:\cex\cex\\`i'\membq`i'`j'.dta,clear
			gen intno=round(((NEWID/10)-int(NEWID/10))*10)
			replace NEWID=(NEWID-intno)/10
			keep NEWID intno CU_CODE SALARYX
			gen salary1=SALARYX if CU_CODE=="1"			/*Ref pers*/
			gen salary2=SALARYX if CU_CODE=="2"			/*Spouse*/
			collapse salary1 salary2,by(NEWID intno)
			compress
			sort NEWID intno
			save c:\papers\partial\revise2\mq`i'`j',replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

			u c:\cex\cex\85\membq861.dta,clear
			gen intno=round(((NEWID/10)-int(NEWID/10))*10)
			replace NEWID=(NEWID-intno)/10
			keep NEWID intno CU_CODE SALARYX
			gen salary1=SALARYX if CU_CODE=="1"
			gen salary2=SALARYX if CU_CODE=="2"
			collapse salary1 salary2,by(NEWID intno)
			compress
			sort NEWID intno
			save c:\papers\partial\revise2\mq861old,replace

			u c:\cex\cex\95\membq961.dta,clear
			gen intno=round(((NEWID/10)-int(NEWID/10))*10)
			replace NEWID=(NEWID-intno)/10
			keep NEWID intno CU_CODE SALARYX
			gen salary1=SALARYX if CU_CODE=="1"
			gen salary2=SALARYX if CU_CODE=="2"
			collapse salary1 salary2,by(NEWID intno)
			compress
			sort NEWID intno
			save c:\papers\partial\revise2\mq961old,replace

			u c:\cex\cex\103\membq1041.dta,clear
			gen intno=round(((NEWID/10)-int(NEWID/10))*10)
			replace NEWID=(NEWID-intno)/10
			keep NEWID intno CU_CODE SALARYX
			gen salary1=SALARYX if CU_CODE=="1"
			gen salary2=SALARYX if CU_CODE=="2"
			collapse salary1 salary2,by(NEWID intno)
			compress
			sort NEWID intno
			save c:\papers\partial\revise2\mq1041,replace

clear
