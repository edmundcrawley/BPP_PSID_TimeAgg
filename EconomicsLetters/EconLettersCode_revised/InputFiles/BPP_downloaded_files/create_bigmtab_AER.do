clear
set mem 250m
cd c:\papers\partial\revise2


cap program drop doit
program def doit
     local i=80
	   while `i' < 86 {
		   local j=1	
		   while `j' < 5 {
				u c:\papers\partial\revise2\expq`i'`j'.dta,clear
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
				u c:\papers\partial\revise2\expq`i'`j'.dta,clear
				gen family=NEWID+70000
				sort family intno
				save,replace
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit

				u c:\papers\partial\revise2\expq961.dta,clear
				keep family intno
				sort family intno
				save c:\papers\partial\revise2\expq961_f.dta,replace

				u c:\papers\partial\revise2\expq861old.dta,clear
				gen family=NEWID
				sort family intno
				save,replace
				u c:\papers\partial\revise2\expq961old.dta,clear
				gen family=NEWID+70000
				sort family intno
				merge family intno using c:\papers\partial\revise2\expq961_f
				keep if _merge==1
				drop _merge
				sort family intno
				save,replace
				u c:\papers\partial\revise2\expq1041.dta,clear
				gen family=NEWID+70000
				sort family intno
				save,replace

clear
set mem 350m

u            c:\papers\partial\revise2\expq801,clear
append using c:\papers\partial\revise2\expq802
append using c:\papers\partial\revise2\expq803
append using c:\papers\partial\revise2\expq804
cap program drop doit
program def doit
     local i=81
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				append using c:\papers\partial\revise2\expq`i'`j'
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit
append using c:\papers\partial\revise2\expq861old
append using c:\papers\partial\revise2\expq961old
append using c:\papers\partial\revise2\expq1041
compress
save c:\papers\partial\revise2\mtabdata,replace


cap program drop doit
program def doit
     local i=80
	   while `i' < 104 {
		   local j=1	
		   while `j' < 5 {
				erase c:\papers\partial\revise2\expq`i'`j'.dta
	         local j=`j'+1
                 }
         local i=`i'+1
                 }
end
qui doit
				erase c:\papers\partial\revise2\expq861old.dta
				erase c:\papers\partial\revise2\expq961old.dta
				erase c:\papers\partial\revise2\expq1041.dta

				erase c:\papers\partial\revise2\expq961_f.dta

clear


