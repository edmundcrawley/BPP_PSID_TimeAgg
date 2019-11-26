****EXTRACTS FROM THE FAMILY FILES THE VARIABLES OF INTEREST *****

cd C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files
clear
set mem 30m
set more off

cap log close
log using C:\Users\edmun\OneDrive\Documents\Research\BPP_PSID_TimeAgg\Code\InputFiles\BPP_downloaded_files\create1_AER,replace t

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f68,clear
rename *, lower
#delimit;
keep v93 v95 v47 v3 v117 v118 v119 v398 v196 v198 v181 v81 v74 v75 v76 v79 
     v337 v239 v439 v313 v115 v37 v83 v53 v45 v315 v10 v5 v14 v143 v251 v99;
#delimit cr    
ren v99	intdate 
ren v251  wages
ren v14    util          
ren v143  ncars
ren v5    house
ren v10    rent
gen       disab=.
gen      outkid=.
ren v315    vet
ren v45   fstmp
ren v53   hourw
ren v93   state
ren v95    smsa
ren v47   hours
ren v3       id
ren v117    age
ren v118   agew
ren v119    sex
ren v398   kids
ren v196  empst
ren v198   self
ren v181   race
ren v81       y
ren v74      ly
ren v75     wly
ren v76    tyhw
ren v79   tyoth
gen       troth=0
gen        trhw=y-tyhw-tyoth
ren v337   avhy
ren v239  marit
ren v439 weight
ren v313   educ
ren v115  fsize
gen       split=0
gen        fchg=0
gen       newhd=0
ren v37    food
gen        fout=.
ren v83    ftax
gen 	    asset=.
gen        trhf=.
gen year=1967
save        fam67,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f69,clear
rename *, lower
#delimit;
keep v537 v539 v909 v542 v465 v493 v442 v1008 v1010 v1011 v550 v639 v641 v791 v801 
     v529 v514 v516 v518 v521 v527 
     v871 v607 v1014 v794 v500 v506 v532 v525 v475 v510 v796 v732 v743 v455 v449 v460 v581 v699 v553;
#delimit cr    
ren v553 intdate 
ren v699   wages
ren v460    util 
ren v581   ncars
ren v449   house
ren v455    rent
ren v743   disab
ren v732  outkid
ren v796     vet
ren v510   fstmp
ren v475   hourw
ren v537   state
ren v539    smsa
ren v909   split
ren v542    fchg
ren v465   hours
ren v493   fsize
ren v442      id
ren v1008    age
ren v1011   agew
ren v1010    sex
ren v550    kids
ren v639   empst
ren v641    self
ren v791   newhd
ren v801    race
ren v529       y
ren v514      ly
ren v516     wly
ren v518    tyhw
ren v521   tyoth
ren v527   troth
ren v525    trhw
ren v871    avhy
ren v607   marit
ren v1014 weight
ren v794    educ
ren v500    food
ren v506    fout
ren v532    ftax
gen        asset=.
gen         trhf=.
gen year=1968
save         fam68,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f70,clear
rename *, lower
#delimit;
keep v1103 v1506 v1106 v1109 v1138 v1167 v1102 v1239 v1240 v1241 v1242 v1278 v1280 v1461 v1490 
     v1514 v1196 v1198 v1205 v1222 v1225 v1220 
     v1567 v1365 v1609 v1485 v1175 v1185 v1122 v1208 v1224 v1148 v1183 v1487 v1399 v1409 v1128 v1122 v1133 v1254 v1191 v1236;
#delimit cr    
ren v1236 intdate 
ren v1191   wages
ren v1133    util 
ren v1254   ncars
ren v1122   house
ren v1128    rent
ren v1409   disab
ren v1399  outkid
ren v1487     vet
ren v1183   fstmp
ren v1148   hourw
ren v1103   state
ren v1506    smsa
ren v1106   split
ren v1109    fchg
ren v1138   hours
ren v1167   fsize
ren v1102      id
ren v1239     age
ren v1241    agew
ren v1240     sex
ren v1242    kids
ren v1278   empst
ren v1280    self
ren v1461   newhd
ren v1490    race
ren v1514       y
ren v1196      ly
ren v1198     wly
ren v1205    tyhw
ren v1222   tyoth
ren v1225   troth
ren v1220    trhw
ren v1567    avhy
ren v1365   marit
ren v1609  weight
ren v1485    educ
ren v1175    food
ren v1185    fout
gen          ftax=v1208+v1224 if v1208!=. & v1224!=.
replace      ftax=v1208       if v1208!=. & v1224==.
replace      ftax=v1224       if v1208==. & v1224!=.
drop v1208 v1224
gen 	      asset=.
gen          trhf=.
gen year=1969
save          fam69,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f71,clear
rename *, lower
#delimit;
keep v1803 v1816 v1806 v1809 v1839 v1868 v1802 v1942 v1943 v1944 v1945 v1983 v1986 v2165 v2202 
     v2226 v1897 v1899 v1906 v1924 v1927 v1922 
     v2279 v2072 v2321 v2197 v1876 v1886 v1910 v1926 v1849 v1884 v2199 v2111 v2121 v1829 v1823 v1834 v1957 v1892 v1939;
#delimit cr    
ren v1939 intdate                              
ren v1892   wages
ren v1834    util 
ren v1957   ncars
ren v1823   house
ren v1829    rent
ren v2121   disab
ren v2111  outkid
ren v2199     vet
ren v1884   fstmp
ren v1849   hourw
ren v1803   state
ren v1816    smsa
ren v1806   split
ren v1809    fchg
ren v1839   hours
ren v1868   fsize
ren v1802      id
ren v1942     age
ren v1944    agew
ren v1943     sex
ren v1945    kids
ren v1983   empst
ren v1986    self
ren v2165   newhd
ren v2202    race
ren v2226       y
ren v1897      ly
ren v1899     wly
ren v1906    tyhw
ren v1924   tyoth
ren v1927   troth
ren v1922    trhw
ren v2279    avhy
ren v2072   marit
ren v2321  weight
ren v2197    educ
ren v1876    food
ren v1886    fout
gen          ftax=v1910+v1926 if v1910!=. & v1926!=.
replace 	 ftax=v1910       if v1910!=. & v1926==.
replace      ftax=v1926       if v1910==. & v1926!=.
drop v1910 v1926
gen         asset=.
gen          trhf=.
gen year=1970
save          fam70,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f72,clear
rename *, lower
#delimit;
keep v2403 v2406 v2407 v2410 v2439 v2468 v2402 v2542 v2543 v2544 v2545 v2581 v2584 v2791 v2828 
     v2852 v2498 v2500 v2507 v2525 v2528 v2523 
     v2906 v2670 v2968 v2823 v2476 v2480 v2511 v2527 v2449 v2478 v2825 v2708 v2718 v2429 v2423 v2434 v2557 v2493 v2539;
#delimit cr    
ren v2539 intdate 
ren v2493   wages
ren v2434    util 
ren v2557   ncars
ren v2423   house 
ren v2429    rent
ren v2718   disab
ren v2708  outkid
ren v2825     vet
ren v2478   fstmp
ren v2449   hourw
ren v2403   state 
ren v2406    smsa
ren v2407   split
ren v2410    fchg
ren v2439   hours
ren v2468   fsize
ren v2402      id
ren v2542     age
ren v2544    agew
ren v2543     sex
ren v2545    kids
ren v2581   empst
ren v2584    self
ren v2791   newhd
ren v2828    race
ren v2852       y
ren v2498      ly
ren v2500     wly
ren v2507    tyhw
ren v2525   tyoth
ren v2528   troth
ren v2523    trhw
ren v2906    avhy
ren v2670   marit
ren v2968  weight
ren v2823    educ
ren v2476    food
ren v2480    fout

gen     ftax=v2511+v2527 if v2511!=. & v2527!=.
replace ftax=v2511       if v2511!=. & v2527==.
replace ftax=v2527       if v2511==. & v2527!=.
drop v2511 v2527
gen asset=.
gen trhf=.
gen year=1971
save          fam71,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f73,clear
rename *, lower
#delimit;
keep v3003 v3006 v3007 v3010 v3027 v3017 v3002 v3095 v3096 v3097 v3098 v3114 v3117 v3217 v3300 
     v3256 v3051 v3053 v3060 v3078 v3081 v3076
     v3275 v3181 v3301 v3241 v3064 v3080 v3035 v3243 v3211 v3244 v3023 v3021 v3046 v3092;
#delimit cr    
ren v3092 intdate     
ren v3046   wages
gen          util=.
gen         ncars=.
ren v3021   house
ren v3023    rent
ren v3244   disab
ren v3211  outkid
ren v3243     vet
gen         fstmp=.
ren v3035   hourw
ren v3003   state
ren v3006    smsa
ren v3007   split
ren v3010    fchg
ren v3027   hours
ren v3017   fsize
ren v3002      id
ren v3095     age
ren v3097    agew
ren v3096     sex
ren v3098    kids
ren v3114   empst
ren v3117    self
ren v3217   newhd
ren v3300    race
ren v3256       y
ren v3051      ly
ren v3053     wly
ren v3060    tyhw
ren v3078   tyoth
ren v3081   troth
ren v3076    trhw
ren v3275    avhy
ren v3181   marit
ren v3301  weight
ren v3241    educ
gen food=.
gen fout=.

gen     ftax=v3064+v3080 if v3064!=. & v3080!=.
replace ftax=v3064       if v3064!=. & v3080==.
replace ftax=v3080       if v3064==. & v3080!=.
drop v3064 v3080
gen asset=.
gen trhf=.
gen year=1972
save          fam72,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f74,clear
rename *, lower
#delimit;
keep v3403 v3406 v3407 v3410 v3423 v3437 v3402 v3508 v3509 v3510 v3511 v3528 v3532 v3639 v3720 
     v3676 v3463 v3465 v3472 v3490 v3493 v3488
     v3695 v3598 v3721 v3663 v3441 v3445 v3476 v3492 v3431 v3443 v3665 v3630 v3666 v3419 v3417 v3458 v3505;
#delimit cr    
ren v3505 intdate     
ren v3458 wages
gen util=.
gen ncars=.
ren v3417 house
ren v3419 rent
ren v3666 disab
ren v3630 outkid
ren v3665 vet
ren v3443   fstmp
ren v3431   hourw
ren v3403   state
ren v3406    smsa
ren v3407   split
ren v3410    fchg
ren v3423   hours
ren v3437   fsize
ren v3402      id
ren v3508     age
ren v3510    agew
ren v3509     sex
ren v3511    kids
ren v3528   empst
ren v3532    self
ren v3639   newhd
ren v3720    race
ren v3676       y
ren v3463      ly
ren v3465     wly
ren v3472    tyhw
ren v3490   tyoth
ren v3493   troth
ren v3488    trhw
ren v3695    avhy
ren v3598   marit
ren v3721  weight
ren v3663    educ
ren v3441    food
ren v3445    fout

gen         ftax=v3476+v3492 if v3476!=. & v3492!=.
replace     ftax=v3476       if v3476!=. & v3492==.
replace     ftax=v3492       if v3476==. & v3492!=.
drop v3476 v3492
gen 	     asset=.
gen         trhf=.
gen year=1973
save          fam73,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f75,clear
rename *, lower
#delimit;
keep v3803 v3806 v3807 v3810 v3823 v3837 v3802 v3921 v3922 v3923 v3924 v3967 v3970 v4114 v4204 
     v4154 v3863 v3865 v3872 v3891 v3905 v3889
     v4174 v4053 v4224 v4198 v3841 v3853 v3876 v3893 v3831 v3851 v4140 v4076 v4145 v3819 v3817 v3931 v3858 v3918;
#delimit cr    
ren v3918 intdate     
ren v3858   wages
gen          util=.    
ren v3931   ncars
ren v3817   house
ren v3819    rent
ren v4145   disab
ren v4076  outkid
ren v4140     vet
ren v3851   fstmp
ren v3831   hourw
ren v3803   state
ren v3806    smsa
ren v3807   split
ren v3810    fchg
ren v3823   hours
ren v3837   fsize
ren v3802      id
ren v3921     age
ren v3923    agew
ren v3922     sex
ren v3924    kids
ren v3967   empst
ren v3970    self
ren v4114   newhd
ren v4204    race
ren v4154       y
ren v3863      ly
ren v3865     wly
ren v3872    tyhw
ren v3891   tyoth
ren v3905   troth
ren v3889    trhw
ren v4174    avhy
ren v4053   marit
ren v4224  weight
ren v4198    educ
ren v3841    food
ren v3853    fout

gen          ftax=v3876+v3893 if v3876!=. & v3893!=.
replace 	 ftax=v3876       if v3876!=. & v3893==.
replace 	 ftax=v3893       if v3876==. & v3893!=.
drop v3876 v3893
gen 		asset=.
gen 		 trhf=.
gen year=1974
save          fam74,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f76,clear
rename *, lower
#delimit;
keep v4303 v4306 v4307 v4310 v4332 v4346 v4302 v4436 v4437 v4438 v4439 v4458 v4461 v4658 v5096 
     v5029 v5031 v4379 v4386 v4406 v4419 v4404
     v5050 v4603 v5099 v5074 v4354 v4368 v4390 v4409 v5076 v4344 v4364 v4683 v4620 v4625 v4326 v4318 v4446 v4373 v4433;
#delimit cr    
ren v4433 intdate          
ren v4373 wages
gen 		 util=.   
ren v4446   ncars
ren v4318   house
ren v4326    rent
ren v4625   disab
ren v4620  outkid
ren v4683     vet
ren v4364   fstmp
ren v4344   hourw
ren v4303   state
ren v4306    smsa
ren v4307   split
ren v4310    fchg
ren v4332   hours
ren v4346   fsize
ren v4302      id
ren v4436     age
ren v4438    agew
ren v4437     sex
ren v4439    kids
ren v4458   empst
ren v4461    self
ren v4658   newhd
ren v5096    race
ren v5029       y
ren v5031      ly
ren v4379     wly
ren v4386    tyhw
ren v4406   tyoth
ren v4419   troth
ren v4404    trhw
ren v5050    avhy
ren v4603   marit
ren v5099  weight
ren v5074    educ
ren v4354    food
ren v4368    fout

gen          ftax=v4390+v4409  if v4390!=. & v4409!=.
replace 	 ftax=v4390        if v4390!=. & v4409==.
replace 	 ftax=v4409        if v4390==. & v4409!=.
replace 	 ftax=ftax-v5076 if v5076!=.
gen 		asset=.
gen 		 trhf=.
drop v4390 v4409 v5076
gen year=1975
save          fam75,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f77,clear
rename *, lower
#delimit;
keep v5203 v5206 v5207 v5210 v5232 v5254 v5202 v5350 v5351 v5352 v5353 v5373 v5376 v5578 v5662 
     v5626 v5627 v5289 v5297 v5318 v5332 v5316
     v5631 v5650 v5665 v5647 v5271 v5273 
     v5301 v5321 v5649 v5244 v5277 v5603 v5555 v5560 v5225 v5217 v5231 v5360 v5283 v5347;
#delimit cr    
ren v5347 intdate 
ren v5283   wages
ren v5231    util 
ren v5360   ncars
ren v5217   house
ren v5225    rent
ren v5560   disab
ren v5555  outkid
ren v5603     vet
ren v5277   fstmp
ren v5244   hourw
ren v5203   state
ren v5206    smsa
ren v5207   split
ren v5210    fchg
ren v5232   hours
ren v5254   fsize
ren v5202      id
ren v5350     age
ren v5352    agew
ren v5351     sex
ren v5353    kids
ren v5373   empst
ren v5376    self
ren v5578   newhd
ren v5662    race
ren v5626       y
ren v5627      ly
ren v5289     wly
ren v5297    tyhw
ren v5318   tyoth
ren v5332   troth
ren v5316    trhw
ren v5631    avhy
ren v5650   marit
ren v5665  weight
ren v5271    food
ren v5273    fout
ren v5647    educ

gen     ftax=v5301+v5321  if v5301!=. & v5321!=.
replace ftax=v5301        if v5301!=. & v5321==.
replace ftax=v5321        if v5301==. & v5321!=.
replace ftax=ftax-v5649 if v5649!=.
drop v5301 v5321 v5649
gen asset=.
gen trhf=.
gen year=1976
save          fam76,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f78,clear
rename *, lower
#delimit;
keep v5703 v5706 v5707 v5710 v5731 v5755 v5702 v5850 v5851 v5852 v5853 v5782 v5875 v6127 v6209 
     v6173 v6174 v5788 v5796 v5817 v5831 v5815
     v6178 v6197 v6212 v6194 v5770 v5772 v5800 v5820 v6196 v5743 v5815 v5831 v5776 v6152 v6097 v6102 v5723 v5717 v5729 v5860 v5872 v5847;
#delimit cr    
ren v5847 intdate         
ren v5782   wages
ren v5729    util 
ren v5860   ncars
ren v5717   house
ren v5723    rent
ren v6102   disab
ren v6097  outkid
ren v6152     vet
ren v5776   fstmp
ren v5815    trhw 
ren v5743   hourw
ren v5703   state
ren v5706    smsa
ren v5707   split
ren v5710    fchg
ren v5731   hours
ren v5755   fsize
ren v5702      id
ren v5850     age
ren v5852    agew
ren v5851     sex
ren v5853    kids
ren v5872   empst
ren v5875    self
ren v6127   newhd
ren v6209    race
ren v6173       y
ren v6174      ly
ren v5788     wly
ren v5796    tyhw
ren v5817   tyoth
ren v5831   troth
ren v6178    avhy
ren v6197   marit
ren v6212  weight
ren v6194    educ
ren v5770    food
ren v5772    fout
gen     ftax=v5800+v5820  if v5800!=. & v5820!=.
replace ftax=v5800        if v5800!=. & v5820==.
replace ftax=v5820        if v5800==. & v5820!=.
replace ftax=ftax-v6196 if v6196!=.
drop v5800 v5820 v6196
gen asset=.
gen year=1977
save          fam77,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f79,clear
rename *, lower
gen truncy	=v6766==99999												/*From now on, define {0,1} for topcoded income*/
gen trunca	=v6400==99999|v6401==99999|v6402==9999 |v6404==99999|v6400==-9999|v6402==-999		/*From now on, define {0,1} for topcoded financial income*/
#delimit;
keep v6303 v6306 v6307 v6310 v6336 v6361 v6302 v6462 v6463 v6464 v6465 v6492 v6493 v6724 v6802 
     v6766 v6767 v6398 v6408 v6428 v6442 v6426
     v6771 v6790 v6805 v6787 v6376 v6378 v6411 v6431 v6789 v6348
     v6400-v6406 v6430 v6426 v6442 v6382 v6749 v6690 v6710 v6326 v6319 v6332 v6473 v6391 v6495 v6496 v6422 v6440 v6459 trunc*;
#delimit cr
ren v6459 intdate
gen     aoth=0            if v6430==0
replace aoth=v6428*0.01   if v6430==1
replace aoth=v6428*0.0375 if v6430==2
replace aoth=v6428*0.085  if v6430==3
replace aoth=v6428*0.1125 if v6430==4
replace aoth=v6428*0.2    if v6430==5
replace aoth=v6428*0.3    if v6430==6
replace aoth=v6428*0.35   if v6430==7
replace aoth=v6428*0.5    if v6430==8
        
ren v6495 unionj  
ren v6496 unioni  
ren v6422  hhelp  
gen        whelp=.  
ren v6440  ohelp 

gen asset=v6400+v6401+v6402+v6403+v6404+v6405+v6406+aoth
drop      v6400-v6406 aoth
ren v6391   wages
ren v6332    util
ren v6473   ncars
ren v6319   house
ren v6326    rent
ren v6710   disab
ren v6690  outkid
ren v6749     vet
ren v6382   fstmp
ren v6348   hourw
ren v6303   state
ren v6306    smsa
ren v6307   split
ren v6310    fchg
ren v6336   hours
ren v6361   fsize
ren v6302      id
ren v6462     age
ren v6464    agew
ren v6463     sex
ren v6465    kids
ren v6492   empst
ren v6493    self
ren v6724   newhd
ren v6802    race
ren v6766       y
ren v6767      ly
ren v6398     wly
ren v6408    tyhw
ren v6428   tyoth
ren v6442   troth
ren v6426    trhw
ren v6771    avhy
ren v6790   marit
ren v6805  weight
ren v6787    educ
ren v6376    food
ren v6378    fout
gen     	 ftax=v6411+v6431  if v6411!=. & v6431!=.
replace 	 ftax=v6411        if v6411!=. & v6431==.
replace 	 ftax=v6431        if v6411==. & v6431!=.
replace 	 ftax=ftax-v6789 if v6789!=.
drop v6411 v6431 v6789 v6430
gen year=1978
save      fam78,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f80,clear
rename *, lower
gen truncy	=v7412==999999
gen trunca	=v6992==9999 |v6994==99999|v6990==-9999|v6991==-9999|v6992==-999
#delimit;
keep v6903 v6906 v6907 v6910 v6934 v6959 v6902 v7067 v7068 v7069 v7070 v7095 v7096 v7357 v7447 
     v7412 v7413 v6988 v6998 v7033 v7046 v7016
     v7417 v7435 v7451 v7433 v6972 v6974 v7399 v7411 v6946
     v6990-v6996 v7035 v7016 v7046 v6976 v7382 v7293 v7343 v6925 v6917 v6931 v7078 v6981 trunc* v7098 v7099 v7012 v7044 v7064;
#delimit cr
ren v7064 intdate        
gen     aoth=0            if v7035==0
replace aoth=v7033*0.01   if v7035==1
replace aoth=v7033*0.0375 if v7035==2
replace aoth=v7033*0.085  if v7035==3
replace aoth=v7033*0.1125 if v7035==4
replace aoth=v7033*0.2    if v7035==5
replace aoth=v7033*0.3    if v7035==6
replace aoth=v7033*0.35   if v7035==7
replace aoth=v7033*0.5    if v7035==8
gen asset=v6990+v6991+v6992+v6993+v6994+v6995+v6996+aoth
drop      v6990-v6996 aoth
ren v7098 unionj  
ren v7099 unioni  
ren v7012  hhelp  
gen        whelp=.  
ren v7044  ohelp 
ren v6981 	wages
ren v6931    util 
ren v7078   ncars
ren v6917 	house
ren v6925 	 rent
ren v7343 	disab
ren v7293  outkid
ren v7382     vet
ren v6976   fstmp
ren v6946   hourw
ren v6903   state
ren v6906    smsa
ren v6907   split
ren v6910    fchg
ren v6934   hours
ren v6959   fsize
ren v6902      id
ren v7067     age
ren v7069    agew
ren v7068     sex
ren v7070    kids
ren v7095   empst
ren v7096    self
ren v7357   newhd
ren v7447    race
ren v7412       y
ren v7413      ly
ren v6988     wly
ren v6998    tyhw
ren v7033   tyoth
ren v7046   troth
ren v7016    trhw
ren v7417    avhy
ren v7435   marit
ren v7451  weight
ren v7433    educ
ren v6972    food
ren v6974    fout
gen     	 ftax=v7399+v7411 if v7399!=. & v7411!=.
replace 	 ftax=v7399       if v7399!=. & v7411==.
replace 	 ftax=v7411       if v7399==. & v7411!=.
drop v7399 v7411 v7035
gen year=1979
save          fam79,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f81,clear
rename *, lower
gen truncy	=v8065==9999999
gen trunca	=v7586==99999|v7587==9999|v7584==-9999
#delimit;
keep v7503 v7506 v7507 v7510 v7530 v7551 v7502 v7658 v7659 v7660 v7661 v7706 v7707 v8009 v8099 
     v8065 v8066 v7580 v7590 v7625 v7638 v7608
     v8069 v8087 v8103 v8085 v7564 v7566 v8052 v8064 v7540
     v7582-v7588 v7627 v7608 v7638 v7568 v8034 v7993 v7974 v7525 v7517 v7699 v7669 v7573 trunc* v7709 v7710 v7604 v7636 v7655;
#delimit cr
ren v7655 intdate
ren v7709 unionj  
ren v7710 unioni  
ren v7604  hhelp  
gen        whelp=.  
ren v7636  ohelp 
gen     aoth=0            if v7627==0
replace aoth=v7625*0.01   if v7627==1
replace aoth=v7625*0.0375 if v7627==2
replace aoth=v7625*0.085  if v7627==3
replace aoth=v7625*0.1125 if v7627==4
replace aoth=v7625*0.2    if v7627==5
replace aoth=v7625*0.3    if v7627==6
replace aoth=v7625*0.35   if v7627==7
replace aoth=v7625*0.5    if v7627==8
gen asset=v7582+v7583+v7584+v7585+v7586+v7587+v7588+aoth
drop v7582-v7588 aoth    
ren v7573   wages
ren v7699    util 
ren v7669   ncars
ren v7517   house
ren v7525    rent
ren v7974   disab
ren v7993  outkid
ren v8034 	  vet
ren v7568   fstmp
ren v7540   hourw
ren v7503   state
ren v7506    smsa
ren v7507   split
ren v7510    fchg
ren v7530   hours
ren v7551   fsize
ren v7502      id
ren v7658     age
ren v7660    agew
ren v7659     sex
ren v7661    kids
ren v7706   empst
ren v7707    self
ren v8009   newhd
ren v8099    race
ren v8065       y
ren v8066      ly
ren v7580     wly
ren v7590    tyhw
ren v7625   tyoth
ren v7638   troth
ren v7608    trhw
ren v8069    avhy
ren v8087   marit
ren v8103  weight
ren v8085    educ
ren v7564    food
ren v7566    fout
gen     	 ftax=v8052+v8064 if v8052!=. & v8064!=.
replace 	 ftax=v8052       if v8052!=. & v8064==.
replace 	 ftax=v8064       if v8052==. & v8064!=.
drop v8052 v8064 v7627
gen year=1980
save      fam80,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f82,clear
rename *, lower
gen truncy	=v8689==999999
gen trunca	=v8275==99999|v8276==99999|v8277==99999|v8277==-9999|v8278==99999|v8279==99999|v8280==9999 |v8281==99999
#delimit;
keep v8203 v8206 v8207 v8210 v8228 v8249 v8202 v8352 v8353 v8354 v8355 v8374 v8375 v8633 v8723 
     v8689 v8690 v8273 v8283 v8318 v8331 v8301
     v8693 v8711 v8727 v8709 v8256 v8258 v8676 v8688 v8238
     v8275-v8281 v8320 v8301 v8331 v8260 v8658 v8612 v8616 v8221 v8217 v8225 v8358 v8265 trunc* v8377 v8378 v8297 v8329 v8349;
#delimit cr
ren v8349 intdate          
  
gen     aoth=0            if v8320==0
replace aoth=v8318*0.01   if v8320==1
replace aoth=v8318*0.0375 if v8320==2
replace aoth=v8318*0.085  if v8320==3
replace aoth=v8318*0.1125 if v8320==4
replace aoth=v8318*0.2    if v8320==5
replace aoth=v8318*0.3    if v8320==6
replace aoth=v8318*0.35   if v8320==7
replace aoth=v8318*0.5    if v8320==8
gen asset=v8275+v8276+v8277+v8278+v8279+v8280+v8281+aoth
drop v8275-v8281 aoth

ren v8377 unionj  
ren v8378 unioni  
ren v8297  hhelp  
gen        whelp=.  
ren v8329  ohelp 
ren v8265   wages
ren v8225    util 
ren v8358   ncars
ren v8217   house
ren v8221    rent
ren v8616   disab
ren v8612  outkid
ren v8658     vet
ren v8260   fstmp
ren v8238   hourw
ren v8203   state
ren v8206    smsa
ren v8207   split
ren v8210    fchg
ren v8228   hours
ren v8249   fsize
ren v8202      id
ren v8352     age
ren v8354    agew
ren v8353     sex
ren v8355    kids
ren v8374   empst
ren v8375    self
ren v8633   newhd
ren v8723    race
ren v8689       y
ren v8690      ly
ren v8273     wly
ren v8283    tyhw
ren v8318   tyoth
ren v8331   troth
ren v8301    trhw
ren v8693    avhy
ren v8711   marit
ren v8727  weight
ren v8709    educ
ren v8256    food
ren v8258    fout
gen     	 ftax=v8676+v8688 if v8676!=. & v8688!=.
replace 	 ftax=v8676       if v8676!=. & v8688==.
replace 	 ftax=v8688       if v8676==. & v8688!=.
drop v8676 v8688 v8320
gen year=1981
save          fam81,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f83,clear
rename *, lower
gen truncy	=v9375==9999999
gen trunca	=v8883==999999|v8884==999999|v8885==99999 |v8885==-9999 |v8886==99999 |v8887==999999|v8888==99999 |v8889==99999
#delimit;
keep v8803 v8806 v8807 v8810 v8830 v8851 v8802 v8961 v8962 v8963 v8964 v9005 v9006 v9319 v9408 
     v9375 v9376 v8881 v8891 v8926 v8939 v8909
     v9379 v9419 v9433 v9395 v8864 v8866 v9362 v9374 v8840
     v8883-v8889 v8928 v8909 v8939 v8868 v9344 v9303 v9290 v8825 v8817 v8998 v8967 v8873 trunc* v9008 v9009 v8905 v8937 v8958;
#delimit cr
ren v8958 intdate     

ren v9008 unionj  
ren v9009 unioni  
ren v8905  hhelp  
gen        whelp=.  
ren v8937  ohelp 
gen asset=v8883+v8884+v8885+v8886+v8887+v8888+v8889+v8928
drop v8883-v8889 v8928  
ren v8873   wages
ren v8998    util
ren v8967   ncars
ren v8817   house
ren v8825    rent
ren v9290   disab
ren v9303  outkid
ren v9344     vet
ren v8868   fstmp
ren v8840   hourw
ren v8803   state
ren v8806    smsa
ren v8807   split
ren v8810    fchg
ren v8830   hours
ren v8851   fsize
ren v8802      id
ren v8961     age
ren v8963    agew
ren v8962     sex
ren v8964    kids
ren v9005   empst
ren v9006    self
ren v9319   newhd
ren v9408    race
ren v9375       y
ren v9376      ly
ren v8881     wly
ren v8891    tyhw
ren v8926   tyoth
ren v8939   troth
ren v8909    trhw
ren v9379    avhy
ren v9419   marit
ren v9433  weight
ren v9395    educ
ren v8864    food
ren v8866    fout
gen          ftax=v9362+v9374 if v9362!=. & v9374!=.
replace      ftax=v9362       if v9362!=. & v9374==.
replace      ftax=v9374       if v9362==. & v9374!=.
drop v9362 v9374
gen year=1982
save      fam82,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f84,clear
rename *, lower
gen truncy	=v11022==999999
gen trunca	=v10265==999999|v10266==999999|v10267==99999 |v10268==99999 |v10270==999999|v10272==999999|v10274==99999 |v10275==999999|v10276==-9999
#delimit;
keep v10003 v10006 v10007 v10010 v10037 v10222 v10002 v10419 v10420 v10421 v10422 v10453 v10456 v10966 v11055 
     v11022 v11023 v10263 v10277 v10382 v10396 v10305
     v11026 v11065 v11079 v11042 v10235 v10237 v11009 v11021 v10131
     v10265-v10268 v10270 v10272 v10274 v10275 v10384 v10305 v10396 v10239 v10991 v10894 v10879 v10026 v10018 v10030 v10430 v10256 trunc*  
     v10458 v10459 v10300 v10394 v10416;
#delimit cr
ren v10416 intdate

ren v10458 unionj 
ren v10459 unioni 
ren v10300  hhelp 
gen         whelp=. 
ren v10394  ohelp

gen asset=v10265+v10266+v10267+v10268+v10270+v10272+v10274+v10275+v10384
drop v10265-v10268 v10270 v10272 v10274 v10275 v10384
ren v10256   wages
ren v10030    util 
ren v10430   ncars
ren v10018   house
ren v10026    rent
ren v10879   disab
ren v10894  outkid
ren v10991     vet
ren v10239   fstmp
ren v10131   hourw
ren v10003   state
ren v10006    smsa
ren v10007   split
ren v10010    fchg
ren v10037   hours
ren v10222   fsize
ren v10002      id
ren v10419     age
ren v10421    agew
ren v10420     sex
ren v10422    kids
ren v10453   empst
ren v10456    self
ren v10966   newhd
ren v11055    race
ren v11022       y
ren v11023      ly
ren v10263     wly
ren v10277    tyhw
ren v10382   tyoth
ren v10396   troth
ren v10305    trhw
ren v11026    avhy
ren v11065   marit
ren v11079  weight
ren v11042    educ
ren v10235    food
ren v10237    fout
gen           ftax=v11009+v11021 if v11009!=. & v11021!=.
replace       ftax=v11009       if v11009!=. & v11021==.
replace       ftax=v11021       if v11009==. & v11021!=.
drop v11009 v11021
gen year=1983
save      fam83,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f85,clear
rename *, lower
gen truncy	=v12371==999999
gen trunca	=v11406==999999|v11407==999999|v11408==99999 |v11406==-99999|v11407==-99999|v11408==-9999 |v11409==99999 |v11412==999999|v11414==999999| ///		
		 v11416==999999|v11417==999999
#delimit;
keep v11103 v11106 v11107 v11112 v11146 v11364 v11102 v11606 v11607 v11608 v11609 v11637 v11640 v11906 v11938 
     v12371 v12372 v11404 v11419 v11561 v11577 v11461
     v12377 v12426 v12446 v12400 v11375 v11377 v12358 v12370 v11258
     v11406-v11409 v11412  v11414  v11416 v11417 v11563 v11461 v11577 v11379 v11940 v11896 v11993 v11133 v11125 v11137 v11635 v11397 trunc* 
     v11649 v11650 v11445 v11457 v11574 v11600;
#delimit cr
ren v11600 intdate

ren v11649 unionj 
ren v11650 unioni 
ren v11445  hhelp 
ren v11457  whelp 
ren v11574  ohelp

gen asset=v11406+v11407+v11408+v11409+v11412+v11414+v11416+v11417+v11563
drop 	v11406-v11409 v11412  v11414  v11416 v11417 v11563
ren v11397   wages
ren v11137    util 
ren v11635   ncars
ren v11125   house
ren v11133    rent
ren v11993   disab
ren v11896  outkid
ren v11940     vet
ren v11379   fstmp
ren v11258   hourw
ren v11103   state
ren v11106    smsa
ren v11107   split
ren v11112    fchg
ren v11146   hours
ren v11364   fsize
ren v11102      id
ren v11606     age
ren v11608    agew
ren v11607     sex
ren v11609    kids
ren v11637   empst
ren v11640    self
ren v11906   newhd
ren v11938    race
ren v12371       y
ren v12372      ly
ren v11404     wly
ren v11419    tyhw
ren v11561   tyoth
ren v11577   troth
ren v11461    trhw
ren v12377    avhy
ren v12426   marit
ren v12446  weight
ren v12400    educ
ren v11375    food
ren v11377    fout
gen           ftax=v12358+v12370 if v12358!=. & v12370!=.
replace       ftax=v12358       if v12358!=. & v12370==.
replace       ftax=v12370       if v12358==. & v12370!=.
drop v12358 v12370
gen year=1984
save           fam84,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f86,clear
rename *, lower
gen truncy	=v13623==999999
gen trunca	=v12806==-99999|v12807==-9999
#delimit;
keep v12503 v12506 v12507 v12510 v12545 v12763 v12502 v13011 v13012 v13013 v13014 v13046 v13049 v13533 v13565 
     v13623 v13624 v12803 v12818 v12968 v12984 v12868
     v13629 v13665 v13687 v13640 v12774 v12776 v13610 v13622 v12657
     v12805-v12808  v12811  v12813  v12815 v12816 v12970 v12868 v12984 v12778 v13567 v13407 v13427 v12532 v12524 v12536 v13044 v12796 trunc*
	 v13052 v13053 v12844 v12864 v12981 v13008;
#delimit cr
ren v13008 intdate

ren v13052 unionj 
ren v13053 unioni 
ren v12844  hhelp 
ren v12864  whelp 
ren v12981  ohelp
 
gen asset=v12805+v12806+v12807+v12808+v12811+v12813+v12815+v12816+v12970
drop v12805-v12808  v12811  v12813  v12815 v12816 v12970
ren v12796   wages
ren v12536    util 
ren v13044   ncars
ren v12524   house
ren v12532    rent
ren v13427   disab
ren v13407  outkid
ren v13567     vet
ren v12778   fstmp
ren v12657   hourw
ren v12503   state
ren v12506    smsa
ren v12507   split
ren v12510    fchg
ren v12545   hours
ren v12763   fsize
ren v12502      id
ren v13011     age
ren v13013    agew
ren v13012     sex
ren v13014    kids
ren v13046   empst
ren v13049    self
ren v13533   newhd
ren v13565    race
ren v13623       y
ren v13624      ly
ren v12803     wly
ren v12818    tyhw
ren v12968   tyoth
ren v12984   troth
ren v12868    trhw
ren v13629    avhy
ren v13665   marit
ren v13687  weight
ren v13640    educ
ren v12774    food
ren v12776    fout
gen           ftax=v13610+v13622 if v13610!=. & v13622!=.
replace       ftax=v13610       if v13610!=. & v13622==.
replace       ftax=v13622       if v13610==. & v13622!=.
drop v13610 v13622
gen year=1985
save           fam85,replace



u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f87,clear
rename *, lower
gen truncy=v14670==999999
gen trunca=v13909==-9999
#delimit;
keep v13703 v13706 v13707 v13710 v13745 v13867 v13702 v14114 v14115 v14116 v14117 v14146 v14149 v14580 v14612 
     v14670 v14671 v13905 v13920 v14070 v14086 v13970
     v14676 v14712 v14737 v14687 v13876 v13878 v14657 v14669 v13809
     v13907-v13910 v13913 v13915 v13917 v13918 v14072 v13970 v14086 v13880 v14614 v14504 v14515 v13732 v13724 v13736 v13898 trunc*
     v14152 v14153 v13946 v13966 v14083 v14111;
#delimit cr
ren v14111 intdate

ren v14152 unionj 
ren v14153 unioni 
ren v13946  hhelp 
ren v13966  whelp 
ren v14083  ohelp

gen asset=v13907+v13908+v13909+v13910+v13913+v13915+v13917+v13918+v14072
drop v13907 v13908 v13909 v13910 v13913 v13915 v13917 v13918 v14072
gen  ncars=.
ren v13898   wages
ren v13736    util
ren v13724   house
ren v13732    rent
ren v14515   disab
ren v14504  outkid
ren v14614     vet
ren v13880   fstmp
ren v13809   hourw
ren v13703   state
ren v13706    smsa
ren v13707   split
ren v13710    fchg
ren v13745   hours
ren v13867   fsize
ren v13702      id
ren v14114     age
ren v14116    agew
ren v14115     sex
ren v14117    kids
ren v14146   empst
ren v14149    self
ren v14580   newhd
ren v14612    race
ren v14670       y
ren v14671      ly
ren v13905     wly
ren v13920    tyhw
ren v14070   tyoth
ren v14086   troth
ren v13970    trhw
ren v14676    avhy
ren v14712   marit
ren v14737  weight
ren v14687    educ
ren v13876    food
ren v13878    fout
gen           ftax=v14657+v14669 if v14657!=. & v14669!=.
replace       ftax=v14657        if v14657!=. & v14669==.
replace       ftax=v14669        if v14657==. & v14669!=.
drop v14657 v14669
gen year=1986
save           fam86,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f88,clear
rename *, lower
gen truncy=v16144==9999999
gen trunca=v14924==-9999
#delimit;
keep v14803 v14806 v14807 v14810 v14835 v14889 v14802 v15130 v15131 v15132 v15133 v15154 v15157 v16054 v16086 
     v16144 v16145 v14920 v14935 v15085 v15101 v14985
     v16150 v16187 v16208 v16161 v16131 v16143 v14865
     v14922-v14925 v14928 v14930 v14932 v14933 v15087 v14985 v15101 v14895 v16088 v15779 v15994 v14824 v14913 trunc*
     v15160 v15161 v14961 v14981 v15098 v15127;
#delimit cr
ren v15127 intdate

ren v15160 unionj
ren v15161 unioni 
ren v14961  hhelp 
ren v14981  whelp 
ren v15098  ohelp

gen asset=v14922+v14923+v14924+v14925+v14928+v14930+v14932+v14933+v15087
drop v14922-v14925 v14928 v14930 v14932 v14933 v15087
gen  ncars=.
gen   util=.
gen    rent=.
ren v14913   wages
ren v14824   house
ren v15994   disab
ren v15779  outkid
ren v16088     vet
ren v14895   fstmp
ren v14865   hourw
ren v14803   state
ren v14806    smsa
ren v14807   split
ren v14810    fchg
ren v14835   hours
ren v14889   fsize
ren v14802      id
ren v15130     age
ren v15132    agew
ren v15131     sex
ren v15133    kids
ren v15154   empst
ren v15157    self
ren v16054   newhd
ren v16086    race
ren v16144       y
ren v16145      ly
ren v14920     wly
ren v14935    tyhw
ren v15085   tyoth
ren v15101   troth
ren v14985    trhw
ren v16150    avhy
ren v16187   marit
ren v16208  weight
ren v16161    educ
gen 		  food=.
gen 		  fout=.
gen           ftax=v16131+v16143 if v16131!=. & v16143!=.
replace       ftax=v16131        if v16131!=. & v16143==.
replace       ftax=v16143        if v16131==. & v16143!=.
drop v16131 v16143
gen year=1987
save           fam87,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f89,clear
rename *, lower
gen truncy=v17533==9999999
gen trunca=v16424==-99999
#delimit;
keep v16303 v16306 v16307 v16310 v16335 v16389 v16302 v16631 v16632 v16633 v16634 v16655 v16658 v17451 v17483 
     v17533 v17534 v16420 v16435 v16585 v16601 v16485 
     v17536 v17565 v17612 v17545 v17528 v17532 v16365
     v16422-v16425 v16428 v16430 v16432 v16433 v16587 v16485 v16601 v16395 v17485 v17308 v17391 v16324 v16413 trunc* 
	 v16661 v16662 v16461 v16481 v16598 v16628;
#delimit cr
ren v16628 intdate

ren v16661 unionj 
ren v16662 unioni 
ren v16461  hhelp 
ren v16481  whelp 
ren v16598  ohelp

gen asset=v16422+v16423+v16424+v16425+v16428+v16430+v16432+v16433+v16587
drop v16422-v16425 v16428 v16430 v16432 v16433 v16587
gen    rent=.
gen  ncars=.
gen   util=.
ren v16413   wages
ren v16324   house
ren v17391   disab
ren v17308  outkid
ren v17485     vet
ren v16395   fstmp
ren v16365   hourw
ren v16303   state
ren v16306    smsa
ren v16307   split
ren v16310    fchg
ren v16335   hours
ren v16389   fsize
ren v16302      id
ren v16631     age
ren v16633    agew
ren v16632     sex
ren v16634    kids
ren v16655   empst
ren v16658    self
ren v17451   newhd
ren v17483    race
ren v17533       y
ren v17534      ly
ren v16420     wly
ren v16435    tyhw
ren v16585   tyoth
ren v16601   troth
ren v16485    trhw
ren v17536    avhy
ren v17565   marit
ren v17612  weight
ren v17545    educ
gen food=.
gen fout=.
gen          ftax=v17528+v17532 if v17528!=. & v17532!=.
replace      ftax=v17528       if v17528!=. & v17532==.
replace      ftax=v17532       if v17528==. & v17532!=.
drop v17528 v17532
gen year=1988
save           fam88,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f90,clear
rename *, lower
gen truncy=v18875==9999999
gen trunca=v17838==999999
#delimit;
keep v17703 v17706 v17707 v17710 v17744 v17798 v17702 v18049 v18050 v18051 v18052 v18093 v18096 v18782 v18814 
     v18875 v18878 v17836 v17851 v18001 v18017 v17901 
     v18887 v18916 v18943 v18898 v17807 v17809 v18862 v18874 v17774
     v17838-v17841 v17844 v17846 v17848 v17849 v18003 v17901 v18017 v17811 v18816 v18712 v18722 v17733 v17724 v17829 trunc*
	 v18099 v18100 v17877 v17897 v18014 v18046;
#delimit cr
ren v18046 intdate

ren v18099 unionj 
ren v18100 unioni 
ren v17877  hhelp 
ren v17897  whelp 
ren v18014  ohelp

gen asset=v17838+v17839+v17840+v17841+v17844+v17846+v17848+v17849+v18003 
drop v17838-v17841 v17844 v17846 v17848 v17849 v18003
gen  ncars=.
gen   util=.
ren v17829   wages
ren v17724   house
ren v17733    rent
ren v18722   disab
ren v18712  outkid
ren v18816     vet
ren v17811   fstmp
ren v17774   hourw
ren v17703   state
ren v17706    smsa
ren v17707   split
ren v17710    fchg
ren v17744   hours
ren v17798   fsize
ren v17702      id
ren v18049     age
ren v18051    agew
ren v18050     sex
ren v18052    kids
ren v18093   empst
ren v18096    self
ren v18782   newhd
ren v18814    race
ren v18875       y
ren v18878      ly
ren v17836     wly
ren v17851    tyhw
ren v18001   tyoth
ren v18017   troth
ren v17901    trhw
ren v18887    avhy
ren v18916   marit
ren v18943  weight
ren v18898    educ
ren v17807    food
ren v17809    fout
replace v18862=. if v18862==999999
replace v18874=. if v18874==999999
gen          ftax=v18862+v18874 if v18862!=. & v18874!=.
replace      ftax=v18862        if v18862!=. & v18874==.
replace      ftax=v18874        if v18862==. & v18874!=.
drop v18862 v18874
gen year=1989
save           fam89,replace

u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f91,clear
rename *, lower
gen truncy=v20175==999999
gen trunca=v19139==-99999|v19140==-9999
#delimit;
keep v19003 v19006 v19007 v19010 v19044 v19098 v19002 v19349 v19350 v19351 v19352 v19393 v19396 v20082 v20114 
     v20175 v20178 v19136 v19151 v19301 v19317 v19201 
     v20187 v20216 v20243 v20198 v19107 v19109 v20162 v20174 v19074
     v19138-v19141 v19144 v19146 v19148 v19149 v19303 v19201 v19317 v19111 v20116 v20012 v20022 v19033 v19024 v19129 trunc* 
	 v19399 v19400 v19177 v19197 v19314 v19346;
#delimit cr
ren v19346 intdate

ren v19399 unionj 
ren v19400 unioni 
ren v19177  hhelp 
ren v19197  whelp 
ren v19314  ohelp

gen asset=v19138+v19139+v19140+v19141+v19144+v19146+v19148+v19149+v19303
drop v19138-v19141 v19144 v19146 v19148 v19149 v19303
gen  ncars=.
gen   util=.
ren v19129   wages
ren v19024   house
ren v19033    rent
ren v20022   disab
ren v20012  outkid
ren v20116     vet
ren v19111   fstmp
ren v19074   hourw
ren v19003   state
ren v19006    smsa
ren v19007   split
ren v19010    fchg
ren v19044   hours
ren v19098   fsize
ren v19002      id
ren v19349     age
ren v19351    agew
ren v19350     sex
ren v19352    kids
ren v19393   empst
ren v19396    self
ren v20082   newhd
ren v20114    race
ren v20175       y
ren v20178      ly
ren v19136     wly
ren v19151    tyhw
ren v19301   tyoth
ren v19317   troth
ren v19201    trhw
ren v20187    avhy
ren v20216   marit
ren v20243  weight
ren v20198    educ
ren v19107    food
ren v19109    fout
gen           ftax=v20162+v20174 if v20162!=. & v20174!=.
replace       ftax=v20162       if v20162!=. & v20174==.
replace       ftax=v20174       if v20162==. & v20174!=.
drop v20162 v20174
gen year=1990
save           fam90,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f92,clear
rename *, lower
gen truncy=v21481==9999999
gen trunca=v20440==-9999
#delimit;
keep v20303 v20306 v20307 v20310 v20344 v20398 v20302 v20651 v20652 v20653 v20654 v20693 v20696 v21388 v21420 
     v21481 v21484 v20436 v20451 v20601 v20617 v20501 
     v21493 v21522 v21547 v21504 v20407 v20409 v20374
     v20438-v20441 v20444 v20446 v20448 v20449 v20603 v20501 v20617 v20411 v21422 v21312 v21322 v20333 v20324 v20429 trunc* 
     v20699 v20700 v20477 v20497 v20614 v20648;
#delimit cr
ren v20648 intdate
ren v20699 unionj 
ren v20700 unioni 
ren v20477  hhelp 
ren v20497  whelp 
ren v20614  ohelp

gen asset=v20438+v20439+v20440+v20441+v20444+v20446+v20448+v20449+v20603
drop v20438-v20441 v20444 v20446 v20448 v20449 v20603
gen  ncars=.
gen   util=.
ren v20429   wages
ren v20324   house
ren v20333    rent
ren v21322   disab
ren v21312  outkid
ren v21422     vet
ren v20411   fstmp
ren v20374   hourw
ren v20303   state
ren v20306    smsa
ren v20307   split
ren v20310    fchg
ren v20344   hours
ren v20398   fsize
ren v20302      id
ren v20651     age
ren v20653    agew
ren v20652     sex
ren v20654    kids
ren v20693   empst
ren v20696    self
ren v21388   newhd
ren v21420    race
ren v21481       y
ren v21484      ly
ren v20436     wly
ren v20451    tyhw
ren v20601   tyoth
ren v20617   troth
ren v20501    trhw
ren v21493    avhy
ren v21522   marit
ren v21547  weight
ren v21504    educ
ren v20407    food
ren v20409    fout
gen           ftax=0
gen year=1991
save           fam91,replace


u C:\Users\edmun\OneDrive\Documents\Research\PSIDdata\DataForBPP2008\f93,clear
rename *, lower
gen truncy=0
gen trunca=0
#delimit;
keep v21603 v21605 v21606 v21608 v21634 v22405 v21602 v22406 v22407 v22408 v22409 v22448 v22451 v23245 v23276 	 
     v23322 v23323 v23324 v21959 v22373 v22397 v22366 
     v23336 v23361 v23333 v21707 v21711 v21670
     v21809 v21813 v21810 v21814 v21811 v21812 v21829 v21845 v21861 v21893 v21877 v21909 v21925 v21957 v21941 v22371
     v22366 v22397 v21727 v23278 v23171 v23181 v21622 v21610 v21739 trunc*
	 v22454 v22455 v22156 v22333 v22392 v22403;
#delimit cr
ren v22403 intdate
ren v22454 unionj 
ren v22455 unioni 
ren v22156  hhelp 
ren v22333  whelp 
ren v22392  ohelp

gen asset=v21809+v21813+v21810+v21814+v21811+v21812+v21829+v21845+v21861+v21893+v21877+v21909+v21925+v21957+v21941+v22371
drop v21809 v21813 v21810 v21814 v21811 v21812 v21829 v21845 v21861 v21893 v21877 v21909 v21925 v21957 v21941 v22371
gen  ncars=.
gen   util=.
ren v21739   wages
ren v21610   house
ren v21622    rent
ren v23181   disab
ren v23171  outkid
ren v23278     vet
ren v21727   fstmp
ren v21670   hourw
ren v21603   state
ren v21605    smsa
ren v21606   split
ren v21608    fchg
ren v21634   hours
ren v22405   fsize
ren v21602      id
ren v22406     age
ren v22408    agew
ren v22407     sex
ren v22409    kids
ren v22448   empst
ren v22451    self
ren v23245   newhd
ren v23276    race
ren v23322       y
ren v23323      ly
ren v23324     wly
ren v21959    tyhw
ren v22373   tyoth
ren v22397   troth
ren v22366    trhw
gen avhy=ly/hours
ren v23336   marit
ren v23361  weight
ren v23333    educ
ren v21707    food
ren v21711    fout
gen          ftax=0
gen year=1992
save           fam92,replace


use fam67, clear
cap program drop doit
program def doit
            local i=68
            while `i' < 93 {
            append using fam`i'
            local i=`i'+1
                  }
end
qui doit
// Restructure intdate to nine date categories
replace intdate=1 if intdate>100 & intdate<315
replace intdate=2 if intdate>=315 & intdate<329
replace intdate=3 if intdate>=329 & intdate<419
replace intdate=4 if intdate>=419 & intdate<503
replace intdate=5 if intdate>=503 & intdate<518
replace intdate=6 if intdate>=518 & intdate<531
replace intdate=7 if intdate>=531 & intdate<701
replace intdate=8 if intdate>=701 & intdate<1231
replace intdate=9 if intdate>=1231
save fam, replace

cap program drop doit
program def doit
            local i=67
            while `i' < 93 {
            erase fam`i'.dta
            local i=`i'+1
                  }
end
qui doit

log close
clear
