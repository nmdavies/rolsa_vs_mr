//Neil Davies
//This creates the descriptive data for Table 1

//Neil Davies 23/05/17
//This cleans the education data in UK Biobank to derived ISCED levels
//This estimates the associations between the Okbay SNPs and each of the outcomes in Biobank study

use  "workingdata/full_merged_dataset",clear
tabstat  cov_male

//Get age at recruitment
joinby n_eid using  "rawdata/biobank_phenotypes_nmd_150417.dta",
tabstat cov_age

//Get proportion with degree
gen out_degree=0 
gen out_post_16=0
gen out_any=0

forvalues i=0(1)5{
	replace out_degree=1 if n_6138_0_0==1|n_6138_0_0==1|n_6138_0_0==6|n_6138_0_0==6|n_6138_0_0==5|n_6138_0_0==5
	replace out_degree=1 if n_6138_1_0==1|n_6138_1_0==1|n_6138_1_0==6|n_6138_1_0==6|n_6138_0_0==5|n_6138_0_0==5
	replace out_post_16=1 if n_6138_0_`i'==2| n_6138_0_`i'==5| n_6138_0_`i'==6|n_6138_0_`i'==1
	replace out_post_16=1 if n_6138_1_`i'==2| n_6138_1_`i'==5| n_6138_1_`i'==6|n_6138_0_`i'==1
	replace out_any=1 if n_6138_0_`i'==1|n_6138_0_`i'==2| n_6138_0_`i'==5| n_6138_0_`i'==6|n_6138_0_`i'==4|n_6138_0_`i'==3
	replace out_any=1 if n_6138_0_`i'==1|n_6138_1_`i'==2| n_6138_1_`i'==5| n_6138_1_`i'==6|n_6138_1_`i'==4|n_6138_1_`i'==3
	}
egen max=rowmax(n_6138*)	
replace out_degree=. if max==-3 & n_6138_0_0!=-7 &  n_6138_1_0!=-7
replace out_post_16=. if  max==-3 & n_6138_0_0!=-7 &  n_6138_1_0!=-7

drop  if eduyears2==.

tab out_degree 
tab out_post_16
tab out_any

tabstat out_highbloodpressure out_diabetes out_stroke out_heartattack out_depression out_cancer out_dead out_smoker out_exsmoker out_income_under_18k out_income_over_31k out_income_over_52k out_income_over_100k if eduyears!=., stats(n mean) c(s)
tabstat out_gripstrength out_arterial_stiffness out_height out_bmi out_dia_bp out_sys_bp out_intell out_happiness out_alcohol out_sedentary out_phys_m_act out_phys_v_act if eduyears!=. , stats(n mean sd) c(s)
tabstat  cov_male if eduyears2!=., stats(N mean count)
tabstat yob eduyears2 if eduyears2!=., stats(N mean sd) c(s)
