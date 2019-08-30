//Neil Davies 20/10/17
//This estimates the effects of staying in school till 18 and 21.


use  "workingdata/full_merged_dataset_SNPs",clear

gen eduyear_16=(eduyears2 >=16) if eduyears2 !=.
gen eduyear_18=(eduyears2 >=18) if eduyears2 !=.
gen eduyear_20=(eduyears2 >=20) if eduyears2 !=.
gen eduyear_21=(eduyears2 >=21) if eduyears2 !=.

//IVREG2 analysis
ivreg2 cov_male (eduyear_18 eduyear_21=rs9537821-rs10061788) cov_male imob_* yob_* pc_* [pweight=weight],cluster(mobi) endog(eduyear_18 eduyear_21) ro partial(cov_male imob_* yob_* pc_*) first
regsave eduyear_18 eduyear_21 using "results/ivreg2_bivar_outcome",detail(all) pval replace
ds out_*
foreach i in `r(varlist)'{
	ivreg2 `i' (eduyear_18 eduyear_21=rs9537821-rs10061788) cov_male imob_* yob_* pc_* [pweight=weight],cluster(mobi) endog(eduyear_18 eduyear_21) partial(cov_male imob_* yob_* pc_*)
	regsave eduyear_18 eduyear_21 using "results/ivreg2_bivar_outcome",detail(all) pval append
	}

/*	
//IVREG2 analysis
ivreg2 cov_male (eduyear_16 eduyear_18=rs9537821-rs10061788) cov_male imob_* yob_* pc_* [pweight=weight],cluster(mobi) endog(eduyear_16 eduyear_18) ro partial(cov_male imob_* yob_* pc_*) first
regsave eduyear_18 eduyear_21 using "results/ivreg2_bivar_16_outcome",detail(all) pval replace
ds out_exsmoke
foreach i in `r(varlist)'{
	ivreg2 `i' (eduyear_16 eduyear_18=rs9537821-rs10061788) cov_male imob_* yob_* pc_* [pweight=weight],cluster(mobi) endog(eduyear_16 eduyear_18) partial(cov_male imob_* yob_* pc_*)
	regsave eduyear_18 eduyear_21 using "results/ivreg2_bivar_16_outcome",detail(all) pval append
	}
*/

use "results/ivreg2_bivar_outcome",clear

//Reorder
gen order=.
#delimit ;
foreach i in 
out_highbloodpressure
out_diabetes
out_stroke
out_heartattack
out_depression
out_cancer
out_dead
out_exsmoker
out_smoker
out_income_under_18k
out_income_over_31k
out_income_over_52k
out_income_over_100k
out_gripstrength
out_arterial_stiffness
out_height
out_bmi
out_dia_bp
out_sys_bp
out_intell
out_happiness
out_alcohol
out_sedentary
out_phys_m_act
out_phys_v_act{;
	local j =`j'+1;
	replace order =`j' if depvar=="`i'";
	};
#delimit cr
sort var cmd order 	

gen cont=0
#delimit ;
foreach i in 
out_gripstrength
out_arterial_stiffness
out_height
out_bmi
out_dia_bp
out_sys_bp
out_intell
out_happiness
out_alcohol
out_sedentary
out_phys_m_act
out_phys_v_act{;
	replace cont =1 if depvar=="`i'";
	};
#delimit cr
drop if order==.
//Rescale 


//Plot results using coefplot
//Genetic scores
mkmat coef stderr if var=="eduyear_18" & cont==0,matrix(results) rownames(order)
matrix results_edu_18=results'

mkmat coef stderr if var=="eduyear_21"& cont==0 ,matrix(results) rownames(order)
matrix results_edu_21=results'

#delimit ;
coefplot (matrix(results_edu_18) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(+0.1))
		 (matrix(results_edu_21) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100",size(small))
	xlabel(-100(20)100,noticks) rescale(100)  legend(off)
	headings(1 = "{bf:Morbidity}"
			 7 = "{bf:Mortality}"
			 8 = "{bf:Health behaviours}"
			 10 = "{bf:Income}"
			 , labsize(small)) 
	coeflabels(1="Hypertension"
2="Diabetes"
3="Stroke"
4="Heart attack"
5="Depression"
6="Cancer"
7="Died"
8="Ever smoked"
9="Currently smoke"
10="Income over £18k"
11="Income over £31k"
12="Income over £52k"
13=".                                 Income over £100k", wrap(51));

graph save "results/figure_X_bin_individual_level_ivreg_effect_het", replace;
#delimit cr

//Cont outcomes
mkmat coef stderr if var=="eduyear_18" & cont==1,matrix(results) rownames(order)
matrix results_edu_18=results'

mkmat coef stderr if var=="eduyear_21"& cont==1 ,matrix(results) rownames(order)
matrix results_edu_21=results'

#delimit ;
coefplot (matrix(results_edu_18) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(+0.1))
		 (matrix(results_edu_21) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
	legend(order(2 "Educated to age 18" 4 "Educated to age 21" ) row(1) size(vsmall))
	/*xlabel(-3(1)3,noticks)*/  xtick(none) ytick(none) 
	headings(14 = "{bf:Indicators of aging}"
			 16 = "{bf:Anthropometry}"
			 18 = "{bf:Blood pressure}"
			 20 = "{bf:Neurocognitive}"
			 22 = "{bf:Health behaviours}"
			 , labsize(small)) 
	coeflabels(14="Gripstrength (kg)"
15="Arterial stiffness"
16="Height (cm)"
17="BMI (kg/m2)"
18="Diastolic (mmHg)"
19="Systolic (mmHg)"
20="Intelligence (0 to 13)"
21="Happiness (0 to 5 Likert)"
22="Alcohol consumption (0 low to 5 high)"
23="Hours watching television per day"
24="Moderate exercise (days/week)"
25="Vigorous exercise (days/week)", wrap(40));	
#delimit cr


graph save "results/figure_X_cont_individual_level_ivreg_effect_het", replace 

graph combine "results/figure_X_bin_individual_level_ivreg_effect_het" "results/figure_X_cont_individual_level_ivreg_effect_het", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/figure_X_bin_individual_level_ivreg_effect_het.eps", as(pdf) replace fontface("Calibri")
