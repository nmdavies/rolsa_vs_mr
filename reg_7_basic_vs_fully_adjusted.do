//Neil Davies 19/06/17
//This sensitivity analysis runs 

use "workingdata/full_merged_dataset",clear

//Drop current location (not a formal confounder as could be an effect of education)
drop cov_current_*

//IVREG2 analysis
ivreg2 cov_male (eduyears2=allele_score_21) imob_* yob_* pc_* cov_matsmoking cov_birthweight cov_dist_lon cov_birth_location_imd ///
	cov_birth_location_urban cov_birth_location_northing cov_birth_location_easting [pweight=weight_hughes],cluster(mobi) endog(eduyears2) ///
	partial(imob_* yob_* pc_* cov_matsmoking cov_birthweight cov_dist_lon cov_birth_location_imd ///
	cov_birth_location_urban cov_birth_location_northing cov_birth_location_easting)
regsave eduyears2 using "results/ivreg2_outcome_full_adj",detail(all) pval replace
ivreg2 cov_male (eduyears2=allele_score_21) imob_* yob_* pc_*  [pweight=weight_hughes] if e(sample),cluster(mobi) endog(eduyears2) ///
	partial(imob_* yob_* pc_*)
regsave eduyears2 using "results/ivreg2_outcome_full_adj_restrict",detail(all) pval replace

ds out_*
foreach i in `r(varlist)'{
	ivreg2 `i' (eduyears2=allele_score_21) imob_* yob_* pc_* cov_matsmoking cov_birthweight cov_dist_lon cov_birth_location_imd ///
		cov_birth_location_urban cov_birth_location_northing cov_birth_location_easting [pweight=weight_hughes], ///
		cluster(mobi) endog(eduyears2) partial(imob_* yob_* pc_* cov_matsmoking cov_birthweight cov_dist_lon ///
		cov_birth_location_imd cov_birth_location_urban cov_birth_location_northing cov_birth_location_easting )
	regsave eduyears2 using "results/ivreg2_outcome_full_adj",detail(all) pval append
	ivreg2 `i' (eduyears2=allele_score_21) imob_* yob_* pc_* [pweight=weight_hughes] if e(sample),cluster(mobi) endog(eduyears2) ///
		partial(imob_* yob_* pc_* )
	regsave eduyears2 using "results/ivreg2_outcome_full_adj_restrict",detail(all) pval append
	}
	
*************************

use "results/ivreg2_outcome_full_adj_restrict",clear
gen adjusted=0
append using "results/ivreg2_outcome_full_adj"
replace adjusted=1 if adjusted==.


drop if depvar=="cov_male"|depvar=="out_effect_af"|depvar=="out_effect"|depvar=="out_se"

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

//Plot results using coefplot
//Genetic scores
mkmat coef stderr if adjusted==0 & cont==0,matrix(results) rownames(order)
matrix results_basic_adjusted=results'

mkmat coef stderr if adjusted==1 & cont==0,matrix(results) rownames(order)
matrix results_full_adjusted=results'

#delimit ;
coefplot (matrix(results_basic_adjusted) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
		 (matrix(results_full_adjusted) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(-0.1)) 
		, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100",size(small))
	xlabel(-20(5)20,noticks) rescale(100)  legend(off)
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
graph save "results/supplementary_figure_5_bin", replace;
#delimit cr

//Continious variables

mkmat coef stderr if adjusted==0 & cont==1,matrix(results) rownames(order)
matrix results_basic_adjusted=results'

mkmat coef stderr if adjusted==1 & cont==1,matrix(results) rownames(order)
matrix results_full_adjusted=results'

#delimit ;
coefplot (matrix(results_basic_adjusted) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
		 (matrix(results_full_adjusted) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(-0.1)) 
		 , graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
		xlabel(-3(1)3,noticks)  xtick(none) ytick(none)   legend(off)
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

graph save "results/supplementary_figure_5_cont", replace 

graph combine "results/supplementary_figure_5_bin" "results/supplementary_figure_5_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/supplementary_figure_5_basic_vs_fully_adjusted.eps", as(pdf) replace fontface("Calibri")


//Create table for spreadsheet
save "workingdata/temp",replace

use "workingdata/temp",clear
keep if adjusted==1

rename coef adj_coef
gen adj_lci=adj_coef-1.96*stderr
gen adj_uci=adj_coef+1.96*stderr

format %9.3f adj_* 
keep depvar var adj_* order
sort  order
order  depvar  adj_coef adj_lci adj_uci 
save "results/fully_adj_results",replace

use "workingdata/temp",clear
keep if adjusted==0

rename coef unadj_coef
gen unadj_lci=unadj_coef-1.96*stderr
gen unadj_uci=unadj_coef+1.96*stderr

format %9.3f unadj_* 
keep depvar var unadj_* order
sort  order
order  depvar  unadj_coef unadj_lci unadj_uci 
save "results/basicadj_results",replace






joinby depvar using  "results/fully_adj_results",
drop var order
