//Neil Davies 23/05/17
//This cleans the education data in UK Biobank to derived ISCED levels
//This estimates the associations between the Okbay SNPs and each of the outcomes in Biobank study

use  "workingdata/full_merged_dataset_SNPs",clear

//Gen updated weights from Hughes et al.
gen weight_hughes=34.29 if weight!=1
replace weight_hughes=12.37 if weight==1

ds out_*
foreach i in `r(varlist)'{
	reg `i' cov_male  imob_* yob_* pc_*  [pweight=weight_hughes], cluster(mobi)
	regsave imob_1 using "results/outcome_`i'",detail(all) pval replace
	ds rs9537821-rs10061788
	foreach j in `r(varlist)'{
		cap: reg `i' `j'  cov_male imob_* yob_* pc_*  [pweight=weight_hughes],  cluster(mobi)
		cap: regsave `j' using "results/outcome_`i'",detail(all) pval append
		}
	}

//Open OLS results and merge with Okbay discovery SNPs
use  "workingdata/full_merged_dataset_SNPs",clear
ds out_*

local i out_alcohol
di "`i'"
use "results/outcome_`i'",clear
	
foreach i in `r(varlist)'{ 
	cap:append using "results/outcome_`i'"
	}
duplicates drop
drop if var=="imob_1"

gen rsid=substr(var, 1,50)
replace rsid=usubinstr(rsid ,"_"," ",2)
gen a1=word(rsid,2)
gen a2=word(rsid,3)
replace rsid=word(rsid,1)

gen out_bb_effect=coef
gen out_bb_se=stderr

//Create weight file from the Okbay discovery SNPs
preserve
use "workingdata/weights",clear
keep if trait=="Eduyears (Okbay discovery) p<5e-08"
compress
save "workingdata/eduyears_disc_p_5e-08",replace
restore

joinby rsid using  "workingdata/eduyears_disc_p_5e-08",

rename out_effect exp_effect
rename out_se exp_se

keep depvar out_bb_effect out_bb_se rsid a1 a2 exp_effect exp_se
gen aw=1/out_bb_se^2

//Create plot titles
gen title="Hypertension" if depvar=="out_highbloodpressure"
replace title="Diabetes" if depvar=="out_diabetes"
replace title="Stroke" if depvar=="out_stroke"
replace title="Heart attack" if depvar=="out_heartattack"
replace title="Depression" if depvar=="out_depression"
replace title="Cancer" if depvar=="out_cancer"
replace title="Died" if depvar=="out_died"
replace title="Ever smoked" if depvar=="out_exsmoker"
replace title="Currently smoke" if depvar=="out_smoker"
replace title="Income over £18k" if depvar=="out_income_under_18k"
replace title="Income over £31k" if depvar=="out_income_over_31k"
replace title="Income over £52k" if depvar=="out_income_over_52k"
replace title="Income over £100k" if depvar=="out_income_over_100k"
replace title="Gripstrength (kg)" if depvar=="out_gripstrength"
replace title="Arterial stiffness" if depvar=="out_arterial_stiffness"
replace title="Height (cm)" if depvar=="out_height"
replace title="BMI (kg/m2)" if depvar=="out_bmi"
replace title="Diastolic blood pressure (mmHg)" if depvar=="out_dia_bp"
replace title="Systolic blood pressure (mmHg)" if depvar=="out_sys_bp"
replace title="Intelligence (0 to 13)" if depvar=="out_intell"
replace title="Happiness (0 to 5 Likert)" if depvar=="out_happiness"
replace title="Alcohol consumption (1 low to 5 high)" if depvar=="out_alcohol"
replace title="Hours watching television per day" if depvar=="out_sedentary"
replace title="Moderate exercise (days/week)" if depvar=="out_phys_m_act"
replace title="Vigorous exercise (days/week)" if depvar=="out_phys_v_act"
replace title="Mortality" if depvar=="out_dead"

//Reorder + generate indicators for continious or binary
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

//Run summary data regressions:
mregger out_bb_effect exp_effect [aweight=aw] if depvar =="out_alcohol",gxse(out_bb_se) ivw  heter
regsave exp_effect using "results/mregger_ivw_outcome",detail(all) pval replace addlabel(depvar, "out_alcohol")

levels depvar 
foreach i in `r(levels)'{
	di "`i'"
	//IVW
	mregger out_bb_effect exp_effect [aweight=aw] if depvar =="`i'",gxse(exp_se) ivw  heter
	regsave exp_effect using "results/mregger_ivw_outcome",detail(all) pval append addlabel(depvar, "`i'") 
	
	//MR-Egger
	mregger out_bb_effect exp_effect [aweight=aw] if depvar =="`i'",gxse(exp_se) 
	regsave using "results/mregger_ivw_outcome",detail(all) pval append addlabel(depvar, "`i'") 	
	
	//Weighted median
	mrmedian out_bb_effect out_bb_se exp_effect exp_se if depvar =="`i'",  w
	regsave using "results/mregger_ivw_outcome",detail(all) pval append addlabel(depvar, "`i'") 	
	
	//Modal
	mrmodal out_bb_effect out_bb_se exp_effect exp_se if depvar =="`i'",  weight
	regsave using "results/mregger_ivw_outcome",detail(all) pval append addlabel(depvar, "`i'") 
	}

levels out_intell
foreach i in `r(levels)'{
	levels title if depvar=="`i'"
	local title=`r(levels)'
	levels cont if depvar=="`i'"
	if `r(levels)'==0{
		local effect="Risk difference"
		}
	else{
		local effect="Mean difference"
		}
	di "XXXXXXXXXXX"
	di "`title'"
	mreggerplot out_bb_effect out_bb_se exp_effect exp_se if depvar =="`i'", gpci
	mrmedian out_bb_effect out_bb_se exp_effect exp_se if depvar =="`i'", weighted
	addplot : function _b[beta]*x if depvar =="`i'", range(0.01 0.05) lc(gs0) lp(shortdash) lw(vthin)
	mrmodal out_bb_effect out_bb_se exp_effect exp_se if depvar =="`i'", phi(.25)
	addplot : function _b[beta]*x if depvar =="`i'", range(0.01 0.05) lc(gs0) lp(longdash) ///
		legend(order(5 "Instruments" 4 "95% CIs" 3 "MR-Egger" 2 "MR-Egger 95% CI" 6 "Weighted median" 7 "Modal") rows(1) si(vsmall) symx(*.5)) ///
	graphregion(color(white))  plotregion(lc(white)) title("`title'")
	graph export "results/supp_figure_X_mreggerplot_`i'.eps", as(pdf) replace fontface("Calibri")	
			
	mrforest out_bb_effect out_bb_se exp_effect exp_se if depvar =="`i'", ivid(rsid)  models(4) modelslabel(All genotypes) ///
		graphregion(color(white))  plotregion(lc(white)) title("`title'") xtitle("`effect'")
	graph export "results/supp_figure_X_mrforest_`i'.eps", as(pdf) replace fontface("Calibri")	
	}	 	 

//Create summary plot for MR-estimates
use "results/mregger_ivw_outcome",clear
duplicates drop

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


mkmat coef stderr if cmd=="mrmedian" & cont==0,matrix(results) rownames(order)
matrix results_median_bin=results'

mkmat coef stderr if cmd=="mrmodal" & cont==0,matrix(results) rownames(order)
matrix results_modal_bin=results'

mkmat coef stderr if cmd=="mregger" & cont==0 & substr(var,-5,5)=="slope",matrix(results) rownames(order)
matrix results_mregger_bin=results'

mkmat coef stderr if cmd=="mregger" & cont==0 & substr(var,-5,5)=="ffect",matrix(results) rownames(order)
matrix results_ivw_bin=results'

mkmat coef stderr if cmd=="mrmedian" & cont==1,matrix(results) rownames(order)
matrix results_median_cont=results'

mkmat coef stderr if cmd=="mrmodal" & cont==1,matrix(results) rownames(order)
matrix results_modal_cont=results'

mkmat coef stderr if cmd=="mregger" & cont==1 & substr(var,-5,5)=="slope",matrix(results) rownames(order)
matrix results_mregger_cont=results'

mkmat coef stderr if cmd=="mregger" & cont==1 & substr(var,-5,5)=="ffect",matrix(results) rownames(order)
matrix results_ivw_cont=results'

#delimit ;
coefplot (matrix(results_ivw_bin) , se(2) ms(T) msize(vsmall) mc(black) ciopts(lc(black) lwidth(vthin)) offset(0.2))
		 (matrix(results_mregger_bin) , se(2) ms(C) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0.1)) 
		 (matrix(results_median_bin) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0)) 
		 (matrix(results_modal_bin) , se(2) ms(D) msize(vsmall) mc(blue) ciopts(lc(blue) lwidth(vthin)) offset(-0.1)) 
	, graphregion(color(white)) plotregion(lc(white))  grid(none) xline(0) ylabel(,labsize(small)) xlabel(,labsize(small)) xtitle("Risk difference*100")
		legend(off)
	xlabel(-50(10)50,noticks) rescale(100)
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
13=".                        Income over £100k", wrap(51));
graph save "results/figure_6_bin_mregger", replace;

//Continious variables
#delimit ;
coefplot (matrix(results_ivw_cont) , se(2) ms(T) msize(vsmall) mc(black) ciopts(lc(black) lwidth(vthin)) offset(0.2))
		 (matrix(results_mregger_cont) , se(2) ms(C) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0.1)) 
		 (matrix(results_median_cont) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0)) 
		 (matrix(results_modal_cont) , se(2) ms(D) msize(vsmall) mc(blue) ciopts(lc(blue) lwidth(vthin)) offset(-0.1)) 
	, graphregion(color(white)) plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(small)) xlabel(,labsize(small)) xtitle("Mean difference")
	legend(order(2 "IVW" 4 "MR-Egger" 6 "Weighted median" 8 "Weighted mode") row(1) size(small))
	xlabel(-12(2)12,noticks)
	coeflabels(14="Gripstrength (kg)"
15="Arterial stiffness"
16="Height (cm)"
17="BMI (kg/m2)"
18="Diastolic blood pressure (mmHg)"
19="Systolic blood pressure (mmHg)"
20="Intelligence (0 to 13)"
21="Happiness (0 to 5 Likert)"
22="Alcohol consumption (1 low to 5 high)"
23="Hours watching television per day"
24="Moderate exercise (days/week)"
25="Vigorous exercise (days/week)");	

graph save "results/figure_6_cont_mregger", replace;
#delimit cr
graph combine "results/figure_6_bin_mregger" "results/figure_6_cont_mregger", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/figure_6_mregger.eps", as(pdf) replace fontface("Calibri")

//Create table
compress
save "workingdata/temp",replace
use "workingdata/temp",clear

preserve
keep if substr(var,-5,5)=="ffect"

rename coef ivw_coef
gen ivw_lci=ivw_coef-1.96*stderr
gen ivw_uci=ivw_coef+1.96*stderr
rename pval ivw_pval
format %9.3f ivw_* 
keep depvar ivw_* order 
order depvar ivw_coef ivw_lci ivw_uci ivw_pval
save "results/ivw",replace
restore

preserve
keep if substr(var,-5,5)=="slope"
rename coef mr_egger_coef
gen mr_egger_lci=mr_egger_coef-1.96*stderr
gen mr_egger_uci=mr_egger_coef+1.96*stderr
rename pval mr_egger_pval
format %9.3f mr_egger_* 
keep depvar mr_egger_* order 
order depvar mr_egger_coef mr_egger_lci mr_egger_uci mr_egger_pval
save "results/mr_egger_slope",replace
restore

preserve
keep if substr(var,-5,5)=="_cons"
rename coef mr_egger_cons_coef
gen mr_egger_cons_lci=mr_egger_cons_coef-1.96*stderr
gen mr_egger_cons_uci=mr_egger_cons_coef+1.96*stderr
rename pval mr_egger_cons_pval
format %9.3f mr_egger_* 
keep depvar mr_egger_* order 
order depvar mr_egger_cons_coef mr_egger_cons_lci mr_egger_cons_uci mr_egger_cons_pval
save "results/mr_egger_cons",replace
restore

preserve
keep if cmd=="mrmedian"
rename coef mr_egger_med_coef
gen mr_egger_med_lci=mr_egger_med_coef-1.96*stderr
gen mr_egger_med_uci=mr_egger_med_coef+1.96*stderr
rename pval mr_egger_med_pval
format %9.3f mr_egger_* 
keep depvar mr_egger_* order 
order depvar mr_egger_med_coef mr_egger_med_lci mr_egger_med_uci mr_egger_med_pval
save "results/mr_egger_med",replace
restore

preserve
keep if cmd=="mrmodal"
rename coef mr_egger_mod_coef
gen mr_egger_mod_lci=mr_egger_mod_coef-1.96*stderr
gen mr_egger_mod_uci=mr_egger_mod_coef+1.96*stderr
rename pval mr_egger_mod_pval
format %9.3f mr_egger_* 
keep depvar mr_egger_* order 
order depvar mr_egger_mod_coef mr_egger_mod_lci mr_egger_mod_uci mr_egger_mod_pval
save "results/mr_egger_mod",replace
restore


keep if I2!=.
rename I2 I2_coef
rename lb_I I2_lci
rename ub_I I2_uci
format %9.3f I2_* 
keep depvar I2_* order 
order depvar I2_coef I2_lci I2_uci 
save "results/I2",replace


use "results/ivw",clear
foreach i in mr_egger_slope mr_egger_cons mr_egger_med mr_egger_mod I2{
	joinby depvar using "results/`i'"
	}

//Create summary plot of the I^2 stats for each trait:
use results/mregger_ivw_outcome.dta,clear

keep if I2!=.

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

keep depvar I2 ub_I2_M1 lb_I2_M1 Q
