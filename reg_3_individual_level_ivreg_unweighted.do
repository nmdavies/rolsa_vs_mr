//Neil Davies 23/05/17
//This cleans the education data in UK Biobank to derived ISCED levels
//Using weighted allele scores

use "workingdata/full_merged_dataset",clear

//IVREG2 analysis
ivreg2 cov_male (eduyears2=allele_score_21) imob_* yob_* pc_*,cluster(mobi) endog(eduyears2) partial(imob_* yob_* pc_*)
regsave eduyears2 using "results/ivreg2_outcome_unweighted",detail(all) pval replace
ds out_*
foreach i in `r(varlist)'{
	ivreg2 `i' (eduyears2=allele_score_21) cov_male imob_* yob_* pc_*,cluster(mobi) endog(eduyears2) partial(cov_male imob_* yob_* pc_*)
	regsave eduyears2 using "results/ivreg2_outcome_unweighted",detail(all) pval append
	}
	
//ROSLA
//Extract from previous paper accounting for age effects
*************************

use "results/ivreg2_outcome_unweighted",clear
gen weight=0
append using "results/ivreg2_outcome"
replace weight=1 if weight==.

drop if depvar=="cov_male"|depvar=="out_effect_af"|depvar=="out_effect"|depvar=="out_se" 
keep if var=="eduyears2"

order depvar
gen order=.

//Reorder

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
//Unweighted
mkmat coef stderr if weight==0 & cont==0,matrix(results) rownames(order)
matrix results_unweighted=results'

mkmat coef stderr if weight==1 & cont==0,matrix(results) rownames(order)
matrix results_weighted=results'


#delimit ;
coefplot (matrix(results_unweighted) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.1))
		 (matrix(results_weighted) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
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
graph save "results/supp_figure_2_weighted_vs_unweighted_bin", replace;
#delimit cr

//Continious variables
mkmat coef stderr if weight==0 & cont==1,matrix(results) rownames(order)
matrix results_unweighted=results'

mkmat coef stderr if weight==1 & cont==1,matrix(results) rownames(order)
matrix results_weighted=results'


#delimit ;
coefplot  (matrix(results_unweighted) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.1))
		 (matrix(results_weighted) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(0)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
	legend(order(2 "Unweighted" 4 "Weighted" ) row(1) size(vsmall))
	xlabel(-3(1)3,noticks)  xtick(none) ytick(none)  legend(off)
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

graph save "results/supp_figure_2_weighted_vs_unweighted_cont", replace 


graph combine "results/supp_figure_2_weighted_vs_unweighted_bin" "results/supp_figure_2_weighted_vs_unweighted_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/supp_figure_2_weighted_vs_unweighted.eps", as(pdf) replace fontface("Calibri")

//Create table for spreadsheet
save "workingdata/temp",replace
use "workingdata/temp",clear
keep if weight==0

rename coef unweighted_coef
gen unweighted_lci=unweighted_coef-1.96*stderr
gen unweighted_uci=unweighted_coef+1.96*stderr

format %9.3f unweighted_* 
keep depvar var unweighted_* order
sort  order
order  depvar  unweighted_coef unweighted_lci unweighted_uci 
save "results/unweighted_results",replace

use "workingdata/temp",clear
keep if weight==1

rename coef weighted_coef
gen weighted_lci=weighted_coef-1.96*stderr
gen weighted_uci=weighted_coef+1.96*stderr

format %9.3f weighted_* 
keep depvar var weighted_* order
sort  order
order  depvar  weighted_coef weighted_lci weighted_uci 
save "results/weighted_results",replace


joinby depvar using  "results/unweighted_results",
drop var order
