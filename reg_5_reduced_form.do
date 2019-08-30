//Neil Davies 19/06/17
//This estimates the association between the Okbay score and the outcomes (the reduced form)


use "workingdata/full_merged_dataset",clear

//Regression
reg cov_male allele_score_21  imob_* yob_* pc_* [pweight=weight_hughes],cluster(mobi)
regsave allele_score_21 using "results/reg_outcome_reduced_form",detail(all) pval replace
ds out_*
foreach i in `r(varlist)'{
	reg `i'  allele_score_21 cov_male imob_*   imob_* yob_* pc_* [pweight=weight_hughes],cluster(mobi)
	regsave allele_score_21 using "results/reg_outcome_reduced_form",detail(all) pval append
	}
	
use "results/reg_outcome_reduced_form",clear

drop if depvar=="cov_male"|depvar=="out_effect_af"|depvar=="out_effect"

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
//Genetic scores
mkmat coef stderr if cont==0,matrix(results) rownames(order)
matrix results=results'

#delimit ;
coefplot (matrix(results) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(+0.1))
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100",size(small))
	xlabel(-30(5)30,noticks) rescale(100)  legend(off)
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
graph save "results/supplementary_figure_3_reduced_form_bin", replace;
#delimit cr

//Continious variables
mkmat coef stderr if cont==1,matrix(results) rownames(order)
matrix results=results'

#delimit ;
coefplot (matrix(results) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(+0.1))
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
	xlabel(-4(1)4,noticks)  xtick(none) ytick(none) legend(off)
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
graph save "results/supplementary_figure_3_reduced_form_cont", replace 

graph combine "results/supplementary_figure_3_reduced_form_bin" "results/supplementary_figure_3_reduced_form_cont", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/supplementary_figure_3_reduced_form.eps", as(pdf) replace fontface("Calibri")

//Create table for spreadsheet
save "workingdata/temp",replace
use "workingdata/temp",clear

rename coef rf_coef
gen rf_lci=rf_coef-1.96*stderr
gen rf_uci=rf_coef+1.96*stderr

format %9.3f rf_* 
keep depvar var rf_* order
sort  order
order  depvar rf_coef rf_lci rf_uci 

