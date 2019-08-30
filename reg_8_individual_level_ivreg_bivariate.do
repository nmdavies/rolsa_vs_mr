//Neil Davies 23/05/17
//This cleans the education data in UK Biobank to derived ISCED levels
//Using weighted allele scores

use "workingdata/full_merged_dataset",clear


joinby n_eid using "workingdata/clean_snps",unmatched(master)
drop rs4819996_A_G-rs9977276_T_G
gen interim=(_m==3)
drop _m

//IVREG2 analysis
ivreg2 out_income_over_31k (eduyears2 out_intell=allele_score_21 allele_score_18)   imob_* yob_* pc_*   if interim!=1,cluster(mobi) endog(eduyears2) first partial(  imob_* yob_* pc_*  )
local sfw_educ=el(e(first),8,1)
local sfw_intell=el(e(first),8,2)
regsave eduyears2 out_intell using "results/ivreg2_bivariate_outcome",detail(all) pval replace addvar(SWF_educ,`sfw_educ',0,SWF_cog,`sfw_intell',0)
ds out_*
foreach i in `r(varlist)'{
	ivreg2 `i' (eduyears2 out_intell=allele_score_21 allele_score_18) cov_male  imob_* yob_* pc_*   [pweight=weight]  if interim!=1,cluster(mobi) endog(eduyears2) first partial(cov_male imob_* yob_* pc_*  )
	local sfw_educ=el(e(first),8,1)
	local sfw_intell=el(e(first),8,2)
	regsave eduyears2 out_intell using "results/ivreg2_bivariate_outcome",detail(all) pval append addvar(SWF_educ,`sfw_educ',0,SWF_cog,`sfw_intell',0)
	}


//ROSLA
//Extract from previous paper accounting for age effects
*************************

use "results/ivreg2_bivariate_outcome",clear

drop if depvar=="cov_male"|depvar=="out_se"


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

//Rescale 

//Plot results using coefplot
//Genetic scores


mkmat coef stderr if var=="eduyears2" & cont==0,matrix(results) rownames(order)
matrix results_eduyears=results'

mkmat coef stderr if var=="out_intell" & cont==0,matrix(results) rownames(order)
matrix results_intell=results'

#delimit ;
coefplot (matrix(results_eduyears) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.1))
		 (matrix(results_intell) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0)) 
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

graph save "results/figure_6_bin_individual_level_ivreg_bivariate", replace;
#delimit cr

//Continious variables

drop if depvar=="out_intell"

mkmat coef stderr if var=="eduyears2" & cont==1,matrix(results) rownames(order)
matrix results_eduyears=results'

mkmat coef stderr if var=="out_intell" & cont==1,matrix(results) rownames(order)
matrix results_intell=results'


#delimit ;
coefplot (matrix(results_eduyears) , se(2) ms(T) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(+0.1))
		 (matrix(results_intell) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
	legend(order(2 "Education" 4 "Intelligence" ) row(1) size(vsmall))
	xlabel(-3(1)3,noticks)  xtick(none) ytick(none) 
	headings(14 = "{bf:Indicators of aging}"
			 16 = "{bf:Anthropometry}"
			 18 = "{bf:Blood pressure}"
			 21 = "{bf:Neurocognitive}"
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

graph save "results/figure_6_cont_individual_level_ivreg_bivariate",  replace 

graph combine "results/figure_6_bin_individual_level_ivreg_bivariate" "results/figure_6_cont_individual_level_ivreg_bivariate", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/figure_6_individual_level_ivreg_bivariate.eps", as(pdf) replace fontface("Calibri")

//Output for the supplementary materials
save "workingdata/temp",replace
use "workingdata/temp",clear
keep if var=="eduyears2" | var=="SWF_educ"
sort depvar var

gen  SWF_educ=coef[_n-1] if var=="eduyears2"

rename coef MR_educ_coef
gen MR_educ_lci=MR_educ_coef-1.96*stderr
gen MR_educ_uci=MR_educ_coef+1.96*stderr
rename estatp MR_educ_estatp

format %9.3f MR_*  SWF_educ
format %9.1e MR_educ_estatp
keep var MR_*  SWF_educ order
sort order
order  var MR_educ_coef MR_educ_lci MR_educ_uci SWF_educ
drop if var=="SWF_educ"
drop if S<10
save "results/MR_bivar_educ_results",replace

use "workingdata/temp",clear
keep if !(var=="eduyears2" | var=="SWF_educ")
sort depvar var

gen  SWF_cogc=coef[_n-1] if var=="out_intell"

rename coef MR_intell_coef
gen MR_intell_lci=MR_intell_coef-1.96*stderr
gen MR_intell_uci=MR_intell_coef+1.96*stderr
rename depvar outcome
format %9.3f MR_intell_*  SWF_cogc
keep outcome var MR_intell_* SWF_cogc order
sort order
order  outcome MR_intell_coef MR_intell_lci MR_intell_uci 
keep if var=="out_intell"
drop if S<10
save "results/MR_bivar_intell_results",replace

use "results/MR_bivar_educ_results",clear
joinby order using "results/MR_bivar_intell_results",
drop var order MR_educ_estatp
order outcome  MR_educ_coef MR_educ_lci MR_educ_uci  MR_intell_coef MR_intell_lci MR_intell_uci SWF_educ SWF_cogc
