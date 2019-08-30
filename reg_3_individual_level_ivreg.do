//Neil Davies 23/05/17
//This cleans the education data in UK Biobank to derived ISCED levels
//Using weighted allele scores

*! 1.0.0 Tom Palmer 8jan2010
program tsci, rclass 
syntax anything, [eform]
* syntax gd segd gp segp cov (if poss)
local n = wordcount("`anything'")
tokenize `anything'
tempname gd segd gp segp cov ratio seratio z p
sca `gd' = `1'
sca `segd' = `2'
sca `gp' = `3'
sca `segp' = `4'
if "`n'" == "4" {
	sca `cov' = 0
}
else {
	sca `cov' = `5'
}
sca `ratio' = `gd'/`gp'
sca `seratio' = sqrt((`segd'^2/`gp'^2) + (`gd'^2/`gp'^4)*`segp'^2 - 2*(`gd'/`gp'^3)*`cov')
sca `z' = abs(`ratio'/`seratio')
sca `p' = 2*normal(-1*`z')
di as res `ratio', `seratio', "(" `ratio' - invnormal(0.975)*`seratio' ", " `ratio' + invnormal(0.975)*`seratio' ")", "Z=" `z', "P=" `p'
if "`eform'" != "" {
di as res exp(`ratio'), "(" exp(`ratio' - invnormal(0.975)*`seratio') ", " exp(`ratio' + invnormal(0.975)*`seratio') ")", "P=" `p'
}
ret sca ratio = `ratio'
ret sca seratio = `seratio'
end


use "workingdata/full_merged_dataset",clear


//IVREG2 analysis
ivreg2 cov_male (eduyears2=allele_score_21) imob_* yob_* pc_* [pweight=weight_hughes],cluster(mobi) endog(eduyears2) ro partial(imob_* yob_* pc_*) first
regsave eduyears2 using "results/ivreg2_outcome",detail(all) pval replace
ds out_*
foreach i in `r(varlist)'{
	ivreg2 `i' (eduyears2=allele_score_21) cov_male imob_* yob_* pc_* [pweight=weight_hughes],cluster(mobi) endog(eduyears2) partial(cov_male imob_* yob_* pc_*)
	regsave eduyears2 using "results/ivreg2_outcome",detail(all) pval append
	}

//Adjusted regression
reg cov_male eduyears2 imob_* yob_* pc_*   [pweight=weight_hughes],cluster(mobi)
regsave eduyears2 using "results/reg_outcome",detail(all) pval replace
ds out_*
foreach i in `r(varlist)'{
	reg `i' eduyears2 cov_male imob_* yob_* pc_1-pc_10 [pweight=weight_hughes],cluster(mobi)
	regsave eduyears2 using "results/reg_outcome",detail(all) pval append
	}
	
//ROSLA
//Extract from previous paper accounting for age effects
*************************

use "results/ivreg2_outcome",clear
drop if _n==_N
append using "results/reg_outcome"
append using "uk_biobank/results.dta",
drop if depvar=="cov_male"|depvar=="out_effect_af"|depvar=="out_effect"|depvar=="out_se"
keep if var=="eduyears2"|var=="ROSLA cohort"
order depvar
//Calculate the effect of remaining in school from the ROSLA estimates using the delta method
//Clark and Royer report the following estimates coef=0.261 se=0.016
//Replace coef (1 year bandwidth results) with diff (diff-in-diff results)
replace coef=diff if var=="ROSLA cohort"
replace stderr=diff_se if var=="ROSLA cohort"

forvalues i=51(1)74{
	local coef=coef[`i']
	local se=stderr[`i']
	tsci `coef' `se' 0.261 0.016
	replace coef =`r(ratio)' in `i'
	replace stderr =`r(seratio)' in `i'
	}



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

//Rescale 

//Plot results using coefplot
//Genetic scores
mkmat coef stderr if cmd=="ivreg2" & cont==0,matrix(results) rownames(order)
matrix results_iv=results'

mkmat coef stderr if cmd=="regress"& cont==0 &var!="ROSLA cohort",matrix(results) rownames(order)
matrix results_ols=results'

mkmat coef stderr if cmd=="regress" & cont==0&var=="ROSLA cohort",matrix(results) rownames(order)
matrix results_rosla=results'


#delimit ;
coefplot (matrix(results_ols) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(+0.1))
		 (matrix(results_iv) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.1)) 
		 (matrix(results_rosla) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0))
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
graph export "results/figure_3_bin_individual_level_ivreg.pdf", as(pdf) replace fontface("Arial");
graph export "results/figure_3_bin_individual_level_ivreg.eps", as(pdf) replace fontface("Calibri (Body)");
graph save "results/figure_3_bin_individual_level_ivreg", replace;
#delimit cr

//Continious variables
mkmat coef stderr if cmd=="ivreg2" & cont==1,matrix(results) rownames(order)
matrix results_iv=results'

mkmat coef stderr if cmd=="regress"& cont==1 &var!="ROSLA cohort",matrix(results) rownames(order)
matrix results_ols=results'

mkmat coef stderr if cmd=="regress" & cont==1&var=="ROSLA cohort",matrix(results) rownames(order)
matrix results_rosla=results'

#delimit ;
coefplot (matrix(results_ols) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(+0.1))
		 (matrix(results_iv) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.1))  	
		 (matrix(results_rosla) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0))
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference",size(vsmall))
	legend(order(2 "ISCED education" 6 "ROSLA" 4 "Mendelian randomization" ) row(1) size(vsmall))
	xlabel(-3(1)3,noticks)  xtick(none) ytick(none) 
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

graph save "results/figure_3_cont_individual_level_ivreg", replace 

graph combine "results/figure_3_bin_individual_level_ivreg" "results/figure_3_cont_individual_level_ivreg", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/figure_3_combined_individual_level_ivreg.eps", as(pdf) replace fontface("Calibri")

//Output for the supplementary materials
save "workingdata/temp",replace
use "workingdata/temp",clear
keep if var=="ROSLA cohort"

rename coef ROSLA_coef
gen ROSLA_lci=ROSLA_coef-1.96*stderr
gen ROSLA_uci=ROSLA_coef+1.96*stderr

format %9.3f ROSLA_* 
keep var ROSLA_* order
sort  order
order  var ROSLA_coef ROSLA_lci ROSLA_uci 
save "results/ROSLA_results",replace

use "workingdata/temp",clear
keep if var=="eduyears2" & endog_ct!=.

rename coef MR_coef
gen MR_lci=MR_coef-1.96*stderr
gen MR_uci=MR_coef+1.96*stderr
rename estatp MR_estatp
rename cdf MR_cdf

format %9.3f MR_* 
format %9.1e MR_estatp
keep var MR_* order
sort order
order  var MR_coef MR_lci MR_uci 
save "results/MR_results",replace

use "workingdata/temp",clear
keep if var=="eduyears2" & endog_ct==.

rename coef OLS_coef
gen OLS_lci=OLS_coef-1.96*stderr
gen OLS_uci=OLS_coef+1.96*stderr
rename depvar outcome
format %9.3f OLS_* 
keep outcome var OLS_* order
sort order
order  outcome OLS_coef OLS_lci OLS_uci 
save "results/OLS_results",replace
joinby order using "results/ROSLA_results",
joinby order using "results/MR_results",
drop var order


gen OLS=string(OLS_coef*100,"%9.2f") +" (95%CI: "+string(OLS_lci*100,"%9.2f")+" to " + string(OLS_uci*100,"%9.2f") +")" if _n<14
replace OLS=string(OLS_coef,"%9.2f") +" (95%CI: "+string(OLS_lci,"%9.2f")+" to " + string(OLS_uci,"%9.2f") +")" if _n>=14

gen MR=string(MR_coef*100,"%9.2f") +" (95%CI: "+string(MR_lci*100,"%9.2f")+" to " + string(MR_uci*100,"%9.2f") +")" if _n<14
replace MR=string(MR_coef,"%9.2f") +" (95%CI: "+string(MR_lci,"%9.2f")+" to " + string(MR_uci,"%9.2f") +")" if _n>=14

gen RD=string(ROSLA_coef*100,"%9.2f") +" (95%CI: "+string(ROSLA_lci*100,"%9.2f")+" to " + string(ROSLA_uci*100,"%9.2f") +")" if _n<14
replace RD=string(ROSLA_coef,"%9.2f") +" (95%CI: "+string(ROSLA_lci,"%9.2f")+" to " + string(ROSLA_uci,"%9.2f") +")" if _n>=14
