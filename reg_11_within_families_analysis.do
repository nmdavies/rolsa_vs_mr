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


use "workingdata/full_merged_dataset_siblings",clear

drop if eduyears2==.
drop fam_n
bys famid :gen fam_n=_N
drop if fam_n==1

//IVREG2 analysis

preserve
drop fam_n
drop if cov_male==.
bys famid :gen fam_n=_N
drop if fam_n==1
xtivreg cov_male (eduyears2=allele_score_21) imob_* yob_* pc_* ,vce(cluster famid ) i(famid) fe
regsave eduyears2 using "results/ivreg2_outcomes_sibs",detail(all) pval replace
restore

ds out_*
foreach i in `r(varlist)'{
	preserve
	drop fam_n
	drop if `i'==.
	bys famid :gen fam_n=_N
	drop if fam_n==1
	xtivreg `i'  (eduyears2=allele_score_21) imob_* yob_* pc_* ,vce(cluster famid ) i(famid) fe
	regsave eduyears2 using "results/ivreg2_outcomes_sibs",detail(all) pval append
	restore
	}


*************************

use "results/ivreg2_outcomes_sibs",clear
drop if _n==_N
append using "results/ivreg2_outcome"

drop if depvar=="cov_male"|depvar=="out_effect_af"|depvar=="out_effect"|depvar=="out_se"
order depvar

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

sort order

//Plot results using coefplot
//Genetic scores
mkmat coef stderr if cmd=="ivreg2" & cont==0,matrix(results) rownames(order)
matrix results_iv=results'

mkmat coef stderr if cmd=="xtivreg"& cont==0 ,matrix(results) rownames(order)
matrix results_xtivreg=results'

#delimit ;
coefplot (matrix(results_iv) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(+0.05))
		 (matrix(results_xtivreg) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.05)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Risk difference*100 (95% CI)",size(small))
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
graph export "results/supp_figure_8_bin_individual_level_xtivreg.pdf", as(pdf) replace fontface("Arial");
graph export "results/supp_figure_8_bin_individual_level_xtivreg.eps", as(pdf) replace fontface("Calibri (Body)");
graph save "results/supp_figure_8_bin_individual_level_xtivreg", replace;
#delimit cr

//Continious variables

mkmat coef stderr if cmd=="ivreg2" & cont==1,matrix(results) rownames(order)
matrix results_iv=results'

mkmat coef stderr if cmd=="xtivreg"& cont==1 ,matrix(results) rownames(order)
matrix results_xtivreg=results'


#delimit ;
coefplot (matrix(results_iv) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(+0.05))
		 (matrix(results_xtivreg) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.05)) 
	, graphregion(color(white))  plotregion(lc(white)) grid(none) xline(0) ylabel(,labsize(vsmall)) xlabel(,labsize(vsmall)) xtitle("Mean difference (95% CI)",size(vsmall))
	legend(order(2 "Unrelated participants" 4 "Sibling fixed effects" ) row(1) size(vsmall))
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

graph save "results/supp_figure_8_cont_individual_level_xtivreg", replace 

graph combine "results/supp_figure_8_bin_individual_level_xtivreg" "results/supp_figure_8_cont_individual_level_xtivreg", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/supp_figure_8_combined_individual_level_xtivreg.eps", as(pdf) replace fontface("Calibri")

//Output for the supplementary materials
save "workingdata/temp",replace
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
