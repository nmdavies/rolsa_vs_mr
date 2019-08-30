//Neil Davies 15/05/17
//This runs the covariate balance tests on the baseline covariates


//Next estimate association of Education, ROSLA, Okbay and Rietveld scores and the covariates
//This programme runs the balance tests

cap prog drop hetero_test
prog def hetero_test

args outcome exposure iv 

macro shift 3
local covar="`*'"

cap drop _const
cap gen _const  = 1

di "outcome=`outcome'"
di "exposure=`exposure'"
di "instrumen=`iv'"
di "covariates=`covar'"

gmm (`outcome' - {xb1:`exposure' `covar' _const})  ///
	(`outcome' - {xb2:`exposure' `covar' _const}) [pweight=weight_hughes], ///
	instruments(1:`exposure' `covar') ///
	instruments(2:`iv' `covar') ///
	winit(unadjusted,independent) onestep  ///
	vce(cluster mob) ///
	deriv(1/xb1 = -1) ///
	deriv(2/xb2 = -1)
drop _const

local outcome2=substr("`outcome'",1,16)
est sto results_`outcome2'

lincom _b[xb1:`exposure']-_b[xb2:`exposure']
local het_p=2*(1-normal(abs(r(estimate)/r(se))))

regsave `exposure' using "results/bias_plots_basic_adjusted_`outcome'_`iv'", detail(all) pval ci replace addvar(het_p,`het_p')

end

use "workingdata/full_merged_dataset",clear

//Drop agression which did not have any valid SNPs in MR-Base
cap:drop allele_score_6

rename allele_score_21 allele_Score_21
rename allele_score_50 allele_Score_50
rename allele_score_51 allele_Score_51

//Run bias calculations for the genome-wide scores
ds allele_score_*
foreach i in `r(varlist)'{ 
	hetero_test `i' eduyears2 allele_Score_21 imob_* cov_male pc_*
	hetero_test `i' more_edu_15 bw12 imob_* cov_male pc_*
	}

//Set negative values of easting and northing to missing
replace cov_birth_location_northing=. if cov_birth_location_northing ==-10
replace  cov_birth_location_easting=. if 	cov_birth_location_easting==-10
		
//Run the bias calculations for the measured baseline phenotypes
//Normalise continuous covariates to allow comparison
foreach i in cov_comp_bodysize10 cov_comp_height10 cov_birthweight cov_dist_lon cov_birth_location_imd cov_birth_location_easting cov_birth_location_northing {
	egen X=std(`i')
	replace `i'=X
	drop X
	}

ds cov_comp_bodysize10 cov_comp_height10 cov_birthweight cov_dist_lon cov_birth_location_imd cov_birth_location_easting cov_birth_location_northing  ///
		cov_breastfed cov_father_alive cov_mother_alive cov_num_sisters cov_num_brothers cov_matsmoking
foreach i in `r(varlist)'{
	hetero_test `i' eduyears2 allele_Score_21 imob_* cov_male pc_*
	hetero_test `i'  more_edu_15 bw12 imob_* cov_male pc_* 
	}
ds allele_score_*
ds cov_*	, varwi(32)

//Clean the results and construct coefficient plots

use "results/bias_plots_basic_adjusted_allele_score_1_bw12",clear
#delimit ;
forvalues i =2(1)52{;
	cap:append using "results/bias_plots_basic_adjusted_allele_score_`i'_bw12";
	};

foreach i in 
cov_mother_alive           cov_comp_bodysize10         cov_GW_EA2_score           cov_birth_location_northing
cov_father_alive           cov_comp_height10           cov_Z_GW_EA2_score         cov_birth_location_easting
cov_num_sisters            cov_matsmoking             cov_dist_lon               cov_current_location_easting
cov_num_brothers           cov_birthweight            cov_birth_location_imd     cov_current_location_northing
cov_breastfed              cov_male                   cov_birth_location_urban{;
	append using "results/bias_plots_basic_adjusted_`i'_bw12";
	};
#delimit cr

gen outcome=substr(word(cmdline ,2),2,.)
gen instrument="bw12" if substr(var,1,3)=="xb2"
replace instrument="more_educ_15" if substr(var,1,3)=="xb1"

gen het_test=coef[_n+1] if substr(var,1,3)=="xb2"
order outcome instrument coef stderr het_test
drop if var=="het_p"

keep outcome instrument coef stderr het_test cmdline
compress
save "results/merged_bw12",replace



use "results/bias_plots_basic_adjusted_allele_score_1_allele_score_21",clear
#delimit ;

forvalues i =2(1)52{;
	cap:append using "results/bias_plots_basic_adjusted_allele_score_`i'_allele_score_21",;
	};
foreach i in 
cov_mother_alive           cov_comp_bodysize10                   cov_birth_location_northing
cov_father_alive           cov_comp_height10                    cov_birth_location_easting
cov_num_sisters            cov_matsmoking             cov_dist_lon               cov_current_location_easting
cov_num_brothers           cov_birthweight            cov_birth_location_imd     cov_current_location_northing
cov_breastfed              cov_male                   cov_birth_location_urban{;
	append using "results/bias_plots_basic_adjusted_`i'_allele_score_21";
	};
#delimit cr

gen outcome=substr(word(cmdline ,2),2,.)
gen instrument="allele_score_21" if substr(var,1,3)=="xb2"
replace instrument="more_educ_15" if substr(var,1,3)=="xb1"

gen het_test=coef[_n+1] if substr(var,1,3)=="xb2"
order outcome instrument coef stderr het_test
drop if var=="het_p"

keep outcome instrument coef stderr het_test cmdline
compress
save "results/merged_allele_score_21",replace
gen analysis="allele_score_21"
append using "results/merged_bw12"
replace analysis="bw12" if analysis==""
save "results/merged",replace
use "results/merged",clear
//Drop results relating to the EA score covariate
drop if outcome=="allele_score_50"|outcome=="allele_score_51"|outcome=="allele_score_21"

//Join in the variable labels
joinby outcome using "workingdata/labels.dta",unmatched(master)
sort order

//Plot results using coefplot
//Genetic scores
mkmat coef stderr if instrument=="bw12" & substr(outcome,1,2)=="al",matrix(results) rownames(order)
matrix results_bw12=results'

mkmat coef stderr if instrument=="allele_score_21"& substr(outcome,1,2)=="al",matrix(results) rownames(order)
matrix results_allele_score_50=results'

mkmat coef stderr if instrument=="more_educ_15" & analysis=="allele_score_21"& substr(outcome,1,2)=="al",matrix(results) rownames(order)
matrix results_more_educ=results'

#delimit ;
coefplot (matrix(results_bw12) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0))
		 (matrix(results_allele_score_50) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.1)) 
		 (matrix(results_more_educ) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(0.1)) 
	, graphregion(color(white)) plotregion(color(white)) grid(none) xline(0) ysize(10.69) xsize(7.27) ylabel(,labsize(tiny)) xlabel(,labsize(tiny)) 
	legend(order(6 "Remained at age 15" 2 "ROSLA" 4 "Mendelian randomization" ) row(1) size(tiny))
	xlabel(-0.5(0.25)0.5,noticks)
	headings(1 = "{bf:Anthropometry}"
			 4 = "{bf:Substance misuse}"
			 9 = "{bf:Childhood phenotypes}"
			 12 = "{bf:Neuropsychatric conditions}"
			 22 = "{bf:Socioeconomic characteristics}"
			 31 = "{bf:Cognition}"
			 42 = "{bf:Nutrition}", labsize(tiny)) 
	coeflabels(1="Height"
2="Body mass index"
3="Body fat percentage"
4="Cigarettes smoked per day"
5="Ever vs never smoked"
6="Age of smoking initiation"
7="Alcohol dependence"
8="Birth weight"
9="Birth length"
10="Infant head circumference"
11="Age at menarche"
12="Depressive symptoms"
13="Major depressive disorder"
14="Autism"
15="Schizophrenia"
16="Bipolar disorder"
17="Migraine in bipolar disorder"
18="PGC cross-disorder traits"
19="Alzheimer's disease"
20="Father’s age at death"
21="Mother’s age at death"
22="Agreeableness"
23="Conscientiousness"
24="Extraversion"
25="Openness to experience"
26="Neuroticism"
27="Internalizing problems"
28="Subjective well being"
29="Chronotype"
30="Sleep duration"
31="G speed factor"
32="Symbol search"
33="Digit symbol"
34="Inspection time"
35="2-choice reaction time"
36="8-choice reaction time"
37="Simple reaction time"
38="Childhood intelligence"
39="Cognition Sniekers et al."
40="Cognition Trampush p<5E-07"
41="High IQ Zabaneh"
42="Omega-3 fatty acids"
43="Omega-6 fatty acids"
44="Omega-9 and saturated fatty acids"
45="Other PUFA"
46="Linoleic acid (LA)"
47="Mono-unsaturated fatty acids"
48="Zinc");
#delimit cr
graph export "results/figure_1_balance_genetic_scores.eps", as(pdf) replace 
save "workingdata/figure_1_bias_plot",replace

//Non-genetic covariates
//Merge in labels
drop _m
joinby outcome using "workingdata/labels2.dta",unmatched(master) update
sort order
drop if order==.

drop if outcome=="cov_current_location_easting"
drop if outcome=="cov_current_location_northing"

mkmat coef stderr if instrument=="bw12" & substr(outcome,1,2)!="al" & outcome!="cov_male",matrix(results) rownames(order)
matrix results_bw12=results'

mkmat coef stderr if instrument=="allele_score_21"& substr(outcome,1,2)!="al"& outcome!="cov_male",matrix(results) rownames(order)
matrix results_allele_score_50=results'

mkmat coef stderr if instrument=="more_educ_15" & analysis=="allele_score_21"& substr(outcome,1,2)!="al"& outcome!="cov_male",matrix(results) rownames(order)
matrix results_more_educ=results'

#delimit ;
coefplot (matrix(results_bw12) , se(2) ms(S) msize(vsmall) mc(edkblue) ciopts(lc(edkblue) lwidth(vthin)) offset(0))  
		 (matrix(results_allele_score_50) , se(2) ms(T) msize(vsmall) mc(eltblue) ciopts(lc(eltblue) lwidth(vthin)) offset(-0.1)) 
		 (matrix(results_more_educ) , se(2) ms(C) msize(vsmall) mc(erose) ciopts(lc(erose) lwidth(vthin)) offset(0.1)) 
	     , graphregion(color(white)) plotregion(color(white)) grid(none) xline(0) ysize(5.345) xsize(7.27) ylabel(,labsize(tiny)) xlabel(,labsize(tiny)) 
	     xlabel(-1(0.5)1,noticks)
		 legend(order(6 "Remained at age 15" 2 "ROSLA" 4 "Mendelian randomization" ) row(1) size(tiny))
		 headings(1 = "{bf:Birth location}"
			 6 = "{bf:Current location}"
			 9 = "{bf:Early life}"
			 14 = "{bf:Family characteristics}"
			 , labsize(tiny)) 
	coeflabels(1="Easting"
2="Northing"
3="Index of Multiple Deprivation"
4="Urban vs. rural"
5="Distance from London"
6="Easting"
7="Northing"
8="Male"
9="Birthweight"
10="Breastfed"
11="Mother smoked in pregnancy"
12="Comparative body size age 10"
13="Comparative height age 10"
14="Father alive"
15="Mother alive"
16="Number of brothers"
17="Number of sisters");

graph export "results/figure_2_balance_phenotypes.eps", as(pdf) replace ;
#delimit cr
//Create Excel spreadsheets of results
save "workingdata/temp",replace
use "workingdata/temp",clear
drop if outcome=="cov_male"
keep if instrument=="bw12"
gen genetic=(substr(outcome,1,2)!="al")
rename coef bw12_coef
gen bw12_lci=bw12_coef-1.96*stderr
gen bw12_uci=bw12_coef+1.96*stderr
rename het_test bw12_het_pval
format %9.3f bw12_* 
keep outcome bw12_* order genetic
sort genetic order
order  outcome bw12_coef bw12_lci bw12_uci bw12_het_pval
save "results/bias_test_bw12",replace

use "workingdata/temp",clear
drop if outcome=="cov_male"
keep if instrument=="allele_score_21"
gen genetic=(substr(outcome,1,2)=="al")

rename coef prs_coef
gen prs_lci=prs_coef-1.96*stderr
gen prs_uci=prs_coef+1.96*stderr
rename het_test prs_het_pval
format %9.3f prs_* 
keep outcome prs_* order genetic
sort genetic order
order  outcome prs_coef prs_lci prs_uci prs_het_pval
save "results/bias_test_prs",replace

use "workingdata/temp",clear
drop if outcome=="cov_male"
keep if instrument=="more_educ_15" & analysis=="allele_score_21"

gen genetic=(substr(outcome,1,2)!="al")

rename coef educ_coef
gen educ_lci=educ_coef-1.96*stderr
gen educ_uci=educ_coef+1.96*stderr
format %9.3f educ_* 
keep outcome educ_* order genetic
sort genetic order
order  outcome educ_coef educ_lci educ_uci 
save "results/bias_test_educ",replace	
	
joinby outcome using "results/bias_test_bw12"
joinby outcome using "results/bias_test_prs"
drop order genetic

format %9.2e *het_pval

//Finally calculating the association between each of the covariates and the outcomes:

use "workingdata/full_merged_dataset",clear

//Gen updated weights from Hughes et al.
gen weight2=34.29 if weight!=1
replace weight2=12.37 if weight==1

drop allele_score_6
reg out_alcohol allele_score_1 imob_* cov_male pc_*[pweight=weight2],ro
regsave allele_score_1 using "results/covariate_outcome",replace detail(all)
ds out_*
foreach j in `r(varlist)'{
	ds allele_score_* cov_*
	foreach i in `r(varlist)'{
		di "`i' `j'"
		reg `j' `i' imob_* cov_male pc_*[pweight=weight2],ro
		regsave `i' using "results/covariate_outcome",append detail(all)
		}
	}
use "results/covariate_outcome",clear

egen pval=2*(1-norm(abs(coef/stderr)))

keep var depvar coef pval
levels depvar 
foreach i in `r(levels)'{
	preserve
	keep if depvar=="`i'"
	gen n=_n
	local j=`j'+1
	save "workingdata/temp`j'"
	restore
	}
use "workingdata/temp1",clear
rename coef coef_1
rename pval pval_1
forvalues i=2(1)25{
	joinby n using "workingdata/temp`i'",
	rename coef coef_`i'
	rename pval pval_`i'
	}
rename var outcome	
joinby outcome using "workingdata/labels.dta",unmatched(master)	
drop _m
joinby outcome using "workingdata/labels2.dta",unmatched(master) update

forvalues i=1(1)25{
	replace coef_`i'=. if pval_`i'>0.05
	}

drop n pval*
replace order=order+48 if _m==3
sort order
