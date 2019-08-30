//Neil Davies 06/11/17
//This estimates the LATE margin for the ROSLA vs MR paper

use "workingdata/full_merged_dataset",clear


/*
gen eduyears3=n_845_0_0 
replace eduyears3=21 if eduyears2>=21 & eduyears2!=.
replace eduyears3=. if eduyears3<0
*/



gen eduyears3=eduyears2
cumul eduyears3, generate(eduyears_cum3)

cumul eduyears3  if post_reform ==0 & bw12!=.,generate(eduyears_cum_pre3)
cumul eduyears3  if post_reform ==1 & bw12!=.,generate(eduyears_cum_post3)

sort eduyears_cum3

//Select the top point 
keep eduyears3 eduyears_cum_pre3 eduyears_cum_post3
drop if eduyears_cum_post3 ==. & eduyears_cum_pre3==.
replace eduyears_cum_pre3=eduyears_cum_pre3[_n-1] if eduyears_cum_pre3==.
replace eduyears_cum_post3=eduyears_cum_post3[_n-1] if eduyears_cum_post3==.
replace eduyears_cum_post3=0 if eduyears_cum_post3==.
bys eduyears3:egen max_eduyears_cum_pre3=max(eduyears_cum_pre3)
bys eduyears3:egen max_eduyears_cum_post3=max(eduyears_cum_post3)
drop eduyears_cum_pre3 eduyears_cum_post3
duplicates drop

sort eduyears3
twoway (line max_eduyears_cum_pre3 eduyears3 ) (line  max_eduyears_cum_post3 eduyears3 )

//Generate difference
gen difference=max_eduyears_cum_post3-max_eduyears_cum_pre3

//Generate the difference in the CDFs
replace difference=max_eduyears_cum_post3-max_eduyears_cum_pre3
line difference eduyears3 if eduyears3>13 & eduyears3<=21 , title("Pre vs. post ROSLA") ///
	 xtitle("Age left full time education") ytitle("CDF difference") graphregion(color(white)) plotregion(lc(white)) ///
	 
graph save "results/figure_X_LATE_ROSLA", replace 





//Repeat for allele score
use "workingdata/full_merged_dataset",clear
xtile edu_AS_split = allele_score_21 , nq(5)

///*
gen eduyears3=n_845_0_0 
replace eduyears3=21 if eduyears2>=21 & eduyears2!=.
replace eduyears3=. if eduyears3<0
//*/



//gen eduyears3=eduyears2

cumul eduyears3, generate(eduyears_cum3)

forvalues i=1(1)5{
	cumul eduyears3  if edu_AS_split==`i',generate(eduyears_cum_AS_`i')
	}
sort eduyears3	
keep eduyears3 eduyears_cum_AS_*
forvalues i=1(1)5{
	replace eduyears_cum_AS_`i'=eduyears_cum_AS_`i'[_n-1] if eduyears_cum_AS_`i'==.
	replace eduyears_cum_AS_`i'=0 if eduyears_cum_AS_`i'==.
	bys eduyears3:egen max_eduyears_cum_AS_`i'=median(eduyears_cum_AS_`i')
	}	

drop eduyears_cum_AS_*
duplicates drop


//Generate the difference
forvalues i=1(1)4{
	gen diff_cum_`i'=max_eduyears_cum_AS_`i'-max_eduyears_cum_AS_5
	}


twoway 	(line max_eduyears_cum_AS_1 eduyears3 if eduyears3>12 & eduyears3<26 , legend(label(1 "1st quintile"))) ///
		(line max_eduyears_cum_AS_2 eduyears3 if eduyears3>13 & eduyears3<25 , legend(label(2 "2nd quintile"))) ///
		(line max_eduyears_cum_AS_3 eduyears3 if eduyears3>13 & eduyears3<25 , legend(label(3 "3rd quintile"))) ///
		(line max_eduyears_cum_AS_4 eduyears3 if eduyears3>13 & eduyears3<25 , legend(label(4 "4th quintile"))) ///
		(line max_eduyears_cum_AS_5 eduyears3 if eduyears3>13 & eduyears3<25 , legend(label(5 "5th quintile"))), xtitle("Age left full time education") 
		
twoway 	(line diff_cum_1 eduyears3 if eduyears3>13 & eduyears3<=21 , lp(solid) legend(label(1 "1st vs 5th quintile"))) ///
		(line diff_cum_2 eduyears3 if eduyears3>13 & eduyears3<=21 , lp(dash_dot) legend(label(2 "2nd vs 5th quintile"))) ///
		(line diff_cum_3 eduyears3 if eduyears3>13 & eduyears3<=21 , lp(dot) legend(label(3 "3rd vs 5th quintile"))) ///
		(line diff_cum_4 eduyears3 if eduyears3>13 & eduyears3<=21 , lp(dash) legend(label(4 "4th vs 5th quintile"))) ///
			, xtitle("Age left full time education") ytitle("CDF difference") graphregion(color(white)) plotregion(lc(white)) title("Mendelian randomiztion")
graph save "results/figure_X_LATE_MR", replace 
graph combine "results/figure_X_LATE_ROSLA" "results/figure_X_LATE_MR", col(1) imargin(0 0 0 0) ysize(8) xsize(6.26) graphregion(color(white))
graph export "results/figure_X_LATE_ROSLA.pdf", as(pdf) replace fontface("Calibri")
	
