//Neil Davies 05/09/17
//This generates the east-west north-south plot of the UK Biobank birth location


use n_eid n_129_?_0 n_130_?_0 using "stata/data.8434.dta",clear

replace n_129_0_0=n_129_1_0 if  n_129_0_0==.
replace n_129_0_0=n_129_2_0 if  n_129_0_0==.

replace n_130_0_0=n_130_1_0 if  n_130_0_0==.
replace n_130_0_0=n_130_2_0 if  n_130_0_0==.

twoway scatter  n_129_0_0 n_130_0_0 , msy(smcircle) mco(black) graphregion(color(white)) plotregion(color(white)) msize(vtiny)
