//Neil Davies 20/02/17
//This cleans and extracts the covariate data from the 8434 extract place of birth

use n_eid n_129_?_0 n_130_?_0 using "stata/data.8434.dta",clear

replace n_129_0_0=n_129_1_0 if  n_129_0_0==.
replace n_129_0_0=n_129_2_0 if  n_129_0_0==.

replace n_130_0_0=n_130_1_0 if  n_130_0_0==.
replace n_130_0_0=n_130_2_0 if  n_130_0_0==.

keep n_eid n_129_0_0 n_130_0_0
compress
//There is a data error and 2,716 of birth easting and northings have been divided by 10.
gen data_error=(n_129_0_0<200000 & n_130_0_0<90000 ) & n_129_0_0!=. & n_130_0_0!=.
replace n_129_0_0=n_129_0_0*10 if data_error==1
replace n_130_0_0=n_130_0_0*10 if data_error==1

save "workingdata/birth_location",replace

use n_eid n_20074_?_0 n_20075_?_0 n_20118* using "stata/data.4263.dta",clear

replace n_20074_0_0=n_20074_1_0 if  n_20074_0_0==.
replace n_20075_0_0=n_20075_1_0 if  n_20075_0_0==.

drop n_20074_1_0  n_20075_1_0

save "workingdata/current_location",replace

//Clean geographic data to get urban/rural indicator
//Raw data downloaded from:
//http://ons.maps.arcgis.com/home/item.html?id=dfa0ff74981b4b228d2030d852f0b14a

import delimited "rawdata/ONSPD_FEB_2017_UK/Data/ONSPD_FEB_2017_UK.csv", delimiter(comma) encoding(ISO-8859-1)
keep oseast1m osnrth1m ru11ind imd
rename oseast1m easting
rename osnrth1m northing
duplicates drop 
//5 duplicated northings and eastings keep easting and northing records with urban/rural indicators:
duplicates tag easting northing ,gen(s)
drop if ru11ind=="" & s==1
drop s
gen x=1
save "workingdata/east_north_rural_urban_imd",replace

//Merge onto UK Biobank data
use "workingdata/birth_location",clear
rename n_129_0_0 northing
rename n_130_0_0 easting
drop if northing ==.|easting ==-1 |northing ==-1 |easting ==-1

//There is a data error and 2,716 of birth easting and northings have been divided by 10.
gen data_error=(northing<200000 & easting<90000 ) & northing!=. & easting!=.
replace northing=northing*10 if data_error==1
replace easting=easting*10 if data_error==1
drop data_error

joinby easting northing using "workingdata/east_north_rural_urban_imd", unmatched(master)
keep if _m==3
drop _m
//Successfully merged most
save "workingdata/birth_location_imd_merged",replace

//Many eastings and northings are not at the centre of a postcode as they've been rounded to the nearest 100.
//Match in nearest IMD/rural/urban indicator

use "workingdata/east_north_rural_urban_imd",clear
gen easting_round=round(easting,100)
gen northing_round=round(northing,100)
compress
save "east_north_rural_urban_imd_round100",replace

use "workingdata/east_north_rural_urban_imd",clear
gen easting_round=round(easting,1000)
gen northing_round=round(northing,1000)
compress
save "east_north_rural_urban_imd_round1000",replace

use "workingdata/east_north_rural_urban_imd",clear
gen easting_round=round(easting,10000)
gen northing_round=round(northing,10000)
compress
save "east_north_rural_urban_imd_round10000",replace

use "workingdata/east_north_rural_urban_imd",clear
gen easting_round=round(easting,100000)
gen northing_round=round(northing,100000)
compress
save "east_north_rural_urban_imd_round100000",replace


use "workingdata/birth_location",clear
rename n_129_0_0 northing
rename n_130_0_0 easting

//There is a data error and 2,716 of birth easting and northings have been divided by 10.
gen data_error=(northing<200000 & easting<90000 ) & northing!=. & easting!=.
replace northing=northing*10 if data_error==1
replace easting=easting*10 if data_error==1
drop data_error

drop if northing ==.|easting ==-1 |northing ==-1 |easting ==-1

joinby easting northing using "workingdata/east_north_rural_urban_imd", unmatched(master)
keep if _m==1

gen easting_round=round(easting,100)
gen northing_round=round(northing,100)
drop _m

rename northing northing_bb
rename easting easting_bb
drop r imd
joinby easting_round northing_round using "east_north_rural_urban_imd_round100",unmatched(master)
preserve
keep if _m==1
save "workingdata/unmatched_imd_rural_urban_round100",replace
restore
keep if _m==3
gen dist=((northing_bb-northing)^2+(easting_bb-easting)^2)^0.5
sort n_eid dist
bys n_eid (dist):keep if _n==1

keep n_eid r imd dist
save "workingdata/birth_location_imd_merged_sixfigure",replace

//Match the remaining individuals who were not matched with 100m
use "workingdata/unmatched_imd_rural_urban_round100",clear

replace easting_round=round(easting_bb,1000)
replace northing_round=round(northing_bb,1000)

drop _m-imd
joinby easting_round northing_round using "east_north_rural_urban_imd_round1000",unmatched(master)

//All bar 8,782 participants were matched
preserve
keep if _m==1
save "workingdata/unmatched_imd_rural_urban",replace
restore
keep if _m==3
gen dist=((northing_bb-northing)^2+(easting_bb-easting)^2)^0.5
sort n_eid dist
bys n_eid (dist):keep if _n==1

keep n_eid r imd dist
save "workingdata/birth_location_imd_merged_fourfigure",replace

//Match remaining at 10000m
use "workingdata/unmatched_imd_rural_urban",clear
replace easting_round=round(easting_bb,10000)
replace northing_round=round(northing_bb,10000)

drop _m-imd
joinby easting_round northing_round using "east_north_rural_urban_imd_round10000",unmatched(master)

//All bar 2,808 participants were matched
preserve
keep if _m==1
save "workingdata/unmatched_imd_rural_urban_10000",replace
restore
keep if _m==3
gen dist=((northing_bb-northing)^2+(easting_bb-easting)^2)^0.5
sort n_eid dist
bys n_eid (dist):keep if _n==1

keep n_eid r imd dist
save "workingdata/birth_location_imd_merged_threefigure",replace

append using "workingdata/birth_location_imd_merged_fourfigure",
append using "workingdata/birth_location_imd_merged_sixfigure",
compress

//Replace values of rural/urban and IMD equal to missing if they were more than 1km from the centre of the relevant postcode
drop if dist>1000
rename ru ruralurban

gen cov_urban=inlist(ruralurban,"A1","A2","B1","B2","C1","C2")
replace cov_urban=1 if inlist(ruralurban,"1","2","3","4","5")
drop ruralurban
compress
save "workingdata/birth_location_imd_rural_urban",replace

//Finally create measure for distance from London (Easting 530310, Northing 179536)
use "workingdata/birth_location",clear

//There is a data error and 2,716 of birth easting and northings have been divided by 10.
rename n_129_0_0 northing
rename n_130_0_0 easting
gen data_error=(northing<200000 & easting<90000 ) & northing!=. & easting!=.
replace northing=northing*10 if data_error==1
replace easting=easting*10 if data_error==1
drop data_error

gen cov_dist_lon=((northing-179536)^2+(easting-530310)^2)^0.5/1000 if northing>0 & easting>0
keep n_eid	cov
compress
save "workingdata/cov_dist_long",replace
