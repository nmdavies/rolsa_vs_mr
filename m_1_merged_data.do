//Neil Davies 08/05/17
//This merges the allele scores to the ROSLA dataset and defines the exclusions for cryptic relatedness and ethnicity

//Clean exclusions list non-europeans
import delimited "rawdata/data.non_europeans_exclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/non-europeans",replace

//Europeans
import delimited "rawdata/data.europeans_inclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/europeans",replace

//Sex mismatches, sex chromosome aneuploidy, and excess heterogeneity
import delimited "rawdata/meta.recommended_exclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/exclusions",replace

//Relateds
import delimited "rawdata/data.relateds_exclusions.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
save "workingdata/relateds",replace


//PCs
import delimited "rawdata/data.pca1-10.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2
forvalues i=3(1)12{
	local j=`i'-2
	rename v`i' pc_`j'
	}
save "workingdata/PCs",replace

//Standard covariates
import delimited "rawdata/data.covariates.txt", delimiter(space) varnames(1) encoding(ISO-8859-1)clear
rename v1 n_eid
drop v2 m
save "workingdata/sample",replace

//Load up the allele scores
use "workingdata/full_sample_allele_scores",clear

//Recomended exclusions
joinby n_eid using  "workingdata/exclusions",unmatched(master)
drop if _m==3
drop _m

//Non-europeans
joinby n_eid using  "workingdata/non-europeans",unmatched(master)
drop if _m==3
drop _m

//Europeans
joinby n_eid using  "workingdata/europeans",unmatched(master)
drop if _m!=3
drop _m

//Relateds
joinby n_eid using   "workingdata/relateds",unmatched(master)
drop if _m==3
drop _m

//PCs
joinby n_eid using "workingdata/PCs",unmatched(master)
drop _m

//Join phenotype data
joinby n_eid using "workingdata/cleaned_biobank_outcomes.dta",unmatched(master)
drop _m

compress

//Merge in additional phenotypic baseline variables

joinby n_eid using  "workingdata/cov_dist_long",unmatched(master)
drop _m
joinby n_eid using   "workingdata/birth_location_imd_rural_urban",unmatched(master)
rename imd cov_birth_location_imd
drop dist
rename cov_urban cov_birth_location_urban
drop _merge
joinby n_eid using   "workingdata/birth_location",unmatched(master)
drop _m
rename n_129_0_0 cov_birth_location_northing
rename n_130_0_0 cov_birth_location_easting

joinby n_eid using  "workingdata/current_location",unmatched(master)
drop _m
rename n_20074_0_0 cov_current_location_easting
rename n_20075_0_0 cov_current_location_northing
replace cov_current_location_easting=n_20074_1_0 if cov_current_location_easting==.
replace cov_current_location_northing=n_20075_1_0 if cov_current_location_northing==.

drop data_error
compress

//Define weights
gen weight=1.8857 if more_edu_15==0
replace weight=1 if weight==.

//Gen updated weights from Hughes et al.
gen weight_hughes=34.29 if weight!=1
replace weight_hughes=12.37 if weight==1

//UK specific eduyears
gen eduyears2=21 if eduyears==20
replace eduyears2=20 if eduyears==19
replace eduyears2=18 if eduyears==15
replace eduyears2=16 if eduyears==13
replace eduyears2=15 if eduyears==10|eduyears==7
compress

//Final exclusion list
joinby n_eid using "rawdata/exclusions_170726.dta", unmatched(master)
drop if _m==3
drop _m
compress

tab yob, gen(yob_)
gen yob_sex=cov_male*yob
tab yob_sex, gen(yob_sex_)
compress
save "workingdata/full_merged_dataset",replace

use "workingdata/full_merged_dataset",clear

//Merge on SNPs
joinby n_eid using "uk_biobank/full_release/final_full_snps_clean.dta",unmatched(both) _merge(XXX)

//Keep the 74 Okbay SNPs, and the cognition SNPs:
#delimit ;
foreach i in  
rs10061788
rs1008078
rs1043209
rs10496091
rs11191193
rs11210860
rs112634398
rs113520408
rs17538393
rs11588857
rs11689269
rs11690172
rs11712056
rs11768238
rs12531458
rs12646808
rs12671937
rs12682297
rs12772375
rs12969294
rs12987662
rs13294439
rs13402908
rs1402025
rs9878943
rs1606974
rs165633
rs16845580
rs17119973
rs17167170
rs1777827
rs17824247
rs1871109
rs2245901
rs2431108
rs2456973
rs2457660
rs2568955
rs2610986
rs2615691
rs2837992
rs2964197
rs2992632
rs301800
rs3101246
rs324886
rs34072092
rs34305371
rs35761247
rs4493682
rs4500960
rs4851251
rs4863692
rs55830725
rs55943044
rs56231335
rs572016
rs61160187
rs62259535
rs62263923
rs62379838
rs6739979
rs6799130
rs7131944
rs7306755
rs76076331
rs7767938
rs7854982
rs7945718
rs7955289
rs8008779
rs895606
rs1487445
rs9537821
rs17393468
rs6669072
rs2297797
rs76114856
rs12500482
rs189977224
rs7640
rs916888
rs10191758
rs10236197
rs10733787
rs11138902
rs12744310
rs12928404
rs13010010
rs1626122
rs16954078
rs2251499
rs2490272
rs2743462
rs41352752
rs4728302
rs6746731
rs6779302
rs7646501
rs78164635{;
	order `i'_;
	};

drop rs185762739_G_A-_merge

save "workingdata/full_merged_dataset_SNPs",replace
