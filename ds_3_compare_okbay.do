//Neil Davies 05/02/18
//This compares the effect of the 162 education variants in Okbay in the full sample compared to results reported in the interim release.


use  "workingdata/full_merged_dataset_SNPs",clear

//Join indicator for interim release
joinby n_eid using "rawdata/biobank_genotype_supp_NMD_150417",unmatched(both)
//Add in phenotype data

forval i = 0/1 {
	forval j = 0/5 {
		di "i=`i', j=`j'"
		//First difference: 568 people with degrees missed here:
		g EA_`i'_`j' = 20 if n_6138_0_0 == 1
		replace EA_`i'_`j' = 13 if n_6138_`i'_`j' == 2
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 3
		replace EA_`i'_`j' = 10 if n_6138_`i'_`j' == 4
		replace EA_`i'_`j' = 19 if n_6138_`i'_`j' == 5
		replace EA_`i'_`j' = 15 if n_6138_`i'_`j' == 6
		replace EA_`i'_`j' = 7 if n_6138_`i'_`j' == -7
		replace EA_`i'_`j' = . if n_6138_`i'_`j' == -3
		}
	}

// take max 
egen EA3 = rmax(EA_*_*)
egen z_eduyears=std(EA3)

gen o_eduyears=1 if EA3 ==7
replace o_eduyears=2 if EA3 ==10
replace o_eduyears=3 if EA3 ==13
replace o_eduyears=4 if EA3 ==15
replace o_eduyears=5 if EA3 ==19
replace o_eduyears=5 if EA3 ==20

//Ordinal scale standardised
egen z_o_eduyears=std(o_eduyears)

//Generate indicator for interim
gen interim=(n_22000_0_0>-12 & n_22000_0_0<23)

reg EA3 pc_*  yob_1-yob_sex_36 cov_male
regsave pc_1 using "results/okbay_replication",replace detail(all) pval ci
forvalues k=0(2)2{
	foreach j in o_eduyears z_eduyears z_o_eduyears EA3{

		#delimit ;

		foreach i in 

		rs301800
		rs11210860
		rs34305371
		rs2568955
		rs1008078
		rs11588857
		rs1777827
		rs2992632
		rs76076331
		rs11689269
		rs1606974
		rs11690172
		rs2457660
		rs10496091
		rs13402908
		rs4851251
		rs12987662
		rs17824247
		rs16845580
		rs4500960
		rs6739979
		rs2245901
		rs55830725
		rs35761247
		rs62259535
		rs11712056
		rs112634398
		rs62263923
		rs6799130
		rs12646808
		rs2610986
		rs34072092
		rs3101246
		rs4863692
		rs4493682
		rs2964197
		rs61160187
		rs324886
		rs10061788
		rs2431108
		rs1402025
		rs62379838
		rs56231335
		rs7767938
		rs2615691
		rs12531458
		rs12671937
		rs113520408
		rs17167170
		rs11768238
		rs12682297
		rs1871109
		rs13294439
		rs895606
		rs7854982
		rs11191193
		rs12772375
		rs7945718
		rs7955289
		rs2456973
		rs7131944
		rs572016
		rs7306755
		rs9537821
		rs1043209
		rs8008779
		rs17119973
		rs55943044
		rs12969294
		rs2837992
		rs165633
		
		//Proxies

		rs17538393 //rs114598875 
		rs9878943 //rs148734725 
		rs1487445  //rs9320913 

		{;
			reg `j' `i' pc_* yob_1-yob_sex_36 cov_male if interim!=`k';
			regsave `i'* using "results/okbay_replication",append detail(all) pval ci;
			};
		};
	};
use "results/okbay_replication",clear	
	
gen rsid=word(subinstr(var,"_"," ",5),1)
gen effect=word(subinstr(var,"_"," ",5),2)
gen other=word(subinstr(var,"_"," ",5),3)

rename coef coef_ukbb_raw_eduyears
rename stderr stderr_ukbb_raw_eduyears

keep rsid depvar effect other coef stderr N
drop if rsid=="pc"
duplicates drop
compress
save "workingdata/UKBB_eduyears",replace

use "workingdata/UKBB_eduyears",clear
duplicates drop

//RSids updated
joinby rsid using  "15_cvd_interaction/rawdata/okbay_results_table_1.15.dta",unmatched(both)

//Harmonize effect alleles
gen flip=(allele1==other)

gen beta_ukbb_interim_harm=beta_ukbb_interim if flip==0
replace beta_ukbb_interim_harm=-1*beta_ukbb_interim if flip==1

gen beta_disc_harm=beta_disc if flip==0
replace beta_disc_harm=-1*beta_disc if flip==1

//Harmonize all variants to be education increasing
gen neg=(beta_disc_harm<0)

replace coef_ukbb_raw_eduyears=coef_ukbb_raw_eduyear*-1 if neg==1
replace beta_disc_harm= beta_disc_harm*-1 if neg==1
replace beta_ukbb_interim_harm= beta_ukbb_interim_harm*-1 if neg==1

log using "logs/comparison_okbay_regressions.txt",text

//Regression of the SSGAC UKBB interim estimates on the full sample
//Raw scale
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="EA3"  & N<100000
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="EA3"  & N>100000

//Standardised scale
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_eduyears" & N<100000
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_eduyears" & N>100000

//Ordinal scale
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="o_eduyears" & N<100000
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="o_eduyears" & N>100000
//Standardised Ordinal scale
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_o_eduyears" & N<100000
reg coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_o_eduyears" & N>100000




//Compared to SSGAC discovery sample
//Raw scale
reg coef_ukbb_raw_eduyears beta_disc_harm if depvar=="EA3"
//Standardised Ordinal scale
reg coef_ukbb_raw_eduyears beta_disc_harm  if depvar=="z_eduyears"
//Ordinal scale
reg coef_ukbb_raw_eduyears beta_disc_harm  if depvar=="o_eduyears"
//Standardised Ordinal scale
reg coef_ukbb_raw_eduyears beta_disc_harm  if depvar=="z_o_eduyears"
log close
//Plot
//Raw scale
twoway scatter coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="EA3"  || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="EA3" 
//Standardised eduyears scale
twoway scatter coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_eduyears" || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_eduyears" 
//Ordinal scale
twoway scatter coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="o_eduyears" || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="o_eduyears"
//Standardised Ordinal scale
twoway scatter coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_o_eduyears" || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_o_eduyears"

//Compared to SSGAC discovery sample
//Raw scale
twoway scatter  coef_ukbb_raw_eduyears beta_disc_harm if depvar=="EA3" || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="EA3" 
//Standardised eduyears scale
twoway scatter  coef_ukbb_raw_eduyears beta_disc_harm  if depvar=="z_eduyears" || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_eduyears" 
//Ordinal scale
twoway scatter  coef_ukbb_raw_eduyears beta_disc_harm  if depvar=="o_eduyears" || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="o_eduyears"
//Standardised Ordinal scale
twoway scatter  coef_ukbb_raw_eduyears beta_disc_harm  if depvar=="z_o_eduyears" || lfit  coef_ukbb_raw_eduyears beta_ukbb_interim_harm if depvar=="z_o_eduyears"

