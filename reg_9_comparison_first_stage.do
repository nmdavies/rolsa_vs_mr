//Neil Davies 26/07/17
//This estimates the association of each of the Okbay SNPs in the full UKBB release

use "full_release/final_full_snps_clean.dta",clear
drop _m
joinby n_eid using  "workingdata/full_merged_dataset",unmatched(master)

keep if _m==3
drop allele_score*
drop _m
joinby n_eid using  "workingdata/full_sample_allele_scores",unmatched(master)

//Estimate association of Okbay lead SNPs and the eduyears
//Include: Array dummy, year-of-birth dummies, year-of-birth by sex interactions, 15 PCs. 		
tab yob, gen(yob_)
gen yob_sex=cov_male*yob
tab yob_sex, gen(yob_sex_)

/*
currently missing SNPs

*/


#delimit ;
foreach i in 
rs2568955
rs114598875
rs13402908
rs4851251
rs6739979
rs35761247
rs62259535
rs148734725
rs112634398
rs324886
rs56231335
rs12772375
rs8008779
rs55943044
rs301800
rs11210860
rs34305371
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
rs12987662
rs17824247
rs16845580
rs4500960
rs2245901
rs55830725
rs11712056
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
rs10061788
rs2431108
rs1402025
rs62379838
rs9320913
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
rs7945718
rs7955289
rs2456973
rs7131944
rs572016
rs7306755
rs9537821
rs1043209
rs17119973
rs12969294
rs2837992
rs165633
rs17538393
rs9878943
rs1487445
rs8008779
rs55943044
rs916888
rs6669072
rs17393468
rs12500482
rs189977224
rs7640
rs2297797
rs76114856
rs4962322
rs2490272 
rs10236197
rs2251499
rs2743462 rs7646501 rs4728302 rs10191758 rs12744310 rs12928404 rs41352752 rs13010010 rs16954078 rs11138902 rs6746731 rs6779302 rs10733787 rs1626122 rs78164635 {;
	order `i'_;
	};
#delimit cr
drop rs185762739_G_A-rs35666724_A_T

//Drop the three non-HRC SSGAC SNPs
drop rs114598875 rs148734725 rs9320913

gen eduyears3=1 if eduyears==7
replace eduyears3=2 if eduyears==10
replace eduyears3=3 if eduyears==13
replace eduyears3=4 if eduyears==15
replace eduyears3=5 if eduyears==19
replace eduyears3=5 if eduyears==20

tabstat rs55943044_G_A -rs2568955_T_C if yob !=. & pc_1!=. & eduyears!=.,c(s) label(40) var(32) stats(mean n)

tabstat rs4962322 if yob !=. & pc_1!=. & eduyears!=.,c(s) label(40) var(32) stats(mean n)
tabstat rs916888 rs6669072 rs17393468 rs12500482 rs189977224 rs7640 rs2297797 rs76114856 if yob !=. & pc_1!=. & eduyears!=.,c(s) label(40) var(32) stats(mean n)

tabstat rs2490272 rs10236197 rs2251499 rs2743462 rs7646501 rs4728302 rs10191758 rs12744310 rs12928404 rs41352752 rs13010010 rs16954078 rs11138902 rs6746731 rs6779302 rs10733787 rs1626122 rs78164635 if yob !=. & pc_1!=. & eduyears!=.,c(s) label(40) var(32) stats(mean n)


reg eduyears rs55943044_G_A yob_* pc_*
regsave rs55943044_G_A using "results/first_stage",replace detail(all) 

ds rs55943044_G_A- rs2568955_T_C

foreach i in `r(varlist)'{
	reg eduyears `i' yob_* pc_*
	regsave `i' using "results/first_stage",replace detail(all) append
	reg eduyears2 `i' yob_* pc_*
	regsave `i' using "results/first_stage",replace detail(all) append
	reg eduyears3 `i' yob_* pc_*
	regsave `i' using "results/first_stage",replace detail(all) append
	}

ds rs916888 rs6669072 rs17393468 rs12500482 rs189977224 rs7640 rs2297797 rs76114856	
foreach i in `r(varlist)'{
	reg out_intell `i' yob_* pc_*
	regsave `i' using "results/first_stage",replace detail(all) append
	}
ds rs2490272 rs10236197 rs2251499 rs2743462 rs7646501 rs4728302 rs10191758 rs12744310 rs12928404 rs41352752 rs13010010 rs16954078 rs11138902 rs6746731 rs6779302 rs10733787 rs1626122 rs78164635
foreach i in `r(varlist)'{
	reg out_intell `i' yob_* pc_*
	regsave `i' using "results/first_stage",replace detail(all) append
	}
		
use "results/first_stage",clear
drop if strpos(var,"rs114598875")!=0 | strpos(var,"rs148734725")!=0 | strpos(var,"rs9320913")!=0

gen SSGAC_SNPs=var
replace SSGAC_SNPs="rs114598875_A_G" if var=="rs17538393_G_T"
replace SSGAC_SNPs="rs148734725_G_A" if var=="rs9878943_G_A"
replace SSGAC_SNPs="rs9320913_C_A" if var=="rs1487445_C_T"




rename var var1 
rename SSGAC_SNPs var

joinby var using "workingdata/ssgac_results", unmatched(both)

//Harmonize effect allele to positive SSGAC
replace coef=coef*-1 if substr(var,-3,1)!=effect_allele

//Plot MR-Egger vs SSGAC
mreggerplot   coef stderr ssgac_beta ssgac_se   if depvar =="eduyears3" , gpci ivw nolci graphregion(color(white))  plotregion(lc(white)) 
mreggerplot   coef stderr ukbi_beta ukbi_se   if depvar =="eduyears" , gpci ivw nolci graphregion(color(white))  plotregion(lc(white)) leg(off)
mreggerplot   coef stderr ukbi_beta ukbi_se   if depvar =="eduyears3" , gpci ivw nolci graphregion(color(white))  plotregion(lc(white)) leg(off)

mregger   coef  ukbi_beta    if depvar =="eduyears3" [aweight=stderr],gxse(ssgac_se) ivw
