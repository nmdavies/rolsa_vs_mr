//Neil Davies 08/05/17
//This constructs the scores using the coefficients from MR-Base

set maxvar 32000
use "$path2/rawdata/final_full_snps_clean.dta",clear
drop _m

joinby n using "workingdata/weights",unmatched(master)

//There are a number of triallelic SNPs
#delimit ;
foreach i in 
rs9907241
rs1659127
rs715040
rs2015436
rs42460
rs142581044
rs16960696
rs13045716
rs1659127
rs12986437
rs11177669
rs16865791
rs4867149
rs1601934
rs13290890
rs9907241
rs415742
rs2843744
rs281288
rs11043207
rs16931824
rs11687170
rs1010587
rs2772572{;
	sum `i'_*;
	};
#delimit cr	
//Drop the rarest allele
drop rs2772572_C_A  rs1010587_G_A rs16931824_T_A rs11043207_C_T rs281288_G_A rs2843744_T_G rs415742_C_A ///
	rs13290890_C_A rs1601934_G_T rs4867149_G_T rs11177669_G_C rs12986437_A_G rs13045716_C_T rs16960696_T_C ///
	rs142581044_G_A rs42460_G_C rs2015436_T_A rs715040_G_A rs11687170_T_A rs16865791_C_T


//Check the allele frequencies
gen ukb_eaf=.
gen ukb_effect=""
gen ukb_other=""

local rsid="X"
cap{
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	cap:ds `rsid'_
	di "`i'"
	di "`rsid'"
	if _rc!=111{
		replace ukb_effect=word(subinstr("`r(varlist)'", "_", " ",2) ,2) in `i'
		replace ukb_other=word(subinstr("`r(varlist)'", "_", " ",2) ,3) in `i'
		cap:tabstat `rsid'_, save
		replace ukb_eaf=el(r(StatTotal),1,1)/2 in `i'
		}
	}
}
	
corr out_effectaf ukb_eaf 
corr out_effectaf ukb_eaf if trait !="Cognition p<5e-05"
//Same as pre-extraction QC	
	
//Program for creating allele scores
levels trait
foreach i in `r(levels)'{
	local j =`j'+1
	gen allele_score_`j'=0
	}
local j=0
local rsid="X"
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local weight=out_effect[`i']
	if "`trait'"!=trait[`i']{
		local trait=trait[`i']
		local j=`j'+1
		di "Trait = `trait'"
		di "Trait ID = `j'"
		}
	cap:replace allele_score_`j'=allele_score_`j'+`weight'*`rsid'_
	if _rc==111{
		di "Error SNP `rsid' not found, trait `trait'"
		di "`j'" 
		di "`i'"
		}
	}

//Label allele scores
local rsid="X"
local j=0
while "`rsid'"!=""{
	local i=`i'+1
	local rsid=rsid[`i']

	if "`trait'"!=trait[`i']{
		local trait=trait[`i']
		local j=`j'+1
		di "Trait = `trait'"
		di "Trait ID = `j'"
		label var allele_score_`j' "`trait'"
		}
	
	}	

preserve
keep n_eid allele_score_*
sum allele_score_*
compress

save "workingdata/full_sample_allele_scores",replace
