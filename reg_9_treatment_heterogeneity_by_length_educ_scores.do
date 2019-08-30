//Neil Davies 20/10/17
//This creates 2 allele scores for educational attainment

use  "workingdata/full_merged_dataset_SNPs",clear

gen eduyear_16=(eduyears2 >=16) if eduyears2 !=.
gen eduyear_18=(eduyears2 >=18) if eduyears2 !=.
gen eduyear_20=(eduyears2 >=20) if eduyears2 !=.
gen eduyear_21=(eduyears2 >=21) if eduyears2 !=.

joinby n using "workingdata/weights",unmatched(master)

/*
//Okbay SNPs start and end at:
2019
2092
*/

gen score_AO_1=0
gen score_AO_2=0
sort n
forvalues i=2019(1)2056{
	local rsid=rsid[`i']
	local weight=out_effect[`i']
	replace score_AO_1=score_AO_1+`weight'*`rsid'_
	}

forvalues i=2057(1)2092{
	local rsid=rsid[`i']
	local weight=out_effect[`i']
	replace score_AO_2=score_AO_2+`weight'*`rsid'_
	}	
	
	

local i=2018
local j=2018
local rsid="X"
while `j'!=2056{
	local i=`i'+1
	local rsid=rsid[`i']
	local weight=out_effect[`i']
	if "`trait'"!=trait[`i']{
		local trait=trait[`i']
		local j=`j'+1
		di "Trait = `trait'"
		di "Trait ID = `j'"
		}
	cap:replace score_AO_1=score_AO_1+`weight'*`rsid'_
	if _rc==111{
		di "Error SNP `rsid' not found, trait `trait'"
		di "`j'" 
		di "`i'"
		}
	}
