//Neil Davies 24/07/17
//This imports the individual SNP data into Stata and cleans it into additive format

cap prog drop import_snps
prog def import_snps
args file
import delimited "`file'.", delimiter(space, collapse) encoding(ISO-8859-1) clear

local snp=v[3]
local effect=v[5]
local other=v[6]

drop if _n<7 
gen id=round(_n/3+1/6)
bys id: gen n=_n
drop if n==3

destring v, replace
replace v=v*2 if n==1 

bys id: egen `snp'_`effect'_`other'=total(v)
drop if n==2
drop v

drop n

save "`file'_clean.dta",replace

end



fs *_trans
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}

//Merge all of the files back togeather
forvalues k=1(1)19{
	use final_snps_`k'.gen_split0000_trans_clean.dta,clear
	fs *snps_`k'.*clean.dta
	foreach i in `r(files)'{
		local j=`j'+1
		di "`j'"
		joinby id using `i' 
		}
	save 	final_snps_`k'.gen_trans_clean.dta,replace
	}

//Check which SNPs have failed to be extracted:
forvalues k=1(1)19{
	use final_snps_`k'.gen_trans_clean.dta,clear
	display c(k)
	}	

forvalues k=1(1)19{
	use final_snps_`k'.gen_trans_clean.dta,clear
	gen rsid=""
	ds *
	foreach i in `r(varlist)'{
		local j=`j'+1
		replace rsid="`i'" in `j'
		}
	keep rsid
	drop if rsid==""
	save rsid_`k',replace
	}		
use rsid_1,clear
gen chr=1
forvalues k=2(1)19{
	append using rsid_`k'
	replace chr=`k' if chr==.
	}
drop if rsid =="id"
gen rsid2=subinstr(rsid, "_", " ",.) 
gen rsid3=word(rsid2,1)
replace rsid=rsid3
drop rsid2 rsid3

save temp2
joinby rsid using "workingdata/weights.dta", unmatched(both)

//Have extracted 369 additional SNPs 
//All other SNPs extracted successfully
//Merge togeather

use final_snps_19.gen_trans_clean.dta,clear
compress
forvalues i=16(1)18{
	joinby id using final_snps_`i'.gen_trans_clean.dta,
	compress
	}
	
save temp,replace

//Clean the ID file
import delimited "rawdata/ukb878_imp_chr10_v2_s487406.sample", delimiter(space) encoding(ISO-8859-1)clear
drop in 1
rename id_1 n_eid
drop id_2
drop missing 
gen n=_n
save "workingdata/sample",replace

use temp,clear
joinby n using "workingdata/sample",unmatched(both)

save "final_full_snps_clean.dta",replace
