//Neil Davies 23/06/17
//This extracts the 18 lead IQ SNPs from the Sniekers et al GWAS results

import delimited "gwas_results/iq - sniekers - full gwas.txt", delimiter(space) encoding(ISO-8859-1)clear

joinby rsid using "gwas_results/sniekers_lead_snps.dta", unmatched(both)
keep if _m==3
compress
save "gwas_results/sniekers_lead_snps_full.dta",replace
use "gwas_results/sniekers_lead_snps_full.dta",clear
drop n
gen n=_n
save "gwas_results/sniekers_lead_snps_full.dta",replace
