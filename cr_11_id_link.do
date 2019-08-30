//Neil Davies 25/07/17
//This cleans the link file

import delimited "rawdata/ukb878_imp_chr10_v2_s487406.sample", delimiter(space) varnames(1) encoding(ISO-8859-1)
drop in 1
gen n=_n
compress
save "workingdata/full_release_ids.dta"

