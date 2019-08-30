//Neil Davies 02/08/17
//This takes all of the SNPs extracted from MR-Base and checks whether they are in the HRC

//First load HRC data
//Consists of 40,405,505 SNPs in the reference panel
import delimited "rawdata/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf", delimiter(tab) varnames(1) rowrange(50) encoding(ISO-8859-1)
rename v2 POS 
rename fileformatvcfv41 CHR
rename v3 rsid
rename v4 REF
rename v5 ALT
rename v6 QUAL
rename v7 FILTER
rename v8 INFO
drop in 1
compress
repace CHR="23" if chr=="X"
replace CHR="23" if chr=="X"
replace CHR="23" if CHR=="X"
destring CHR, replace
compress
drop QUAL FILTER INFO

//Drop the duplicate SNPs
bys CHR POS : keep if _n==1

//Save to RDSF
save "workingdata/hrc_reference_panel.dta",replace
//Save to desktop (quick access)
save "uk_biobank/hrc_reference_panel.dta",replace

//Next load up the list of SNPs extracted from MR-base
//Consists of 5,626 SNPs
use "workingdata/coefficients.dta", clear
rename SNP rsid
joinby rsid using "uk_biobank/hrc_reference_panel.dta" ,unmatched(master)

gen in_hrc=(_m==3)
drop _merge POS-INFO
duplicates drop

//343 SNPs are not in the HRC.
//These will be dropped from the analysis
drop if in_hrc!=1
compress
sort rsid

save "workingdata/coefficients_hrc",replace

//Next check whether these SNPs are in the UK Biobank data full release release
//Clean UKB SNP-STATS files
forvalues i=13(1)13{
	if `i'<10{
		local i="0`i'"
		}
	import delimited "rawdata/data.chr`i'.snp-stats", delimiter(tab) varnames(11) rowrange(11) encoding(ISO-8859-1)clear
	sort rsid
	save "uk_biobank/ukb_snp_stats_chr`i'.dta",replace
	}
	
//Merge on the UKBB SNPSTATS data	
use "workingdata/coefficients_hrc",clear
//Total SNPs 5,288

forvalues i=1(1)22{
	if `i'<10{
		local i="0`i'"
		}
	preserve
	keep if CHR==`i'
	joinby rsid using "uk_biobank/ukb_snp_stats_chr`i'.dta",unmatched(none) 
	save "uk_biobank/temp_`i'.dta",replace
	restore
	}
//Fix INFO variable to be a string
forvalues i=1(1)22{
	if `i'<10{
		local i="0`i'"
		}	
	use "uk_biobank/temp_`i'.dta",clear
	tostring info impute_info hw_lrt_p_value, replace force
	save "uk_biobank/temp_`i'.dta",replace
	}
	
//Append files
use "uk_biobank/temp_01.dta",clear
forvalues i=2(1)22{
	if `i'<10{
		local i="0`i'"
		}	
	tostring hw_lrt_p_value,replace
	append using "uk_biobank/temp_`i'.dta",force
	}
	
//There are 25 triallelic SNPs
//Drop the rarest allele
bys rsid exposure (minor_allele_frequency): keep if _N==_n
 
//One cognition SNP was not in UKBB total SNPs =5,282

//Harmonize effect allele

//Check for effects with the correct direction	
gen correct=0
replace correct=1 if effect_allele_exposure==allelea &  other_allele_exposure==alleleb 
 
//Check for SNPs which are coded in the opposite direction
gen flipped=0
replace flipped=1 if   other_allele_exposure==allelea &  effect_allele_exposure==alleleb 

count if correct==0 & flipped==0

//No SNPs are not correct or flipped

//Create indicator for UK Biobank effect allele frequency (for consistency with MR-Base)
gen ukb_maf_ea=minor_allele_frequency if minor_allele==allelea
replace ukb_maf_ea=1-minor_allele_frequency if minor_allele==alleleb

//Check ambigous SNPs
//ambiguous if has the following alleles A/T and C/G 
gen ambiguous=((allelea=="A"&alleleb=="T")|(allelea=="T"&alleleb=="A")|(allelea=="C"&alleleb=="G")|(allelea=="G"&alleleb=="C"))

//256 ambiguous alleles with MAF>0.3 in biobank

//Create new summary data harmonized to the UK Biobank effect alleles
gen out_effect=beta_exposure if correct==1
replace out_effect=beta_exposure*-1 if flipped==1

gen out_effectaf=eaf_exposure if correct==1
replace  out_effectaf=1-eaf_exposure if flipped==1

gen out_effect_allele=allelea
gen out_other_allele=alleleb

gen out_se=se_exposure

//Effect allele frequency correlated 0.9913
corr out_effectaf ukb_maf_ea 
gen diff=out_effectaf-ukb_maf_ea 
sum diff 

//The cognition SNPs allele frequency is a MAF, not a EAF.
//This is not a problem because the alleles are correct.
corr out_effectaf ukb_maf_ea if exposure!="Cognition p<5e-05"
corr out_effectaf ukb_maf_ea if exposure=="Cognition p<5e-05"

twoway scatter out_effectaf ukb_maf_ea  if exposure=="Cognition p<5e-05"
twoway scatter out_effectaf ukb_maf_ea  if exposure!="Cognition p<5e-05"
	
//The next section removes any SNPs which are within a 500kb region of the 74 Okbay education SNPs.
gen okbay_snps=(id_exposure=="1001" & pval_exposure <5*10^-8)
gsort - okbay_snps pval_exposure

gen snp_proximal_okbay=0

//Loop through each of the Okbay SNPs to check the remaining hits on all the other traits to see if they're proximal to the Okbay SNPs
forvalues i=1(1)74{
	local okbay_pos_min=position[`i']-500000
	local okbay_pos_max=position[`i']+500000
	local okbay_chr=chromosome[`i']
	replace snp_proximal_okbay=1 if `okbay_chr'==chromosome & (`okbay_pos_min'<position & `okbay_pos_max'>position)
	}

//Drop the 365 SNPs that are within 500kb of the Okbay SNPs
drop if snp_proximal_okbay==1 & exposure!="Eduyears (Okbay discovery) p<5e-08"& exposure!="Years of schooling || SSGAC || 2016 || SD (years)" ///
			& exposure!="Years of schooling || SSGAC || 2013 || SD (years)"

//Create the score weight file
//This is file of all the RSids, their weight and the effect allele
//Space delimited

compress


//Keep the RSID coef SE effect and other allele and save the coefficient dataset
keep rsid out_effect out_se pval_exposure out_effectaf out_effect_allele out_other_allele exposure
rename pval_exposure  pvalue
rename exposure trait
sort trait
gen n=_n
compress

//Note of the 74 Okbay SNPs, 5 were not in the HRC. These have been replaced by proxies in perfect LD.
/*
rs114598875   > Potential proxy rs17538393
rs148734725   > Potential proxy rs9878943
rs192818565   > Potential proxy rs8008779
  rs8005528   > Potential proxy rs55943044
  rs9320913   > Potential proxy rs1487445
*/

save "workingdata/weights",replace
