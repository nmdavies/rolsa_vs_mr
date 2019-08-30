#Neil Davies 15/03/17
#This runs the MR-Base scripts for the education MR paper

#Load MR base packages
library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(foreign)
#Load TwoSampleMR package
library(TwoSampleMR)

# List available GWASs
ao <- available_outcomes()#Neil Davies 15/03/17
#This runs the MR-Base scripts for the education MR paper

#Install MR base instruments
devtools::install_github("MRCIEU/MRInstruments")

library(MRInstruments)

#Get GWAS catalog
data(gwas_catalog)
head(gwas_catalog)

#List available traits in the GWAS catalog
subset(ao, select=c(trait, id))

#Create list of IDs for all phenotypes
a<-subset(ao, select=c(id))

#Extract IDs, coefficients and effect alleles for the covariate allele scores
#Because we're not using these scores as instruments, but to assess bias, we're using a lower threshold.
#This is particularly important for phenotypes that we do not have any genomewide significant hits for.

phenotypes <-c(89,835,999,961,962,964,813,1083,29,28,16,1095,1000,805,806,22,298,801,803,1019,1092,1093,113,114,115,117,118,1028,1029,1001,1009,1087,1088,1061,1062,1063,1064,1065,1066,1067,855,856,857,893,916,917,1079,837)
instruments<-extract_instruments(outcomes=phenotypes, p1 = 5e-05)

#Save instruments file
write.table(instruments,file="workingdata/mr_base_extract.rda")
instruments<-read.table(file="workingdata/mr_base_extract.rda")

#Read in Okbay 540 SNPs
file<-"rawdata/EduYears_Discovery_5000_okbay_ex_UK_Biobank.txt"
okbay_disc_exp_dat<-read.table(file, header = TRUE)
okbay_coefficients<-clump_data(format_data(subset(okbay_disc_exp_dat, Pval<5e-05), type="exposure", snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", eaf_col = "EAF", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "Pval"))
okbay_coefficients2<-clump_data(format_data(subset(okbay_disc_exp_dat, Pval<5e-08), type="exposure", snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", eaf_col = "EAF", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "Pval"))
okbay_coefficients3<-format_data(subset(okbay_disc_exp_dat, Pval<5e-08), type="exposure", snp_col = "MarkerName", beta_col = "Beta", se_col = "SE", eaf_col = "EAF", effect_allele_col = "A1", other_allele_col = "A2", pval_col = "Pval")

okbay_coefficients3$exposure<-"Eduyears (Okbay discovery) p<5e-08"

instruments_combine <- rbind(subset(instruments, select=c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure" )), okbay_coefficients3)

#Read in the 373 lead IQ SNPs from the Sniekers et al GWAS results
file<-"gwas_results/iq - sniekers - full gwas.txt"
snieker_disc_exp_dat<-read.table(file, header = TRUE)
snieker_coefficients <-clump_data(format_data(subset(snieker_disc_exp_dat, p_value<5e-08), type="exposure", snp_col = "rsid", beta_col = "Beta", se_col = "SE", eaf_col = "MAF", effect_allele_col = "ref", other_allele_col = "alt", pval_col = "p_value"))
snieker_coefficients$exposure<-"Cognition p<5e-05"

#This results in 373 IQ SNPs p<5e-05", there are 17 at p<5e-08

#Combining with the other SNPs from MR-Base
instruments_combine <- rbind(subset(instruments_combine, select=c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure" )), snieker_coefficients)

#Read in intelligence SNPs from the Zabaneh and Trampush papers
zabaneh_disc_exp_dat<-head(read.table("gwas_results/zabaneh.txt",header=TRUE),122)
zabaneh_formated<-format_data(zabaneh_disc_exp_dat,type="exposure", snp_col = "SNP", beta_col = "Odds.ratio", se_col = "S.E.", eaf_col = "SSGACallelefreq", effect_allele_col = "Allele1", other_allele_col = "Allele2", pval_col = "P.value")
zabaneh_formated$beta.exposure=log(zabaneh_formated$beta.exposure)
zabaneh_formated$exposure="High IQ Zabaneh"

instruments_combine <- rbind(subset(instruments_combine, select=c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure" )), zabaneh_formated)

tram_disc_exp_dat<-head(read.table("gwas_results/trampush.txt",header=TRUE),122)
tram_coefficients <-clump_data(format_data(tram_disc_exp_dat,type="exposure", snp_col = "SNP", beta_col = "Beta", se_col = "S.E.", eaf_col = "SSGACallelefreq", effect_allele_col = "Allele1", other_allele_col = "Allele2", pval_col = "P.value"))
tram_coefficients$exposure="IQ Trampush p<5E-07"
instruments_combine <- rbind(subset(instruments_combine, select=c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure" )), tram_coefficients)

proxies<-read.table("gwas_results/okbay_proxies.txt", header=T)
proxies$exposure<-"Eduyears (Okbay discovery) p<5e-08"
instruments_combine <- rbind(subset(instruments_combine, select=c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","exposure","mr_keep.exposure","pval_origin.exposure","id.exposure" )), proxies)

#Create list of SNPs to extract from Biobank
#Select RS IDs
write.table(subset(instruments_combine,select=c(SNP)), "workingdata/snps.txt", sep="\t" , quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(instruments_combine, "workingdata/coefficients.txt", sep="\t" , quote=TRUE, row.names = FALSE, col.names = TRUE)
write.dta(as.data.frame(instruments_combine), "workingdata/coefficients.dta", convert.factors = "string")

#Copy this file across to BlueCyrstal and extract all 5631 SNPs from the UK Biobank imputed data
#Uses the following scripts 
# ./uk_biobank/scripts/qctool_submission.sh
# ./uk_biobank/scripts/qctool_all_snps_qsub.sh   