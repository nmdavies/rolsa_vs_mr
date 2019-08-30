//Neil Davies 03/08/17
//Selection of education proxies for the 3 Okbay SNPs not in the HRC

use "workingdata/hrc_reference_panel.dta",clear

//rs114598875

gen proxy_rs114598875=0
gen proxy_rs148734725=0
gen proxy_rs9320913=0

//Proxies for rs114598875 as reported by https://analysistools.nci.nih.gov/LDlink/
//Limited to proxies with an r-squared of 1
#delimit ;
foreach i in rs114598875
rs141577663
rs116096252
rs115559740
rs114836262
rs199590928
rs149117522
rs35028795
rs6761132
rs10201111
rs58291450
rs10196947
rs35252652
rs10172614
rs11887292
rs13030367
rs13398041
rs13411655
rs7562375
rs35905145
rs67895117
rs34977426
rs13011212
rs13003538
rs68175272
rs13003349
rs75246165
rs12996706
rs12104783
rs35836454
rs55771126
rs70959857
rs13389728
rs13027741
rs10183870
rs57705389
rs56949301
rs13417624
rs67779326
rs10206990
rs13028000
rs10186916
rs13416936
rs17538393{;
	count if rsid=="`i'";
	replace proxy_rs114598875=1 if rsid=="`i'";
	};
	
#delimit ;
foreach i in 
rs148734725
rs34588335
rs67485053
rs13087851
rs6809216
rs1987628
rs7623659
rs7648841
rs3811699
rs6446264
rs1800668
rs1050450
rs11711536
rs13090388
rs35169793
rs17080528
rs6793308
rs199933423
rs111903592
rs13086611
rs71627384
rs11716974
rs71627385
rs9878943
rs71080505
rs9871380
rs11706370
rs6779524
rs6997
rs9814873
rs10640
rs11715915
rs9859556
rs35024166
rs11922013
rs6446272
rs7646366
rs6446277
rs35115732
rs13079643
rs35200461
rs71324962
rs60482815
rs6767355{;
	count if rsid=="`i'";
	replace proxy_rs148734725=1 if rsid=="`i'";
	};
	
#delimit ;
foreach i in rs9320913
rs368928614
rs71721600
rs201168005
rs138027849
rs12206087
rs2388334
rs9372734
rs12202969
rs968050
rs77910749
rs1487445
rs71806471
rs9375188
rs1487441
rs1906252
rs9401593{;
	count if rsid=="`i'";
	replace proxy_rs9320913=1 if rsid=="`i'";
	};
	
//Keep the potential proxies
keep if proxy_rs114598875==1| proxy_rs148734725==1| proxy_rs9320913==1;

//Check whether these SNPs are available in UKB:
joinby rsid using  "uk_biobank/ukb_snp_stats_chr06.dta",unmatched(master)
joinby rsid using  "uk_biobank/ukb_snp_stats_chr02.dta",unmatched(master) _merge(_merge2) update
joinby rsid using  "uk_biobank/ukb_snp_stats_chr03.dta",unmatched(master) _merge(_merge3) update

//Keep proxy with the highest INFO score
/*
Okbay SNP > Proxy
rs114598875 > rs17538393
rs148734725 > 3:49434654_G_A
rs9320913 > rs1487445

In each case the Major allele is the education increasing variant

*/

