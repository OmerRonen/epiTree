# epiTree


Requirements: Plink, Python (as in PrediXcan), R, all command line comments are assumed to be executed on a Mac or Unix machine


1. Get list of SNPs for PrediXcan weights

File: can be done with “utilities_snpFromDB.R”

2. Generate Plink files with subset of SNPs from PrediXcan database

File: submit_plinkSubset.sh

Notes: fam file can just be fam file from plink when you don’t want to remove any subjects, or it could be a subset of subjects. This was done on Sherlock, one would need to modify potentially on other servers (ml biology, plink, etc.). Plink files are expected to be one per chromosome, remember that this should include as many SNPs as available including imputed SNPs, as PrediXcan will only consider a subset of those. You have to run this script for each chromosome separately, so on command line you would do the following: for…

3. Create dosage files out of Plink files

File: submit_plink2dos.sh

Note: Again, you have to run this for each chromosome separately, it is a wrapper of the script from “convert_plink_to_dosage.py” from (https://github.com/hakyimlab/PrediXcan/blob/master/Software/convert_plink_to_dosage.py) see comment here (https://github.com/hakyimlab/PrediXcan/tree/master/Software)


4. Run PrediXcan on dosage files to get imputed gene expression

File: submit_dos2predix.sh

Note: Make sure that there is a fam file in the dosage folder! This is a wrapper from PrediXcan.py from PrediXcan. This is the verision we used for the analysis, you might download latest version directly from.

5. Run iRF on PrediXcan

File: analysis_iRF_gene.R

Note: One has to specify several path: PrediXcan file, phenotype file, name of data field in phenotype file, encoding of cases and NAs in the data field, mapping file between phone (1st column) and geno (2nd column) IDs, for different phenotypes one has to adapt the load_pheno.R file accordingly


6. Compute p-values

File analysis_pcsPvalues_gene.R

Notes: Same paths as above need to be specified. The four tests functions are implemented in utitilities_test.R can can be used independently. Mutli is a wrapper using glm R package together with LR function. You need to specify paths.


7. Extract SNP coordinates for genes of step 4. and then respective SNPs from Plink 

File 1. utilities_extractSNPsfromGenes.R 
       2. submit_plinkSubset.sh

Note: Plink needs to be loaded for this, as this includes wrapper for plink. First step generated list of SNPs from Plink files which should be included in the analysis. Second script generates subset plink files with these SNPs (as explained above). This should be the input plink file for the next stip.

8. Run iRF on SNPs

File: analysis_iRF_snp.R

Note: Specify correct plink file location and phenotype location 

9. Compute p-values for SNPs

File analysis_pcsPvalues_snp.R

Note: you need to specify paths, PCA
