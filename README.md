# epiTree -- Learning epistatic polygenic phenotypes with boolean interactions

## General Information

In this GitHub repository we provide code and a step-by-step guideline for the ***epiTree*** pipeline, to extract epistatic higher-order interactions from genotype - phenotype data.

A detailed explanaition of the pipeline can be found at the following manuscript

Behr M, Kumbier K, Cordova-Palomera A, Matthew Aguirre M, Ashley E, Butte E, Arnaout R, Brown B, Priest J†, Yu B† (2020)
**Learning epistatic polygenic phenotypes with boolean interactions**


## Software requirements

To run the pipeline, as it is presented in the following, the following software requirements are needed

- Linux or Mac OS
- R (version 3.6.0 or higher)
- Python 2.7
  - numpy package
- Plink 1.9
- slurm to run the .sh scripts


## Data input

The ***epiTree*** pipeline requires the following data input:

1. Genotype files in Plink format (bim/fam/bam), one file per chromosome named as *name_plink_file_chr1*, *name_plink_file_chr2*, etc.. It is ok when some chromosomes are missing.
2. A corresponding phenotype file, where a single, binary phenotype is considered. Real valued phenotypes are also possible, with slight modifications to the scripts.
3. A mapping file that maps the subject identification numbers of the phenotype file (first column) to the subject identification numbers of the genotype file (second column).
5. A file with PCAs for individual subjects, which is assumed to have the same subject IDs as the genotype file.
6. A PrediXcan [1] data base file, which is used to impute tissue specific gene expression data and can be downloaded at http://predictdb.org.



## Step-by-step guidline to run the epiTree pipeline

### Biologically inspired dimension reduction via PrediXcan

The first step of the ***epiTree*** pipeline is to perform a biologically inspired dimension reduction step via imputing tissue specific gene expression levels from SNP data, which can be done with the PrediXcan software [1]. Depending on the PrediXcan database that is used, this reduces the feature dimension from several million SNPs to a few thousand gene level feautures. 

A detailed description for how to obtain imputed gene expression levels, is provided by the authors of PrediXcan, see https://github.com/hakyimlab/PrediXcan and manuscript [1]. For sake of completeness, we provide a step-by-step description below.

#### Extract SNPs from the PrediXcan database
To speed up computation time and avoid memory issues, in a first step we will create new Plink files that only contain the SNPs that actually enter in the particular PrediXcan model. To this end, we first obtain a list of all SNPs that enter the particular PrediXcan data base. An R script, which does this is provided in 

`scripts/utilities_snpFromDB.R`

Note that you might need to adjust `name.db` and `path.db` to provide the correct name and path for your database file.

#### Generate Plink files which only contain PrediXcan SNPs

Given the list of SNPs from the previous step, we now generate new Plink files that only contain those SNPs. The following shell script submitts an slurm job which does this

`scripts/submit_plinkSubset.sh`

Note that you might need to modify the slurm arguments in the script accordingly, as well as module load (ml) commands.
Further, you will need to specify 

1. path and name of your input and output Plink files, 
2. path and name of the SNP list that you generated in the previous step, 
3. path and name of a .fam file with the subject IDs that you want to consider.

When you only want to consider a subset of subjects from your original Plink file, than this can be specified in the .fam file. Otherwise, you can take the .fam file from your input Plink files. Recall that your input Plink files are assumed to be one (bim/bam/fam) file per chromosome, named as *name_plink_file_chr1*, *name_plink_file_chr2*, etc..

You will need to run this script for each chromosome separately, so for example, when you consider chromosome 1 - 22, you can run

`for i in {1..22}; do ./submit_plinkSubset.sh $i; done`


#### Generate dosage files

PrediXcan requires genotype input files as dosage files. The following script submits slurm jobs which transform the Plink files from the previous step into dosage files, using the python script “convert_plink_to_dosage.py” from https://github.com/hakyimlab/PrediXcan/blob/master/Software/convert_plink_to_dosage.py.

`scripts/submit_plink2dos.sh`

Note that you will need to specify 

1. path and name of the plink files from the previous step,
2. path and name of the output dosage files.


#### Run PrediXcan on dosage files to get imputed gene expression

Having the dosage files from the previous step, we can finally run PrediXcan. The following script submits a slurm job for this. 

`scripts/submit_dos2predix.sh`

Note that you will need to run the script only once and not for each chromosome separately (it is assumed that there is one dosage file per chromosome, as was generated by the previous script). You will need to specify 

1. path and name of the dosage files, 
2. path and name of the PrediXcan data base as in previous steps,
3. name of a .fam file, which is assumed to be in the same folder as the dosage files. 

The `PrediXcan.py` file that we provide in the `scripts` folder here, is the one that we used for our analysis on the red-hair phenotpye, see [2], you may download the latest version at https://github.com/hakyimlab/MetaXcan (see also https://github.com/hakyimlab/PrediXcan).


### Gene level analysis

Once we obtained the imputed gene expression data from the previous step, we can run the gene level analysis of the ***epiTree***  pipeline. To this end, we need to do a data split into training and test data. As an example, in the following scripts we select a balanced sample of 26K training and 4K test samples, as in [2]. For different sample size, the scripts can be easily adapted.

#### Run iRF candidate interaction selection for imputed gene expression features 

As a first step, using the *training data* only, we extract candidate interactions on the gene level using the iRF pipeline [3], [4]. The following R script uses the *iRF* R package to do this. It also runs two competing prediction methods, namely penalized logistic regression with L1 penality using the *glm* R package and random forest as in the *ranger* R package implementation with default parameters.

`scripts/analysis_iRF_gene.R`

Note that, as before, you will need specify

1. path and name for imputed gene expression files from the previous steps as *path.predx*, 
2. path and name for the PrediXcan data basefile *path.db*, 
3. path and name for the phenotype files *path.pheno*, 
4. path and name for the mapping files between genotype and phenotype subject IDs *path.key*.

Note that there were a couple of genes which, by default, were excluded by the analysis and are specified in the file `/data/geneExclusionList.txt` and passed over to the script via the *path.excl* argument. However, none of those genes was present in the PrediXcan database `/data/gtex_v7_Skin_Sun_Exposed_Lower_leg_imputed_europeans_tw_0.5_signif.db` that was used for our red hair analysis. In practice, one can modify this file accordingly.

#### Compute CART based PCS p-values on gene level

For the candidate interactions from the previous step, we can now compute CART based PCS p-values, using the *test data* to evaluate their significane. Note that the this step also requires the *training data* from the previous step, as it is used to fit the individual CART components. This is done in the following R script, which also computes standard p-values from logistic regression with a multiplicative interaction term for comparison, using the *glm* and *lmtest* R packages. As before, you will need to adapt the paths accordingly.

`scripts/analysis_pcsPvalues_gene.R`

We stress that PCS p-values can also be computed independently from the iRF model selection step of the previous section. An implementation of CART based PCS p-values can be found in the script `scripts/utilities_tests.R`.

### SNP level analysis

From the candidate interactions between imputed gene expression features, we can go back to the SNP level, to seach for interactions between SNPs that correspond to genes that were extracted from the previous step.

#### Extract SNP coordinates and subset Plink files

The following R script uses the results from the previous step to extract for each gene that appears in an iRF interaction the respective coordinates (start/end location +/- 1K base pairs). It also includes the top 50 genes in terms of Gini importance. Then it extracts a list of all SNPs within those regions from the input Plink files. As before, one has to specify the path of the plink files and PrediXcan database accordingly.

`scripts/utilities_extractSNPsfromGenes.R`

Given this list of SNPs, one can use again the script `submit_plinkSubset.sh` to generate Plink files that only contain those SNPs.


#### Run iRF on SNP level 

Given the Plink files from the previous step, with reduced number of SNPs, we can now again use iRF to search for stable interactions on the SNP level, as done by the following script (with paths specified accordingly).

`scripts/analysis_iRF_snp.R`


#### Compute PCS p-values on SNP level

For the candidate interactions from the previous step, we can now compute PCS p-values on the SNP level, as in the following script (with paths specified accordingly).


`scripts/analysis_pcsPvalues_snp.R`



### Evaluate results

All results and output files are stored in the `results` folder. 

We have added the results that we obtained in our red hair analysis [2] to this github repository and demonstrate some visualization and reproduction of figures from [2] in the R Jupyter Notebook *red_hair_analysis.ipynb*.



## References

[1] Gamazon ER†, Wheeler HE†, Shah KP†, Mozaffari SV, Aquino-Michaels K,
Carroll RJ, Eyler AE, Denny JC, Nicolae DL, Cox NJ, Im HK. (2015)
**A gene-based association method for mapping traits using reference
transcriptome data**. Nat Genet. doi:10.1038/ng.3367.
([Link to paper](http://www.nature.com/ng/journal/v47/n9/full/ng.3367.html),
[Link to Preprint on BioRxiv](http://biorxiv.org/content/early/2015/06/17/020164))

[2] Behr M, Kumbier K, Cordova-Palomera A, Matthew Aguirre M, Ashley E, Butte E, Arnaout R, Brown B, Priest J†, Yu B† (2020)
**Learning epistatic polygenic phenotypes with boolean interactions**

[3] Basu S, Kumbier K, Brown B, Yu B (2018)
**Iterative random forests to discover predictive and stable high-order interactions**. PNAS. https://doi.org/10.1073/pnas.1711236115

[4] Kumbier K, Basu S, Brown B, Celniker S, Yu B (2018).
***Refining interaction search through signed iterative Random Forests***
https://arxiv.org/abs/1810.07287 

