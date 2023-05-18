# DriverGenePathway
DriverGenePathway is an R package developed to identify driver genes and pathways in cancer using MutSigCV and statistical methods. To provide a user-friendly interface for those who may not be familiar with R programming, we have also created the DriverGenePathway web server ([http://www.hello-ai.cloud:3838/DriverGenePathwayAPP/](http://www.hello-ai.cloud:3838/DriverGenePathwayAPP/)). The web server allows users to upload their data using a simple interface and perform analyses using the same algorithms as those used in the R package.

## Access
DriverGenePathway is free for non-commerical use only.

## Installation

### Github Installation
```
library(devtools)
install_github("bioinformatics-xu/DriverGenePathway")
```

## Computational procedure
DriverGene and DriverPathway are two main functions in DriverGenePathway to identify driver genes and driver pathways.

### DriverGene

#### Run DriverGene with user-input preprocessing data:

Users can download and input the necessary coverage, covariates, mutation dictionary and reference genome data applicable to their mutation data for preprocessing. For instance, users can download these files from MutSigCV and ensure that they are stored in the current working directory of R as shown below:

```
Folder Structure of current working directory of R
|__ brca_maf.txt (example mutation data)
|__ exome_full192.coverage.txt
|__ gene.covariates.txt
|__ mutation_type_dictionary_file.txt
|__ chr_files_hg19
```

Identify driver genes using DriverGene() function

```
BRCAmutation <- as.data.frame(fread("brca_maf.txt"))
C <- as.data.frame(data.table::fread(file = "exome_full192.coverage.txt"))
dict <- as.data.frame(data.table::fread(file = "mutation_type_dictionary_file.txt"))
V <- as.data.frame(data.table::fread(file = "gene.covariates.txt"))
DriverGenes(Mutation=BRCAmutation,Coverage=C,Covariate=V,MutationDict=dict, chr_files_directory="chr_files_hg19",categ_flag=NaN, bmr=1.2e-6, p_class="binomial", sigThreshold = 0.05)
```

#### Run DriverGene with default preprocessing data:


If the default files for preprocessing are NULL, they will be automatically downloaded. However, note that this process may take some time. Once the necessary preprocessing data is ready, identifying driver genes can be easily conducted using the DriverGene() function as shown below:

```
BRCAmutation <- fread("brca_maf.txt")
DriverGenes(Mutation=BRCAmutation)
```                              

#### Description of parameters:

|Parameters|Description|
|----|----|
|Mutation|Mutation maf data, mandatory data.|
|Coverage|Coverage raw data which can be input by users. By default, the function will automatically download exome_full192.coverage.txt from MutSigCV as Coverage.|
|Covariate|Covariate data which can be input by users. By default, the function will automatically download "gene.covariates.txt" from MutSigCV as Covariate.|
|MutationDict|Mutation dictionary to map Variant_Classification to mutation effect in mutation data. Be default, the function will automatically download "mutation_type_dictionary_file.txt" from MutSigCV as MutationDict.|
|chr_files_directory|Chromosome files directory, hg19 or hg18, which can be input by users. By default, the function will automatically download "chr_files_hg19.zip" from MutSigCV as chr_files_directory.|
|categ_flag|Mutation category number, should be either NaN or numeric, defaulted to NaN.(category number is set to 4, and Method3 is adopted).|
|bmr|The default background mutation rate is 1.2e-6, and the value alters when function ends.|
|output_filestem|The parameters to name the output files, defaulted to "output".|
|p_class|Hypothesis test methods. "BB" represents beta binomial distribution test; "FCPT" represents Fisher combined P-value test; "LRT" represents likelihood ratio test; "CT" represents convolution test; "projection" represents projection test method; "allTest" represents the mutual results of all methods.|
|sigThreshold|The threshhold of q-value to judge if the gene is significant.|


#### Description of output:
|Output file|Description|
|----|----|
|output_filestem$+"_covaraite.txt"(output_covaraite.txt by default)| Filled covaraite.|
|output_filestem$+"_mutations.txt"(output_mutations.txt by default)| Preprocessed mutation data corresponding to the selected categories.|
|output_filestem$+"_coverage"(output_coverage.txt by default)| Preprocessed coverage data corresponding to the selected categories|
|output_filestem$+"_mutcateg_discovery.txt"(output_mutcateg_dicovery.txt by default)| Mutation categories for k=1,2,...,ncat.|

### DriverPathway

#### Default data preparation:

Mutation matrix or MAF.

```
Folder Structure of your working directory of Rstudio
|__ mutation_matrix.txt (example mutation matrix data)
```

#### Running command:
```
mutation_matrix <- fread("mutation_matrix.txt")
DriverPathway(mutation_matrix,
                       driver_size=3,
                       pop_size=200,
                       iters=500,
                       permut_time=1000,outfile="denovoDriverPathway.txt")
```
#### Description of parameters:

|Parameters|Description|
|----|----|
|mutation_matrix|0-1 mutation matrix where rows represent patients, columns represent genes or MAF file. If the mutation data is a MAF file, then the preprocessing procedure will be performed using coverage data, covariate data, mutation dictionary, and chr files directory, which can be specified by the user. The user should name these files "coverage.txt.$" for coverage data, "covariates.txt$" for covariate data, "dictionary file.txt$" for mutation dictionary, and "^chr files hg" for "chr files directory". If the specified named files can be found, they will be read and utilized to generate mutation matrix. If the files cannot be detected, the function will automatically download them from the MutSigCV website, which may take some time.|
|driver_size|The size of identified driver gene set, defaulted to 3.|
|pop_size|The population size of GA, defaulted to 200.|
|iters|The iteration time of GA, defaulted to 500.|
|permut_time|The times of permutation test, defaulted to 1000.|
|preprocess_bmr|The background mutation rate to preprocess MAF file to mutation matrix using binomial hypothesis testing.|

#### Description of output:
|Output|Description|
|----|----|
|driverlist|A list including the mutation matrix, identified driver gene set, and p-value.|


### Developer: 

Xiaolu Xu 

lu.xu@lnnu.edu.cn

Liaoning Normal University.
