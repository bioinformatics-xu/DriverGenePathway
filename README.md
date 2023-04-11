# DriverGenePathway
DriverGenePathway is to identify driver genes and driver pathways in cancer based on MutSigCV and statistical methods.

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
|Mutation|Mutation maf data (mandatory).|
|Coverage|Coverage raw data, read from exome_full192.coverage.txt in MutSigCV by default.|
|Covariate|Covariate data, read from "gene.covariates.txt" in MutSigCV by default.|
|MutationDict|Mutation dictionary to map Variant_Classification to mutation effect in mutation data, read from "mutation_type_dictionary_file.txt" in MutSigCV by default.|
|chr_files_directory|Chromosome files directory, hg19 or hg38.|
|categ_flag|Mutation category number, should be either NaN or numeric, default to be NaN(category number is set to 4, and Method3 is adopted).|
|bmr|The default background mutation rate, 1.2e-6 by default.|
|p_class|Hypothesis test methods. "binomial" represents binomial distribution test, which is default;"betabinomial" represents beta binomial distribution test;"fisher" represents Fisher combined P-value test;"LRT" represents likelihood ratio test; "CT" represents convolution test; "projection" represents 2D-projection test; "allTest" represents the mutual results of all test methods.|
|sigThreshold|The significance level for q-value to idenfity driver genes.|
|output_filestem|The parameters to name the output files, defaulted to "output"|

#### Description of output files:
|Output file|Description|
|----|----|
|MutationCategories.txt| Information about selected mutation categories.|
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
|mutation_matrix|0-1 mutation matrix where rows represent patients, columns represent genes.|
|driver_size|The size of identified driver gene set, defaulted to 3.|
|pop_size|The population size of GA, defaulted to 200.|
|iters|The iteration time of GA, defaulted to 500.|
|permut_time|The times of permutation test, defaulted to 1000.|
|outfile|The output file of driver gene set, defaulted to "denovoDriverPathway.txt".|

#### Description of output files:
|Output file|Description|
|----|----|
|denovoDriverPathway.txt|The identified driver gene set and statistical significance of permutation test.|

### The Default files are available in the following URL.

.......

### Developer: 

Xiaolu Xu 

lu.xu@lnnu.edu.cn

Liaoning Normal University.
