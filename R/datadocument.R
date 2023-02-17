#'Sample mutation data
#'
#'A Maf data set containing raw information of mutations
#'
#'@usage data(SampleMutationMaf)
#'
#'@format A data frame with 816 rows and 12 variables:
#'\describe{
#'  \item{Hugo_Symbol}{The gene where the mutation happens}
#'  \item{Chromosome}{The Chromosome where the mutation happens}
#'  \item{Start_position}{The start position of the mutation}
#'  \item{End_position}{The end position of the mutation}
#'  \item{Variant_Classification}{Variant type of mutation}
#'  \item{Reference_Allele}{Reference allele}
#'  \item{Tumor_Seq_Allele1}{The first tumor allele}
#'  \item{Tumor_Seq_Allele2}{The second tumor allele}
#'  \item{Tumor_Sample_Barcode}{Code name of the patient}
#'  \item{is_coding}{If the mutation is in the coding area}
#'  \item{is_silent}{If the mutation is in the silent area}
#'  \item{categ}{The initial category of the mutation}
#'  }
"M"

#'Sample coverage data
#'
#'A data set containing raw information of coverages
#'
#'@usage data(SampleCoverage)
#'
#'@format A data frame with 31104 rows and 4 variables:
#'\describe{
#'  \item{gene}{The gene of mutation}
#'  \item{effect}{The effect of mutation}
#'  \item{categ}{The specific category}
#'  \item{coverage}{The number of coverage}
#'  }
"C"


#'Mutation type dictionary
#'
#'A data set containing variant mutation types and their correspondings for use
#'
#'@usage data(SampleDict)
#'
#'@format A data frame with 50 rows and 2 variables:
#'\describe{
#'  \item{Variant_Classification}{Variant types of effect}
#'  \item{effect}{Specified types for use}
#'  }
"dict"

#'Preprocessed mutation data
#'
#'The output mutation data of preprocessing
#'
#'@usage data(SampleOutMutation)
#'
#'@format A data frame with 816 rows and 16 variables:
#'\describe{
#'  \item{Hugo_Symbol}{The gene where the mutation happens}
#'  \item{Chromosome}{The Chromosome where the mutation happens}
#'  \item{Start_position}{The start position of the mutation}
#'  \item{End_position}{The end position of the mutation}
#'  \item{Variant_Classification}{Variant type of mutation}
#'  \item{Reference_Allele}{Reference allele}
#'  \item{Tumor_Seq_Allele1}{The first tumor allele}
#'  \item{Tumor_Seq_Allele2}{The second tumor allele}
#'  \item{Tumor_Sample_Barcode}{Code name of the patient}
#'  \item{gene}{The gene where the mutation happens}
#'  \item{patient}{Code name of the patient}
#'  \item{effect}{The effect of mutation}
#'  \item{chr}{The Chromosome where the mutation happens}
#'  \item{start}{The start position of the mutation}
#'  \item{ref_allele}{Reference allele}
#'  \item{categ}{The available category of the mutation}
#'  }
"preOutM"

#'Preprocessed coverage data
#'
#'The output coverage data of preprocessing
#'
#'@usage data(SampleOutCoverage)
#'
#'@format A data frame with 810 rows and 4 variables:
#'\describe{
#'  \item{gene}{The gene of mutation}
#'  \item{effect}{The effect of mutation}
#'  \item{categ}{The specific category}
#'  \item{coverage}{The number of coverage}
#'  }
"preOutC"

#'Covariate data
#'
#'The covariate data used for background mutation rate discovery
#'
#'@usage data(SampleCovariate)
#'
#'@format A data frame with 54 rows and 4 variables:
#'\describe{
#'  \item{gene}{Specific gene}
#'  \item{expr}{Expression level}
#'  \item{reptime}{Replication time}
#'  \item{hic}{Chromatin compartment}
#'  }
"V"

#'Output of BMR function
#'
#'A list containing result of BMR method
#'
#'@usage data(SampleBMRout)
#'
#'@format A list of 11 elements:
#'\describe{
#'  \item{G}{Specific gene}
#'  \item{N_silent}{The number of silent coverages}
#'  \item{n_silent}{The number of silent mutations}
#'  \item{N_nonsilent}{The number of nonsilent coverages}
#'  \item{n_nonsilent}{The number of nonsilent mutations}
#'  \item{N_noncoding}{The number of noncoding coverages}
#'  \item{n_noncoding}{The number of noncoding mutations}
#'  \item{X_gcp}{Used in gene detecting methods}
#'  \item{x_gcp}{Used in gene detecting methods}
#'  \item{null_categ}{The number of categories}
#'  \item{nanalysis}{Number of rows}
#'  }
"BMR_out"

#'Mutation matrix
#'
#'A matrix of mutations occur in which patients and genes.
#'The values are either 1 or 0, representing mutated or not.
#'
#'@usage data(SampleMutationMatrix)
#'
"mutation_matrix"

