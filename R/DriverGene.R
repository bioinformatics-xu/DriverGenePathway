#' Identify driver genes
#'
#' This function depends on a chromosome files folder, hg19 or hg18, which includes the txt files of each chromosome.
#' This function first preprocesses input data via the preprocessing function.
#' Then the BMR function is used to calculate the background mutation rate. In order to identify driver genes, there are five hypothesis test methods to select
#' from according to the parameter p_class, i.e. binomial distribution test, represents beta binomial distribution test,
#' Fisher combined P-value test, likelihood ratio test, convolution test and 2D-projection method, which are all encapsulated in the sigGenes function.
#' @param Mutation Mutation maf data, mandatory data.
#' @param Coverage Coverage raw data, read from exome_full192.coverage.txt in MutSigCV by default.
#' @param Covariate Covariate data, read from "gene.covariates.txt" in MutSigCV by default.
#' @param MutationDict Mutation dictionary to map Variant_Classification to mutation effect in mutation data, read from "mutation_type_dictionary_file.txt" in MutSigCV by default.
#' @param chr_files_directory Chromosome files directory, hg19 or hg38.
#' @param categ_flag Mutation category number, should be either NaN or numeric, defaulted to NaN.
#' @param bmr The default background mutation rate is 1.2e-6, and the value alters when
#' function ends.
#' @param p_class Hypothesis test methods. "BB" represents beta binomial distribution test; "FCPT" represents Fisher combined P-value test;
#' "LRT" represents likelihood ratio test; "CT" represents convolution test; "projection" represents projection test
#' method; "allTest" represents the mutual results of all methods.
#' @param sigThreshold The threshhold of q-value to judge if the gene is significant.
#' @param output_filestem The parameters to name the output files, defaulted to "output".
#' @return Return of DriverGene function is driver gene table. Besides, output files contains txt files including the preprocessed mutation and coverage data, mutation categories, and significant driver genes. There are several plots output as
#' pdf files as well.
#' @author Xiaolu Xu <lu.xu@@lnnu.edu.cn>
#' @export
DriverGene <- function(Mutation=NULL,Coverage=NULL,Covariate=NULL,MutationDict=NULL,
                       chr_files_directory=NULL,categ_flag=NaN, bmr=1.2e-6,
                       p_class="BB", sigThreshold = 0.05,output_filestem="output")
{
  if(!is.null(Mutation)){
    preOut <- preprocessing(M=Mutation, C=Coverage, dict=MutationDict,V=Covariate,
                            chr_files_directory=chr_files_directory,categ_flag=categ_flag,output_filestem=output_filestem)
    plotCategory(preOut$M)
    plotEffect(preOut$M)

    BMR_out <- BMR(preOut$M, preOut$C, preOut$V, bmr)
    sigGenes <- sigGenes(BMR_out, p_class,output_filestem = output_filestem,sigThreshold)
    return(sigGenes)
  }else{
    cat("Error:invalid input mutation data")
  }
}
