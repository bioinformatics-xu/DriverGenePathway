#' Identify driver genes based on MutSigCV and statistical methods
#'
#' This function depends on a chromosome files folder, hg19 or hg18, which includes the txt files of each chromosome.
#' This function first preprocesses input data via the preprocessing function.
#' Then the BMR function is used to calculate the background mutation rate. In order to identify driver genes, there are five hypothesis test methods to select
#' from according to the parameter p_class, i.e. binomial distribution test, represents beta binomial distribution test,
#' Fisher combined P-value test, likelihood ratio test, convolution test and 2D-projection method, which are all encapsulated in the sigGenes function.
#' @param Mutation Mutation maf data, mandatory data.
#' @param Coverage Coverage raw data which can be input by users. By default, the function will automatically download exome_full192.coverage.txt from MutSigCV as Coverage.
#' @param Covariate Covariate data which can be input by users. By default, the function will automatically download "gene.covariates.txt" from MutSigCV as Covariate.
#' @param MutationDict Mutation dictionary to map Variant_Classification to mutation effect in mutation data. Be default, the function will automatically download "mutation_type_dictionary_file.txt" from MutSigCV as MutationDict.
#' @param chr_files_directory Chromosome files directory, hg19 or hg18, which can be input by users. By default, the function will automatically download "chr_files_hg19.zip" from MutSigCV as chr_files_directory.
#' @param categ_flag Mutation category number, should be either NaN or numeric, defaulted to NaN.
#' @param bmr The default background mutation rate is 1.2e-6, and the value alters when function ends.
#' @param p_class Hypothesis test methods. "BB" represents beta binomial distribution test; "FCPT" represents Fisher combined P-value test;
#' "LRT" represents likelihood ratio test; "CT" represents convolution test; "PJ" represents projection test
#' method; "allTest" represents the mutual results of all methods.
#' @param sigThreshold The threshhold of q-value to judge if the gene is significant.
#' @param output_filestem The parameters to name the output files, defaulted to "output".
#' @return Driver gene table.
#' @author Xiaolu Xu <lu.xu@@lnnu.edu.cn>
#' @export
DriverGene <- function(Mutation=NULL,Coverage=NULL,Covariate=NULL,MutationDict=NULL,
                       chr_files_directory=NULL,categ_flag=NaN, bmr=1.2e-6,output_filestem="output",
                       p_class="allTest", sigThreshold = 0.05)
{
  if(!is.null(Mutation)){
    preOut <- preprocessing(M=Mutation, C=Coverage, dict=MutationDict,V=Covariate,
                            chr_files_directory=chr_files_directory,categ_flag=categ_flag, output_filestem=output_filestem)

    BMR_out <- BMR(preOut$M, preOut$C, preOut$V, bmr)
    sigGenes <- sigGenes(BMR_out, p_class,sigThreshold, output_filestem=output_filestem)
    return(sigGenes)
  }else{
    stop("Invalid input mutation data")
  }
}
