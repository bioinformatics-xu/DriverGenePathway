% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datadocument.R
\docType{data}
\name{preOutM}
\alias{preOutM}
\title{Preprocessed mutation data}
\format{
A data frame with 816 rows and 16 variables:
\describe{
 \item{Hugo_Symbol}{The gene where the mutation happens}
 \item{Chromosome}{The Chromosome where the mutation happens}
 \item{Start_position}{The start position of the mutation}
 \item{End_position}{The end position of the mutation}
 \item{Variant_Classification}{Variant type of mutation}
 \item{Reference_Allele}{Reference allele}
 \item{Tumor_Seq_Allele1}{The first tumor allele}
 \item{Tumor_Seq_Allele2}{The second tumor allele}
 \item{Tumor_Sample_Barcode}{Code name of the patient}
 \item{gene}{The gene where the mutation happens}
 \item{patient}{Code name of the patient}
 \item{effect}{The effect of mutation}
 \item{chr}{The Chromosome where the mutation happens}
 \item{start}{The start position of the mutation}
 \item{ref_allele}{Reference allele}
 \item{categ}{The available category of the mutation}
 }
}
\usage{
data(SampleOutMutation)
}
\description{
The output mutation data of preprocessing
}
\keyword{datasets}
