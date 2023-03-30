#DriverGenePathway R package developed by xulu 20220213

#' Preprocess input data
#'
#' This function depends on a chromosome files folder, hg19 or hg38, which includes the txt files of each chromosome.
#' The hg19 folder is able to download at the url of our Github https://github.com/bioinformatics-xu/DriverGenePathway.
#' This function preprocesses the mutation data and coverage data and selects important mutation categories.
#' It generates intermediate results including preprocessed mutation file, coverage file, categ file and mutcateg_discovery file.
#'
#' @param M Mutation maf data.
#' @param C Coverage raw data.
#' @param V Covariate raw data.
#' @param dict Mapping Variant_Classification to mutation effect in mutation data.
#' @param chr_files_directory Chromosome files directory, hg19 or hg18.
#' @param categ_flag Mutation category number, should be either NaN or numeric.
#' @param output_filestem Unified prefix of output data, could be modified as well.
#' @param preprocessedOutput If output the preprocessed files, if preprocessedOutput is TRUE, then output preprocessed mutation, coverage, and covariate files; otherwise not.
#' @importFrom utils download.file unzip
#' @importFrom stats p.adjust
#' @return The output is a list includes the preprocessed mutation and coverage data, which are the inputs of BMR function.
preprocessing <- function(M,C,dict,V,
                          chr_files_directory,categ_flag=NaN,
                          output_filestem = "Preprocessed",preprocessedOutput=TRUE)
{is_coding = is_silent = categ = gene = effect = categ_idx = NULL

if(is.null(C)){
  # download and read coverage data
  cat("Download coverage file")
  download.file("http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/exome_full192.coverage.zip",
                destfile = "exome_full192.coverage.txt.zip")
  unzip("exome_full192.coverage.txt.zip", exdir = ".")
  file.remove("exome_full192.coverage.txt.zip")
  C <- as.data.frame(data.table::fread(file = "exome_full192.coverage.txt"))
}

if(is.null(dict)){
  cat("Download dictionary file")
  download.file("http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/mutation_type_dictionary_file.txt",
                destfile = "mutation_type_dictionary_file.txt")
  dict <- as.data.frame(data.table::fread(file = "mutation_type_dictionary_file.txt"))
}

if(is.null(V)){
  cat("Download covariate file")
  download.file("http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/gene.covariates.txt",
                destfile = "gene.covariates.txt")
  V <- as.data.frame(data.table::fread(file = "gene.covariates.txt"))
}

if(is.null(chr_files_directory)){
  cat("Download chromosome files (about 900M), which will take some time")
  download.file("http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutsig/reference_files/chr_files_hg19.zip",
                destfile = "chr_files_hg19.zip")
  unzip("chr_files_hg19.zip", exdir = ".")
  file.remove("chr_files_hg19.zip")
  chr_files_directory = "chr_files_hg19"
}

############### begin preprocess ###################
# make sure mutation data has gene+patient
#GENE
if (("gene" %in% colnames(M))&("Hugo_Symbol" %in% colnames(M))){
  # cat("NOTE: Both 'gene' and 'Hugo_Symbol' are present in mutation_file.
  #     Using 'gene'. \n")
}else if("Hugo_Symbol" %in% colnames(M)){
  M$gene <- M$Hugo_Symbol
}else if("gene" %in% colnames(M)){

}else{
  stop("mutation_file lacks 'gene' or 'Hugo_Symbol' column.")
}
#M <- M[order(M$gene),]

#PATIENT
if (("patient" %in% colnames(M))&("Tumor_Sample_Barcode" %in% colnames(M))){
  # cat("NOTE: Both 'patient' and 'Tumor_Sample_Barcode' are present
  #     in mutation_file. Using 'patient'. \n")
}else if("patient" %in% colnames(M)){

}else if("Tumor_Sample_Barcode" %in% colnames(M)){
  M$patient <- M$Tumor_Sample_Barcode
}else{
  stop("mutation_file lacks 'patient' or 'Tumor_Sample_Barcode' column.")
}

if (length(unique(M$patient)) < 2){
  stop("DriverGenePathway is not applicable to single patients. \n")
}

# ensure coverage data has gene+effect+categ
if (!("gene" %in% colnames(C))){
  stop("no 'gene' column in coverage_file")
}
if (!("effect" %in% colnames(C))&("zone" %in% colnames(C))){
  zone_col <- which(colnames(C) %in% "zone")
  colnames(C)[zone_col] <- "effect"
}
if (!("effect" %in% colnames(C))){
  stop("no 'effect' column in coverage_file")
}
C$effect <- gsub("^flank.*","noncoding",C$effect)
C_effectnames <- c("noncoding","silent","nonsilent")
if (!any(unique(C$effect) %in% C_effectnames)){
  stop("in coverage_file, 'effect' must be one of noncoding/silent/nonsilent")
}
if (!("categ" %in% colnames(C))){
  stop("no 'categ' column in coverage_file")
}
f <- colnames(C)
coverage_patient_names <- f[4:ncol(C)]
# if there are many patients, then coverage_patient_names is a list of patients


##############
#   EFFECT   #
##############
# take the effect of M, by dictionary.Variant_Classification
# and M.Variant_Classification

cat("Preprocessing mutation data ... \n")
if (("is_coding" %in% colnames(M)) || ("is_silent" %in% colnames(M))){
  # cat("NOTE: This version now ignores 'is_coding' and 'is_silent' \n")
  # cat("Requires Variant_Classification/type column and
  #     mutation_type_dictionary so we can assign nulls.\n")
  M <- subset(M,select = -is_coding)
  M <- subset(M,select = -is_silent)
  # cat("having run M subset")
}
if ("effect" %in% colnames(M)){
  # cat("Will use the pre-existing effect column.\n")
  M$effect <- gsub("^flank.*","noncoding",M$effect)
  M_effectnames <- c("noncoding","silent","nonsilent","null")
  if (!any(unique(M$effect) %in% M_effectnames)){
    stop("in mutation_file, 'effect' must be one of
           noncoding/silent/nonsilent/null")
  }
}else{
  if(!("Variant_Classification" %in% colnames(M))&("type" %in% colnames(M))){
    # cat("after if variant_classification in colnames ")
    M$Variant_Classification <- M$type
  }
  if(!("Variant_Classification" %in% colnames(M))){
    stop("mutation_file is missing Variant_Classification")
  }
  flag_num <- match(toupper(M$Variant_Classification),
                    toupper(dict$Variant_Classification),nomatch =nrow(dict)+1)
  dict <- rbind(dict,dict[nrow(dict),])
  dict[nrow(dict)+1,] <- "unknown"

  M$effect <- dict$effect[flag_num]
  bad <- which(M$effect == "unknown")
  if(length(bad)>0){
    cat(sprintf("WARNING: %d/%d mutations could not be mapped to
                    effect using mutation_type_dictionary_file: \n",
                length(bad),length(M$effect)))
    #table(bad_variant <- M$Variant_Classification[bad])
    cat("They will be removed from the analysis. \n")
    M <- M[-bad]
  }
  if(nrow(M) == 0){
    stop("No mutations left!")
  }
}

##############
#   CATEG   #
##############
ucc <- unique(C$categ)

# generate_192_categ_names
bases <- c("A","C","G","T")
names192 <- data.frame(names=0)
i=1
for (from in 1:4) {
  for (left in 1:4) {
    for (right in 1:4) {
      for (to in 1:4) {
        if(from==to){
          next
        }
        names192[i,1] <- paste(bases[left],"(",bases[from],"->",bases[to],")",
                               bases[right],sep = "")
        i <- i+1
      }
    }
  }
}

coverage_is_on_full192 <- (length(ucc) == 192 &
                             length(intersect(ucc,as.matrix(names192))) == 192)
categs_already_present <- FALSE

if("categ" %in% colnames(M)){
  mcc <- unique(M$categ)
  ucc <- unique(C$categ)
  if(length(ucc) != length(mcc) || sum(!(ucc %in% mcc)) > 0 ||
     sum(!(mcc %in% ucc)) > 0){
    # cat("categ of mutation_file does not match coverage_file. Ignoring it.\n")
    M <- subset(M,select = -categ)
  }else{
    categs_already_present <- TRUE
  }
}

# find out if we can do category discovery
can_do_category_discovery <- TRUE
if(!coverage_is_on_full192){
  can_do_category_discovery <- FALSE
  reason <- "coverage_file not on full192"
}
if(!("chr" %in% tolower(colnames(M))) & "chromosome" %in%
   tolower(colnames(M))){
  colnum <- which(tolower(colnames(M)) %in% "chromosome")
  M$chr <- M$Chromosome
}
if(!("start" %in% tolower(colnames(M))) & "start_position" %in%
   tolower(colnames(M))){
  colnum <- which(tolower(colnames(M)) %in% "start_position")
  M$start <- M[[colnum]]
}
if(!("ref_allele" %in% colnames(M)) & "reference_allele" %in%
   tolower(colnames(M))){
  colnum <- which(tolower(colnames(M)) %in% "reference_allele")
  M$ref_allele <- M[[colnum]]
}
if(!("newbase" %in% tolower(colnames(M)))){
  if("tumor_seq_allele1" %in% tolower(colnames(M))){
    colnum <- which(tolower(colnames(M)) %in% "tumor_seq_allele1")
    M$newbase <- M[[colnum]]
    if("tumor_seq_allele2" %in% tolower(colnames(M))){
      colnum_allele1 <- which(tolower(colnames(M)) %in% "tumor_seq_allele1")
      idx <- which(M$ref_allele == M[[colnum_allele1]])
      colnum_allele2 <- which(tolower(colnames(M)) %in% "tumor_seq_allele2")
      M$newbase[idx] <- M[[colnum_allele2]][idx]
    }
  }
}
if(!("chr" %in% colnames(M)) || !("start" %in% colnames(M)) ||
   !("ref_allele" %in% colnames(M)) || !("newbase" %in% colnames(M))){
  can_do_category_discovery <- FALSE
  reason <- "Chromosome/Start_position/Reference_Allele/Tumor_Seq_Allele1/
    Tumor_Seq_Allele2 missing from mutation data"
}else{
  M$start <- as.numeric(M$start)
  bad <- which(is.nan(M$start))
  if(length(bad) > 0){
    cat(sprintf("WARNING: %d/%d mutations had non-numeric Start_position.
                  Excluding them from analysis. \n",
                length(bad),length(M$start)))
    M <- M[-bad]
  }
  if(nrow(M) == 0){
    stop("No mutations left!\n")
  }
  if(!(length(list.files(chr_files_directory)))){
    can_do_category_discovery <- FALSE
    reason <- "no chr_files_directory available"
  }else{
    f1_flag <- grep("^chr.*txt$",list.files(chr_files_directory))
    f1 <- list.files(chr_files_directory)[f1_flag]
    uchr <- unique(M$chr)
    uchr <- uchr[order(uchr)]
    M$chr_idx <- match(M$chr,uchr)
    uchr <- gsub("chr","",uchr)
    f2 <- paste("chr",uchr,".txt",sep = "")
    chr_file_available <- f2 %in% f1
    if(sum(chr_file_available) == 0){
      can_do_category_discovery <- FALSE
      reason <- "no chr files available"
    }else{
      bad <- which(!chr_file_available[M$chr_idx])
      if(length(bad) > 0){
        cat(sprintf("WARNING: %d/%d mutations are on chromosomes not found in
                      chr_files_directory Excluding them from analysis. \n",
                    length(bad),nrow(M)))
        M <- M[-bad]
      }
      if(nrow(M) == 0){
        stop("No mutations left!\n")
      }
    }
  }
}

# DECIDE WHAT TO DO ABOUT CATEGORIES
# METHODS:     1. use the existing categories
#              2. have only one category for non-nulls
#              3. discover k categories for non-nulls
#categ_flag <- NaN
if(is.nan(categ_flag)){
  if(can_do_category_discovery){
    method <- 3
    ncategs <- 4
  }else if(categs_already_present){
    method <- 1
  }else{
    cat(sprintf("NOTE: unable to perform category discovery, because %s. ",
                reason))
    method<- 2
  }
}else if(categ_flag == 0){
  if(!categs_already_present){
    stop("when setting categ_flag==0, categ column must be already present in
           mutation data")
  }
  method <- 1
}else if(categ_flag==1){
  method <-2
}else if(categ_flag > 1){
  if(categ_flag>6){
    cat("NOTE: maximum categories that can be discovered is 6. \n")
    categ_flag <- 6
  }
  if(!can_do_category_discovery){
    stop(sprintf("unable to perform category discovery, because %s", reason))
  }else{
    method <- 3
    ncategs <- categ_flag
  }
  if(ncategs > 4){
    cat("NOTE: It may take dozens of hours to finish the process when the number of categories is set larger than 4.")
  }
}

# find categs by method
if(method==1){
  # cat(sprintf("Using the categories already presented. \n"))
  K <- data.frame(names = rep(1,each=length(unique(M$categ))))
  K$names <- sort(unique(M$categ))
}else if(method==2){
  # cat(sprintf("Will use two categories: missense and null+indel. \n"))

  K1 <- data.frame(left=0)
  K1$left <- "ACGT"
  K1$from <-"AC"
  K1$change <- "in"
  K1$right <- "ACGT"
  K1$autoname <- "missense"
  K1$name <- "missense"
  K1$type <- "point"

  K2 <- data.frame(left=0)
  K2$left <- "ACGT"
  K2$from <-"AC"
  K2$change <- "in"
  K2$right <- "ACGT"
  K2$autoname <- "null+indel"
  K2$name <- "null+indel"
  K2$type <- "non-point"

  K <- data.frame(left=0,from=0,change=0,right=0,autoname=0,name=0,type=0)
  K[1,] <- K1[1,]
  K[2,] <- K2[1,]

  # assign categories
  M$categ <- rep("---",times=nrow(M))
  flag_null <- which(M$effect == "null")
  M$categ[flag_null] <- rep(K$name[2],each=length(flag_null))
  flag_nonull <- which(M$effect != "null")
  M$categ[flag_nonull] <- K$name[1]

  # collapse coverage
  cat(sprintf("Preprocessing Coverage data... \n"))
  temp <- unique(C$categ)
  C$categ_idx <- match(C$categ,temp)

  C <- dplyr::arrange(C,gene,effect,categ_idx)
  #order as gene, effect, categ_idx, if decrease, then desc(gene)

  ug <- unique(C$gene)
  ng <- length(ug)
  ue <- unique(C$effect)
  ne <- length(ue)
  nk <- nrow(K)

  idx <- which(C$categ_idx <= nk)
  C2 <- C[idx,]
  C2$categ[which(C2$categ_idx == 1)] <- K$name[1]
  C2$categ[which(C2$categ_idx == 2)] <- K$name[2]
  C2 <- subset(C2,select = c(gene,effect,categ))
  np <- length(coverage_patient_names)

  for (p in 1:np) {
    oldcov <- array(as.integer(unlist(C[[coverage_patient_names[p]]])),
                    c(192,ne,ng))
    newcov <- rep(colSums(oldcov),each=2)
    C2[[coverage_patient_names[p]]] <- newcov
  }
  C <- C2
  remove(C2)
}else if(method == 3){
  cat(sprintf("Generating mutation categories...\n"))
  C_coverage <- C[,4:ncol(C)]
  if(length(dim(C_coverage)) == 0){
    C$totcov <- C_coverage
  }else{
    C$totcov <- apply(C_coverage,1,sum)
  }

  npm <- length(unique(M$patient))
  if((length(coverage_patient_names)==1) & (npm > 1)){
    C$totcov <- C$totcov * npm
  }

  # will use only the coding mutations+coverage to do this
  C$is_coding[which((C$effect == "nonsilent")|(C$effect == "silent"))] <- 1
  C$is_coding[which(!((C$effect == "nonsilent")|(C$effect == "silent")))] <- 0

  # collapse coverage to 192
  X <- data.frame(categ=rep(1,each=length(unique(C$categ))))
  X$categ <- unique(C$categ)
  C$categ_idx <- match(C$categ,unique(C$categ),nomatch = 0)
  X$left <- X$categ
  X$split <- strsplit(X$categ,split = "")
  X$left <- substr(X$categ,start=1,stop = 1)
  X$from <- substr(X$categ,start = 3,stop = 3)
  X$to <- substr(X$categ,start = 6,stop = 6)
  X$right <- substr(X$categ,start = 8,stop = 8)

  X$yname <- paste(X$from,"in",X$left,sep = " ")
  X$yname <- paste(X$yname,"_",X$right,sep = "")
  C_iscoding <- C[which(C$is_coding == 1),]
  X$N <- tapply(C_iscoding$totcov,C_iscoding$categ_idx,sum)
  X$newbase_idx <-  match(X$to,c("A","C","G","T"))

  Y_name <- preprocessGenerateCategContext65()
  Y <- data.frame(num=1:65,name=Y_name)
  X$context65 <- match(X$yname,Y_name,nomatch = 65)
  N <- tapply(X$N,X$context65,sum)
  Y_N <- matrix(data = 0,nrow = nrow(Y))
  Y_N[which(as.numeric(rownames(N)) %in% as.numeric(Y[,1])),1] <-  N
  Y$N <- Y_N
  Y$N <- round(Y$N/3)

  # STEPE 2
  #MUTATIONS：get context65 by looking up from reference genome
  cat(sprintf("Looking up trinucleotide contexts from chr files...\n"))

  #f2 <- paste("chr",uchr,".txt",sep = "")
  f2 <- paste(chr_files_directory,f2,sep = "/")
  triplet <- data.frame(name = rep(1,each=nrow(M)))
  for (ci in 1:length(uchr)){
    cat(sprintf("%d/%d ", ci,length(uchr)))
    midx <- which(M$chr_idx == ci)
    chrfile <- f2[ci]
    d <- file.info(chrfile)
    if(is.na(d$size)){
      next
    }
    filesize <- d$size
    ff <- data.table::fread(f2[ci],header = F)

    triplet$name[midx] <- stringr::str_sub(ff[1],M$start[midx]-1,M$start[midx]+1)
  }
  flag_non_triplet <- which(triplet == 1)
  if(length(flag_non_triplet)){
    triplet$name[flag_non_triplet] <- "---"
  }
  M$triplet <- toupper(triplet$name)
  M$triplet_middle <- stringr::str_sub(M$triplet,2, 2)
  midx <- which((M$ref_allele != "-") & (M$newbase != "-"))
  matchfrac <- sum(M$ref_allele[midx] == M$triplet_middle[midx])/length(midx)
  matchfrac
  if(matchfrac < 0.9){
    adj <- "possible"
    if(matchfrac < 0.7){
      adj <- "probable"
    }
    cat(sprintf("WARNING: %s build mismatch between mutation_file and chr_files", adj))
  }
  M$yname <- paste(substr(M$triplet,2,2),"in",substr(M$triplet,1,1),sep = " ")
  M$yname <- paste(M$yname,substr(M$triplet,3,3),sep = "_")
  M$context65 <- match(M$yname,Y$name,nomatch = 65)
  M$newbase_idx <-  match(substr(M$newbase,1,1),c("A","C","G","T"),
                          nomatch = NaN)

  midx <- which((M$ref_allele != "-") & (M$newbase != "-") &
                  (M$context65 >= 1) & (M$context65 <= 65) &
                  (M$newbase_idx >= 1) & (M$newbase_idx <= 4))
  M_Na <- subset(M,select = c("context65","newbase_idx"))[midx,]
  M_N <- plyr::count(M_Na, names(M_Na))
  freq <- M_N$freq
  bases <- c("A","C","G","T")
  Y$A <- 0
  Y$C <- 0
  Y$G <- 0
  Y$T <- 0

  for (i in 1:4) {
    M_N_i <- M_N[which(M_N$newbase_idx == i),]
    Y[[bases[i]]][M_N_i$context65] <- M_N_i$freq
  }

  #STEP3 Category Discovery
  #Category Discovery
  preprocessed_file_stem <- paste(output_filestem,"outputs",sep = "_")
  if(!dir.exists(preprocessed_file_stem)){
    dir.create(preprocessed_file_stem)
  }
  setwd(preprocessed_file_stem)

  Nn <- preprocessCollapseNn65to32(Y)
  P <- data.frame(max_k = ncategs,mutcategs_report_filename =
                    paste(output_filestem,"_mutcateg_discovery.txt",sep = ""))
  PP <- P
  Ks <- preprocessFindMutCateg(Nn,P)
  K <- Ks[[ncategs]]
  c <- preprocessAssignCateg(K)
  X$kidx <- matrix(data = NaN,nrow = nrow(X),ncol = 1)
  for (i in 1:nrow(X)) {
    X$kidx[i] <- which(c[X$context65[i],,X$newbase[i]]==1)
  }

  #STEP4
  #assign mutation categories
  cat(sprintf("Assigning mutation categories...\n"))
  M$categ <- rep("---",times=nrow(M))
  for (i in 1:nrow(X)) {
    idx <- which((M$context65==X$context65[i]) &
                   (M$newbase_idx==X$newbase_idx[i]))
    M$categ[idx] <- rep(as.character(K$name[X$kidx[i]]),each=length(idx))
  }

  # add null+indel category
  K2 <- as.data.frame(matrix(data = character(0),nrow = 1))
  K2$left <- "ACGT"
  K2$from <- "AC"
  K2$change <- "in"
  K2$right <- "ACGT"
  K2$autoname <-"null+indel"
  K2$name <- "null+indel"
  K2$type <- "non-point"
  K <- preprocessConcatStructsKeepFields(K,K2)
  K$N[nrow(K)] <- sum(K$N[1:(nrow(K)-1)])
  midx <- which(M$effect == "null")
  M$categ[midx] <- rep(K2$name[nrow(K2)],each=length(midx))
  K$n[nrow(K)] <- length(midx)
  K$rate[nrow(K)] <- K$n[nrow(K)]/K$N[nrow(K)]
  K$relrate[nrow(K)] <- K$rate[nrow(K)]/K$rate[1]*K$relrate[1]

  # STEP5
  # collapse coverage
  cat(sprintf("Collapsing coverages...\n"))
  C <- dplyr::arrange(C,gene,effect,categ_idx)
  #order as gene, effect, categ_idx, if decrease, then desc(gene)

  ug <- unique(C$gene)
  ng <- length(ug)
  ue <- unique(C$effect)
  ne <- length(ue)
  nk <- nrow(K)

  idx <- which(C$categ_idx <= nk)

  K_name <- unlist(unique(K$name))
  C2 <- C[idx,]
  # find categ of Coverage
  #C2_categ <- rep(NaN,times=nrow(C2))
  idx1 <- which((!is.nan(C2$categ_idx))&(C2$categ_idx>=1) &
                  (C2$categ_idx<=length(K_name)))
  C2$categ <- K_name[C2$categ_idx[idx1]]
  C2 <- subset(C2,select = c(gene,effect,categ))

  np <- length(coverage_patient_names)

  for (p in 1:np) {
    oldcov <- array(as.integer(unlist(C[[coverage_patient_names[p]]])),
                    c(192,ne,ng))
    newcov <- array(data = NaN,c(nk,ne,ng))
    for (ki in 1:nk) {
      if(ki==nk){
        cidx <- c(1:192)
      }else{
        cidx <- which(X$kidx==ki)
      }
      newcov[ki,,] <- apply(oldcov[cidx,,], c(2,3), sum)
    }

    C2[[coverage_patient_names[p]]] <- as.vector(newcov)
  }

  C <- C2
  remove(C2)

} # else if(method == 3)

#SAVE OUTPUT FILES
# (1) mutation file
rem_col_M <- c("newbase","chr_idx","triplet","yname","context65",
               "newbase_idx","context65_right","triplet_middle")
#col_M <- c("gene","patient","effect","chr","start","ref_allele","categ")
M <- subset(M,select = colnames(M)[which(!(colnames(M) %in% rem_col_M))])
V <- covariate_preprocessing(M,V)

if(preprocessedOutput){
  utils::write.table(M,file = paste(output_filestem,"_mutations.txt",sep = ""),
                     sep = "\t",quote = F,row.names = F)

  #(2) coverage file
  utils::write.table(C,file = paste(output_filestem,"_coverage.txt",sep = ""),
                     sep = "\t",quote = F,row.names = F)

  utils::write.table(V,file = paste(output_filestem,"_covariate.txt",sep = ""),
                     sep = "\t",quote = F,row.names = F)

  #(3) categories file
  utils::write.table(as.matrix(K),file = "MutationCategories.txt",
                     sep = "\t",quote = F,row.names = F)
}

setwd("..")
cat(sprintf("Preprocessing finished. \n"))

return(list(M=M,C=C,V=V))


} # of preprocess

