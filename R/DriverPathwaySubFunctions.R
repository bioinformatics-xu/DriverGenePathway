select_operator <- function(evalVals_vector){
  len_eval <- length(evalVals_vector)
  order_eval <- order(evalVals_vector,decreasing = T)
  p <- rep(0,len_eval)
  p[order_eval] <- 2*c(1:len_eval)/(len_eval*(len_eval+1))
  cump <- cumsum(p)
  randD <- stats::runif(1,0,1)
  select1 <- min(which(cump>randD))
  randD <- stats::runif(1,0,1)
  select2 <- min(which(cump>randD))
  select_index <- c(select1,select2)
  return(select_index)
}

cross_operator <- function(parent1,parent2,genevars){
  total_num <- sum(parent1)
  new_pop <- rep(0,genevars)
  index_optim <- which(parent1+parent2 == 2)
  new_pop[index_optim] <- 1
  parent1[index_optim] <- 0
  parent2[index_optim] <- 0
  index_or <- which(parent1 + parent2 == 1)
  index_rand <- sample(which(parent1+parent2==1),total_num-sum(new_pop))
  new_pop[index_rand] <- 1
  return(new_pop)
}

mutation_operator_main <- function(pop_x,nx){
  nd <- sum(pop_x)
  index_x0 <- sample(c(1:nx),1)
  while (pop_x[index_x0]){
    index_x0 <- sample(c(1:nx),1)
  }
  pop_x_nonzero <- which(pop_x == 1)
  index_x1 <- sample(pop_x_nonzero,1)
  mut_pop_x <- pop_x
  mut_pop_x[index_x0] <- 1
  mut_pop_x[index_x1] <- 0
  return(mut_pop_x)
}

mutation_operator <- function(mutation_matrix,population_i,mutn,mutN){
  for (mi in 1:mutN) {
    population_j <- mutation_operator_main(population_i,mutn)
    if(evalFunc(mutation_matrix,population_j) < evalFunc(mutation_matrix,population_i)){
      population_i <- population_j
    }
  }
  return(population_i)
}

evalFunc <- function(mutation_matrix,chromosome){
  returnVal = 0
  matrix <- mutation_matrix[,which(chromosome == 1)]
  sumcol <- apply(matrix,2,sum)
  sumrow <- apply(matrix,1,sum)
  sumg <- sumcol %*% (exp(-sumcol/nrow(mutation_matrix))/sum(exp(-sumcol/nrow(mutation_matrix))))
  #print(2*sum(Ai.sumrow>0))
  #print(Ai.sumg)
  returnVal <- as.numeric(-(2*sum(sumrow>0)-sumg))
  return(returnVal)
}

mulExclusive_significance_subSet <- function(mutation_matrix, driver_geneset, permut_time = 1000){
  m <- nrow(mutation_matrix)
  weight_score <- rep(0,permut_time)
  n <- length(driver_geneset)
  flag <- which(colnames(mutation_matrix) %in% driver_geneset)
  chromosome_data <- rep(0,ncol(mutation_matrix))
  chromosome_data[flag] <- 1

  for (j in 1:permut_time) {
    mutMatrix_temp <- mutation_matrix
    mutMatrix_temp[,flag] <- 0
    for (i in 1:n) {
      temp <- sum(mutation_matrix[,flag[i]])
      index <- sample(1:m,temp,replace = F)
      mutMatrix_temp[index,flag[i]] <- 1
    }
    weight_score[j] <- evalFunc(mutMatrix_temp,chromosome_data)
  }

  p_value <- sum(evalFunc(mutation_matrix,chromosome_data) >= weight_score)/permut_time
  return(p_value)
}

mulExclusive_significance <- function(mutation_matrix,driver_geneset,permut_time=1000){
  if(any(class(driver_geneset) == "matrix")){
    p_value_list <- c()
    for (gs in 1:nrow(driver_geneset)) {
      driver_geneset_sub <- driver_geneset[gs,]
      p_value <- mulExclusive_significance_subSet(mutation_matrix, driver_geneset_sub, permut_time)
      p_value_list <- c(p_value_list, p_value)
    }
    return(p_value_list)
  }else if(class(driver_geneset) == "character"){
    p_value <- mulExclusive_significance_subSet(mutation_matrix, driver_geneset, permut_time)
    return(p_value)
  }
}

# generate input mutation matrix for AWEMP
generate_mutation_matrix <- function(M,gene.name){
  patient.mut <- base::unique(subset(M,select = c("patient","gene")))

  gene.flag <- which(patient.mut$gene %in% gene.name)
  pat.mut.col <- as.data.frame(patient.mut[gene.flag,])
  gene <- sort(unique(pat.mut.col$gene))
  patient <- sort(unique(pat.mut.col$patient))
  pat.mut.matrix <- as.data.frame(matrix(data = NA,nrow = length(patient),ncol = length(gene)))
  dimnames(pat.mut.matrix) <- list(patient,gene)

  for(i in 1:(length(patient))){
    maf.patient <- pat.mut.col[which(pat.mut.col$patient == patient[i]),]
    flag.g <- which(gene %in% maf.patient$gene)
    pat.mut.matrix[i,flag.g] <- 1
  }

  pat.mut.matrix[is.na(pat.mut.matrix)] <- 0
  pat.mut.matrix <- pat.mut.matrix[,!(colnames(pat.mut.matrix) %in% "TTN")]
  return(pat.mut.matrix)
}

preprocessing_mutation_data <- function(mutation_data,bmr=1.2e-6){
  gene_idx = NULL
  if(length(which(mutation_data == 0)) + length(which(mutation_data == 1)) == nrow(mutation_data)*ncol(mutation_data)){
    # the mutation data contains only 0 and 1, then identify driver pathway directly
    mutation_matrix <- mutation_data
    return(mutation_matrix)
  }else if(all(c("Hugo_Symbol", "Tumor_Sample_Barcode","Variant_Classification") %in% colnames(mutation_data))){
    # The mutation data is MAF file, then preprocess it to 0/1 mutation matrix
    message("The input mutation data is MAF file, preprocess it to 0/1 mutation matrix")
    current.files <- list.files()
    if(any(grepl("coverage.txt",current.files))){
      C <- as.data.frame(data.table::fread(file = current.files[grepl("coverage.txt",current.files)]))
    }else{
      C <- NULL
    }

    if(any(grepl("dictionary_file.txt",current.files))){
      dict <- as.data.frame(data.table::fread(file = current.files[grepl("dictionary_file.txt",current.files)]))
    }else{
      dict <- NULL
    }

    if(any(grepl("covariates.txt",current.files))){
      V <- as.data.frame(data.table::fread(file = current.files[grepl("covariates.txt",current.files)]))
    }else{
      V <- NULL
    }

    if(any(grepl("chr_files_hg",current.files))){
      chr_files_directory <- current.files[grepl("chr_files_hg",current.files)]
    }else{
      chr_files_directory <- NULL
    }

    preOut <- preprocessing(M=mutation_data, C=C, dict=dict,V=V,chr_files_directory=chr_files_directory,preprocessedOutput=FALSE)
    M <- preOut$M
    C <- preOut$C
    V <- preOut$V

    G <- data.frame(gene=as.character(unique(C$gene)))

    message("Loading covariate data")
    #V <- fread(covariate_file,header = TRUE,sep = '\t')
    f <- colnames(V)
    cvnames <- f[2:ncol(V)]
    nv <- length(cvnames)
    gidx <- match(G$gene,V$gene)
    for (i in 1:nv){
      covariate <- V[[cvnames[i]]]
      G[[cvnames[i]]] <- covariate[gidx]
    }

    f <- colnames(C)
    coverage_patient_names <- f[4:ncol(C)]

    #remove any genes that we don't have coverage for
    bad_gene  <- setdiff(M$gene,C$gene)
    if (length(bad_gene)!=0){
      #message("%d/%d gene names could not be mapped to coverage information. Excluding them. \n",length(bad_gene),
      #            length(unique(M$gene)))
      flag_remove <- which(!(M$gene %in% bad_gene))
      M <- M[flag_remove,]
    }

    #map categories
    K <- data.frame(name = stringr::str_sort(unique(C$categ),locale = "C"))
    C$categ_idx <- match(C$categ,K$name)
    ncat <- length(K$name)
    M$categ_idx <- match(M$categ,K$name)

    #make sure there is a null+indel category
    if (is.numeric(unique(C$categ))){
      null_categ <- length(K$name)
    }else{
      null_categ <- grep('null|indel',K$name)
      if (length(null_categ)==0){
        stop("ERROOR: no null/indel category. \n")
      }
      if (length(null_categ)>1){
        stop("ERROR: multiple null/indel categories. \n")
      }
    }

    # make sure C is sorted by the same gene order as in G
    C$gene_idx <- match(C$gene,G$gene)
    C <- dplyr::arrange(C,gene_idx)

    #map genes
    M$gene_idx <- match(M$gene,G$gene)

    #regularize the sample name in the mutation table
    name_before <- M$patient
    M$patient <- gsub("-Tumor$","",M$patient)
    if (any(name_before != M$patient)){
      #message("Trimming '-Tumor' from patient names.")
    }
    name_before <- M$patient
    M$patient <- gsub("-","_",M$patient)
    if (any(name_before != M$patient)){
      #message("Converting '-' to '_' in patient names.")
    }

    patient <- stringr::str_sort(unique(M$patient),locale = "C")
    pat <- data.frame(name = patient)
    pat$cov_idx <- match(pat$name,coverage_patient_names)
    M$patient_idx <- match(M$patient,patient)
    np <- nrow(pat)
    if (np < 2){
      stop("DriverPathway is not applicable to single patients.")
    }

    #is generic coverage data given?
    generic_column_name <- "coverage"
    if ((length(coverage_patient_names)>1) || (length(coverage_patient_names)==1 &
                                               (coverage_patient_names[1]!=generic_column_name))){
      if (any(coverage_patient_names %in% generic_column_name)){
        stop(gettextf("reserved name '%s' cannot appear in list of patient names",
                      generic_column_name))
      }
      if (length(coverage_patient_names)!=length(unique(coverage_patient_names))){
        stop("patient names in coverage_file must be unique")
      }
      #make sure all patients are accounted for in coverage file
      if (any(is.na(pat$cov_idx))){
        stop("some patients in mutation_file are not accounted for in
                 coverage_file")
      }
      generic_coverage_flag <- FALSE
    }else{
      pat$cov_idx <- 1
      generic_coverage_flag <- TRUE
    }

    #BUILD n and N tables
    message("Building n and N tables...\n")

    gene <- sort(as.character(unique(C$gene)))
    categ_name <- as.character(K$name)
    categ_name[length(categ_name)+1] <- "total"
    ng <- length(unique(C$gene))
    n_silent <- array(0,c(ng,ncat+1,np),dimnames=list(gene,categ_name,patient))
    n_nonsilent <- array(0,c(ng,ncat+1,np),dimnames=list(gene,categ_name,patient))
    n_noncoding <- array(0,c(ng,ncat+1,np),dimnames=list(gene,categ_name,patient))

    for (i in 1:ncat){
      for (j in 1:np){
        silent_hist <- graphics::hist(M$gene_idx[(M$effect %in% "silent") &
                                                   (M$categ_idx %in% i) & (M$patient %in% patient[j])],
                                      c(0:length(gene)), plot =FALSE)
        n_silent[,i,j] <-  silent_hist$counts
      }
    }

    for (i in 1:ncat){
      for (j in 1:np){
        nonsilent_hist <- graphics::hist(M$gene_idx[((M$effect %in% "nonsilent")|
                                                       (M$effect %in% "null")) & (M$categ_idx %in% i) &
                                                      (M$patient %in% patient[j])],c(0:length(gene)),
                                         plot =FALSE)
        n_nonsilent[,i,j] <-  nonsilent_hist$counts
      }
    }


    for (i in 1:ncat){
      for (j in 1:np){
        noncoding_hist <- graphics::hist(M$gene_idx[(M$effect %in% "noncoding") &
                                                      (M$categ_idx %in% i) & (M$patient %in% patient[j])],
                                         c(0:length(gene)), plot =FALSE)
        n_noncoding[,i,j] <-  noncoding_hist$counts
      }
    }


    N_silent <- array(0,c(ng,ncat+1,np),dimnames=list(gene,categ_name,patient))
    N_nonsilent <- array(0,c(ng,ncat+1,np),dimnames=list(gene,categ_name,patient))
    N_noncoding <- array(0,c(ng,ncat+1,np),dimnames=list(gene,categ_name,patient))

    for (i in 1:ncat)
      for (j in 1:np)
      {
        N_silent[,i,j] <-  C$coverage[(C$categ %in% categ_name[i]) & (C$effect %in% "silent")]
        N_nonsilent[,i,j] <-  C$coverage[(C$categ %in% categ_name[i]) & (C$effect %in% "nonsilent")]
        N_noncoding[,i,j] <-  C$coverage[(C$categ %in% categ_name[i]) & (C$effect %in% "noncoding")]
      }


    #MAKE SURE ALL NUMBERS ARE INTEGERS
    n_silent <- round(n_silent)
    n_nonsilent  <- round(n_nonsilent)
    n_noncoding  <- round(n_noncoding)
    N_silent  <- round(N_silent)
    N_nonsilent  <- round(N_nonsilent)
    N_noncoding  <- round(N_noncoding)

    #REMOVE MUTATIONS IN BINS WITH EXTREMELY LOW COVERAGE
    n_silent[n_silent>N_silent] <- 0
    n_nonsilent[n_nonsilent>N_nonsilent] <- 0
    n_noncoding[n_noncoding>N_noncoding] <- 0

    #SANITY CHECKS ON TOTALS
    tot_n_nonsilent <- sum(n_nonsilent)
    tot_N_nonsilent <- sum(N_nonsilent)
    tot_n_silent <- sum(n_silent)
    tot_N_silent <- sum(N_silent)
    tot_n_noncoding <- sum(n_noncoding)
    tot_N_noncoding <- sum(N_noncoding)
    tot_rate_nonsilent <- tot_n_nonsilent/tot_N_nonsilent
    tot_rate_silent <- tot_n_silent/tot_N_silent
    tot_rate_noncoding <- tot_n_noncoding/tot_N_noncoding
    tot_rate_coding <- (tot_n_nonsilent+tot_n_silent)/
      (tot_N_nonsilent+tot_N_silent)

    min_tot_n_nonsilent <- 50
    min_tot_n_silent <- 50
    min_tot_n_noncoding <- 50
    min_rate_nonsilent <- 1e-9
    max_rate_nonsilent <- 1e-3
    min_rate_silent <- 1e-9
    max_rate_silent <- 1e-3
    min_rate_noncoding <- 1e-9
    max_rate_noncoding <- 1e-3
    max_abs_log2_difference_nonsilent_silent <- 1.0
    max_abs_log2_difference_noncoding_coding <- 1.0

    #see if silent and nonsilent are OK: if not, give warning
    # if (tot_n_nonsilent<min_tot_n_nonsilent || tot_n_silent<min_tot_n_silent)
    #   stop("not enough mutations to analyze")
    if (tot_rate_nonsilent<min_rate_nonsilent ||
        tot_rate_nonsilent>max_rate_nonsilent)
      stop("nonsilent mutation rate out of range")
    if (tot_rate_silent<min_rate_silent || tot_rate_silent>max_rate_silent)
      stop("silent mutation rate out of range")
    abs_log2_difference_nonsilent_silent = abs(log2(tot_rate_nonsilent/
                                                      tot_rate_silent));
    if (abs_log2_difference_nonsilent_silent>
        max_abs_log2_difference_nonsilent_silent)
      warning('silent and nonsilent rates are too different')

    ## see if noncoding is OK: if not, give warning and zero it all out
    ok = FALSE
    {if (tot_n_noncoding==0)
      message('no noncoding mutations.')
      else
      { if (tot_n_noncoding<min_tot_n_noncoding)
        warning('not enough noncoding mutations to analyze')
        else
        { if (tot_rate_noncoding<min_rate_noncoding ||
              tot_rate_noncoding>max_rate_noncoding)
          warning('noncoding mutation rate out of range')
          else
          { abs_log2_difference_noncoding_coding <-
            abs(log2(tot_rate_noncoding/tot_rate_coding));
          if (abs_log2_difference_noncoding_coding >
              max_abs_log2_difference_noncoding_coding)
            warning("coding and noncoding rates are too different")
          else
            ok = TRUE
          }
        }
      }
    }

    if (!ok){
      message('Zeroing out all noncoding mutations and coverage for the rest of the calculation.');
      n_noncoding[,,]= 0;
      N_noncoding[,,] = 0;
    }

    #add total columns
    N_silent[,ncat+1,] <- N_silent[,null_categ,]
    N_nonsilent[,ncat+1,] <- N_nonsilent[,null_categ,]
    N_noncoding[,ncat+1,] <- N_noncoding[,null_categ,]
    n_silent[,ncat+1,] <- apply(n_silent,c(1,3),sum)
    n_nonsilent[,ncat+1,] <- apply(n_nonsilent,c(1,3),sum)
    n_noncoding[,ncat+1,] <- apply(n_noncoding,c(1,3),sum)

    #total across patients, save in G
    G$N_nonsilent <- apply(N_nonsilent[,ncat+1,],1,sum)
    G$N_silent <- apply(N_silent[,ncat+1,],1,sum)
    G$N_noncoding <- apply(N_noncoding[,ncat+1,],1,sum)
    G$n_nonsilent <- apply(n_nonsilent[,ncat+1,],1,sum)
    G$n_silent <- apply(n_silent[,ncat+1,],1,sum)
    G$n_noncoding <- apply(n_noncoding[,ncat+1,],1,sum)

    G1 <- G
    G1$pb <- NaN
    for (i in 1:nrow(G)) {
      G1$pb[i] <- stats::binom.test(G1$n_nonsilent[i],G1$N_nonsilent[i],bmr,alternative = "greater")$p.value
    }
    G1$qb <- p.adjust(G1$pb,method = "BH",n=length(G1$pb))
    ord_pb <- order(G1$qb)
    G1 <- G1[ord_pb,]
    analysis_gene <- G1$gene[which(G1$pb <= 0.1)]
    if(length(analysis_gene) < 50){
      analysis_gene <- G1$gene[order(G1$pb)[1:50]]
    }

    if(length(analysis_gene) > 500){
      analysis_gene <- G1$gene[which(G1$pb <= 0.05)]
      if(length(analysis_gene) > 500){
        analysis_gene <- G1$gene[which(G1$pb <= 0.01)]
        if(length(analysis_gene) > 500){
          analysis_gene <- G1$gene[which(G1$qb <= 0.05)]
          if(length(analysis_gene) > 500){
            analysis_gene <- G1$gene[which(G1$qb <= 0.01)]
            if(length(analysis_gene) > 500){
              analysis_gene <- G1$gene[which(G1$qb <= 0.005)]
              if(length(analysis_gene) > 500){
                analysis_gene <- G1$gene[which(G1$qb <= 0.001)]
              }
            }
          }
        }
      }
    }
    select.M <- M[M$gene %in% analysis_gene,]
    mutation_matrix <- generate_mutation_matrix(select.M, analysis_gene)
    message("Return preprocessed mutation matrix")
    return(mutation_matrix)
  }else{
    stop("Input data of DriverPathway function should be 0/1 mutation matrix or MAF file.")
  }
}
