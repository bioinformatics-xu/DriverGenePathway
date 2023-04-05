#' Find background mutation rate
#'
#' It uses the data output by the preprocessing function in order to get
#' background mutation rate. It uses the method of bagel. The result will be then
#' used in varified specific gene detecting methods.
#'
#' @param preOutM Mutation data output by preprocessing function or designated by users.
#' @param preOutC Coverage data output by preprocessing function or designated by users.
#' @param preOutV Covariate data which includes expression level, replication time and chromatin compartment defaultly, and could be altered by users.
#' @param bmr The default background mutation rate is 1.2e-6, and the value alters when function ends.
#' @return The output BMR_out is an input to the sigGenes function.
#' @author Xiaolu Xu <lu.xu@@lnnu.edu.cn>
BMR <- function(preOutM, preOutC, preOutV, bmr=1.2e-6)
{
  message("Building BMR")
  gene_idx = NULL
  options(digits=15)
  M<-preOutM
  C<-preOutC
  V<-preOutV

  G <- data.frame(gene=as.character(unique(C$gene)))
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
    #message(sprintf("%d/%d gene names could not be mapped to coverage information. Excluding them.",length(bad_gene),length(unique(M$gene))))
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
      stop("no null/indel category.")
    }
    if (length(null_categ)>1){
      stop("multiple null/indel categories.")
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
    stop("Single patients is not applicable.")
  }

  #is generic coverage data given?
  generic_column_name <- "coverage"
  if ((length(coverage_patient_names)>1) || (length(coverage_patient_names)==1 &
                                             (coverage_patient_names[1]!=generic_column_name))){
    if (any(coverage_patient_names %in% generic_column_name)){
      stop(gettextf("Reserved name '%s' cannot appear in list of patient names",
                    generic_column_name))
    }
    if (length(coverage_patient_names)!=length(unique(coverage_patient_names))){
      stop("Patient names in coverage_file must be unique")
    }
    #make sure all patients are accounted for in coverage file
    if (any(is.na(pat$cov_idx))){
      stop("Some patients in mutation_file are not accounted for in coverage_file")
    }
    generic_coverage_flag <- FALSE
  }else{
    pat$cov_idx <- 1
    generic_coverage_flag <- TRUE
  }

  #BUILD n and N tables
  message("Building n and N tables.")

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
    message('silent and nonsilent rates are too different.')

  ## see if noncoding is OK: if not, give warning and zero it all out
  ok = FALSE
  {if (tot_n_noncoding==0)
    message('no noncoding mutations')
    else
    { if (tot_n_noncoding<min_tot_n_noncoding)
      warning('not enough noncoding mutations to analyze')
      else
      { if (tot_rate_noncoding<min_rate_noncoding ||
            tot_rate_noncoding>max_rate_noncoding)
        warning('noncoding mutation rate out of range')
        else
        {abs_log2_difference_noncoding_coding <-
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
    message('Zeroing out all noncoding mutations and coverage for the rest of the calculation');
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

  # PROCESS COVARIATES
  #message("Processing covariates")

  V <- matrix(data = NaN,nrow = ng,ncol = nv)
  for (i in 1:nv) {
    V[,i] <- G[[cvnames[i]]]
  }
  colnames(V) <- cvnames

  # Find Bagels
  message("Finding bagels")
  max_neighbors <- 50
  qual_min <- 0.05

  Z <- scale(V)

  G$nnei <- NaN
  G$x <- NaN
  G$X <- NaN

  G1 <- G
  G1$pb <- NaN
  for (i in 1:nrow(G)) {
    G1$pb[i] <- stats::binom.test(G1$n_nonsilent[i],G1$N_nonsilent[i],bmr,alternative = "greater")$p.value
  }
  G1$qb <- stats::p.adjust(G1$pb,method = "BH",n=length(G1$pb))
  ord_pb <- order(G1$qb)
  G1 <- G1[ord_pb,]
  analysis_gene <- G1$gene[which(G1$qb <= 0.05)]

  for (g in analysis_gene){
    gi <- which(G$gene %in% g)
    if (which(analysis_gene %in% g) %% 1000 == 0) print(gi)
    df2= (Z-matrix(rep(Z[gi,],ng),nrow=ng,ncol=nv,byrow = TRUE))^2
    dfz  = df2
    dfz[is.nan(dfz)] <-0
    dist2 <- apply(dfz,1,sum)/apply(!is.nan(df2),1,sum)
    dist2 <- round(dist2,15)
    ord <- order(dist2,decreasing = FALSE,na.last = TRUE)
    ord <- c(gi,ord[ord!=gi])
    nfit <-0
    Nfit <-0
    for (ni in 0:max_neighbors)
    {
      gidx = ord[ni+1]
      ngene = G$n_silent[gidx]+ G$n_noncoding[gidx]
      Ngene = G$N_silent[gidx] + G$N_noncoding[gidx]
      if (ni==0)
      {
        ngene0=ngene
        Ngene0=Ngene
      }

      hc <- fun_hc(ngene,Ngene,ngene0,Ngene0)
      qual_left = min(hc,1-hc)
      qual <- 2*qual_left

      nfit=nfit+ngene
      Nfit=Nfit+Ngene
      if (ni>0 && qual<qual_min)
        break
      G$nnei[gi] <- ni
      G$x[gi] <- nfit
      G$X[gi] <- Nfit
    }

  } # of for (gi in 1:ng)

  message("Expanding to (x,X).gcp")
  n_gcp = n_nonsilent + n_silent + n_noncoding
  N_gcp = N_nonsilent + N_silent + N_noncoding

  n_cp = apply(n_gcp,c(2,3),sum)
  N_cp = apply(N_gcp,c(2,3),sum)

  n_c = apply(n_cp,1,sum)
  N_c = apply(N_cp,1,sum)
  mu_c = n_c/N_c

  n_tot = n_c[length(n_c)]
  N_tot = N_c[length(N_c)]
  mu_tot = n_tot/N_tot
  f_c = mu_c/mu_tot
  f_Nc = N_c/N_tot

  n_p = n_cp[nrow(n_cp),]
  N_p = N_cp[nrow(N_cp),]
  mu_p = n_p/N_p
  f_p = mu_p/mu_tot
  f_Np = N_p/mean(N_p)

  x_gcp = array(rep(G$x,(ncat+1)*np),c(nrow(G),ncat+1,np))
  X_gcp = array(rep(G$X,(ncat+1)*np),c(nrow(G),ncat+1,np))

  for (ci in 1:(ncat+1)) {
    x_gcp[,ci,] <- x_gcp[,ci,]*(f_c[ci]*f_Nc[ci])
    X_gcp[,ci,] <- X_gcp[,ci,]*f_Nc[ci]
  }
  for (pi in 1:np) {
    x_gcp[,,pi] <- x_gcp[,,pi]*(f_p[pi]*f_Np[pi])
    X_gcp[,,pi] <- X_gcp[,,pi]*f_Np[pi]
  }
  message("Building BMR finished")
  return(list(G=G,N_silent=N_silent,n_silent=n_silent,
              N_nonsilent=N_nonsilent,n_nonsilent=n_nonsilent,
              N_noncoding=N_noncoding,n_noncoding=n_noncoding,
              X_gcp=X_gcp,x_gcp=x_gcp,null_categ=null_categ,
              analysis_gene=analysis_gene))

} # of MutSig_runCV function

