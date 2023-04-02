# Run different hypothesis test method to identify significant genes using the
# files from preprocess function.
# p_class The hypothesis test method which could be designated by users.
# p_class options include binomial(binomial test method),
# beta.binomial(beta binomial test method),fisher(Fisher combined binomial method),
# LRT(likelihood ratio test method), CT(Convolution test method) and allTest


#' Identify significante genes
#'
#' Applying different hypothesis test methods to identify significant genes on output
#' of the BMR function
#'
#' @param BMR_out Data output by the BMR function.
#' @param p_class "BB" represents beta binomial distribution test; "FCPT" represents Fisher combined P-value test;
#' "LRT" represents likelihood ratio test; "CT" represents convolution test; "projection" represents projection test
#' method; "allTest" represents the mutual results of all methods.
#' @param output_filestem Unified prefix of output data, could be modified as well
#' @param sigThreshold The threshhold of q-value to judge if the gene is significant
#' @importFrom ggsci pal_npg
#' @importFrom pheatmap pheatmap
#' @importFrom magrittr %>%
#' @importFrom grDevices dev.off pdf
#' @importFrom grid grid.draw
#' @importFrom ggplot2 ggplot geom_point theme_minimal theme labs ggsave aes element_rect element_blank
#' @importFrom ggrepel geom_text_repel
#' @details This function searches the significant genes using different hypothesis test methods,
#' including binomial distribution test, beta binomial distribution test, Fisher combined P-value test,
#' likelihood ratio test and convolution test.
#' @return Signifiant genes

sigGenes <- function(BMR_out, p_class = "allTest",output_filestem = "output",sigThreshold = 0.05)
{
  q.btBinom = p.btBinom = q.fisher = q.lrt = q.ct = q.projection = NULL
  sig_filestem <- paste(output_filestem,"sigGenes",sep = "_")
  if(!dir.exists(sig_filestem)){
    dir.create(sig_filestem)
  }
  setwd(sig_filestem)
  options(digits = 15)
  G <- BMR_out$G
  N_silent_gc <- apply(BMR_out$N_silent,c(1,2),sum)
  N_nonsilent_gc <- apply(BMR_out$N_nonsilent,c(1,2),sum)
  N_noncoding_gc <- apply(BMR_out$N_noncoding,c(1,2),sum)
  n_silent_gc <- apply(BMR_out$n_silent,c(1,2),sum)
  n_nonsilent_gc <- apply(BMR_out$n_nonsilent,c(1,2),sum)
  n_noncoding_gc <- apply(BMR_out$n_noncoding,c(1,2),sum)

  alln_gc <- n_silent_gc + n_nonsilent_gc + n_noncoding_gc
  allN_gc <- N_silent_gc + N_nonsilent_gc + N_noncoding_gc

  X_gc <- apply(BMR_out$X_gcp,c(1,2),sum)
  x_gc <- apply(BMR_out$x_gcp,c(1,2),sum)

  # use different BMR for different genes
  cat_num <- ncol(x_gc)-1

  BMR_N <- X_gc
  BMR_n <- x_gc
  BMR <- x_gc/X_gc

  N_nonsilent_g <- G$N_nonsilent
  n_nonsilent_g <- G$n_nonsilent
  N_silent_g <- G$N_silent
  n_silent_g <- G$n_silent
  N_noncoding_g <- G$N_noncoding
  n_noncoding_g <- G$n_noncoding
  X_g <- G$X
  x_g <- G$x
  analysis_gene = BMR_out$analysis_gene


  if(!p_class %in% c("BB","FCPT","LRT","CT","PJ","allTest")){
    stop("p_class must be one of 'BB','FCPT','LRT','CT','PJ', or 'allTest'")
  }

  if(p_class == "BB"){

    sig_betabinomial <- as.data.frame(matrix(numeric(0),nrow = length(x_g),ncol = 3))
    colnames(sig_betabinomial) <- c("gene","p.btBinom","q.btBinom")
    sig_betabinomial$gene <- G$gene
    sig_betabinomial$p.btBinom <- NA
    for (g in analysis_gene) {
      i <- which(G$gene %in% g)
      sig_betabinomial$p.btBinom[i] <- 1-fun_hc(n_nonsilent_g[i],
                                                N_nonsilent_g[i],
                                                x_g[i],X_g[i])

    }
    sig_betabinomial <- sig_betabinomial[which(!is.na(sig_betabinomial$p.btBinom)),]
    sig_betabinomial$q.btBinom <- stats::p.adjust(sig_betabinomial$p.btBinom,method="BH",nrow(G))
    sig_betabinomial <- sig_betabinomial[order(sig_betabinomial$q.btBinom,decreasing = F),]
    output_file <-  paste(output_filestem,"_",p_class,"_sigGenes.csv",sep = "")
    utils::write.csv(sig_betabinomial,file = output_file,row.names = F)
    setwd("..")
    return(sig_betabinomial)
  }

  if(p_class == "FCPT"){

    # Fisher combined binomial p-value
    sig_fisher <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_fisher) <- c("gene","p.fisher","q.fisher")
    sig_fisher$gene <- G$gene
    sig_fisher$p.fisher <- NA
    FCPT_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_fisher[g,]$p.fisher <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          FCPT_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        q = ( -2 ) * sum( log( FCPT_binom[g,flag] ))
        df = 2 * cat_num
        sig_fisher[g,]$p.fisher <- 1 - stats::pchisq( q, df )
      }
    }
    sig_fisher <- sig_fisher[which(!is.na(sig_fisher$p.fisher)),]
    stats::p.adjust( sig_fisher$p.fisher, method="BH" ,nrow(G)) -> sig_fisher$q.fisher
    sig_fisher <- sig_fisher[order(sig_fisher$q.fisher,decreasing = F),]
    output_file <-  paste(output_filestem,"_",p_class,"_sigGenes.csv",sep = "")
    utils::write.csv(sig_fisher,file = output_file,row.names = F)
    setwd("..")
    return(sig_fisher)
  } # end of if(p_class == "fisher")

  if(p_class == "LRT"){

    sig_lrt <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_lrt) <- c("gene","p.lrt","q.lrt")
    sig_lrt$gene <- G$gene
    sig_lrt$p.lrt <- NA
    # Likelihood ratio test
    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    LRT_lh1 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_lrt[g,]$p.lrt <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
          BMR_raw_n <- n_nonsilent_gc[g,i] + n_silent_gc[g,i] + n_noncoding_gc[g,i]
          BMR_raw_N <- N_nonsilent_gc[g,i] + N_silent_gc[g,i] + N_noncoding_gc[g,i]
          LRT_lh1[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR_raw_n/BMR_raw_N,log = T)
          #LRT_lh1[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],n_nonsilent_gc[g,i]/N_nonsilent_gc[g,i],log = T)
        }
        q = 2 * ( sum( LRT_lh1[g,flag] ) - sum( LRT_lh0[g,flag]))
        df = sum( LRT_lh1[g,flag] != 0 );
        if( df > 0 ) p_lrt = 1 - stats::pchisq( q, df )
        if( df == 0 ) p_lrt = 1
        sig_lrt$p.lrt[g] <- p_lrt
      }
    }
    sig_lrt <- sig_lrt[which(!is.na(sig_lrt$p.lrt)),]
    stats::p.adjust( sig_lrt$p.lrt, method="BH" , nrow(G)) -> sig_lrt$q.lrt
    sig_lrt <- sig_lrt[order(sig_lrt$q.lrt,decreasing = F),]
    output_file <-  paste(output_filestem,"_",p_class,"_sigGenes.csv",sep = "")
    utils::write.csv(sig_lrt,file = output_file,row.names = F)
    setwd("..")
    return(sig_lrt)
  }

  if(p_class == "CT"){

    # Convolution test
    sig_ct <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = 3))
    colnames(sig_ct) <- c("gene","p.ct","q.ct")
    sig_ct$gene <- G$gene
    sig_ct$p.ct <- NA

    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
        }
      }
    }

    gethist <- function( xmax, n, p, ptype = "positive_log" ) {
      stats::dbinom( 0:xmax, n, p ) -> ps;
      ps = ps[ps > 0];
      lastp = 1 - sum( ps );
      if( lastp > 0 ) ps = c( ps, lastp );
      if( ptype == "positive_log" ) ps = -log( ps );
      return( ps );
    }

    binit <- function( x, hmax, bin, dropbin = T ) {
      bs = as.integer( x / bin );
      bs[bs > hmax/bin] = hmax / bin;
      bs[is.na( bs )] = hmax / bin;
      tapply( exp(-x), as.factor( bs ), sum ) -> bs;
      bs = bs[bs>0];
      bs = -log( bs );
      if( dropbin ) bs = as.numeric( bs );
      return( bs );
    }

    convolute_b <- function( a, b ) {
      tt = NULL;
      for( j in b ) { tt = c( tt, ( a + j )); }
      return( tt );
    }

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          ni <- N_nonsilent_gc[g,i]
          ei <- BMR[g,i]
          xmax = 100
          hmax = 25
          bin = 0.001
          gethist(xmax,ni,ei,ptype = "positive_log") -> bi
          binit(bi,hmax, bin ) -> bi
          if( i == flag[1] ) { hist0 = bi; }
          if( i > flag[1] & i < flag[length(flag)] ) { hist0 = convolute_b( hist0, bi ); binit( hist0, hmax, bin ) -> hist0 }
          if( i ==  flag[length(flag)]) { hist0 = convolute_b( hist0, bi ) }
        }
        # Convolution test
        bx = -sum( LRT_lh0[g,flag] )
        p.ct = sum( exp( -hist0[hist0>=bx] ))
        qc = sum( exp( -hist0 ))
        sig_ct$p.ct[g] <- p.ct
      }
    }
    sig_ct <- sig_ct[which(!is.na(sig_ct$p.ct)),]
    stats::p.adjust( sig_ct$p.ct, method="BH" ,nrow(G)) -> sig_ct$q.ct
    sig_ct <- sig_ct[order(sig_ct$q.ct,decreasing = F),]
    output_file <- paste(output_filestem,"_",p_class,"_sigGenes.csv",sep = "")
    utils::write.csv(sig_ct,file = output_file,row.names = F)
    setwd("..")
    return(sig_ct)
  }

  if(p_class == "PJ"){
    #options(digits=15)

    G=BMR_out$G
    N_nonsilent=BMR_out$N_nonsilent
    n_nonsilent=BMR_out$n_nonsilent
    X_gcp=BMR_out$X_gcp
    x_gcp=BMR_out$x_gcp
    null_categ=BMR_out$null_categ

    ######### 2D projection to calculate p value and q ##########

    #PROJECTION
    cat("Calculating p.value using 2D Projection method...")
    null_score_boost = 3
    min_effect_size =1.25
    convolution_numbins =1000

    ncat <- dim(X_gcp)[2]-1
    np <- dim(X_gcp)[3]

    G$p <-NaN
    for (g in analysis_gene){
      g <- which(G$gene %in% g)
      if (g %% 1000 == 0) print(g)
      P0 <- matrix(0,nrow=np,ncol = ncat)
      P1 <- matrix(0,nrow=np,ncol = ncat)
      priority <- matrix(0,nrow=np,ncol = ncat)
      N <- t(N_nonsilent[g,1:ncat,])
      n <- t(n_nonsilent[g,1:ncat,])
      x <- t(x_gcp[g,1:ncat,])
      X <- t(X_gcp[g,1:ncat,])

      P0 <- fun_hp(0,N,x,X)
      P1 <- fun_hp(1,N,x,X)
      priority <- t(apply(P1,1,order,decreasing=T))

      shft <- priority - matrix(rep(1:ncat,np),nrow=np,ncol=ncat,byrow = TRUE)
      map <- matrix(1:(np*ncat),nrow = np,ncol = ncat,byrow = FALSE)
      newmap <- map +shft*np
      P0 <- matrix(P0[newmap[1:length(newmap)]],nrow = np)
      P1 <- matrix(P1[newmap[1:length(newmap)]],nrow = np)
      P2 <- 1-(P0+P1)
      P2[P2<0] <- 0

      Pdeg <- array(0,c(np,ncat+1,ncat+1))

      for (d1 in 0:ncat){
        for (d2 in 0:d1)
        {
          if (d1 == ncat)
            p <- matrix(1,nrow = np,ncol = 1)
          else
            p <- apply(as.matrix(P0[,(d1+1):ncol(P0)]),1,prod)
          if (d1>0)
          { if (d1==d2)
            p = p*P2[,d1]
          else
          {
            p = p * P1[,d1]
            if (d2<(d1-1))
              p = p * apply(as.matrix(P0[,(d2+1):(d1-1)]),1,prod)
            if (d2>0)
              p = p * (P1[,d2]+P2[,d2])
          }
          }
          Pdeg[,d1+1,d2+1] <- p
        }
      }


      Sdeg = array(0,c(np,ncat+1,ncat+1))
      for (d1 in 1:ncat)
      {
        for (d2 in 0:d1)
        {
          if (d1==d2){
            p = P2[,d1]
          }else
          {
            if (d2>0)
              p = P1[,d1]*P1[,d2]
            else
              p = P1[,d1]
          }

          Sdeg[,d1+1,d2+1] = -log10(p)
        }
      }

      priority2 <- cbind(matrix(0,nrow = np,ncol = 1),priority)

      #Sdeg[priority2==null_categ] = Sdeg[priority2==null_categ] + null_score_boost
      s1 <- Sdeg[,,1]
      Sdeg[,,1][priority2==null_categ] <- s1[priority2==null_categ]+null_score_boost
      degree <-  matrix(0,nrow=np,ncol=2)
      score_obs = 0

      for(p in  1:np)
      {
        i = 1
        for (d in ncat:1)
        {
          c = priority[p,d]
          if (i==1)
          {
            if (n[p,c]>=2)
            {
              degree[p,] = c(d,d)
              i = 3
            }
            else if (n[p,c]==1)
            {
              degree[p,i] = d
              i=i+1
            }
          }
          else if (i==2)
          {
            if (n[p,c]>=1)
            {
              degree[p,i] = d
              i=i+1
            }
            else
              break
          }
        }
        score_sample = Sdeg[p,degree[p,1]+1,degree[p,2]+1]
        if(!is.finite(score_sample)){
          cat(sprintf("g=%d",g))
          cat(sprintf("p=%d",p))
          flag_finite <- is.finite(Sdeg[,degree[p,1]+1,degree[p,2]+1])
          score_sample <- max(Sdeg[,degree[p,1]+1,degree[p,2]+1][flag_finite])
          cat(sprintf("maxscore=%f",score_sample))
        }
        score_obs = score_obs + score_sample
      }

      score_obs = score_obs / min_effect_size

      if (score_obs<=0)
      {
        G$p[g]=1
        next
      }
      numbins = convolution_numbins
      binsize = score_obs / numbins
      H = matrix(0,nrow=numbins,ncol=1)
      H[1] <- 1

      offset = pmin(array(numbins,c(np,ncat+1,ncat+1)), round(Sdeg/binsize))
      ncols = (ncat+1)*(ncat+2)/2
      newH = matrix(0,nrow=numbins,ncol=ncols)
      for (p in 1:np)
      {
        newH[,] <- 0
        col=1
        for (d1 in 0:ncat)
          for (d2 in 0:d1)
          {
            o = offset[p,d1+1,d2+1]
            if (o<length(H))
              newH[(o+1):nrow(newH),col] = Pdeg[p,d1+1,d2+1] * H[1:(length(H)-o)]
            col=col+1
          }
        H = apply(newH,1,sum,drop=FALSE)
      }

      G$p[g] = max(0,1-sum(H))

    }  # of for (g in 1:ng)
    G <- G[which(!is.na(G$p)),]
    G$q <- stats::p.adjust(G$p,method="BH",nrow(G))
    G <- G[order(G$q),]
    sig_projection <- subset(G,select=c("gene","p","q"))
    colnames(sig_projection) <- c("gene","p.projection","q.projection")
    # output_file <- paste(output_filestem,"_",p_class,"_sigGenes.txt",sep = "")
    # utils::write.table(sig_projection,file = output_file,quote = F,sep = "\t",row.names = F)
    output_file <- paste(output_filestem,"_",p_class,"_sigGenes.csv",sep = "")
    utils::write.csv(sig_projection,file = output_file,row.names = F)
    setwd("..")
    return(sig_projection)
  }

  if(p_class == "allTest"){

    ### beta binomial test to calculate p value and q ###
    sig_betabinomial <- as.data.frame(matrix(numeric(0),nrow = length(x_g),ncol = 3))
    colnames(sig_betabinomial) <- c("gene","p.btBinom","q.btBinom")
    sig_betabinomial$gene <- G$gene
    for (i in 1:length(x_g)) {
      sig_betabinomial$p.btBinom[i] <- 1-fun_hc(n_nonsilent_g[i],
                                                N_nonsilent_g[i],x_g[i],X_g[i])
    }
    sig_betabinomial$q.btBinom <- stats::p.adjust(sig_betabinomial$p.btBinom,method="BH",length(sig_betabinomial$p.btBinom))
    sig_betabinomial <- sig_betabinomial[order(sig_betabinomial$q.btBinom,decreasing = F),]
    output_file <-  paste(output_filestem,"_BB_sigGenes.csv",sep = "")
    utils::write.csv(sig_betabinomial,file = output_file,row.names = F)
    ### beta binomial test to calculate p value and q ###

    ### Fisher combined binomial p-value  ###
    # Fisher combined binomial p-value
    sig_fisher <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_fisher) <- c("gene","p.fisher","q.fisher")
    sig_fisher$gene <- G$gene
    sig_fisher$p.fisher <- NA
    FCPT_binom <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_fisher[g,]$p.fisher <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & (N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num]) & BMR[g,1:cat_num]>0)
        for(i in flag){
          FCPT_binom[g,i] <- stats::binom.test(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],alternative = "greater")$p.value
        }
        q = ( -2 ) * sum( log( FCPT_binom[g,flag] ))
        df = 2 * cat_num
        sig_fisher[g,]$p.fisher <- 1 - stats::pchisq( q, df )
      }
    }

    stats::p.adjust( sig_fisher$p.fisher, method="BH" ,length(sig_fisher$p.fisher)) -> sig_fisher$q.fisher
    sig_fisher <- sig_fisher[order(sig_fisher$q.fisher,decreasing = F),]
    output_file <-  paste(output_filestem,"_FCPT_sigGenes.csv",sep = "")
    utils::write.csv(sig_fisher,file = output_file,row.names = F)


    ### Fisher combined binomial p-value end  ###

    ### likelihood ratio test start ###
    # Likelihood ratio test

    sig_lrt <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = 3))
    colnames(sig_lrt) <- c("gene","p.lrt","q.lrt")
    sig_lrt$gene <- G$gene
    sig_lrt$p.lrt <- NA
    # Likelihood ratio test
    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))
    LRT_lh1 <- as.data.frame(matrix(data = NA,nrow = length(x_g),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_lrt[g,]$p.lrt <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
          BMR_raw_n <- n_nonsilent_gc[g,i] + n_silent_gc[g,i] + n_noncoding_gc[g,i]
          BMR_raw_N <- N_nonsilent_gc[g,i] + N_silent_gc[g,i] + N_noncoding_gc[g,i]
          LRT_lh1[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR_raw_n/BMR_raw_N,log = T)
          #LRT_lh1[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],n_nonsilent_gc[g,i]/N_nonsilent_gc[g,i],log = T)
        }
        q = 2 * ( sum( LRT_lh1[g,flag] ) - sum( LRT_lh0[g,flag]))
        df = sum( LRT_lh1[g,flag] != 0 );
        if( df > 0 ) p_lrt = 1 - stats::pchisq( q, df )
        if( df == 0 ) p_lrt = 1
        sig_lrt$p.lrt[g] <- p_lrt
      }
    }


    stats::p.adjust( sig_lrt$p.lrt, method="BH" ,length(sig_lrt$p.lrt)) -> sig_lrt$q.lrt
    sig_lrt <- sig_lrt[order(sig_lrt$q.lrt,decreasing = F),]
    output_file <-  paste(output_filestem,"_LRT_sigGenes.csv",sep = "")
    utils::write.csv(sig_lrt,file = output_file,row.names = F)

    ### likelihood ratio test end ###


    ### Convolution test start ###
    sig_ct <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = 3))
    colnames(sig_ct) <- c("gene","p.ct","q.ct")
    sig_ct$gene <- G$gene
    sig_ct$p.ct <- NA

    LRT_lh0 <- as.data.frame(matrix(data = NA,nrow = nrow(G),ncol = cat_num))

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          LRT_lh0[g,i] <- stats::dbinom(n_nonsilent_gc[g,i],N_nonsilent_gc[g,i],BMR[g,i],log = T)
        }
      }
    }

    gethist <- function( xmax, n, p, ptype = "positive_log" ) {
      stats::dbinom( 0:xmax, n, p ) -> ps;
      ps = ps[ps > 0];
      lastp = 1 - sum( ps );
      if( lastp > 0 ) ps = c( ps, lastp );
      if( ptype == "positive_log" ) ps = -log( ps );
      return( ps );
    }

    binit <- function( x, hmax, bin, dropbin = T ) {
      bs = as.integer( x / bin );
      bs[bs > hmax/bin] = hmax / bin;
      bs[is.na( bs )] = hmax / bin;
      tapply( exp(-x), as.factor( bs ), sum ) -> bs;
      bs = bs[bs>0];
      bs = -log( bs );
      if( dropbin ) bs = as.numeric( bs );
      return( bs );
    }

    convolute_b <- function( a, b ) {
      tt = NULL;
      for( j in b ) { tt = c( tt, ( a + j )); }
      return( tt );
    }

    for(g in analysis_gene){
      g <- which(G$gene %in% g)
      if(n_nonsilent_gc[g,cat_num+1]<=0 | N_nonsilent_gc[g,cat_num+1] <= 0 | BMR[g,cat_num+1] <= 0){
        sig_ct[g,]$p.ct <- 1
      }else{
        flag <- which(N_nonsilent_gc[g,1:cat_num] >0 & N_nonsilent_gc[g,1:cat_num] > n_nonsilent_gc[g,1:cat_num] & BMR[g,1:cat_num]>0)
        for(i in flag){
          ni <- N_nonsilent_gc[g,i]
          ei <- BMR[g,i]
          xmax = 100
          hmax = 25
          bin = 0.001
          gethist(xmax,ni,ei,ptype = "positive_log") -> bi
          binit(bi,hmax, bin ) -> bi
          if( i == flag[1] ) { hist0 = bi; }
          if( i > flag[1] & i < flag[length(flag)] ) { hist0 = convolute_b( hist0, bi ); binit( hist0, hmax, bin ) -> hist0 }
          if( i ==  flag[length(flag)]) { hist0 = convolute_b( hist0, bi ) }
        }
        # Convolution test
        bx = -sum( LRT_lh0[g,flag] )
        p.ct = sum( exp( -hist0[hist0>=bx] ))
        qc = sum( exp( -hist0 ))
        sig_ct$p.ct[g] <- p.ct
      }
    }
    sig_ct <- sig_ct[which(!is.na(sig_ct$p.ct)),]
    stats::p.adjust( sig_ct$p.ct, method="BH" ,nrow(G)) -> sig_ct$q.ct
    sig_ct <- sig_ct[order(sig_ct$q.ct,decreasing = F),]
    output_file <- paste(output_filestem,"_CT_sigGenes.csv",sep = "")
    utils::write.csv(sig_ct,file = output_file,row.names = F)


  G=BMR_out$G
  N_nonsilent=BMR_out$N_nonsilent
  n_nonsilent=BMR_out$n_nonsilent
  X_gcp=BMR_out$X_gcp
  x_gcp=BMR_out$x_gcp
  null_categ=BMR_out$null_categ

  ######### 2D projection to calculate p value and q ##########

  #PROJECTION
  cat("Calculating p.value using 2D Projection method...")
  null_score_boost = 3
  min_effect_size =1.25
  convolution_numbins =1000

  ncat <- dim(X_gcp)[2]-1
  np <- dim(X_gcp)[3]

  G$p <-NaN
  for (g in analysis_gene){
    g <- which(G$gene %in% g)
    if (g %% 1000 == 0) print(g)
    P0 <- matrix(0,nrow=np,ncol = ncat)
    P1 <- matrix(0,nrow=np,ncol = ncat)
    priority <- matrix(0,nrow=np,ncol = ncat)
    N <- t(N_nonsilent[g,1:ncat,])
    n <- t(n_nonsilent[g,1:ncat,])
    x <- t(x_gcp[g,1:ncat,])
    X <- t(X_gcp[g,1:ncat,])

    P0 <- fun_hp(0,N,x,X)
    P1 <- fun_hp(1,N,x,X)
    priority <- t(apply(P1,1,order,decreasing=T))

    shft <- priority - matrix(rep(1:ncat,np),nrow=np,ncol=ncat,byrow = TRUE)
    map <- matrix(1:(np*ncat),nrow = np,ncol = ncat,byrow = FALSE)
    newmap <- map +shft*np
    P0 <- matrix(P0[newmap[1:length(newmap)]],nrow = np)
    P1 <- matrix(P1[newmap[1:length(newmap)]],nrow = np)
    P2 <- 1-(P0+P1)
    P2[P2<0] <- 0

    Pdeg <- array(0,c(np,ncat+1,ncat+1))

    for (d1 in 0:ncat){
      for (d2 in 0:d1)
      {
        if (d1 == ncat)
          p <- matrix(1,nrow = np,ncol = 1)
        else
          p <- apply(as.matrix(P0[,(d1+1):ncol(P0)]),1,prod)
        if (d1>0)
        { if (d1==d2)
          p = p*P2[,d1]
        else
        {
          p = p * P1[,d1]
          if (d2<(d1-1))
            p = p * apply(as.matrix(P0[,(d2+1):(d1-1)]),1,prod)
          if (d2>0)
            p = p * (P1[,d2]+P2[,d2])
        }
        }
        Pdeg[,d1+1,d2+1] <- p
      }
    }


    Sdeg = array(0,c(np,ncat+1,ncat+1))
    for (d1 in 1:ncat)
    {
      for (d2 in 0:d1)
      {
        if (d1==d2){
          p = P2[,d1]
        }else
        {
          if (d2>0)
            p = P1[,d1]*P1[,d2]
          else
            p = P1[,d1]
        }

        Sdeg[,d1+1,d2+1] = -log10(p)
      }
    }

    priority2 <- cbind(matrix(0,nrow = np,ncol = 1),priority)

    #Sdeg[priority2==null_categ] = Sdeg[priority2==null_categ] + null_score_boost
    s1 <- Sdeg[,,1]
    Sdeg[,,1][priority2==null_categ] <- s1[priority2==null_categ]+null_score_boost
    degree <-  matrix(0,nrow=np,ncol=2)
    score_obs = 0

    for(p in  1:np)
    {
      i = 1
      for (d in ncat:1)
      {
        c = priority[p,d]
        if (i==1)
        {
          if (n[p,c]>=2)
          {
            degree[p,] = c(d,d)
            i = 3
          }
          else if (n[p,c]==1)
          {
            degree[p,i] = d
            i=i+1
          }
        }
        else if (i==2)
        {
          if (n[p,c]>=1)
          {
            degree[p,i] = d
            i=i+1
          }
          else
            break
        }
      }
      score_sample = Sdeg[p,degree[p,1]+1,degree[p,2]+1]
      if(!is.finite(score_sample)){
        cat(sprintf("g=%d",g))
        cat(sprintf("p=%d",p))
        flag_finite <- is.finite(Sdeg[,degree[p,1]+1,degree[p,2]+1])
        score_sample <- max(Sdeg[,degree[p,1]+1,degree[p,2]+1][flag_finite])
        cat(sprintf("maxscore=%f",score_sample))
      }
      score_obs = score_obs + score_sample
    }

    score_obs = score_obs / min_effect_size

    if (score_obs<=0)
    {
      G$p[g]=1
      next
    }
    numbins = convolution_numbins
    binsize = score_obs / numbins
    H = matrix(0,nrow=numbins,ncol=1)
    H[1] <- 1

    offset = pmin(array(numbins,c(np,ncat+1,ncat+1)), round(Sdeg/binsize))
    ncols = (ncat+1)*(ncat+2)/2
    newH = matrix(0,nrow=numbins,ncol=ncols)
    for (p in 1:np)
    {
      newH[,] <- 0
      col=1
      for (d1 in 0:ncat)
        for (d2 in 0:d1)
        {
          o = offset[p,d1+1,d2+1]
          if (o<length(H))
            newH[(o+1):nrow(newH),col] = Pdeg[p,d1+1,d2+1] * H[1:(length(H)-o)]
          col=col+1
        }
      H = apply(newH,1,sum,drop=FALSE)
    }

    G$p[g] = max(0,1-sum(H))

  }  # of for (g in 1:ng)
  G <- G[which(!is.na(G$p)),]
  G$q <- stats::p.adjust(G$p,method="BH",nrow(G))
  G <- G[order(G$q),]
  sig_projection <- subset(G,select=c("gene","p","q"))
  colnames(sig_projection) <- c("gene","p.projection","q.projection")
  output_file <- paste(output_filestem,"_PJ_sigGenes.csv",sep = "")
  utils::write.csv(sig_projection,file = output_file,row.names = F)



    sigGenes_combined <- c(
                           as.character(sig_betabinomial$gene[sig_betabinomial$q.btBinom <= sigThreshold]),
                           as.character(sig_fisher$gene[sig_fisher$q.fisher <= sigThreshold]),
                           as.character(sig_ct$gene[sig_ct$q.ct <= sigThreshold]),
                           as.character(sig_lrt$gene[sig_lrt$q.lrt <= sigThreshold]),
                           as.character(sig_projection$gene[sig_projection$q.projection <= sigThreshold]))
    sigGenes_table <- as.data.frame(table(sigGenes_combined))
    sigGenes_table <- sigGenes_table[order(sigGenes_table$Freq,decreasing = T),]

    sigGenes_final <- sigGenes_table[which(sigGenes_table$Freq == max(sigGenes_table$Freq)),1]

    q_result <- sig_betabinomial

    gene_beta <- q_result[,1]

    q_result <- subset(q_result,q.btBinom<=0.1,select=-p.btBinom)

    q_new <- sig_fisher

    gene_fisher <- q_new[,1]

    q_new <- subset(q_new,q.fisher<=0.1)
    q_union <- data.frame(gene=union(q_new$gene,q_result$gene))
    flag_old <- match(q_union$gene,q_result$gene)
    flag_new <- match(q_union$gene,q_new$gene)
    q_union[,2] <- q_result[flag_old,2]
    q_union$q.fisher <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_lrt

    gene_lrt <- q_new[,1]

    q_new <- subset(q_new,q.lrt<=0.1)
    q_union <- data.frame(gene=union(q_new$gene,q_result$gene))
    flag_old <- match(q_union$gene,q_result$gene)
    flag_new <- match(q_union$gene,q_new$gene)
    q_union[,2:3] <- q_result[flag_old,2:3]
    q_union$q.lrt <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_ct

    gene_ct <- q_new[,1]

    q_new <- subset(q_new,q.ct<=0.1)
    q_union <- data.frame(gene=union(q_new$gene,q_result$gene))
    flag_old <- match(q_union$gene,q_result$gene)
    flag_new <- match(q_union$gene,q_new$gene)
    q_union[,2:4] <- q_result[flag_old,2:4]
    q_union$q.ct <- q_new[flag_new,3]
    q_result <- q_union

    q_new <- sig_projection

    gene_projection <- q_new[,1]

    q_new <- subset(q_new,q.projection<=0.1)
    q_union <- data.frame(gene=union(q_new$gene,q_result$gene))
    flag_old <- match(q_union$gene,q_result$gene)
    flag_new <- match(q_union$gene,q_new$gene)
    q_union[,2:5] <- q_result[flag_old,2:5]
    q_union$q.projection <- q_new[flag_new,3]
    q_result <- q_union

    names(q_result) <- c("gene","q.btBinom","q.fisher","q.lrt","q.ct","q.projection")

    #sigGenes_plot(q_result)

    utils::write.csv(q_result,file = "q_values.csv",quote = F,row.names = F)

    # create plots
    if(!dir.exists("sigGenes_plots")){
      dir.create("sigGenes_plots")
    }
    setwd("sigGenes_plots")
    gene = lgq.btBinom = lgq.fisher = lgq.LRT = lgq.CT = lgq.PJ = gene = NULL

    # heatmap for q values
    q2 <- q_result%>%
      dplyr::filter(!is.na(q.ct),!is.na(q.lrt),!is.na(q.fisher),!is.na(q.btBinom),!is.na(q.projection))
    q2 <- as.data.frame(q2)
    row.names(q2) <- q2[,1]
    q2[,-1]%>%
      as.matrix()%>%
      t()%>%
      pheatmap(show_colnames = T, fontsize_col = 7, filename = "q_value_heatmap.pdf", display_numbers = T)

    result_BB <- subset(sig_betabinomial,q.btBinom<0.1)
    if(nrow(result_BB)>100){
      result_BB1 <- result_BB[1:100,]
    }
    else{
      result_BB1 <- result_BB
    }
    gene_betabinomial <- result_BB[,1]
    result_BB1$lgq.btBinom <- -log10(result_BB1$q.btBinom+0.00000001)
    BB <- ggplot(result_BB1, aes(x=gene,y=lgq.btBinom))+
      geom_point(aes(size=2*(lgq.btBinom^2)), color = "#00AFBB") +
      geom_text_repel(aes(label = gene), size = 5)+
      theme_minimal()+
      theme(legend.background = element_rect(fill="lightblue"))+
      theme(axis.text.x = element_blank())+
      labs(title="Plot of result genes tested by beta-binomial", x ="Genes", y = "-log(q.beta-binomial)")
    ggsave(BB, filename="BB.pdf", width = 15)

    result_FCPT <- subset(sig_fisher,q.fisher<0.1)
    if(nrow(result_FCPT)>100){
      result_FCPT1 <- result_FCPT[1:100,]
    }
    else{
      result_FCPT1 <- result_FCPT
    }
    gene_fisher <- result_FCPT[,1]
    result_FCPT1$lgq.fisher <- -log10(result_FCPT1$q.fisher+0.00000001)
    FCPT <- ggplot(result_FCPT1, aes(x=gene,y=lgq.fisher))+
      geom_point(aes(size=2*(lgq.fisher^2)), color = "#00AFBB") +
      geom_text_repel(aes(label = gene), size = 5)+
      theme_minimal()+
      theme(legend.background = element_rect(fill="lightblue"))+
      theme(axis.text.x = element_blank())+
      labs(title="Plot of result genes tested by fisher", x ="Genes", y = "-log(q.fisher)")
    ggsave(FCPT, filename="FCPT.pdf", width = 15)

    result_LRT <- subset(sig_lrt,q.lrt<0.1)
    if(nrow(result_LRT)>100){
      result_LRT1 <- result_LRT[1:100,]
    }
    else{
      result_LRT1 <- result_LRT
    }
    gene_lrt <- result_LRT[,1]
    result_LRT1$lgq.LRT <- -log10(result_LRT1$q.lrt+0.00000001)
    LRT <- ggplot(result_LRT1, aes(x=gene,y=lgq.LRT))+
      geom_point(aes(size=2*(lgq.LRT^2)), color = "#00AFBB") +
      geom_text_repel(aes(label = gene), size = 5)+
      theme_minimal()+
      theme(legend.background = element_rect(fill="lightblue"))+
      theme(axis.text.x = element_blank())+
      labs(title="Plot of result genes tested by lrt", x ="Genes", y = "-log(q.lrt)")
    ggsave(LRT, filename="LRT.pdf", width = 15)

    result_CT <- subset(sig_ct,q.ct<0.1)
    if(nrow(result_CT)>100){
      result_CT1 <- result_CT[1:100,]
    }
    else{
      result_CT1 <- result_CT
    }
    gene_ct <- result_CT[,1]
    result_CT1$lgq.CT <- -log10(result_CT1$q.ct+0.00000001)
    CT <- ggplot(result_CT1,aes(x=gene,y=lgq.CT))+
      geom_point(aes(size=2*(lgq.CT^2)), color = "#00AFBB") +
      geom_text_repel(aes(label = gene), size = 5)+
      theme_minimal()+
      theme(legend.background = element_rect(fill="lightblue"))+
      theme(axis.text.x = element_blank())+
      labs(title="Plot of result genes tested by ct", x ="Genes", y = "-log(q.ct)")
    ggsave(CT, filename="CT.pdf", width = 15)

    result_PJ <- subset(sig_projection,q.projection<0.1)
    if(nrow(result_PJ)>100){
      result_PJ1 <- result_PJ[1:100,]
    }
    else{
      result_PJ1 <- result_PJ
    }
    gene_projection <- result_PJ[,1]
    result_PJ1$lgq.PJ <- -log10(result_PJ1$q.projection+0.00000001)
    PJ <- ggplot(result_PJ1,aes(x=gene,y=lgq.PJ))+
      geom_point(aes(size=2*(lgq.PJ^2)), color = "#00AFBB") +
      geom_text_repel(aes(label = gene), size = 5)+
      theme_minimal()+
      theme(legend.background = element_rect(fill="lightblue"))+
      theme(axis.text.x = element_blank())+
      labs(title="Plot of result genes tested by 2D projection", x ="Genes", y = "-log(q.projection)")
    ggsave(PJ, filename="PJ.pdf", width = 15)

    vennplot <- VennDiagram::venn.diagram(x=list(gene_projection=gene_projection,gene_betabinomial=gene_betabinomial,gene_fisher=gene_fisher,gene_lrt=gene_lrt,gene_ct=gene_ct),
                 filename = NULL,
                 fill =pal_npg()(5),
                 cat.col =pal_npg()(5),cat.cex = 1,cat.pos = 0, cat.dist = 0.07,cat.fontfamily = "serif")
    pdf("vennplot.pdf")
    grid.draw(vennplot)
    dev.off()
    #fill =c("cornflowerblue","green","yellow","darkorchid1","red")
    #cat.col =c("darkblue", "darkgreen", "orange","darkorchid4","black")
    setwd("..")
    setwd("..")
    return(sigGenes_final)
  }  # end of if(p_class == "all")
} # of MutSig_runCV function



