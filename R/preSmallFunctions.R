#' Assigning categories for mutation types
#'
#' For each of the 192 mutation types, assigning a corresponding category from number 1 to 4.
#' @param K Dataframe of previously detected mutation categories
#' @author Xiaolu Xu <lu.xu@@lnnu.edu.cn>
preprocessAssignCateg <- function(K){

  parse1 <- function(S,f){
    tokens <- stringr::str_extract_all(S, "[[:alnum:]]+")
    x <- t(as.data.frame(tokens))[,c(1,3,4)]
    colnames(x) <- f
    rownames(x) <- c(1:64)
    return(x)
  }

  nk <- nrow(K)
  idx <- grep("indel|null",K$type)
  if(length(idx)!=0){
    K$type[idx] <- rep("non-point",times = length(idx))
  }
  if("name" %in% colnames(K)){
    idx <- grep("indel|null",K$name)
    if(length(idx)!=0){
      K$type[idx] <- rep("non-point",times = length(idx))
    }
  }
  X <- as.data.frame(matrix(numeric(0),nrow = 64,ncol = 2))
  colnames(X) <- c("num","name")
  X$name <- preprocessGenerateCategContext65()[1:64]
  X$num <- c(1:64)
  X <- X[order(X$num),]
  X <- as.data.frame(parse1(X$name,c("from","left","right")))

  base <- c("A","C","G","T")
  complement <- matrix(data = c("T","G","C","A"),ncol = 4)
  colnames(complement) <- base

  whatchange <- matrix(data = 0,ncol = 4,nrow = 4)
  dimnames(whatchange) <- list(base,base)
  whatchange[1,] <- c("n","s","t","f")
  whatchange[2,] <- c("s","n","f","t")
  whatchange[3,] <- c("t","f","n","s")
  whatchange[4,] <- c("f","t","s","n")

  c <- array(data = 0,dim = c(64,nk,4))

  for (x in 1:64) {
    from <- X$from[x]
    if(length(from)==0){
      next
    }
    left <- X$left[x]
    right <- X$right[x]
    if(grepl(from,"GT")){
      from <- complement[1,from]
      left <- complement[1,X$right[x]]
      right <- complement[1,X$left[x]]
    }
    for (k in 1:nk) {
      if((grepl(from,K$from[k]))&(grepl(left,K$left[k]))&
         (grepl(right, K$right[k]))){
        if(K$type[k]=="non-point"){
          c[x,k,] = 1
        }else if(K$type[k]=="point"){
          for (n in 1:4) {
            oldbase <- X$from[x]
            newbase <- base[n]
            change <- whatchange[oldbase,newbase]
            if(grepl(change,K$change[k])){
              c[x,k,n] = 1
            }
          } # of for(n in 1:4)
        }else{
          stop(sprintf("Unknown type: %s",K$type[k]))
        }
      }
    }
  } #of for (x in 1:64)

  return(c)

} # of preprocessAssignCateg

#' Collapsing category number to 32
#'
#' Collapsing previous category number from 65 to 32.
#' @param Y Dataframe of numbers of coverage and 4 mutation types of 65 categories
#' @details This function collapses the previous category number from 65 to 32 by first removing the "any N" type,
#' then fusing a specific trinucleotide and its opposite.
preprocessCollapseNn65to32 <- function(Y){
  Y_order <- subset(Y,select = c("name","N","T","G","C","A"))
  Y_32 <- subset(Y,select = c("N","A","C","G","T"))[1:32,]
  Y_order32 <- Y_order[64:33,1:6]
  Y_ordername <- Y_order32$name
  Y_name6 <- substr(Y_ordername,6,6)
  Y_name8 <- substr(Y_ordername,8,8)
  stringr::str_sub(Y_ordername,6,6) <- Y_name8
  stringr::str_sub(Y_ordername,8,8) <- Y_name6
  #替换  A in C_G  中C和G的位置，因为正链三连碱基的顺序和逆链相反
  flag_order <- match(Y_ordername,Y_order32$name)
  Y_order32 <- Y_order32[flag_order,]
  Nn <- Y_32+Y_order32[,2:6]
  return(Nn)
}

#' Appending a new catagory to the previous ones
#'
#' Adding a new category representing 'null+indel' to the previously detected 4 categories.
#' @param K Dataframe of previously detected mutation categories
#' @param K2 Dataframe of the newly defined category
#' @details The number of categories is set 4 defaultly, this function append the newly defined 'null+indel' category
#' to the previous ones.
preprocessConcatStructsKeepFields <- function(K,K2){
  col_K <- colnames(K)
  col_K2 <- colnames(K2)
  union_col <- union(col_K2,col_K)
  nocol_K <- setdiff(union_col,col_K)
  nocol_K2 <- setdiff(union_col,col_K2)

  if(length(nocol_K2) != 0){
    x <- as.data.frame(matrix(data = NaN,nrow = nrow(K2),
                              ncol = length(nocol_K2)))
    colnames(x) <- nocol_K2
    K2 <- cbind(K2,x)
    K2 <- K2[,order(colnames(K2))]
  }

  if(length(nocol_K) != 0){
    x <- as.data.frame(matrix(data = NaN,nrow = nrow(K),ncol = length(nocol_K)))
    colnames(x) <- nocol_K2
    K <- cbind(K,x)
    K <- K[,order(colnames(K))]
  }
  K <- rbind(K,K2)
  return(K)
}

############### find_mut_categs function  begin ##############
#' Finding mutation categories
#'
#' Turning the number of categories from 32 to a reasonable expected value, through analysing the numbers of
#' mutations and coverages.
#' @param Nn The frame containing numbers of coverages and mutations for each of the 32 types
#' @param P The frame containing the expected number of categs and the name of output mutation file
#' @details Through analysing the mutation numbers of the 32 types previously defined, this function finds several
#' reasonable mutation categories as categs. The default number of categs is 4.
preprocessFindMutCateg <- function(Nn,P){

  #####   subfunctions of find_mut_categs  #####

  categ_to_int <- function(c){
    xx <- c(c[[1]],c[[2]]+4,c[[3]]+8,c[[4]]+12)-1
    i <- sum(2^xx)
    return(i)
  }

  int_to_categ <- function(i){
    binary_rev_i <- as.integer(intToBits(i))
    x <- binary_rev_i[1:16]
    b <- c(1:4)
    c <- list(b[which(x[1:4] > 0)],b[which(x[5:8] > 0)],b[which(x[9:12] > 0)],b[which(x[13:16] > 0)])
    return(c)
  }

  sumpart <- function(m,cut){

    x <- m[cut[[1]],cut[[2]],cut[[3]],cut[[4]]]
    x <- sum(x)
    return(x)
  }

  entropy <- function(p){
    if(length(p)==0){
      p1 <- NaN
      p2 <- NaN
    }else if(is.nan(p)){
      p1 <- NaN
      p2 <- NaN
    }else if((p==0)|(p==1)){
      p1 <- 0
      p2 <- 0
    }else if((p<0)|(p>1)){
      stop(sprintf("p must be between zero and one"))
    }else{
      p1 <- p*log2(p)
      p2 <- (1-p)*log2(1-p)
    }
    H <- -(p1+p2)
    return(H)
  }

  entropy_by_parts2 <- function(partsi){
    H <- 0
    for (pp in 1:length(partsi)) {
      part <- int_to_categ(partsi[pp])
      ns <- sumpart(n,part)
      Ns <- sumpart(N,part)
      f <- Ns/Ntot
      H_part <- entropy(ns/Ns)
      H <- H +f*H_part
    }
    return(H)
  }

  convert_parts_to_rule <- function(parts){
    bases <- c("A","C","G","T")
    change <- c("t","f","s")
    np <- length(parts)
    s <- list()
    for (p in 1:np) {
      part <- parts[[p]]
      bases_p1 <- paste(bases[part[[1]]],collapse = "")
      bases_p3 <- paste(bases[part[[3]]],collapse = "")
      change_p4 <- paste(change[part[[4]]],collapse = "")
      bases_p2 <- paste(bases[part[[2]]],collapse = "")
      s[[p]] <- paste(bases_p1,"[",bases_p3,"->",change_p4,"]",bases_p2,sep = "")
    }
    return(s)
  }

  rulestats <- function(parts){
    np <- length(parts)
    s <- as.data.frame(matrix(data=0,nrow = np))
    s <- s[,-1]
    for (p in 1:np) {
      s$N[p] <- sumpart(N,parts[[p]])
      s$n[p] <- sumpart(n,parts[[p]])
    }
    s$rate <- s$n/s$N
    rate_tot <- sum(n)/sum(N)
    s$relrate <- s$rate/rate_tot
    return(s)
  }

  parse <- function(S,f){
    tokens <- unlist((stringr::str_extract_all(S, "[[:alnum:]]+")[[1]]))
    x <- list()
    x[f] <- tokens
    return(x)
  }

  find_good_names_for_mutation_categories <- function(autonames){

    z <- list()
    z$code <- c("At","Af","As","Atf","Afs","Ats","Atfs",
                "Ct","Cf","Cs","Ctf","Cfs","Cts","Ctfs",
                "ACt","ACf","ACs","ACtf","ACfs","ACts","ACtfs")
    z$name = c("A->G","A->T","A->C","A->(G/T)","A->(T/C)","A->(C/G)","A->mut",
               "C->T","C->G","C->A","C->(T/G)","C->(G/A)","C->(A/T)","C->mut",
               "N->transit","N->flip","N->skew","N->nonskew","N->transver","N->nonflip","N->mut")
    if(length(autonames)>0){
      p <- list()
      names <- list()
      for (i in 1:length(autonames)) {
        p[[i]] <- parse(autonames[[i]],c("before","at","change","after"))
        p[[i]]$muttype <- z$name[match(paste(p[[i]]$at,p[[i]]$change,sep = ""),z$code)]
        muttype_split <- unlist(strsplit(p[[i]]$muttype,split = "\\->"))
        names[[i]] <- paste(p[[i]]$before,"p*",muttype_split[1],"p",p[[i]]$after,"->",muttype_split[2],sep = "")
        names[[i]] <- gsub("ACGTp|pACGT","",names[[i]])

        names_xxxp <- regmatches(names[[i]],gregexpr("[ACGT][ACGT][ACGT]p",names[[i]]))
        if(length(unlist(names_xxxp))>0){
          names_xxxpwell <- paste("(",stringr::str_sub(names_xxxp[[1]],1,1),"/",stringr::str_sub(names_xxxp[[1]],2,2),"/",stringr::str_sub(names_xxxp[[1]],3,3),")",stringr::str_sub(names_xxxp[[1]],4,4),sep = "")
          for (j in 1:length(names_xxxp[[1]])) {
            names[[i]] <- sub("[ACGT][ACGT][ACGT]p",names_xxxpwell[j],names[[i]])
          }
        }

        names_xxp <- regmatches(names[[i]],gregexpr("[ACGT][ACGT]p",names[[i]]))
        if(length(unlist(names_xxp))>0){
          names_xxpwell <- paste("(",stringr::str_sub(names_xxp[[1]],1,1),"/",stringr::str_sub(names_xxp[[1]],2,2),")",stringr::str_sub(names_xxp[[1]],3,3),sep = "")
          for (j in 1:length(names_xxp[[1]])) {
            names[[i]] <- sub("[ACGT][ACGT]p",names_xxpwell[j],names[[i]])
          }
        }

        names_pxxx <- regmatches(names[[i]],gregexpr("p[ACGT][ACGT][ACGT]->",names[[i]]))
        if(length(unlist(names_pxxx))>0){
          names_pxxxwell <- paste("p(",stringr::str_sub(names_pxxx[[1]],2,2),"/",stringr::str_sub(names_pxxx[[1]],3,3),"/",stringr::str_sub(names_pxxx[[1]],4,4),")","->",sep = "")
          for (j in 1:length(names_pxxx[[1]])) {
            names[[i]] <- sub("p[ACGT][ACGT][ACGT]->",names_pxxxwell[j],names[[i]])
          }
        }

        names_pxx <- regmatches(names[[i]],gregexpr("p[ACGT][ACGT]->",names[[i]]))
        if(length(unlist(names_pxx))>0){
          names_pxxwell <- paste("p(",stringr::str_sub(names_pxx[[1]],2,2),"/",stringr::str_sub(names_pxx[[1]],3,3),")","->",sep = "")
          for (j in 1:length(names_pxx[[1]])) {
            names[[i]] <- sub("p[ACGT][ACGT]->",names_pxxwell[j],names[[i]])
          }
        }

        names_x <- regmatches(names[[i]],gregexpr("^\\*[ACN]->",names[[i]]))
        if(length(unlist(names_x))>0){
          names_x_well <- paste(stringr::str_sub(names_x[[1]],2,2),"->",sep = "")
          for (j in 1:length(names_x[[1]])) {
            names[[i]] <- sub("^\\*[ACN]->",names_x_well[j],names[[i]])
          }
        }

        names[[i]] <- gsub("^N->","",names[[i]])
      }
    }
    return(names)
  } # of find_good_names_for_mutation_categories <- function(autonames)

  subfprintf <- function(str,...){
    if(!is.null(P$mutcategs_report_filename)){
      # i <- ord[j]
      sink(file = as.character(P$mutcategs_report_filename),append = TRUE,split = FALSE)
      #print(sprintf(str,varargin[1:length(varargin)]))
      cat(sprintf(str,...))
      cat("\n")
      sink()
    }else{
      #print(sprintf(str,varargin[1:length(varargin)]))
      cat(sprintf(str,...))
    }
  }



  reportrule2 <- function(partsi){
    parts <- list()
    for (i in 1:length(partsi)) {
      parts[[i]] <- int_to_categ(partsi[i])
    }
    rules <- convert_parts_to_rule(parts)
    stats <- rulestats(parts)

    if(length(rules)>0){
      tmp <- list()
      for (i in 1:length(rules)) {
        tmp[[i]] <- unlist(parse(rules[[i]],c("left","from","change","right")))
      }
    }
    tmpp <- as.data.frame(matrix(numeric(0),ncol=4))
    for (i in 1:length(tmp)) {
      tmpp <- rbind(tmpp,t(as.data.frame(unlist(tmp[[i]]))))
    }
    rownames(tmpp) <- NULL

    ck <- cbind(tmpp,stats)
    ck$autoname <- rules
    ck$name <- find_good_names_for_mutation_categories(ck$autoname)
    ck$type <- matrix(rep("point",each=nrow(ck)))
    ord <- order(stats$relrate,decreasing = T)
    for (j in 1:length(parts)) {
      i <- ord[j]
      subfprintf("%s  n %.0f  N %.0f   rate  %e (%.4fx)...",ck$name[i],stats$n[i],stats$N[i],stats$rate[i],stats$relrate[i])
      #print(sprintf("%s  n %.0f  N %.0f   rate  %e (%.4fx)...\n",ck$name[i],stats$n[i],stats$N[i],stats$rate[i],stats$relrate[i]))
    }
    ck <- ck[ord,]
    out <-list(ck,stats)
    return(out)
  }  # of reportrule2 <- function(partsi)


  ###  subfunctions of find_mut_categs end   ###



  ######## begin find_mut_categs ########
  if(is.null(P$max_k)){
    P$max_k <- 5
  }
  if(is.null(P$mutcategs_report_filename)){
    P$mutcategs_report_filename <- ""
  }
  if(dim(Nn)[1] != 32){
    stop(sprintf("input must have 32 rows"))
  }
  if(dim(Nn)[2] != 5){
    stop(sprintf("input must have 5 columns (N A C G T)"))
  }

  orig_N <- matrix(data = NA,nrow = 6,ncol = 16)
  orig_N[1:3,] <- matrix(rep(Nn$N[1:16],times=3),nrow = 3,ncol = 16,byrow = T)
  orig_N[4:6,] <- matrix(rep(Nn$N[17:32],times=3),nrow = 3,ncol = 16,byrow = T)
  orig_n <- matrix(data = NA,nrow = 6,ncol = 16)
  orig_n <- rbind(t(Nn[1:16,c(3,4,5)]),t(Nn[17:32,c(2,4,5)]))

  n <- array(0,c(4,4,2,3))
  N <- array(0,c(4,4,2,3))

  for (base5 in 1:4) {
    for (base3 in 1:4) {
      for (oldbase in 1:2) {
        for (muttype in 1:3) {
          col <- 4*(base5-1)+base3
          row <- switch(oldbase,switch(muttype,2,3,1),switch(muttype,6,5,4))
          n[base5,base3,oldbase,muttype] <- orig_n[row,col]
          N[base5,base3,oldbase,muttype] <- orig_N[row,col]
        }
      }
    }
  }
  Ntot <- sum(N)
  unsplit <- list(base5=c(1,2,3,4),base3=c(1,2,3,4),oldbase=c(1,2),muttype=c(1,2,3))

  C <- list()
  Sc <- list()

  # case method = "best"
  unsplit <- categ_to_int(list(base5=c(1,2,3,4),base3=c(1,2,3,4),oldbase=c(1,2),muttype=c(1,2,3)))
  mask <- 15*16^c(0:3)
  nybs <- as.matrix(c(1:14)) %*% t(as.matrix(16^c(0:3)))
  initial <- unsplit
  leaves <- as.matrix(initial)
  for (k in 1:(P$max_k)) {
    subfprintf("k=%d",k)
    #cat(sprintf("k=%d",k))
    if(k>1){
      old.leaves <- leaves
      leaves <- matrix(data = 0,nrow = 34^(k-1),ncol = k)
      lastleaf <- 0

      for (l in 1:nrow(old.leaves)) {
        leaf <- old.leaves[l,]
        for (c in 1:(k-1)) {
          parent <- leaf[c]
          for (d in 1:4) {
            tobreak <- bitwAnd(parent,mask[d])
            #tobreak <- strtoi(paste(as.integer(rev(as.integer(intToBits(parent))) & rev(as.integer(intToBits(mask[d])))),collapse =""),base = 2)
            #parent mask[d] 进行二进制逻辑与，再转成十进制??? 赋值给tobreak
            tokeep <- parent - tobreak
            frags1 <- bitwAnd(tobreak,nybs[,d])
            frags2 <- tobreak - frags1
            frags <- cbind(frags1,frags2)
            frags <- frags[which(!(frags[,1] == 0 | frags[,2] == 0)),]
            children <- tokeep + frags
            nc <- dim(children)[1]

            if(nc==0){
              next
            }else if(lastleaf+nc <= nrow(leaves)){
              leaves[(lastleaf+1):(lastleaf+nc),1:(k-1)] <- matrix(rep(leaf,each=nc),nrow = nc,ncol = k-1)
              leaves[(lastleaf+1):(lastleaf+nc),c(c,k)] <- children
            }else{
              temp <- matrix(data = NA,nrow = nc,ncol = k)
              temp[1:nc,1:(k-1)] <- matrix(rep(leaf,each=nc),nrow = nc,ncol = k-1)
              temp[1:nc,c(c,k)] <- children
              leaves <- rbind(leaves[1:lastleaf,],temp)
            }

            lastleaf <- lastleaf + nc
          } #of for (d in 1:4)
        }
      }
      leaves_orderrow <- unique(t(apply(leaves[1:lastleaf,],1,sort,decreasing=F)))  #把leaves的每行按大小排序,并去???
      leaves <- leaves_orderrow[order(leaves_orderrow[,1]),]
    }  # of k>1

    best_l <- 0
    best_H <- Inf
    if(k<4){
      for (l in 1:nrow(leaves)) {
        H <- entropy_by_parts2(leaves[l,])
        if(H<best_H){
          best_l <- l
          best_H <- H
        }
      }
    }else{
      if(k==4){
        mx <- (2^16)-1
        lookup <- matrix(data = NaN,nrow = mx,1)
        for (i in 1:mx) {
          binary_rev_i <- as.integer(intToBits(i))[1:16]
          if(any(binary_rev_i[c(11,12,16)])){
            next
          }
          #cat(sprintf("%d",i))
          c <- int_to_categ(i)
          ns <- sumpart(n,c)
          Ns <- sumpart(N,c)
          f <- Ns/Ntot
          H_part <- entropy(ns/Ns)
          lookup[i] <- f*H_part
        }
      }  # of if(k==4)
      for (l in 1:nrow(leaves)) {
        H <- 0
        for (i in 1:k) {
          H <- H + lookup[leaves[l,i]]
        }
        if(H<best_H){
          best_l <- l
          best_H <- H
        }
      } # of for (l in 1:nrow(leaves))
    } # of if(k<4) else{

    best_leaf <- leaves[best_l,]
    reportrule2_out <- reportrule2(best_leaf)
    C[[k]] <- reportrule2_out[[1]]
    stats <- reportrule2_out[[2]]
    Sc[[k]] <- unlist(list(H=best_H))
  }  # of k in 1:(P$max_k)

  return(C)

}  #of preprocessFindMutCateg

############### preprocess_find_mut_categs function  end ##############

# collapse coverage to context65
#' Generating 65 mutation categories
#'
#' Generating 65 types of mutaion categories for subsequently classification and analysis.
#' @details Considering the contexts for each of the 4 bases(A,C,G,T), i.e. the total number of categories reaches
#' to 65(64 types of trinucleotides and a supplement called "any N").
preprocessGenerateCategContext65 <- function()
{
  x <- c("A","C","G","T")
  x_in <- rep(x,each=16)
  x_left <- rep(rep(x,each=4),times=4)
  x_right <- rep(x,times=16)
  y <- paste(x_in,"in",x_left,sep = " ")
  y <- paste(y,"_",x_right,sep = "")
  y[65] <- "any N"
  return(y)
}

#preprocessing covariate files
# fill missing data in covariate file
covariate_preprocessing <- function(M,V,k=100){
  cluster=expr=reptime=hic=NULL
  # table gene effects(noncoding, nonsilent, silent, null) for each gene
  Mtable <- tapply(M$effect,M$gene,table)
  gene <- sort(as.character(unique(M$gene)))
  ng <- length(gene)
  effect <- unique(M$effect)
  ne <- length(effect)
  geneEffect <- as.data.frame(matrix(0,nrow = ng,ncol = ne,dimnames = list(gene,effect)))

  for (g in 1:ng) {
    flag <- match(names(Mtable[[gene[g]]]),effect)
    # find effect names in Mtable
    geneEffect[g,flag] <- Mtable[[gene[g]]]
    # fill geneEffect with values in mutation M table
  }

  km_result <- stats::kmeans(geneEffect, k, nstart = 100,iter.max = 2000)
  # print(km_result)
  geneEffect<- cbind(geneEffect, cluster = km_result$cluster)

  geneEffect2<-data.frame(geneEffect,row.names = NULL)


  geneEffect2$gene<-gene

  means <- data.frame(cluster=1:k)
  # MVinteraction <- data_frame(gene=interaction(geneeffect2$gene,V$gene))
  means$expr <- 0
  means$reptime <- 0
  means$hic <- 0
  for(i in 1:k){
    select1 <- subset(geneEffect2,cluster==i)
    V1 <- subset(V,!is.nan(expr))
    flag <- match(select1$gene,V1$gene)
    flag<-flag[stats::complete.cases(flag)]
    means[i,2] <- mean(V1[flag,2])

    V1 <- subset(V,!is.nan(reptime))
    flag <- match(select1$gene,V1$gene)
    flag<-flag[stats::complete.cases(flag)]
    means[i,3] <- mean(V1[flag,3])

    V1 <- subset(V,!is.nan(hic))
    flag <- match(select1$gene,V1$gene)
    flag<-flag[stats::complete.cases(flag)]
    means[i,4] <- mean(V1[flag,4])
  }

  Vnan <- subset(V,is.nan(expr))
  flag1 <- match(Vnan$gene,V$gene)
  flag2 <- match(Vnan$gene,geneEffect2$gene)
  V[flag1,2]<- means[geneEffect2$cluster[flag2],2]

  Vnan <- subset(V,is.nan(reptime))
  flag1 <- match(Vnan$gene,V$gene)
  flag2 <- match(Vnan$gene,geneEffect2$gene)
  V[flag1,3]<- means[geneEffect2$cluster[flag2],3]


  Vnan <- subset(V,is.nan(hic))
  flag1 <- match(Vnan$gene,V$gene)
  flag2 <- match(Vnan$gene,geneEffect2$gene)
  V[flag1,4]<- means[geneEffect2$cluster[flag2],4]

  V<-subset(V,!is.na(expr),!is.na(reptime))
  return(V)

}

