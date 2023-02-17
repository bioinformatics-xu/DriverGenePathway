#' Search driver pathway using de novo method
#'
#' This function searches driver pathway using de novo method based on mutual exclusivity and coverage
#' It outputs the driver pathway and its p-value as a txt file.
#'
#' @param mutation_matrix 0-1 mutation matrix where rows represent patients, columns represent genes.
#' @param driver_size The size of identified driver gene set, defaulted to 3.
#' @param pop_size The population size of GA, defaulted to 200.
#' @param iters The iteration time of GA, defaulted to 500.
#' @param permut_time The times of permutation test, defaulted to 1000.
#' @param outfile The output file of driver gene set.
#' @author Xiaolu Xu <lu.xu@@lnnu.edu.cn>
#' @examples
#' data(SampleMutationMatrix)
#' DriverPathway(mutation_matrix, driver_size=3, pop_size=20, iters=50, permut_time=100)
#' @export
DriverPathway <- function(mutation_matrix,
                       driver_size=3,
                       pop_size=200,
                       iters=500,
                       permut_time=1000,outfile="denovoDriverPathway.txt") {

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

    mulExclusive_significance <- function(mutation_matrix,driver_geneset,permut_time=1000){
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



    obj_weight <- vector()

    # start with an random population
    total_gene_size <- ncol(mutation_matrix)
    population = matrix(nrow=(pop_size*2), ncol=total_gene_size);
    # fill values
    for (child in 1:(pop_size*2)) {
        population[child,] <- rep(0,total_gene_size)
        population[child,sample(1:total_gene_size,driver_size)] <- 1
    }

    # calculate each object
    evalVals = rep(NA, pop_size*2);
    for (object in 1:(pop_size*2)) {
        evalVals[object] = evalFunc(mutation_matrix,population[object,]);
    }

    order_flag <- order(evalVals,decreasing = F)
    population <- population[order_flag,]
    Record <- 0
    iter <- 1
    pj <- 1
    while( (iter < iters) & (Record < 10) ){
        weight_vector <- evalVals[order_flag]
        t_minweight <- min(weight_vector)
        for (pj in 1:pop_size) {
            selected_index <- select_operator(weight_vector)
            population[pop_size+pj,] <- cross_operator(population[selected_index[1],], population[selected_index[2],],total_gene_size)
            population[pop_size+pj,] <- mutation_operator(mutation_matrix,
                                                          population[pop_size+pj,],total_gene_size,1)
            weight_vector[pop_size+pj] <- evalFunc(mutation_matrix,
                                                   population[pop_size+pj,])
        }
        population <- population[order(weight_vector,decreasing = F),]
        weight_vector <- sort(weight_vector,decreasing = F)
        obj_weight[iter] <- weight_vector[1]

        if(Record == 2){
            temp <- sort(unique(weight_vector),decreasing = F)
            index <- which(weight_vector <= temp[1])
            for (j in min(length(index),1)) {
                population[index[j],] <- mutation_operator(mutation_matrix,
                                                           population[index[j],],
                                                           total_gene_size,
                                                           floor(sqrt(total_gene_size)))
                weight_vector[index[j]] <- evalFunc(mutation_matrix,
                                                    population[index[j],])
            }

        }

        if(Record == 5){
            temp <- sort(unique(weight_vector),decreasing = F)
            index <- which(weight_vector <= temp[1])
            for (j in min(length(index),1)) {
                population[index[j],] <- mutation_operator(mutation_matrix,
                                                           population[index[j],],
                                                           total_gene_size,
                                                           total_gene_size)
                weight_vector[index[j]] <- evalFunc(mutation_matrix,
                                                    population[index[j],])
            }

        }


        min_weight <- min(weight_vector)

        if(min_weight == t_minweight){
            Record <- Record + 1
        }else{
            Record <- 0
        }
        iter <- iter + 1

    }

    order_flag <- order(weight_vector)
    population <- population[order_flag,]
    max_pop <- population[1:pop_size,]

    # delete the duplicate rows
    len_maxpop <- dim(max_pop)[1]
    index_del <- rep(0,times=len_maxpop)
    max_pop <- max_pop[!duplicated(max_pop),]

    if(length(nrow(max_pop) > 1) > 0){
        driver_geneset <- matrix(data = NA, ncol = driver_size, nrow = nrow(max_pop))
        for (i in 1:nrow(max_pop)) {
            driver_geneset[i,] <- colnames(mutation_matrix)[which(max_pop[i,] == 1)]
        }
    }else{
        driver_geneset <- colnames(mutation_matrix)[which(max_pop==1)]

    }

    p_value <- mulExclusive_significance(mutation_matrix,driver_geneset,permut_time = permut_time)
    G <- list(driver_geneset=driver_geneset,p_value=p_value)
    utils::write.table(G,file = outfile,quote = F,sep = "\t",row.names = F)
    return(G)
}
