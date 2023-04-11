#' Identify driver pathway using a de novo method, AWRMP(https://www.frontiersin.org/articles/10.3389/fgene.2019.00233/full).
#'
#' @param mutation_data 0-1 mutation matrix where rows represent patients, columns represent genes or MAF file.
#' If the mutation data is a MAF file, then the preprocessing procedure will be performed using coverage data, covariate data, mutation dictionary, and chr files directory, which can be specified by the user.
#' The user should name these files "coverage.txt.$" for coverage data, "covariates.txt$" for covariate data, "dictionary file.txt$" for mutation dictionary, and "^chr files hg" for "chr files directory".
#' If the specified named files can be found, they will be read and utilized to generate mutation matrix.
#' If the files cannot be detected, the function will automatically download them from the MutSigCV website, which may take some time.
#' @param driver_size The size of identified driver gene set, defaulted to 3.
#' @param pop_size The population size of GA which should be adjusted according the number of samples, defaulted to 200.
#' @param iters The iteration time of GA, defaulted to 1000.
#' @param permut_time The times of permutation test, defaulted to 1000.
#' @param process_bmr The background mutation rate to preprocess MAF file to mutation matrix using binomial hypothesis testing.
#' @return A list including the mutation matrix, identified driver gene set, and p-value.
#' @author Xiaolu Xu <lu.xu@@lnnu.edu.cn>
#' @examples
#' data(SampleMutationMatrix)
#' driverSet <- DriverPathway(mutation_matrix, driver_size=3, pop_size=20, iters=50, permut_time=100)
#' @export
DriverPathway <- function(mutation_data,
                          driver_size=3,
                          pop_size=200,
                          iters=1000,
                          permut_time=1000,process_bmr=1.2e-6) {
  obj_weight <- vector()
  mutation_matrix <- preprocessing_mutation_data(mutation_data,process_bmr)
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
  G <- list(mutation_matrix=mutation_matrix,driver_geneset=driver_geneset,p_value=p_value)
  message("Driver pathway identification finished")
  #utils::write.table(G,file = outfile,quote = F,sep = "\t",row.names = F)
  return(G)
}
