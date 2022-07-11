# Select the parents
gena.hybrid <- function(population,
                        fitness,
                        hybrid.n = 1,
                        method,
                        par,
                        opt.par,
                        info = FALSE,
                        ...)          
{      
  # ---
  # The algorithm brief description:
  # 1. Select the chromosomes to become the hybrids
  # 2. Perform local optimization on chromosomes
  #    to make them the hybrids
  # ---
  
  # Prepare some values
  
  n_pop <- nrow(population)                                # number of chromosomes
  n_genes <- ncol(population)                              # number of genes
  
  dots <- list(...)                                        # fn specific parameters
  all.args <- append(opt.par, dots)                        # all parameters to be
                                                           # passed to optim
  hybrids <- matrix(NA, nrow = n_pop, ncol = n_genes)      # matrix to store hybrids
  hybrids_fitness <- rep(NA, n_pop)                        # matrix to store fitness of hybrids
  
  # To select hybrids use the same function
  # as for a parents
  
  hybrids_list <- gena.mating(population = population,     # select the chromosomes
                              fitness = fitness,           # to become a hybrids
                              n_parents = hybrid.n,
                              method = method,
                              par = par,
                              iter = iter,
                              self = TRUE)
  hybrids_ind <- unique(hybrids_list$ind)                  # indexes of chromosomes that
                                                           # will become a hybrids
  # Perform the optimization
  
  counter <- 1
  for (i in hybrids_ind)
  {
    hybrid_message <- paste0("hybrid ", counter," fitness: before = ",
                             fitness[i],
                             ", after = ")
    
    all.args$par <- population[i, ]
    opt.result <- NULL
    tryCatch(
    {
      opt.result <- do.call(optim, args = all.args)
      if (fitness[i] < opt.result$value)
      {
        population[i, ] <- opt.result$par
        fitness[i] <- opt.result$value
      }
    })
    hybrid_message <- paste0(hybrid_message, fitness[i], "\n")
    if (info)
    {
      cat(hybrid_message)
    }
    counter <- counter + 1
  }
  
  return_list <- list(population = population,
                      fitness = fitness)
  
  return(return_list)
}

# References
# Best BHGA: https://www.hindawi.com/journals/jam/2013/103591/