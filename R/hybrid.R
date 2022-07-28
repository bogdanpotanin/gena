#' Hybridization
#' @description Hybridization method (algorithm) to be used in the
#' genetic algorithm.
#' @param population numeric matrix which rows are chromosomes i.e. vectors of 
#' parameters values.
#' @param fitness numeric vector which \code{i}-th element is the value of 
#' \code{fn} at point \code{population[i, ]}.
#' @param hybrid.n positive integer representing the number of hybrids.
#' @param method hybridization method to improve chromosomes via local search.
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param opt.par parameters of the local optimization function
#' to be used for hybridization algorithm.
#' @param info logical; if \code{TRUE} then some optimization related 
#' information will be printed each iteration.
#' @details This function uses \code{\link[gena]{gena.mating}} function to 
#' select hybrids. Therefore \code{method} and \code{par} arguments will
#' be passed to this function.
#' 
#' Currently \code{\link[stats]{optim}} is the only available local
#' optimizer. Therefore \code{opt.par} is a list containing parameters
#' that should be passed to \code{\link[stats]{optim}}.
#' 
#' @return The function returns a list with the following elements:
#' \itemize{
#' \item \code{population} - matrix which rows are chromosomes including hybrids.
#' \item \code{fitness} - vector which i-th element is the fitness of the
#' i-th chromosome.
#' }
#' 
#' @examples 
#' # Consider the following fitness function
#' fn <- function(x)
#' {
#'   val <- x[1] * x[2] - x[1] ^ 2 - x[2] ^ 2
#' }
#' 
#' # Also let's provide it's gradient (optional)
#' gr <- function(x)
#' {
#'   val <- c(x[2] - 2 * x[1],
#'            x[1] - 2 * x[2])
#' }
#' 
#' # Randomly initialize the population
#' set.seed(123)
#' n_population <- 10
#' population <- gena.population(pop.n = n_population,
#'                               lower = c(-5, -5), 
#'                               upper = c(5, 5))
#'
#' # Calculate fitness of each chromosome
#' fitness <- rep(NA, n_population)
#' for(i in 1:n_population)
#' {
#'   fitness[i] <- fn(population[i, ])
#' }
#' 
#' # Perform hybridization
#' hybrids <- gena.hybrid(population = population,
#'                        fitness = fitness,
#'                        opt.par = list(fn = fn,
#'                                       gr = gr,
#'                                       method = "BFGS",
#'                                       control = list(fnscale = -1,
#'                                                      abstol = 1e-10,
#'                                                      reltol = 1e-10,
#'                                                      maxit = 1000)),
#'                        hybrid.n = 2,
#'                        method = "rank",
#'                        par = 0.8)
#' print(hybrids)
#' 

# Make hybrids
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
  genes.n <- ncol(population)                              # number of genes
  
  dots <- list(...)                                        # fn specific parameters
  all.args <- append(opt.par, dots)                        # all parameters to be
                                                           # passed to optim
  hybrids <- matrix(NA, nrow = n_pop, ncol = genes.n)      # matrix to store hybrids
  hybrids_fitness <- rep(NA, n_pop)                        # vector to store fitnesses 
                                                           # of hybrids
  
  # To select hybrids use the same function
  # as for the parents
  
  hybrids_list <- gena.mating(population = population,     # select the chromosomes
                              fitness = fitness,           # to become a hybrids
                              parents.n = hybrid.n,
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