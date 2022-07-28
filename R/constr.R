#' Constraints
#' @description Impose constraints on chromosomes.
#' @param population numeric matrix which rows are chromosomes i.e. vectors of 
#' parameters values.
#' @param method method used to impose constraints.
#' @param par additional parameters to be passed depending on the \code{method}.
#' @param iter iteration number of the genetic algorithm.
#' @details If \code{method = "bounds"} then chromosomes will be bounded
#' between \code{par$lower} and \code{par$upper}.
#' @return The function returns a matrix which rows are chromosomes satisfying
#' the constraints.
#' @examples 
#' # Randomly initialize population
#' set.seed(123)
#' population <- gena.population(pop.n = 10,
#'                               lower = c(-5, -5), 
#'                               upper = c(5, 5))
#'                            
#' # Impose lower and upper bounds constraints
#' pop.constr <- gena.constr(population, 
#'                           method = "bounds",
#'                           par = list(lower = c(-1, 2),
#'                                      upper = c(1, 5)))
#' print(pop.constr)
#' 
gena.constr <- function(population,
                        method = "bounds",
                        par,
                        iter)              
{
  if (method == "bounds")
  {
      genes.n <- ncol(population)
      for (i in 1:genes.n)
      {
        population[population[, i] < lower[i], i] <- par$lower[i]
        population[population[, i] > upper[i], i] <- par$upper[i]
      }
  }
  
  return(population)
}

# Assign default parameters for
# constraint algorithm depending
# on the "method"
gena.constr.validate <- function(method, par)
{
  # Validate the "method"
  
  methods <- c("bounds")                                   # the list of all
                                                           # available methods
  
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect constr.method argument. ",      # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
}