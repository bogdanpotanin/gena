#' Genetic Algorithm
#' @description This function allows to use genetic algorithm for
#' numeric global optimization of real-valued functions.
#' @param fn function to be maximized.
#' @param gr gradient of the \code{fn}.
#' @param lower lower bound of the search space.
#' @param upper upper bound of the search space.
#' @param pop.initial numeric matrix which rows are chromosomes to be
#' included into the initial population.
#' @param pop.n integer representing the size of the population.
#' @param pop.method the algorithm to be applied for a criation of 
#' the initial population.
#' @param mating.method the algorithm to be applied for a mating i.e. selection
#' of parents.
#' @param mating.par parameters of the mating algorithm.
#' @param mating.self logical; if \code{TRUE} then the chromosome may
#' mate with itself i.e. both parents may be the same chromosome.
#' @param crossover.method an algorithm to be applied for crossover i.e.
#' creation of the children.
#' @param crossover.par parameters of the crossover algorithm.
#' @param crossover.prob probability of the crossover for each pair of parents.
#' @param mutation.method algorithm to be applied for mutation i.e. random
#' change in some genes of the children.
#' @param mutation.par parameters of the mutation algorithm.
#' @param mutation.prob mutation probability for the chromosomes.
#' @param mutation.genes.prob mutation probability for the genes.
#' @param elite.n number of the elite children i.e. those which have the
#' highest function value and will be preserved for the next population.
#' @param elite.duplicates logical; if \code{TRUE} then some elite children
#' may have the same genes.
#' @param hybrid.method hybrids selection algorithm i.e. mechanism 
#' determining which chromosomes should be subject to local optimization.
#' @param hybrid.par parameters of the hybridization algorithm.
#' @param hybrid.prob probability of generating the hybrids each iteration.
#' @param hybrid.opt.par parameters of the local optimization function
#' to be used for hybridization algorithm.
#' @param hybrid.n number of hybrids that appear if hybridization
#' should take place during the iteration.
#' @param maxiter maximum number of iterations of the algorithm.
#' @param constr.method the algorithm to be applied for imposing constraints
#' on the chromosomes.
#' @param constr.par parameters of the constraint algorithm.
#' @param info logical; if \code{TRUE} then some optimization related 
#' information will be printed each iteration.
#' @param ... additional parameters to be passed to 
#' \code{fn} and \code{gr} functions.
#' @details To find information on particular methods available via
#' \code{mating.method}, \code{crossover.method}, \code{mutation.method},
#' \code{hybrid.method} and \code{constr.method}
#' arguments please see 'Details' section of 
#' \code{\link[gena]{gena.population}}, \code{\link[gena]{gena.crossover}},
#' \code{\link[gena]{gena.mutation}}, \code{\link[gena]{gena.hybrid}} 
#' and \code{\link[gena]{gena.constr}} correspondingly.
#' 
#' It is possible to provide manually implemented functions for 
#' population initialization, mating, crossover, mutation and hybridization. 
#' For example manual mutation function may be provided through
#' \code{mutation.method} argument. It should have the same signature 
#' (arguments) as \code{\link[gena]{gena.mutation}} function and return 
#' the same object i.e. the matrix of chromosomes of the appropriate size.
#' Manually implemented functions for other operators (crossover, mating
#' and so on) may be provided in a similar way.
#' 
#' By default function does not impose any constraints upon the parameters.
#' If \code{constr.method = "bounds"} then \code{lower} and \code{upper} 
#' constraints will be imposed. Lower bounds should be strictly smaller
#' then upper bounds.
#' 
#' Currently the only available termination condition is \code{maxiter}. We
#' are going to provide some additional termination conditions during
#' future updates.
#' 
#' Infinite values in \code{lower} and \code{upper} are substituted with
#' \code{.Machine$double.xmin * 0.9} and \code{.Machine$double.xmax * 0.9}
#' correspondingly.
#' 
#' By default if \code{gr} is provided then BFGS algorithm will be used inside
#' \code{\link[stats]{optim}} during hybridization.
#' Otherwise \code{Nelder-Mead} will be used. 
#' Manual values for \code{\link[stats]{optim}} arguments may be provided 
#' through \code{hybrid.opt.par} argument.
#' 
#' Arguments \code{pop.n} and \code{elite.n} should be even integers and
#' \code{elite.n} should be greater then 2. If these arguments are odd integers 
#' then they will be coerced to even integers by adding 1. 
#' Also \code{pop.n} should be greater then \code{elite.n} at least by 2.
#' 
#' @return This function returns an object of class \code{gena} that is a list
#' containing the following elements:
#' \itemize{
#' \item \code{par} - chromosome (solution) with the highest fitness
#' (objective function) value.
#' \item \code{value} - value of \code{fn} at \code{par}.
#' \item \code{chromosomes} - matrix of chromosomes (solutions) of the 
#' last iteration of the algorithm.
#' }
#' @examples
#' ## Consider Ackley function
#' \donttest{
#' fn <- function(par, a = 20, b = 0.2)
#' {
#'   val <- a * exp(-b * sqrt(0.5 * (par[1] ^ 2 + par[2] ^ 2))) +
#'          exp(0.5 * (cos(2 * pi * par[1]) + cos(2 * pi * par[2]))) -
#'          exp(1) - a
#'   return(val)
#' }
#' 
#' # Maximize this function using classical
#' # genetic algorithm setup
#' set.seed(123)
#' lower <- c(-5, -100)
#' upper <- c(100, 5)
#' opt <- gena(fn = fn, 
#'             lower = lower, upper = upper,
#'             hybrid.prob = 0,
#'             a = 20, b = 0.2)
#' print(opt$par)
#' 
#' # Replicate optimization using hybridization
#' opt <- gena(fn = fn, 
#'             lower = lower, upper = upper,
#'             hybrid.prob = 0.2,
#'             a = 20, b = 0.2)
#' print(opt$par)
#' }
#' 
#' ## Consider Rosenbrock function
#' fn <- function(par, a = 100)
#' {
#'   val <- -(a * (par[2] - par[1] ^ 2) ^ 2 + (1 - par[1]) ^ 2 +
#'            a * (par[3] - par[2] ^ 2) ^ 2 + (1 - par[2]) ^ 2)
#'   return(val)
#' }
#' 
#' # Apply genetic algorithm
#' lower <- rep(-10, 3)
#' upper <- rep(10, 3)
#' set.seed(123)
#' opt <- gena(fn = fn,
#'             lower = lower, upper = upper,
#'             hybrid.prob = 0,
#'             a = 100,)
#' print(opt$par)
#' 
#' \donttest{
#' # Improve the results by hybridization
#' opt <- gena(fn = fn,
#'             lower = lower, upper = upper,
#'             hybrid.prob = 0.2,
#'             a = 100)
#' print(opt$par)
#' }
#' 
#' # Provide manually implemented mutation function
#' # which simply randomly sorts genes. 
#' # Note that this function should have the same 
#' # arguments as gena.mutation.
#' mutation.my <- function(children, lower, upper, 
#'                         prob, prob.genes, 
#'                         method, par, iter)
#' {
#'   # Get dimensional data
#'   children.n <- nrow(children)
#'   genes.n <- ncol(children)
#'   
#'   # Select chromosomes that should mutate
#'   random_values <- runif(children.n, 0, 1)
#'   mutation_ind <- which(random_values <= prob)
#'    
#'   # Mutate chromosomes by randomly sorting
#'   # their genes
#'   for (i in mutation_ind)
#'   {
#'     children[i, ] <- children[i, sample(1:genes.n)]
#'   }
#'   
#'   # Return mutated chromosomes
#'   return(children)
#' }
#' 
#' opt <- gena(fn = fn,
#'             lower = lower, upper = upper,
#'             hybrid.prob = 0,
#'             mutation.method = mutation.my,
#'             a = 100)
#' print(opt$par)
#' 
gena <- function(
  fn,
  gr = NULL,
  lower,
  upper,
  pop.n = 100,
  pop.initial = NULL,
  pop.method = "uniform",
  mating.method = "rank",
  mating.par = NULL,
  mating.self = FALSE,
  crossover.method = "local",
  crossover.par = NULL,
  crossover.prob = 0.8,
  mutation.method = "constant",
  mutation.par = NULL,
  mutation.prob = 0.2,
  mutation.genes.prob = 1 / length(lower),
  elite.n = min(10, 2 * round(pop.n / 20)),
  elite.duplicates = FALSE,
  hybrid.method = "rank",
  hybrid.par = 2,
  hybrid.prob = 0.1,
  hybrid.opt.par = NULL,
  hybrid.n = 1,
  maxiter = 100,
  constr.method = NULL,
  constr.par = NULL,
  info = TRUE,
  ...)
{
  # Initialize some variables
  
  pop.n <- ceiling(pop.n)
  elite.n <- ceiling(elite.n)
  genes.n <- length(lower)
  
  # Perform initial validation
  
    # fn
  
  if (!is.function(fn))
  {
    stop("Please, insure that 'fn' is a function\n")
  }
  
    # gr
  if(!is.null(gr))
  {
    if (!is.function(gr))
    {
      stop("Please, insure that 'gr' is a function\n")
    }
  }
  
    # lower and upper
  
  if (length(lower) != length(upper))
  {
    stop("Please, insure that 'lower' and 'upper' are of the same length\n")
  }
  
  if (any(lower > upper))
  {
    stop(paste0("Please, insure that all elements of the 'lower' are ",
                "smaller than the corresponding elements of the 'upper'",
                "\n"))
  }
  
  lower[is.infinite(lower)] <- .Machine$double.xmin * 0.9
  upper[is.infinite(upper)] <- .Machine$double.xmin * 0.9
  
    # pop.initial
  
  if (!is.null(pop.initial))
  {
    if (!is.numeric(pop.initial))
    {
      stop("Please, insure that 'pop.initial' is a numeric matrix\n")
    }
    
    if (!is.matrix(pop.initial))
    {
      pop.initial <- matrix(pop.initial, nrow = 1)
    }

    if (ncol(pop.initial) != genes.n)
    {
      stop(paste0("Please, insure that the number of columns of 'pop.initial' ",
                  "equals to the number of estimated parameters i.e. genes",
                  "\n"))
    }
  }
  
    # pop.n
  
  if (pop.n <= 1)
  {
    stop("Please, insure that 'pop.n' is an integer greater than one")
  }
  
  if ((pop.n %% 2) != 0) 
  {
    pop.n <- pop.n + 1
    warning("Since 'pop.n' is not even it has been incremented by one\n")
  }
  
    # elite.n
  
  if (elite.n <= 1)
  {
    stop("Please, insure that 'elite.n' is an integer greater than one\n")
  }
  
  if ((elite.n %% 2) != 0) 
  {
    elite.n <- elite.n + 1
    warning("Since 'elite.n' is not even it has been incremented by one\n")
  }
    
  # Initialize some additional variables  
  
  parents.n <- pop.n - elite.n
  fitness <- rep(NA, pop.n)
  
  # Continue validation
  
    # method.mating and mating.par
  if (!is.function(mating.method))
  {
    mating.par <- gena.mating.validate(method = mating.method,
                                       par = mating.par,
                                       parents.n = parents.n)
  }
  else
  {
    gena.mating <- mating.method
  }
  
    # mating.self
  
  if (!is.logical(mating.self))
  {
    if ((mating.self == 0) | (mating.self == 1))
    {
      mating.self <- as.logical(mating.self)
    } else {
      stop("Please, insure that 'mating.self' is a logical\n")
    }
  }
  
    # crossover.method and crossover.par
  
  if (!is.function(crossover.method))
  {
    crossover.par <- gena.crossover.validate(method = crossover.method,
                                             par = crossover.par)
  }
  else
  {
   gena.crossover <- crossover.method
  }
  
    # crossover.prob  
  
  if (is.numeric(crossover.prob))
  {
    if ((crossover.prob < 0) | (crossover.prob > 1))
    {
      stop("Please, insure that 'crossover.prob' is between 0 and 1\n")
    }
  } else {
    stop("Please, insure that 'crossover.prob' is a numeric value\n")
  }
  
    # mutation.method and mutation.par
  
  if (!is.function(mutation.method))
  {
    mutation.par <- gena.mutation.validate(method = mutation.method,
                                           par = mutation.par,
                                           genes.n = genes.n)
  }
  else
  {
    gena.mutation <- mutation.method
  }
  
    # mutation.prob
  
  if (is.numeric(mutation.prob))
  {
    if ((mutation.prob < 0) | (mutation.prob > 1))
    {
      stop("Please, insure that 'mutation.prob' is between 0 and 1\n")
    }
  } else {
    stop("Please, insure that 'mutation.prob' is a numeric value\n")
  }
  
    # mutation.genes.prob
  
  mutation.genes.prob_n <- length(mutation.genes.prob)
  
  if (mutation.genes.prob_n == 1)
  {
    mutation.genes.prob <- rep(mutation.genes.prob, genes.n)
  }
  
  if (mutation.genes.prob_n > genes.n)
  {
    stop(paste0("Please, insure that the length of 'mutation.genes.prob_n' ",
                "equals 1 or the same as the length of 'lower'",
                "\n"))
  }
  
  if (any(mutation.genes.prob < 0) | any(mutation.genes.prob > 1))
  {
    stop(paste0("Please, insure that values in 'mutation.genes.prob' ",
                "are between 0 and 1'",
                "\n"))
  }
  
  if ((pop.n - elite.n) < 2)
  {
    stop(paste0("Please, insure that ((pop.n - elite.n) >= 2) i.e. ",
                "it seems that 'elite.n' is too large relative to 'pop.n'",
                "\n"))
  }
  
    # elite.duplicates
  
  if (!is.logical(elite.duplicates))
  {
    if ((elite.duplicates == 0) | (elite.duplicates == 1))
    {
      elite.duplicates <- as.logical(elite.duplicates)
    } else {
      stop("Please, insure that 'elite.duplicates' is a logical\n")
    }
  }
  
    # hybrid.method and hybrid.par
  
  if (!is.function(hybrid.method))
  {
    hybrid.par <- gena.mating.validate(method = hybrid.method,
                                       par = hybrid.par,
                                       parents.n = parents.n)
  }
  else
  {
    gena.hybrid <- hybrid.method
  }
  
    # hybrid.prob
  
  if (is.numeric(hybrid.prob))
  {
    if ((hybrid.prob < 0) | (hybrid.prob > 1))
    {
      stop("Please, insure that 'hybrid.prob' is between 0 and 1\n")
    }
  } else {
    stop("Please, insure that 'hybrid.prob' is a numeric value\n")
  }
  
    # hybrid.n
  
  if (hybrid.n < 0)
  {
    stop("Please, insure that 'hybrid.n' is a non-negative integer\n")
  }
  
  if (hybrid.n > pop.n)
  {
    stop("Please, insure that 'hybrid.n' is not greater than 'pop.n'")
  }
  
    # maxiter
  
  if (maxiter < 1)
  {
    stop("Please, insure that 'maxiter' is a positive integer")
  }
  
    # constr.method and constr.par
  
  if(!is.null(constr.method))
  {
    if (!(is.function(constr.method) | is.character(constr.method)))
    {
      stop("Please, insure that 'constr.method' is a function or character\n")
    }
    if (is.character(constr.method))
    {
      if ((constr.method == "bounds") & is.null(constr.par))
      {
        par$lower <- lower
        par$upper <- upper
      }
      gena.constr.validate(method = constr.method, par = constr.par)
    }
    if (is.function(constr.method))
    {
      gena.constr <- constr.method
    }
  }
  
    # info
  
  if (!is.logical(info))
  {
    if ((info == 0) | (info == 1))
    {
      info <- as.logical(info)
    } else {
      stop("Please, insure that 'info' is a logical\n")
    }
  }
  
  # Deal with hybrid optimizer default parameters
  
  hybrid.opt.par.default = list(fn = fn,
                                gr = gr,
                                method = ifelse(is.null(gr), 
                                                "Nelder-Mead", 
                                                "BFGS"),
                                control = list(fnscale = -1,
                                               abstol = 1e-10,
                                               reltol = 1e-10,
                                               maxit = 1000000000))
  
  if (!is.null(hybrid.opt.par))
  {
    if (!is.null(hybrid.opt.par$control))
    {
      hybrid.opt.par.default$control <- modifyList(hybrid.opt.par.default$control,
                                                   hybrid.opt.par$control)
    }
    hybrid.opt.par = modifyList(hybrid.opt.par.default, 
                                hybrid.opt.par)
  } else {
    hybrid.opt.par <- hybrid.opt.par.default
  }

  # Initialize the population
  
  if (is.function(pop.method))
  {
    gena.population <- pop.method
  }
  
  population <- gena.population(pop.n = pop.n,
                                lower = lower, 
                                upper = upper,
                                method = pop.method,
                                pop.initial = pop.initial)
  
  # -------------------------------------
  # Perform the genetic algorithm routine
  # -------------------------------------
  
  for(i in 1:maxiter)
  {
    # Get the parameters to be passed to fitness
    
    dots <- list(...)
    
    # Impose the constraints on the parameters
    
    if (!is.null(constr.method))
    {
      if (is.function(constr.method))
      {
        constr.par$population <- population
        constr.par$iter <- iter
        population <- do.call(what = constr.method, args = constr.par)
      }
      else
      {
        population <- gena.constr(population = population, 
                                  iter = iter,
                                  method = constr.method,
                                  par = constr.par)
      }
    }
    
    # Estimate the fitness of each chromosome

    for(j in 1:pop.n)
    {
      dots$par <- population[j, ] 
      fitness[j] <- do.call(what = fn, args = dots)
    }
    fitness[is.na(fitness) | is.nan(fitness)] <- -Inf
    
    dots$par <- NULL
    
    # Make the hybrids
    
    hybrid_value <- runif(1, 0, 1)
    
    if (((hybrid_value <= hybrid.prob) |
         (i == 1) | (i == maxiter)) & 
        (hybrid.n >= 1) & (hybrid.prob > 0))
    {
      hybrids_list <- gena.hybrid(population = population,
                                  fitness = fitness,
                                  hybrid.n = hybrid.n,
                                  method = hybrid.method,
                                  par = hybrid.par,
                                  opt.par = hybrid.opt.par,
                                  info = info,
                                  ... = ...)
      
      population <- hybrids_list$population
      fitness <- hybrids_list$fitness
    }
    
    # Return the results if need
    if (i == maxiter)
    {
      best_ind <- which.max(fitness)
      return_list <- list(chromosomes = population,
                          par = population[best_ind, ],
                          value = fitness[best_ind])
      class(return_list) <- "gena"
      return(return_list)
    }
    
    # pop.initialize the parents
    
    parents_list <- gena.mating(population = population,
                                parents.n = parents.n,
                                fitness = fitness,
                                method = mating.method,
                                self = mating.self,
                                par = mating.par)
    parents <- parents_list$parents
    parents_fitness <- parents_list$fitness

    # Make the children via the crossover
    
    children_co <- NULL
    
    children_co <- gena.crossover(parents = parents,
                                  fitness = parents_fitness,
                                  prob = crossover.prob,
                                  method = crossover.method,
                                  par = crossover.par,
                                  iter = i)
    
    # Make the elite children (no need for a parents, just the best chromosomes)
    
    fitness_order <- order(fitness, decreasing = TRUE)
    children_elite <- population[fitness_order[1:elite.n], ]
    
    # Remove duplicate elite if need
    
    if(!elite.duplicates)
    {
      elite_counter <- 2
      j_prev <- fitness_order[1]
      for (j in fitness_order[-1])
      {
        if (elite_counter > elite.n)
        {
          break
        }
        
        if (all(population[j, ] != population[j_prev, ]))
        {
          children_elite[elite_counter, ] <- population[j, ]
          elite_counter <- elite_counter + 1
        }
        
        j_prev <- j
      }
    }

    # Mutate the children
    
    children_co <- gena.mutation(children = children_co,
                                 lower = lower, 
                                 upper = upper,
                                 method = mutation.method,
                                 prob = mutation.prob,
                                 prob.genes = mutation.genes.prob,
                                 par = mutation.par,
                                 iter = i)
    
    # Combine crossover and elite children to
    # make a new population

    population <- rbind(children_elite, children_co)
    
    if (info)
    {
      cat(paste0("Iter = ", i, ", ",
                 "mean = ", round(mean(fitness), 5), ", ",
                 "median = ", round(median(fitness), 5), ", ",
                 "best = ", round(max(fitness), 5), "\n"))
    }
  }

  return(NULL)
}