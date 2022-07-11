#' Genetic Algorithm
#' @description This function allows to use genetic algorithm for
#' numeric global optimization purposes.
#' @param fn function to be maximized
#' @param gr gradient of the \code{fn}
#' @param lower lower bound of the search space
#' @param upper upper bound of the search space
#' @param pop.initial numeric matrix which rows are chromosomes to be
#' included into the initial population
#' @param pop.n integer representing the size of the population
#' @param mating.method the algorithm to be applied for a mating i.e. selection
#' of parents
#' @param mating.par parameters of the mating algorithm (see details)
#' @param mating.self logical; if \code{TRUE} then the chromosome may
#' mate with itself i.e. both parents may be the same chromosome
#' @param crossover.method an algorithm to be applied for crossover i.e.
#' creation of the children
#' @param crossover.par parameters of the crossover algorithm (see details)
#' @param crossover.prob probability of the crossover for each pair of parents
#' @param mutation.method algorithm to be applied for mutation i.e. random
#' change in some genes of the children
#' @param mutation.par parameters of the mutation algorithm (see details)
#' @param mutation.prob mutation probability for the chromosomes
#' @param mutation.genes.prob mutation probability for the genes
#' @param elite.n number of the elite children i.e. those which have the
#' highest function value and will be preserved for the next population
#' @param elite.duplicates logical; if \code{TRUE} then some elite children
#' may have the same genes
#' @param hybrid.method hybrids selection algorithm i.e. mechanism 
#' determining which chromosomes should be subject to local optimization
#' @param hybrid.par parameters of the hybridization algorithm (see details)
#' @param hybrid.prob probability of generating the hybrids each iteration
#' @param hybrid.opt.par parameters of the local optimization function
#' to be used for hybridization algorithm
#' @param hybrid.n number of hybrids that appear if hybridization
#' should take place during the iteration
#' @param maxiter maximum number of iterations of the algorithm
#' @param constr a function imposing a constraints on each of the chromosomes
#' @param info logical; if \code{TRUE} then some optimization related 
#' information will be printed each iteration
#' @param ... additional parameters to be passed to 
#' \code{fn}, \code{gr} and \code{constr} functions
#' @return This function returns the list containing the information on the
#' results of the genetic algorithm
gena <- function(
  fn,
  gr = NULL,
  lower,
  upper,
  pop.n = 100,
  pop.initial = NULL,
  mating.method = "rank",
  mating.par = NULL,
  mating.self = FALSE,
  crossover.method = "local",
  crossover.par = NULL,
  crossover.prob = 0.8,
  mutation.method = "percent",
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
  constr = NULL,
  info = TRUE,
  ...)
{   
  # Initialize some variables
  
  pop.n <- ceiling(pop.n)
  elite.n <- ceiling(elite.n)
  n_genes <- length(lower)
  n_parents <- pop.n - elite.n
  
  fitness <- rep(NA, pop.n)
  
  # Perform pop.initial validation
  
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

    if (ncol(pop.initial) != n_genes)
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
    warning("Since 'pop.n' is not even it will be incremented by one\n")
  }
  
    # method.mating and mating.par
  
  mating.par <- gena.mating.validate(method = mating.method,
                                     par = mating.par,
                                     n_parents = n_parents)
  
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
  
  crossover.par <- gena.crossover.validate(method = crossover.method,
                                           par = crossover.par)
  
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
  
  mutation.par <- gena.mutation.validate(method = mutation.method,
                                         par = mutation.par,
                                         n_genes = n_genes)
  
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
    mutation.genes.prob <- rep(mutation.genes.prob, n_genes)
  }
  
  if (mutation.genes.prob_n > n_genes)
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
  
    # elite.n
  
  if (elite.n <= 1)
  {
    stop("Please, insure that 'elite.n' is an integer greater than one\n")
  }
  
  if ((elite.n %% 2) != 0) 
  {
    elite.n <- elite.n + 1
    warning("Since 'elite.n' is not even it will be incremented by one\n")
  }
  
  if ((pop.n - elite.n) < 2)
  {
    stop(paste0("Please, insure that ((pop.n - elite.n) >= 2) i.e. ",
                "it seems that 'elite.n' is too large",
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
  
  hybrid.par <- gena.mating.validate(method = hybrid.method,
                                       par = hybrid.par,
                                       n_parents = n_parents)
  
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
  
    # constr
  
  if(!is.null(constr))
  {
    if (!is.function(constr))
    {
      stop("Please, insure that 'constr' is a function\n")
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
  
  population <- gena.population(pop.n = pop.n,
                                lower = lower, 
                                upper = upper,
                                pop.initial = pop.initial)
  
  # -------------------------------------
  # Perform the genetic algorithm routine
  # -------------------------------------
  
  for(i in 1:maxiter)
  {
    # Get the parameters to be passed to fitness
    
    dots <- list(...)
    
    # Impose the constraint on the parameters
    
    if (!is.null(constr))
    {
      for (j in 1:pop.n)
      {
        dots$par <- population[j, ]
        population[j, ] <- do.call(what = constr, args = dots)
      }
      
      dots$par <- NULL
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
      return_list <- list(par = population[best_ind, ],
                          value = fitness[best_ind],
                          chromosomes = population)
      
      return(return_list)
    }
    
    # pop.initialize the parents
    
    parents_list <- gena.mating(population = population,
                                n_parents = n_parents,
                                fitness = fitness,
                                method = mating.method,
                                self = mating.self,
                                par = mating.par)
    parents <- parents_list$parents
    parents_fitness <- parents_list$fitness

    # Make the children via the crossover
    
    children_co <- NULL
    
    if (!is.function(crossover.method))
    {
      children_co <- gena.crossover(parents = parents,
                                    fitness = parents_fitness,
                                    prob = crossover.prob,
                                    method = crossover.method,
                                    par = crossover.par,
                                    iter = i)
    } else {
      crossover.par$parents <- parents
      crossover.par$fitness <- parents_fitness
      crossover.par$iter <- i
      crossover.par$n_genes <- n_genes
      crossover.par$pop.n <- pop.n
      
      children_co <- do.call(what = crossover.method, 
                             args = crossover.par)
    }
    
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
    
    if (!is.function(mutation.method))
    {
      children_co <- gena.mutation(children = children_co,
                                   lower = lower, 
                                   upper = upper,
                                   method = mutation.method,
                                   prob = mutation.prob,
                                   prob.genes = mutation.genes.prob,
                                   par = mutation.par,
                                   iter = i)
    } else {
      mutation.par$parents <- children_co
      mutation.par$fitness <- parents_fitness
      mutation.par$iter <- i
      children_co <- do.call(what = mutation.method, 
                             args = mutation.par)
    }
    
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