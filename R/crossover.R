# Crossover the parents
gena.crossover <- function(parents,                        # the parents
                           fitness,                        # fitness of the parents
                           prob = 0.8,                     # crossover probability
                           method = "local",               # crossover type
                           par = NULL,                     # crossover parameters
                           iter = NULL,                    # iteration of the
                                                           # genetic algorithm
                           ...)                            # additional parameters                
{
  # ---
  # The algorithm brief description:
  # 1. Randomly pick the parents for the crossover
  # 2. Perform the crossover for each pair of parents
  # ---
  
  # Prepare some variables
  
  n_parents <- nrow(parents)                               # number of parents
  n_genes <- ncol(parents)                                 # number of genes
  
  children <- parents                                      # matrix to store children
                                                           # initially contains
                                                           # the parents

  # Perform the crossover for each chromosome
  
  random_values <- runif(n_parents / 2, 0, 1)              # values to control for
                                                           # weather crossover should
                                                           # take place for a given
                                                           # pair of chromosomes
  crossover_ind <- which(random_values <= prob) * 2        # select the pairs
                                                           # for a crossover
  
  for (i in crossover_ind)
  {
    ind <- c(i - 1, i)                                     # indexes of the fist 
                                                           # and the second parents 
    i_adj <- i / 2                                         # half the index
    
    if (method == "split")
    {
      split_ind <- rbinom(n_genes, 1, 0.5)                 # indexes to get from
                                                           # the other parent
      children[ind[1], split_ind] <- parents[ind[2],       
                                             split_ind] 
      children[ind[2], split_ind] <- parents[ind[1], 
                                             split_ind]
    }
    
    if (method == "arithmetic")
    {
      children[ind[1], ] <- parents[ind[1], ] * par[1] +      # mix the parents with
        parents[ind[2], ] * (1 - par[1])                      # a given weight
      children[ind[2], ] <- parents[ind[2], ] * par[1] +
        parents[ind[1], ] * (1 - par[1])
    }
    
    if (method == "local")                                    # local crossover
    {
      weight <- runif(n = n_genes, min = 0, max = 1)          # random weight
      children[ind[1], ] <- parents[ind[1], ] * weight +      # mix the parents with
                            parents[ind[2], ] * (1 - weight)  # a random weight
      children[ind[2], ] <- parents[ind[2], ] * weight +
                            parents[ind[1], ] * (1 - weight)
    }
    
    if (method == "flat")                                          # flat crossover
    {
      for (j in 1:n_genes)                                         # for each gene
      {                                                            # and for both
        for (t in 1:2)                                             # children
        {
          children[ind[t], j] <- runif(n = 1,                      # randomly pick    
                                       min(parents[ind[1], j],     # the value between   
                                           parents[ind[2], j]),    # the genes of
                                       max(parents[ind[1], j],     # the parents
                                           parents[ind[2], j]))
        }
      }
    }
  }
  
  return(children)                                         # return the children
}

# Assign default parameters for
# crossover algorithm depending
# on the "method"
gena.crossover.validate <- function(method, par)
{
  # Validate the "method"
  
  methods <- c("split", "arithmetic", "local", "flat")     # the list of all
                                                           # available methods
  
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect crossover.method argument. ",   # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
  
  # Assign default parameters
  
  if (method == "arithmetic")
  {
    if (!is.null(par))
    {
      if ((length(par) != 1) | (!is.numeric(par)))
      {
        stop(paste0("Incorrect crossover.par agrument. Please, insure that ",
                    "(length(crossover.par) == 1) and ",
                    "is.numeric(crossover.par)",
                    "\n"))
      }
    } else {
      par <- 0.5
    }
  }
  
  return(par)
}