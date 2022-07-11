# Mutation operator
gena.mutation <- function(children,                        # children
                          lower,                           # lower bound of search space
                          upper,                           # upper bound of search space
                          prob = 0.2,                      # mutation probability for children
                          prob.genes = 1 / nrow(children), # mutation probabilities for genes
                          method = "constant",             # mutation type
                          par = 1,                         # mutation parameters
                          iter = NULL)                     # iteration of the
                                                           # genetic algorithm 
{
  # ---
  # The algorithm brief description:
  # 1. Randomly pick the children for mutation
  # 2. Perform the mutation for each gene of each children
  # ---
  
  # Prepare some variables
  
  n_children <- nrow(children)                             # number of children
  n_genes <- ncol(children)                                # number of genes
  
  # Select children for mutation
  
  random_values <- runif(n_children, 0, 1)                 # values to control for
                                                           # weather mutation should
                                                           # take place for a given
                                                           # child
  mutation_ind <- which(random_values <= prob)             # select the children
                                                           # to mutate
  
  for (i in mutation_ind)
  {
    random_values_genes <- runif(n_genes, 0, 1)
    genes_ind <- which(prob.genes <= random_values_genes)
    
    if (method == "constant")
    {
      for (j in genes_ind)
      {
        children[i, j] <- runif(1, lower[j], upper[j])
      }
    }
    
    if (method == "gaussian")
    {
      for (j in genes_ind)
      {
        children[i, j] <- children[i, j] + rnorm(1, mean = 0, sd = par[j])
        #children[i, j] <- min(max(children[i, j], lower[j]), upper[j])
      }
    }
    
    if (method == "percent")
    {
      par_adj <- par / 100
      children_adj <- abs(children[i, ]) * par_adj
      mutation_lower <- children[i, ] - children_adj
      mutation_upper <- children[i, ] + children_adj
      for (j in genes_ind)
      {
        children[i, j] <- runif(1, mutation_lower[j], mutation_upper[j])
        #children[i, j] <- min(max(children[i, j], lower[j]), upper[j])
      }
    }
  }
  
  return(children)                                         # return the children
}

# Assign default parameters for
# mutation algorithm depending
# on the "method"
gena.mutation.validate <- function(method, par, n_genes)
{
  # Validate the "method"
  
  methods <- c("constant", "gaussian", "percent")
  
  if (!(method %in% methods))
  {
    stop(paste0("Incorrect mutation.method argument. ",    # if the user has provided
                "It should be one of: ",                   # incorrect argument
                paste(methods, collapse = ", "),
                "\n"))
  }
  
  # Assign default parameters
  
  if (method == "gaussian")
  {
    if (is.null(par))
    {
      par <- 1
    } else {
      if (any((par <= 0)))
      {
        stop("par should be a vector with posive values")
      }
    }
    
    if (length(par) != n_genes)
    {
      if (length(par) == 1)
      {
        par <- rep(par, n_genes)
      } else {
        stop("length(par) should be equal to length(lower)")
      }
    }
  }
  
  if (method == "percent")
  {
    if (is.null(par))
    {
      par <- rep(20, n_genes)
    } else {
      if (length(par) == 1)
      {
        par <- rep(par, n_genes)
      } else {
        if (length(par) > n_genes)
        {
          stop(paste0("Please, insure that 'par' has the same size as the ",
                      "number of genes i.e. 'length(par)==ncol(lower)'"))
        }
      }
      
      if (any((par < 0)))
      {
        stop("par should be have non-negative values")
      }
    }
    
    if (length(par) != n_genes)
    {
      if (length(par) == 1)
      {
        par <- rep(par, n_genes)
      } else {
        stop("length(par) should be equal to length(lower)")
      }
    }
  }
  
  return(par)
}