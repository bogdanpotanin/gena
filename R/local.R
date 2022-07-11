# Estimate the gradient
gena.grad <- function(fn, par, 
                      eps = NULL, fn_args)
{
  n_par <- length(par)
  
  if(is.null(eps))
  {
    eps <- sqrt(.Machine$double.eps) * par
  }
  
  if(length(eps) == 1)
  {
    eps <- rep(eps, n_par)
  }
    
  par_plus <- par
  par_minus <- par
    
  val <- rep(0, n_par)

  for(i in 1:n_par)
  {
    par_plus[i] <- par[i] + eps[i]
    par_minus[i] <- par[i] - eps[i]
    
    fn_args$par <- par_plus
    fn_plus <- do.call(fn, fn_args)
 
    fn_args$par <- par_minus
    fn_minus <- do.call(fn, fn_args)
    
    par_plus <- par
    par_minus <- par
      
    val[i] <- (fn_plus - fn_minus) / (2 * eps[i])
  }
  
  return(val)
}

# Estimate the Hessian
gena.hessian <- function(gr, par, 
                         eps = NULL, fn_args, gr_args = NULL)
{
  n_par <- length(par)
  
  if (is.null(eps))
  {
    eps <- sqrt(.Machine$double.eps) * par
  }
  
  par_plus <- par
  par_minus <- par
  
  val <- matrix(0, n_par, n_par)
  
  for(i in 1:n_par)
  {
    par_plus[i] <- par[i] + eps[i]
    par_minus[i] <- par[i] - eps[i]
      
    gr_args$par <- par_plus
    gr_plus <- do.call(gr, gr_args)
      
    gr_args$par <- par_minus
    gr_minus <- do.call(gr, gr_args)
      
    par_plus <- par
    par_minus <- par
      
    val[i, ] <- (gr_plus - gr_minus) / (2 * eps[i])
  }
  
  for (i in 1:n_par)
  {
    for (j in 1:i)
    {
      if (i != j)
      {
        val[i, j] <- (val[i, j] + val[j, i]) / 2
        val[j, i] <- val[i, j]
      }
    }
  }
  
  return(val)
}

gena.lineSearch <- function(fn, gr = NULL, 
                            fn.val, gr.val, 
                            par, 
                            p, a = 1,
                            tau = 0.5, c = 0.5,
                            maxiter = 100, 
                            fn.args = NULL, gr.args = NULL)
{
  # Estimate a special value for Armijo condition
  m <- sum(gr.val * p)
  
  # Initialize variables to store parameters and
  # update parameters information from previous step
  par_old <- par
  a_old <- a
  
  # A list to return
  return.list <- list(a = a,
                      par = par)
  
  # If TRUE then update parameter will be increasing.
  # Otherwise it will be decreasing.
  is_increase = FALSE

  for (i in 1:maxiter)
  {
    return.list$iter <- i
    
    fn.args$par <- par + a * p                             # update parameter
    
    fn.val.new <- do.call(fn, fn.args)                     # estimate function value
    cond.val <- fn.val.new - fn.val - (a * c * m)          # estimate Armijo condition
                                                           # related value

    if (cond.val >= 0)                                     # check whether Armijo
    {                                                      # condition is satisfied
      if (i == 1)                                          # decide whether to increase
      {                                                    # or decrease update parameter
        is_increase <- TRUE
      }

      if (!is_increase)                                    # if the parameter should
      {                                                    # be decreased
        return.list$par <- fn.args$par
        return.list$a <- a
        return(return.list)
      }
    }
    
    if ((cond.val < 0) & is_increase)                      # if Armijo condition is not
    {                                                      # but we need to increase
      return.list$par <- par_old                           # the update parameter
      return.list$a <- a_old
      return(return.list)
    }
    
    par_old <- fn.args$par
    a_old <- a
    
    a <- ifelse(is_increase, a / tau, a * tau)
  }
  
  return(return.list)
}

gena.gd <- function(fn,                                    # function to optimize
                    gr = NULL,                             # gradient of the function
                    x0,                                    # initial point
                    iter = 1000,                           # the number of iterations
                    reltol = 1e-10,                        # relative tolerance
                    lr = 1,                                # learning rate
                    type = "GD",                           # optimization type
                    momentum = 0.9,
                    ...)                                   # additional arguments for fn and gr  
{
  dots <- list(...)
  
  # Get the number of estimated parameters
  n_x <- length(x0)
  
  # Initialize the matrix to store the
  # points for each iteration
  x <- matrix(NA, 
              nrow = iter,
              ncol = length(x0))
  x[1, ] <- x0
  
  # Initialize variable to store gradient
  # information
  gr.val <- matrix(NA, nrow = iter, ncol = n_x)

  # Start the algorithm
  for (i in 2:iter)
  {
    # Estimate the function value
    dots$par <- x[i - 1, ]
    fn.val <- do.call(fn, dots)

    # Estimate the analytical gradient if it has
    # been provided by the user. Otherwise apply
    # numeric optimization routine.
    if(is.null(gr))
    {
      gr.val[i, ] <- gena.grad(fn, x[i - 1, ], fn.args = dots)
    } else {
      dots$par <- x[i - 1, ]
      gr.val[i, ] <- do.call(gr, dots)
    }

    # Determine the direction
    dir <- gr.val[i, ]
    if (i >= 3)
    {
     G <- t(gr.val[3:(i - 1), ]) %*% gr.val[3:(i - 1), ]
     dir <- dir / sqrt(diag(G) + (1e-8))
    }
    # Perform the line search
    ls.result <- gena.lineSearch(fn = fn, gr = gr, 
                                 fn.val = fn.val, gr.val = gr.val[i, ],
                                 par = x[i - 1, ], p = dir, a = lr,
                                 fn.args = dots)
    x[i, ] <- ls.result$par
    lr <- ls.result$a
    print(lr)
    print(ls.result$iter)
    
    # Check whether termination conditions
    # are satisfied
    if (all(abs((x[i - 1, ] - x[i, ]) /
                x[i - 1, ]) < reltol))
    {
      x <- x[1:i, ]
      break
    }
    
    cat(paste0("iter-", i, ": ", fn.val, "\n"))
  }
  
  # Get the final point
  solution <- tail(x, 1)
  
  return_list <- list("points" = x,
                      "solution" = solution,
                      "fn_value" = fn(solution,    
                                      ... = ...),
                      "fn_grad" = gr.val)
  
  return(return_list)
}