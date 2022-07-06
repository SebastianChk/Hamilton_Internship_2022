# setwd(...)

library(dplyr)
library(ggplot2)
library(MASS)
library(data.table)
library(reshape2)
library(tidyverse)
library(stringr)
library(fs)

PD <- 'plots_lda'

file.path(getwd(), PD) %>%
  dir.create()

# Error estimators --------------------------------------------------------


# First I need functions to evaluate explicitly all the errors.

# Anderson's W statistic
W.statistic <- function(Sigma.inv, Xbar0, Xbar1, Xi){
  # This function assumes that Sigma is a diagonal matrix
  a <- (Xi - (Xbar0 + Xbar1)/2)
  b <- diag(Sigma.inv)
  c <- (Xbar0 - Xbar1)
  return(sum(a*b*c))
}

# Resubstitution error estimator
resub <- function(Sigma.inv, Xbar0, Xbar1, X, n, n0=0, n1=0){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  W.func <- function(X){return(W.statistic(Sigma.inv, Xbar0, Xbar1, X))}
  W <- apply(X, 1, W.func)
  
  W0 <- W[1:n0]
  e0 <- mean( W0 <= 0)
  a0 <- n0/n
  
  W1 <- W[round(n0+1):n]
  e1 <- mean( W1 >  0)
  a1 <- n1/n
  return(a0*e0 + a1*e1)
}

# Plug-in error estimator
plugin <- function(Sigma.inv, Xbar0, Xbar1, X){
  a <- (Xbar0 - Xbar1)^2
  b <- diag(Sigma.inv)
  deltahat <- sqrt( sum(a*b) )
  return(pnorm(-deltahat/2))
}

# Smoothed resubstitution error estimator
smoothed.resub <- function(Sigma.inv, Xbar0, Xbar1, X, b, deltahat, n, n0=0, n1=0){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  g <- function(x){pnorm(-x/(b*deltahat))}
  
  W.func <- function(X){return(W.statistic(Sigma.inv, Xbar0, Xbar1, X))}
  W <- apply(X, 1, W.func)

  W0 <- W[1:n0]
  e0 <- mean(   g(W0) )
  a0 <- n0/n
  
  W1 <- W[round(n0+1):n]
  e1 <- mean( 1-g(W1) )
  a1 <- n1/n

  return(a0*e0 + a1*e1)
}

# True error
true.error <- function(Sigma.inv, Xbar0, Xbar1, mu0, mu1, n, n0=0, n1=0){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  a0 <- n0/n
  a1 <- n1/n

  numerator0 <- sum( 
    ( mu0-(Xbar0+Xbar1)/2 ) * diag(Sigma.inv) * (Xbar0-Xbar1) 
    )
  
  numerator1 <- sum(
    ( mu1-(Xbar0+Xbar1)/2 ) * diag(Sigma.inv) * (Xbar0-Xbar1) 
    )
  
  denominator <- sqrt(sum(
    (Xbar0-Xbar1)^2 * diag(Sigma.inv)
    ))

  e0 <- pnorm(-numerator0/denominator)
  e1 <- pnorm(numerator1/denominator)

  return(a0*e0 + a1*e1)
}

# Asymptotic-finite error expectations ------------------------------------
# Note:

# The limit of the expectation of the resub and plug-in estimators are the same.
# Hence I only define a function for the resubsitution error expectation

resub.af <- function(n, p, n0=0, n1=0, delta.squared=4){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  J0 <- p/n0
  J1 <- p/n1
  
  return(pnorm(-sqrt(delta.squared + J0 + J1)/2))
}

smoothed.resub.af <- function(n, p, b, n0=0, n1=0, delta.squared=4){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  J0 <- p/n0
  J1 <- p/n1
  
  return(pnorm( -sqrt(delta.squared + J0 + J1) / (2*sqrt(1+b^2)) ))
}

true.error.af <- function(n, p, n0=0, n1=0, delta.squared=4){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  J0 <- p/n0
  J1 <- p/n1
  
  a0 <- n0/n
  a1 <- n1/n
  
  denom <- sqrt(delta.squared + J0 + J1)
  
  e0 <- pnorm(-0.5*(delta.squared+J1-J0)/denom)
  e1 <- pnorm(-0.5*(delta.squared-J1+J0)/denom)

  return(a0*e0 + a1*e1)
}




# Redundant-finite error expectations -------------------------------------
f <- function(a, h, n0, n1, J0, J1, a0, a1){
  return(a0*pnorm(-0.5*a/sqrt(h(n0, n1, J0, J1))) + a1*pnorm(-0.5*a/sqrt(h(n1, n0, J1, J0))))
  }

true.error.rf <- function(n0, n1, J0, J1, a0, a1, delta.squared=4){
  h <- function(n0, n1, J0, J1){
    return(delta.squared * (1+1/n1)+ J0+ J1+ J1/(2*n1) + J0/(2*n0) )
  }
  a <- delta.squared+J1-J0
  return(f(a, h, n0, n1, J0, J1, a0, a1))
}

# The paper makes the distinction between alpha and alpha hat, but in these
# simulations we take these to be the same
resub.rf <- function(n0, n1, J0, J1, a0, a1, delta.squared=4){
  h <- function(n0, n1, J0, J1){
    return(delta.squared * (1+1/n1)+ J0+ J1+ J1/(2*n1) + J0/(2*n0) )
  }
  a <- delta.squared+J1+J0
  return(f(a, h, n0, n1, J0, J1, a0, a1))
}

plugin.rf <- function(n0, n1, J0, J1, a0, a1, delta.squared=4){
  h <- function(n0, n1, J0, J1){
    return(delta.squared * (1+1/n1+1/n0)+ J0+ J1+ J1/(2*n1) + J0/(2*n0) + J0/n1)
  }
  a <- delta.squared+J1+J0
  return(f(a, h, n0, n1, J0, J1, a0, a1))
}

smoothed.resub.rf <- function(n0, n1, J0, J1, a0, a1, b, delta.squared=4){
  h <- function(n0, n1, J0, J1){
    return(delta.squared * (1+ b^2 +1/n1+1/(2*n0)) + J0/(2*n1)+ J1/(2*n1) + (1+b^2) * (J0 + J1))
  }
  a <- delta.squared+J1+J0
  return(f(a, h, n0, n1, J0, J1, a0, a1))
}

# Monte Carlo Simulations -------------------------------------------------

# Function to carry out one MC simulation
lda.mc <- function(N, p, n, seed=99, n0=n/2, n1=n/2, delta.squared=4){
  start <- Sys.time()
  # For comments, see the previous section, 'Test simulation'
  J0 <- p/n0
  J1 <- p/n1
  
  set.seed(seed)
  mu0 <- rnorm(p, mean=1, sd=2)
  mu1 <- rep(0, p)
  Sigma <- diag( p/4 *mu0^2, nrow=p, ncol=p)
  Sigma.inv <- diag(1/diag(Sigma))
  
  # This mean and covariance works too, I just used the less trivial one that
  # follows because I was experimenting and making sure that it did not affect
  # the results. Both of these setups satisfy the key property: that the
  # mahalanobis distance between the distribution means is 2
  
  # mu0 <- c(1, rep(0,p-1))
  # mu1 <- c(0, rep(0,p-1))
  # Sigma <- diag(c(1/delta.squared, rep(1,p-1)), nrow=p, ncol=p)
  # Sigma.inv <- diag(1/diag(Sigma))
  
  errors <- matrix(rep(-10, 4*N), nrow=N, ncol=4)
  colnames(errors) <- c('resub', 'plugin', 'sresub', 'true')
  for(i in 1:N){
    X <- rbind(mvrnorm(n0, mu0, Sigma), mvrnorm(n1, mu1, Sigma))
    
    Xbar0 <- apply(X[1:n0,], 2, mean)
    Xbar1 <- apply(X[round(n0+1):n,], 2, mean)
    
    deltahat <- ((Xbar0 - Xbar1)^2 * diag(Sigma.inv)) %>% sum %>% sqrt
    
    # Relies on n_0=n_1
    bopt <- 1/delta.squared * sqrt( (J1 + J0)*(J0 + J1 + 2*delta.squared) )
    
    errors[i, 'resub'] <- resub(Sigma.inv, Xbar0, Xbar1, X, n)
    errors[i, 'plugin'] <- plugin(Sigma.inv, Xbar0, Xbar1, X)
    errors[i, 'sresub'] <- smoothed.resub(Sigma.inv, Xbar0, Xbar1, X, bopt, deltahat, n)
    errors[i, 'true'] <- true.error(Sigma.inv, Xbar0, Xbar1, mu0, mu1, n)
  }
  end <- Sys.time()
  print(end-start)
  return(errors)
}

# Asymptotic-finite sample approximation ----------------------------------
lda.af <- function(p, n, n0=0, n1=0, delta.squared=4){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  J0 <- p/n0
  J1 <- p/n1
  
  bopt <- 1/delta.squared * sqrt( (J1 + J0)*(J0 + J1 + 2*delta.squared) )
  
  errors <- c()
  c('resub', 'plugin', 'sresub', 'true')
  errors['resub'] <- resub.af(n, p)
  errors['plugin'] <- resub.af(n, p)
  errors['sresub'] <- smoothed.resub.af(n, p, bopt)
  errors['true'] <- true.error.af(n, p)
  
  return(errors)
}


# Redundant-finite sample approximation -----------------------------------
lda.rf <- function(p, n, n0=0, n1=0, delta.squared=4){
  if(n0==0){
    n0 <- n/2
    n1 <- n/2
  }
  J0 <- p/n0
  J1 <- p/n1
  
  bopt <- 1/delta.squared * sqrt( (J1 + J0)*(J0 + J1 + 2*delta.squared) )
  
  errors <- c()
  c('resub', 'plugin', 'sresub', 'true')
  errors['resub'] <- resub.rf(n0, n1, J0, J1, n0/n, n1/n, delta.squared)
  errors['plugin'] <- resub.rf(n0, n1, J0, J1, n0/n, n1/n, delta.squared)
  errors['sresub'] <- smoothed.resub.rf(n0, n1, J0, J1, n0/n, n1/n, bopt, delta.squared)
  errors['true'] <- true.error.rf(n0, n1, J0, J1, n0/n, n1/n, delta.squared)
  
  return(errors)
}

# Running everything for all n and p --------------------------------------

# If this were python I would use a multi-index pandas dataframe in order to
# plot this. My solution instead here is to basically create nested lists.

# n and p values to find errors for
ns <- c(40, 60, 80, 100, 120, 140, 160, 180, 200)
ps <- c(3, 9, 30, 100)

# estimators to find errors for
est_names <- c('resub', 'plugin', 'sresub', 'true')
# est_names <- c('plugin')


# approximations to use
approx_names <- c('mc', 'af')#, 'rf')

# Number of Monte-Carlo samples
N <- 1000

# df is a dataframe of which copies will be made. results is a list that  will
# have such a dataframe at its lowest level, and the dataframe will correspond
# to the error estimates for a given estimator and approximation. The
# rows/columns of the dataframe correspond to sample size/dimensionality
# respectively
results <- list()
df <- data.frame(matrix(data=-1,
                        ncol=length(ps),
                        nrow=length(ns)),
                 row.names=ns)
colnames(df) <- ps

for(i in approx_names){
  est_list <- list()
  for(j in est_names){
    est_list[[j]] <- copy(df)
  }
  results[[i]] <- est_list
}

# Insert Monte-Carlo Simulation errors into results$mc

for(n in ns){
  for(p in ps){
    mc.errors <- lda.mc(N, p, n)
    mc.errors <- apply(mc.errors, 2, mean)
    n.str <- as.character(n)
    p.str <- as.character(p)
    for(est_name in est_names){
      # only doing the rounding temporarily
      results$mc[[est_name]][n.str, p.str] <- mc.errors[est_name] 
    }
  }
}

# Insert asymptotic-finite sample approximation errors into results$af
for(n in ns){
  for(p in ps){
    af.errors <- lda.af(p, n)
    n.str <- as.character(n)
    p.str <- as.character(p)
    for(est_name in est_names){
      # only doing the rounding temporarily
      results$af[[est_name]][n.str, p.str] <- af.errors[est_name]
    }
  }
}

# redundant-finite (didn't include in report)
# for(n in ns){
#   for(p in ps){
#     rf.errors <- lda.rf(p, n)
#     n.str <- as.character(n)
#     p.str <- as.character(p)
#     for(est_name in est_names){
#       # only doing the rounding temporarily
#       results$rf[[est_name]][n.str, p.str] <- rf.errors[est_name]
#     }
#   }
# }


# Plots -------------------------------------------------------------------


# Need to create an array/df for each individual plot. So, for fixed
# est_name(plot), I define 8 curves, each corresponding to an element of the
# cartesian product of {3,9,30,100} and {'mc', 'af'}

# In hindsight, I should have really just made the original 'results' list to
# have this structure at the lowest level. But at this point this is unnecessary
# to fix since my code works
df.cnames <- c('n', '3mc', '9mc', '30mc', '100mc', '3af', '9af', '30af', '100af')

# actually not entirely sure if this works for the redundant-finite sample approximation, will have to re-check

for(est_name in est_names){
  # Need to create an array/df for each individual plot. So, for fixed
  # est_name(plot), I define 8 curves, each corresponding to an element of the
  # cartesian product of {3,9,30,100} and {'mc', 'af', 'rf'}
  
  # In hindsight, I should have really just made the original 'results' list to
  # have this structure at the lowest level. But at this point this is
  # unnecessary to fix since my code works
  df <- data.frame(matrix(nrow=9, ncol=9))
  
  colnames(df) <- df.cnames
  df$n <- ns
  
  # extract the corresponding results from 'results' list
  for(ni in 1:length(ns)){
    for(p in ps){
      for(app in approx_names){
        n.str <- as.character(ns[ni])
        p.str <- as.character(p)
        df[ni, paste(p, app, sep='')] <- results[[app]][[est_name]][n.str, p.str] 
      }
    }
  }
  
  # Prepare dataframe for plotting
  df.l <- melt(df, id.var='n')
  df.l$approx <- df.l$variable %>% str_sub(-2, -1)
  df.l$p <- df.l$variable %>% str_sub(end=-3) %>% factor(levels=c('3', '9', '30', '100'))
  
  convert <- function(string){
    if(string=='af'){
      return('Asymptotic approximation')
    } else if(string=='mc'){
      return('Monte Carlo simulation')
    } else if(string=='rf'){
      return('Redundant-finite approximation')
    }
  }
  
  df.l$approx <- apply(matrix(df.l$approx), 1, FUN=convert)
  
  # Create plot
  plot.lda <- ggplot(data=df.l, aes(x=n, y=value, group=variable)) + 
    geom_point(aes(colour=approx, shape=approx), size=2) +
    scale_shape(solid = FALSE) +
    geom_line(aes(colour=approx, linetype=p)) + 
    labs(colour='Estimation method', shape='Estimation method') +
    theme(aspect.ratio = 4/8)
  
  # Set y limits to (approximately) match those from the paper
  ylims <- if (est_name=='true' || est_name=='sresub') c(0.15, 0.3) else c(0.025, 0.16)
  plot.lda <- plot.lda + coord_cartesian(ylim=ylims) 
  
  # Bayes Error line
  plot.lda <- plot.lda + geom_hline(yintercept=0.1586)
  
  # Remove y-label, cannot typeset greek letters in R
  plot.lda <- plot.lda + ylab(NULL)
  plot.lda
  # Save plots
  ggsave(path_join(c(PD, paste0(est_name,'.png'))), plot.lda, dpi = 500, height=4.5, width=10)
}

# save(errors, file='data/lda-errors.RData')
# load('data/lda-errors.RData')
# df.l
# 
# df.l$approx <- df.l$variable %>% str_sub(-2, -1)
# 
# convert <- function(string){
#   if(string=='af'){
#     return('Asymptotic-finite approximation')
#   } else if(string=='mc'){
#     return('Monte Carlo simulation')
#   } else if(string=='rf'){
#     return('Redundant-finite approximation')
#   }
# }
# df.l$approx <- apply(matrix(df.l$approx), 1, FUN=convert)
# df.l
