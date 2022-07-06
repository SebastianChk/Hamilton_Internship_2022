# setwd(...)
# Setup -------------------------------------------------------------------

# plots folder/directory
PD <- 'plots'

# data folder/directory
DD <- 'data'

# experiment name (so that the resulting files can be differentiated by file name alone)
exp.name <- ''

library(fs)
library(dplyr)
library(purrr)
library(ggplot2)
# The Fig 4 and 5 plots take a very long time to run
# setwd(...)
set.seed(100)

# set n and p
n <- 4000
p <- 800

# sigmoid function for convenience
sigmoid <- function(xb){return(exp(xb)/(1+exp(xb)))}

# function that simulates covariate matrix - I created this function so that the
# code can easily be generalised to an arbitrary distribution of the covariate
# matrix by simply editing this function
rcovariate <- function(n, p){
  return(rnorm(n*p, 0, 1/sqrt(n)) %>% matrix(n, p))
}

# Fig 2 -------------------------------------------------------------------

# create design matrix and parameter vector
X <- rcovariate(n, p)
beta <- c(rep(10, 1/8 * p),  rep(-10, 1/8 * p), rep(0, 3/4 * p))

# Generate Y sample
probs <- sigmoid(X %*% beta)
Y <- rbinom(n, 1, probs)

# Fit logistic regression to Y sample
model1 <- glm(Y ~ X-1, family=binomial)

# Prepare data for plotting
df1 <- as.data.frame(beta)
df1$index <- seq(1,p)
df1$mle_beta <- model1$coefficients

ggplot(df1) +
  geom_point(aes(index, mle_beta), color='blue', size=1) +
  geom_point(aes(index, beta), size=0.5) +
  xlab('Index') +
  ylab('Coefficients(true & fitted)') +
  theme(aspect.ratio = 4/9)

# Create 'plots' folder if it doesn't already exist
file.path(getwd(), PD) %>%
  dir.create()

# Save figure
ggsave(path_join(c(PD, paste(exp.name, 'fig2.png', sep=''))), dpi=900)


# Fig 3 -------------------------------------------------------------------
set.seed(100)
# replace previous beta with N(3, 16) beta
beta <- rnorm(p, 3, 4)

# Sample Y and fit model
probs <- sigmoid(X %*% beta)
Y <- rbinom(n, 1, probs)
model2 <- glm(Y ~ X-1, family=binomial)

# put together data for fig 2A
df2A <- as.data.frame(beta)
df2A$mle_beta <- model2$coefficients

# Fig 3A ------------------------------------------------------------------

ggplot(df2A) +
  geom_point(aes(beta, mle_beta), color='blue', size=1) +
  geom_abline(intercept=0, slope=1, color='black', size=1) + 
  geom_abline(intercept=0, slope=1.499, color='red', size=1) + 
  xlab('True signal') +
  ylab('MLE estimate') +
  theme(aspect.ratio = 8/8)


ggsave(path_join(c(PD, paste(exp.name, 'fig3A.png', sep=''))), dpi=500, height=5, width=5)

# Fig 3B ------------------------------------------------------------------
pred_probs <- predict(model2, as.data.frame(X), type='response')
Xb <- X %*% beta

df2B <- as.data.frame(cbind(Xb, probs, pred_probs))
colnames(df2B) <- c('Xbeta', 'true_probs', 'pred_probs')
df2B %>% head
df2B <- df2B[order(df2B$Xbeta),]
df2B %>% head
df2B$order <- seq(1,n)

# ggplot(df2B) +
#   geom_point(aes(Xbeta, pred_probs), color='blue', size=0.8, alpha=1) +
#   geom_line(aes(Xbeta, true_probs), size=1.4) +
#   theme(aspect.ratio = 8/8) +
#   ylab('Probabilities (True and predicted)')#+
#   # scale_x_sqrt() #+ 
#   # 
# 
# ggsave(path_join(c(PD, paste(exp.name, 'fig3B.png', sep=''))), dpi=500, height=5, width=5)

ggplot(df2B) +
  geom_point(aes(order, pred_probs), color='blue', size=0.9, alpha=0.4) +
  geom_line(aes(order, true_probs), size=1.4) +
  theme(aspect.ratio = 8/8) +
  ylab('Probabilities (True and predicted)') + 
  xlab('Observation index')

ggsave(path_join(c(PD, paste(exp.name, 'fig3Bscaled.png', sep=''))),  dpi=500, height=5, width=5)


# Fig 4 -------------------------------------------------------------------

### NOTE: the p-value is actually 1 - p-value as I have defined it here


set.seed(100)
# keep beta fixed, with half zero and the others ~ N(7, 1)
beta <- c(rep(0,p/2), rnorm(n=p/2, mean=7, sd=1)) %>% matrix(nrow=p, ncol=1)

# number of monte carlo samples
m <- 1000

# placeholder matrices which will be replaced with the relevant samples
# more efficient than appending on each iteration I think
nullcoeffs <- matrix(rep(-1,p/2), nrow=p/2, ncol=m)
SEbeta1 <- rep(-1, m)

LLs <- matrix(rep(-1,2*m), nrow=m)

start <- Sys.time()
for(j in 1:m){

  # sample X and Y
  Xj <- rcovariate(n, p)
  probsj <- sigmoid(Xj %*% beta)
  Yj <- rbinom(n, 1, probsj) %>% matrix(n, 1)
  
  # fit full model
  modelj <- glm(Yj ~ Xj-1, family=binomial)
  
  # fit null model
  Xjnull <- Xj[,2:p]
  modeljnull <- glm(Yj ~ Xjnull-1, family=binomial)
  
  # Extract log-likelihoods
  LLj <- modelj$deviance
  LLjnull <- modeljnull$deviance
  # LLR <- LLjnull - LLj
  # pval <- pchisq(LLR, df=1)
  LLs[j,] <- c(LLj, LLjnull)
  
  # update matrix with null coefficient estimates
  nullcoeffs[,j] <- modelj$coefficients[1:round(p/2)]
  
  # update matrix with the estimate of beta0's SE
  SEbeta1[j] <- summary(modelj)$coefficients[, 2][1]
  
  # keep track of progress
  print(j)
}
pbeta1 <- 1 - pchisq(LLs[,2]-LLs[,1], df=1)
pbeta1 %>% hist


# create /data folder if it doesn't already exist
file.path(getwd(), "data") %>%
  dir.create()

# Save results of the above for loop (it takes a while to run)
save(nullcoeffs, file=path_join(c(DD, paste(exp.name, 'nullcoeffs-fig4.RData', sep=''))))
save(SEbeta1, file=path_join(c(DD, paste(exp.name, 'SEbeta1-fig4.RData', sep=''))))
save(pbeta1, file=path_join(c(DD, paste(exp.name, 'pbeta1-fig5.RData', sep=''))))

load(path_join(c(DD, paste(exp.name, 'nullcoeffs-fig4.RData', sep=''))))
load(path_join(c(DD, paste(exp.name, 'SEbeta1-fig4.RData', sep=''))))
load(path_join(c(DD, paste(exp.name, 'pbeta1-fig5.RData', sep=''))))

# Fig 4A ------------------------------------------------------------------
set.seed(100)
# Re-load the beta vector if necessary
# beta <- c(rep(0,400), rnorm(400, 7, 1))

mc_ses <- apply(nullcoeffs, 1, sd)
ggplot() + 
  geom_histogram(aes(mc_ses), color='black', fill='#990000', bins=50) +
  geom_vline(xintercept=2.663, linetype='dashed') +
  theme(aspect.ratio = 8/8) +
  ylab('Relative Counts') +
  xlab('SEs of coefficients')
  
colors <- c("y=1.499x, using modern big p theory"="")
ggsave(path_join(c(PD, paste(exp.name, 'fig4A.png', sep=''))), dpi=900, height=5, width=5)


# Fig 4B ------------------------------------------------------------------

ggplot() + 
  geom_histogram(aes(SEbeta1), color='black', fill='#6699FF') +
  theme(aspect.ratio = 8/8) +
  ylab('Relative Counts') +
  xlab('SEs of coefficients')

ggsave(path_join(c(PD, paste(exp.name, 'fig4B.png', sep=''))), dpi=900, height=5, width=5)


# Fig 5 -------------------------------------------------------------------

pbeta1 %>% hist

ggplot() +
  geom_histogram(aes(pbeta1),
                 color='red',
                 fill='#000099',
                 boundary=0, bins=40) +
  theme(aspect.ratio=9/18) +
  xlim(0,1) +
  # coord_cartesian(xlim=c(0,1)) +
  ylab('Relative Counts') +
  xlab('P-Values')

ggsave(path_join(c(PD, paste(exp.name, 'fig5temp.png', sep=''))), dpi=900, height=5, width=10)

pbeta1
  



