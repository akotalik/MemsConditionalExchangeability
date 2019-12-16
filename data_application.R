#application of the method described in "Dynamic Borrowing in the Presence of Treatment Effect Heterogeneity"
#simulate dataset with structure closely resembling the one used in the data application section (not available)

# libraries used
library(mvtnorm)
library(R2jags)


#define parameters:
n=c(200, 300); muz=c(10, 18); sdz=c(3,4); mu0=c(3,5); beta_z=0.3; beta_t=c(5,5); beta_tz=-0.2; sd=c(3,4)

set.seed(1)

#simulate primary study: one treatment effect modifier Z
#treatment indicator
trt<- rep(0, n[1])
trt[sample(1:n[1],n[1]*0.5)] <- 1
#effect modifier
Z <- rnorm(n[1], mean = muz[1], sd=sdz[1])
#outcome
Y <- rnorm(n[1], mean= mu0[1] +Z*(beta_z) + beta_t[1]*trt + beta_tz*trt*Z, sd=sd[1])
primary <- data.frame(Y,Z,trt)
primary$df <- 1
primary$secondary <- 0

#simulate secondary study
trt<- rep(0, n[2])
trt[sample(1:n[2],n[2]*0.5)] <- 1
Z <- rnorm(n[2], mean = muz[2], sd=sdz[2])
Y <- rnorm(n[2], mean= mu0[2] +Z*(beta_z) + beta_t[2]*trt + beta_tz*trt*Z, sd=sd[2])
secondary <- data.frame(Y,Z,trt)
secondary$df <- 2
secondary$secondary <- 1

df <- rbind(primary, secondary)

#apply the proposed method:
#only one supplemental source, 2 exchangeability patterns:

#model 1: no borrowing
m1 <- lm(Y ~ 1 + trt + secondary + I(Z*(1-secondary)) + I(Z*secondary) + I(Z*trt*(1-secondary)) + I(Z*trt*secondary) + I(trt*secondary), data=df)

#model 2: borrowing
m2 <- lm(Y ~ 1 + trt + secondary + I(Z*(1-secondary)) + I(Z*secondary) + I(Z*trt*(1-secondary)) + I(Z*trt*secondary), data=df)

#these are used only to obtain the design matrices, which are then further modified
allModelsResults <- list(m1, m2)

#same models but centered by sample mean- useful only for calculating BIC as sample mean=MLE
sampmean <- mean(df$Z[df$df==1])
allModelsResults2 <- list(lm(Y ~ 1 + trt + secondary + I((Z-sampmean)*(1-secondary)) + I((Z-sampmean)*secondary) + I((Z-sampmean)*trt*(1-secondary)) + I((Z-sampmean)*trt*secondary) + I(trt*secondary), data=df), lm(Y ~ 1 + trt + secondary + I((Z-sampmean)*(1-secondary)) + I((Z-sampmean)*secondary) + I((Z-sampmean)*trt*(1-secondary)) + I((Z-sampmean)*trt*secondary), data=df))

#prior values
var_b=100^2; a=c(0.001, 0.001); b=c(0.001, 0.001); burnin=1000; ndraw=5000; n_sources=2; N_z=1

#type variable: 4 "types" of observations: primary and treated, primary untreated, same for secondary
#this is only relevant for what columns in the design matrix to center by mux
#this variable is used to tell which columns to center by mux for each row
df$type <- 1*(df$df==1 & df$trt==0)+ 2*(df$df==1 & df$trt==1) + 3*(df$df==2 & df$trt==0)+ 4*(df$df==2 & df$trt==1)+ 5*(df$df==3 & df$trt==0)+ 6*(df$df==3 & df$trt==1)


bic <- NULL
betas <- NULL
for(i in 1:length(allModelsResults)) {
  fit <- allModelsResults[[i]]
  #get design matrix for borrowing/not models
  X <- as.matrix(model.matrix(fit))
  
  #prior for the regression parameters: very diffused
  mub <- rep(0, ncol(X))
  VB <- var_b*diag(ncol(X))
  
  #set initial values for jags
  jags.inits <- function(){
    list("inv.var"=rep(1, n_sources))
  }
  
  #JAGS model:
  M <- function() {
    # Likelihood
    for(i in 1:n_tot){
      Y[i]   ~ dnorm(mu[i],prec[i])
      mu[i] <- (X[i, ]-cen[type[i], ]) %*% beta
    }
    
    #random design matrix centering: all the following code does is subtracts the random mux from the appropriate columns of design matrix in each row
    #this needs to be incorporated into jags, as mux is random itself and needs to be estimated jointly with beta
    for(j in 1:(4-1)){
      cen[1, j] <- 0.0
    }
    cen[1, 4] <- mux
    for(j in 5:J){
      cen[1, j] <- 0.0
    }
    
    for(j in 1:(4-1)){
      cen[2, j] <- 0.0
    }
    cen[2, 4] <- mux
    cen[2, 5] <- 0.0
    cen[2, 6] <- mux
    for(j in 7:J){
      cen[2, j] <- 0.0
    }
    
    for(j in 1:(5-1)){
      cen[3, j] <- 0.0
    }
    cen[3, 5] <- mux
    for(j in 6:J){
      cen[3, j] <- 0.0
    }
    
    for(j in 1:(5-1)){
      cen[4, j] <- 0.0
    }
    cen[4, 5] <- mux
    cen[4, 6] <- 0.0
    cen[4, 7] <- mux
    for(j in 8:J){
      cen[4, j] <- 0.0
    }
    
    #there are 4 "cen" values in total: one for each type of observation, all the above code does is it tells JAGS where to subtract mux
    
    #centering values likelihood
    for (j in 1:n_1){
      X[j, 4] ~ dmnorm(mux[], E)
    }
    
    # Prior for mux
    for(j in 1:N_z){
      mux[j] ~ dnorm(0.0, 1.0E-8)
    }
    
    #prior for E
    E   ~ dgamma(0.001, 0.001)
    #######
    
    # Prior for beta
    for(j in 1:J){
      beta[j] ~ dnorm(mub[j],1/(var_b))
    }
    
    # Prior for the inverse variance of outcome
    for(j in 1:nsource){
      inv.var[j]   ~ dgamma(a[j], b[j])
      sigma[j] <- 1/inv.var[j]
    }
    for(i in 1:n_tot){
      prec[i] <- inv.var[df[i]]
    }
  }
  
  
  
  #run MCMC:
  out <- jags(list(Y=df$Y, X=X, n_tot=nrow(df), df=df$df, type=df$type, nsource=n_sources, mub = mub, var_b= var_b, a = a, b =b,
                   J=ncol(X), N_z= N_z, n_1= n[1]),  inits = jags.inits, n.chains=1, n.thin=1, n.iter=ndraw, n.burnin = burnin, DIC=F, parameters.to.save =  c('beta'), M)
  
  #save MCMC output, only interested in single parameter: the main effect on trt
  out_mcmc <- as.mcmc(out)
  temp <- data.frame(unlist(out_mcmc[[1]]))
  names(temp) <- gsub("\\.", "", names(temp))
  names(temp) <- gsub("beta", "", names(temp))
  temp <- temp[ , order(as.numeric(names(temp)))]
  betas <- cbind(betas, temp[,2])
  
  #BIC weights:
  #iterated weighted least squares
  #use regression with centered X (by sample mean since that is MLE already)
  fit2 <- allModelsResults2[[i]]
  X_cen <- as.matrix(model.matrix(fit2))
  
  cf <- coef(fit2)
  s <- rep(NA, n_sources)
  for(k in 1:1000) {
    old <- cf
    for(j in 1:n_sources){
      s[j] <- sum((df$Y[df$df==j] - X_cen[df$df==j,] %*% cf)^2)/(length(df$Y[df$df==j]))
    }

    cf <- coef(lm(df$Y ~ X_cen + 0, weights = sqrt(rep(s, times=n))))
    if (sum(abs(old - cf)) < 0.0000000000001)
    {
      break;
    }
  }
  
  #BIC:
  bic <- c(bic, log(nrow(df))*(ncol(X_cen)+n_sources) - 2*dmvnorm(df$Y, mean= X_cen %*% cf, sigma = diag(rep(s, times=n)), log=TRUE))

}

#bic approximation
wts <- exp(-0.5*(bic-min(bic)))/sum(exp(-0.5*(bic-min(bic))))

#ESS calculation: use empirical estimate of posterior precision
rats <- (1/diag(var(betas)))/(1/var(betas[,1])) -1
esss <- n[1] *(sum(wts*rats))
ess <- n[1] + esss

#function that draws from the mixture for beta with given weights
rmixture = function(wts, betas){
  ##random generation of the indices
  id = sample(1:length(wts),prob=wts,size=nrow(betas),replace=TRUE)  
  id = cbind(1:nrow(betas),id)
  betas[id]
}

#posterior: draw from individual posteriors according to weights
samples <- rmixture(wts=wts, betas=betas)
#credible interval
CI <- c(quantile(samples, 0.025), quantile(samples, 0.975))
#point estimate is posterior mean
pointest <- mean(samples)

#the true value: primary study population treatment effect
#truth <- beta_t[1] + muz[1]*beta_tz

