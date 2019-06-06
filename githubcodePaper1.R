
#manuscript simulation code

# libraries used
library(parallel)
cores=8
library(mvtnorm)
library(ggplot2)
library(invgamma)
library(rjags)
library(compiler)
library(R2jags)


#function that draws from the mixture for beta
rmixture = function(wts, betas){
  ##random generation of the indices
  id = sample(1:length(wts),prob=wts,size=nrow(betas),replace=TRUE)  
  id = cbind(1:nrow(betas),id)
  betas[id]
}


#simulation function

doSim <- function(seed, n, mu0, beta_z, p, beta_t, beta_m, muz, varz, sd, var_b=NULL, a=NULL, b=NULL, method, burnin=NULL, ndraw=NULL, BICiterations=NULL, nint=NULL, marginal=NULL) {
  set.seed(seed) ## Set the random seed for this simulation
  n_tot=sum(n)
  n_sources = length(n)
  n_m = length(beta_m)
  #simulate data:
  #make first n_m covariates effect modifiers
  
  #primary
  trt <- rbinom(n[1], 1, p[1])
  Z <- rmvnorm(n[1], mean = muz[[1]], sigma=as.matrix(varz[[1]]))
  Y <- rnorm(n[1], mean= mu0[1] +Z%*%as.matrix(beta_z) + beta_t[1]*trt + (Z[,1:n_m]%*%as.matrix(beta_m))*trt, sd=sd[1])
  primary <- data.frame(Y,Z,trt)
  primary$df <- 1
  primary$secondary <- 0
  
  #secondary
  secondary <- NULL
  for(i in 2:n_sources) {
    trt <- rbinom(n[i], 1, p[i])
    Z <- rmvnorm(n[i], mean = muz[[i]], sigma=as.matrix(varz[[i]]))
    Y <- rnorm(n[i], mean= mu0[i] +Z%*%as.matrix(beta_z) + beta_t[i]*trt + (Z[,1:n_m]%*%as.matrix(beta_m))*trt, sd=sd[i])
    df <- i
    secondary <- rbind(secondary, data.frame(Y,Z,trt,df))
  }
  secondary$secondary <- 1
  
  #data frame together
  data <- data.frame(rbind(primary, secondary))
  
  Z<- as.matrix(data[,2:(length(beta_z)+1)])
  for(i in 1:n_sources) {
    for(j in 1:length(beta_z)) {
      x <- 1*(data$df==i) * (Z[,j])
      data <- cbind(data, x)
    }
  }
  
  Z<- as.matrix(data[,2:(length(beta_z)+1)])
  for(i in 1:n_sources) {
    for(j in 1:length(beta_z)) {
      x <- 1*(data$df==i) * (Z[,j] - mean(Z[data$secondary==0,j]))
      data <- cbind(data, x)
    }
  }
  
  #dummies for source
  for(i in 2:n_sources) {
    x <- 1*(data$df==i)
    data <- cbind(data, x)
  }
  
  #source*trt*z interactions
  for(i in 1:n_sources) {
    for(j in 1:n_m) {
      x <- 1*(data$df==i) * (Z[,j])*data$trt
      data <- cbind(data, x)
    }
    
  }
  
  for(i in 1:n_sources) {
    for(j in 1:n_m) {
      x <- 1*(data$df==i) * (Z[,j] - mean(Z[data$secondary==0,j]))*data$trt
      data <- cbind(data, x)
    }
    
  }
  
  
  #create interactions for source*trt which will decide whether I borrow or not on the trt effect
  for(i in 2:n_sources) {
    x <- 1*(data$df==i) * data$trt
    data <- cbind(data, x)
  }
  names(data) <- c("Y", paste("Z", 1:length(beta_z), sep=""), "trt", "df", "secondary", paste("Z", rep(1:n_sources, each=length(beta_z)), rep(1:length(beta_z), n_sources), sep=""), paste("Z_centered", rep(1:n_sources, each=length(beta_z)), rep(1:length(beta_z), n_sources), sep="") , paste("S", 2:n_sources, sep=""),paste("Z_trt", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""), paste("Z_centered_trt", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""),paste("trt_S", 2:n_sources, sep=""))
  
  data$type <- 1*(data$df==1 & data$trt==0)+ 2*(data$df==1 & data$trt==1) + 3*(data$df==2 & data$trt==0)+ 4*(data$df==2 & data$trt==1)+ 5*(data$df==3 & data$trt==0)+ 6*(data$df==3 & data$trt==1)
  
  #first marginal MEM: only considers main effects
  if(method == "MEM1") {
    
    ls <- vector("list", n_sources-1)
    for(i in 1:(n_sources-1)) {
      x<-c(TRUE, FALSE)
      ls[[i]] <- x
    }
    
    regMat <- expand.grid(ls)
    regressors <- paste("trt_S", 2:n_sources, sep="")
    
    allModelsList <- apply(regMat, 1, function(x) as.formula(
      paste(c("Y ~ 1 + trt",paste("S", 2:n_sources, sep=""), regressors[x]),
            collapse=" + ")) )
    
    allModelsResults <- lapply(allModelsList,
                               function(x) lm(x, data=data))
    
    bic <- NULL
    aic <- NULL
    betas <- NULL
    for(i in 1:length(allModelsResults)) {
      fit <- allModelsResults[[i]]
      X <- as.matrix(model.matrix(fit))
      mub <- rep(0, ncol(X))
      VB <- var_b*diag(ncol(X))
      jags.inits <- function(){
        list("inv.var"=rep(1, n_sources))
      }
      
      #JAGS model:
      M <- function() {
        # Likelihood
        for(i in 1:n_tot){
          Y[i]   ~ dnorm(mu[i],prec[i])
          mu[i] <- X[i, ] %*% beta
        }
        
        # Prior for beta
        for(j in 1:J){
          beta[j] ~ dnorm(mub[j],1/(var_b))
        }
        
        # Prior for the inverse variance
        for(j in 1:nsource){
          inv.var[j]   ~ dgamma(a[j], b[j])
          sigma[j] <- 1/inv.var[j]
        }
        for(i in 1:n_tot){
          prec[i] <- inv.var[df[i]]
        }
      }
      
      
      out <- jags(list(Y=data$Y, X=X, n_tot=sum(n), df=data$df, nsource=n_sources, mub = mub, var_b= var_b, a = a, b =b,
                       J=ncol(X)),  inits = jags.inits, n.chains=1, n.thin=2, n.iter=ndraw, n.burnin = burnin, DIC=F, parameters.to.save =  c('beta'), M)
      
      out_mcmc <- as.mcmc(out)
      temp <- data.frame(unlist(out_mcmc[[1]]))
      names(temp) <- gsub("\\.", "", names(temp))
      names(temp) <- gsub("beta", "", names(temp))
      temp <- temp[ , order(as.numeric(names(temp)))]
      betas <- cbind(betas, temp[,2])
      
      if(marginal == "BIC" | marginal == "AIC") {
        #BIC weights:
        #iterated weighted least squares
        cf <- coef(fit)
        s <- rep(NA, n_sources)
        for(k in 1:BICiterations) {
          old <- cf
          for(j in 1:n_sources){
            s[j] <- sum((data$Y[data$df==j] - X[data$df==j,] %*% cf)^2)/(length(data$Y[data$df==j]))
          }
          #W <- diag(rep(s, times=n))
          #cf <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% data$Y
          
          cf <- coef(lm(data$Y ~ X + 0, weights = sqrt(rep(s, times=n))))
          if (sum(abs(old - cf)) < 0.0000000000001)
          {
            break;
          }
        }
        
        #BIC:
        bic <- c(bic, - log(nrow(data))*(ncol(X)+n_sources) + 2*dmvnorm(data$Y, mean= X %*% cf, sigma = diag(rep(s, times=n)), log=TRUE))
        aic <- c(aic, - 2*(ncol(X)+n_sources) + 2*dmvnorm(data$Y, mean= X %*% cf, sigma = diag(rep(s, times=n)), log=TRUE))
        
      }
      
    }
    
    if(marginal == "BIC") {
      wts <- exp(0.5*(bic-max(bic)))/sum(exp(0.5*(bic-max(bic))))
    }else if(marginal == "AIC") {
      wts <- exp(0.5*(aic-max(aic)))/sum(exp(0.5*(aic-max(aic))))
    }
    
    #posterior: draw from individual posteriors according to weights
    samples <- rmixture(wts=wts, betas=betas)
    CI <- c(quantile(samples, 0.025), quantile(samples, 0.975))
    pointest <- mean(samples)
    vars <- diag(var(betas))
  } 
  
  ##################
  #proposed method
  else if(method == "MEM2") {
    
    ls <- vector("list", n_sources-1)
    for(i in 1:(n_sources-1)) {
      x<-c(TRUE, FALSE)
      ls[[i]] <- x
    }
    
    regMat <- expand.grid(ls)
    regressors <- paste("trt_S", 2:n_sources, sep="")
    
    allModelsList <- apply(regMat, 1, function(x) as.formula(
      paste(c("Y ~ 1 + trt",paste("S", 2:n_sources, sep=""), paste("Z", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""), paste("Z_trt", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""), regressors[x]),
            collapse=" + ")) )
    
    allModelsResults <- lapply(allModelsList,
                               function(x) lm(x, data=data))
    #results from simply centering at sample mean
    allModelsList2 <- apply(regMat, 1, function(x) as.formula(
      paste(c("Y ~ 1 + trt",paste("S", 2:n_sources, sep=""), paste("Z_centered", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""), paste("Z_centered_trt", rep(1:n_sources, each=n_m), rep(1:n_m, n_sources), sep=""), regressors[x]),
            collapse=" + ")) )
    
    allModelsResults2 <- lapply(allModelsList2,
                                function(x) lm(x, data=data))
    
    lik2 <- NULL
    bic <- NULL
    aic <- NULL
    betas <- NULL
    for(i in 1:length(allModelsResults)) {
      fit <- allModelsResults[[i]]
      X <- as.matrix(model.matrix(fit))
      mub <- rep(0, ncol(X))
      VB <- var_b*diag(ncol(X))
      #run gibbs sampler
      #system.time(res <- gibbs(X, data$Y, data$df, burnin = burnin, ndraw = ndraw, mub=mub, VB=VB, as=a, bs=b, N_z=length(beta_m), n=n))
      #betas <- cbind(betas, res$beta[,1])
      jags.inits <- function(){
        list("inv.var"=rep(1, n_sources))
      }
      
      if(length(beta_m)==1){
        if(n_sources==2){
          #JAGS model:
          M <- function() {
            # Likelihood
            for(i in 1:n_tot){
              Y[i]   ~ dnorm(mu[i],prec[i])
              #mu[i] <- X[i, ] %*% beta
              mu[i] <- (X[i, ]-cen[type[i], ]) %*% beta
            }
            
            
            ###### 1
            
            #the random centering works for 2 or 3 sources only
            
            #primary type 1
            
            for(j in 1:(pos[1]-1)){
              cen[1, j] <- 0.0
            }
            cen[1, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in (pos[1]+N_z):J){
              cen[1, j] <- 0.0
            }
            
            #primary type 2
            
            for(j in 1:(pos[1]-1)){
              cen[2, j] <- 0.0
            }
            cen[2, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in pos[2]:(pos[2]+N_z-1)){
              cen[2, j] <- 0.0
            }
            cen[2, (pos[2]+N_z):(pos[2]+2*N_z-1)] <- mux
            for(j in (pos[2]+2*N_z):J){
              cen[2, j] <- 0.0
            }
            
            #secondary 1 type 1
            for(j in 1:(pos[2]-1)){
              cen[3, j] <- 0.0
            }
            cen[3, pos[2]:(pos[2]+N_z-1)] <- mux
            for(j in (pos[2]+N_z):J){
              cen[3, j] <- 0.0
            }
            
            #secondary 1 type 2
            for(j in 1:(pos[2]-1)){
              cen[4, j] <- 0.0
            }
            cen[4, pos[2]:(pos[2]+N_z-1)] <- mux
            for(j in (pos[2]+N_z):(pos[2]+2*N_z-1)){
              cen[4, j] <- 0.0
            }
            cen[4, (pos[2]+2*N_z):(pos[2]+3*N_z-1)] <- mux
            for(j in (pos[2]+3*N_z):J){
              cen[4, j] <- 0.0
            }
            
            
            #centering values
            for (j in 1:n_1){
              X[j, pos[1]:(N_z+pos[1]-1)] ~ dmnorm(mux[], E)
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
            
            # Prior for the inverse variance
            for(j in 1:nsource){
              inv.var[j]   ~ dgamma(a[j], b[j])
              sigma[j] <- 1/inv.var[j]
            }
            for(i in 1:n_tot){
              prec[i] <- inv.var[df[i]]
            }
          }
          
          pos=c(2+(n_sources), 2+(n_sources)+length(beta_m))
        } else if(n_sources==3){
          
          #JAGS model:
          M <- function() {
            # Likelihood
            for(i in 1:n_tot){
              Y[i]   ~ dnorm(mu[i],prec[i])
              #mu[i] <- X[i, ] %*% beta
              mu[i] <- (X[i, ]-cen[type[i], ]) %*% beta
            }
            
            
            ###### 2
            
            #primary type 1
            for(j in 1:(pos[1]-1)){
              cen[1, j] <- 0.0
            }
            cen[1, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in (pos[1]+N_z):J){
              cen[1, j] <- 0.0
            }
            
            #primary type 2
            
            for(j in 1:(pos[1]-1)){
              cen[2, j] <- 0.0
            }
            cen[2, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in (pos[1]+N_z):(pos[2]-1)){
              cen[2, j] <- 0.0
            }
            cen[2, (pos[2]):(pos[2]+N_z-1)] <- mux
            for(j in (pos[2]+N_z):J){
              cen[2, j] <- 0.0
            }
            
            #secondary 1 type 1
            for(j in 1:(pos[3]-1)){
              cen[3, j] <- 0.0
            }
            cen[3, pos[3]:(pos[3]+N_z-1)] <- mux
            for(j in (pos[3]+N_z):J){
              cen[3, j] <- 0.0
            }
            
            #secondary 1 type 2
            for(j in 1:(pos[3]-1)){
              cen[4, j] <- 0.0
            }
            cen[4, pos[3]:(pos[3]+N_z-1)] <- mux
            for(j in (pos[3]+N_z):(pos[4]-1)){
              cen[4, j] <- 0.0
            }
            cen[4, pos[4]:(pos[4]+N_z-1)] <- mux
            for(j in (pos[4]+N_z):J){
              cen[4, j] <- 0.0
            }
            
            #secondary 2 type 1
            for(j in 1:(pos[5]-1)){
              cen[5, j] <- 0.0
            }
            cen[5, pos[5]:(pos[5]+N_z-1)] <- mux
            for(j in (pos[5]+N_z):J){
              cen[5, j] <- 0.0
            }
            
            #secondary 2 type 2
            for(j in 1:(pos[5]-1)){
              cen[6, j] <- 0.0
            }
            cen[6, pos[5]:(pos[5]+N_z-1)] <- mux
            for(j in (pos[5]+N_z):(pos[6]-1)){
              cen[6, j] <- 0.0
            }
            cen[6, pos[6]:(pos[6]+N_z-1)] <- mux
            for(j in (pos[6]+N_z):J){
              cen[6, j] <- 0.0
            }
            
            
            #centering values
            for (j in 1:n_1){
              X[j, pos[1]:(N_z+pos[1]-1)] ~ dmnorm(mux[], E)
            }
            
            # Prior for mux
            for(j in 1:N_z){
              mux[j] ~ dnorm(0.0, 1.0E-8)
            }
            
            #prior for E
            E   ~ dgamma(0.001, 0.001)
            
            # Prior for beta
            for(j in 1:J){
              beta[j] ~ dnorm(mub[j],1/(var_b))
            }
            
            # Prior for the inverse variance
            for(j in 1:nsource){
              inv.var[j]   ~ dgamma(a[j], b[j])
              sigma[j] <- 1/inv.var[j]
            }
            for(i in 1:n_tot){
              prec[i] <- inv.var[df[i]]
            }
          }
          
          pos=c(2+(n_sources), 2+(n_sources)+(n_sources)*length(beta_m), 2+(n_sources)+length(beta_m), 2+(n_sources)+(n_sources)*length(beta_m) + length(beta_m), 2+(n_sources)+2*length(beta_m), 2+(n_sources)+2*length(beta_m) + (n_sources-1)*length(beta_m)+1)
          
        }
      }else{
        if(n_sources==2){
          #JAGS model:
          M <- function() {
            # Likelihood
            for(i in 1:n_tot){
              Y[i]   ~ dnorm(mu[i],prec[i])
              #mu[i] <- X[i, ] %*% beta
              mu[i] <- (X[i, ]-cen[type[i], ]) %*% beta
            }
            
            
            ###### 
            #primary type 1
            
            for(j in 1:(pos[1]-1)){
              cen[1, j] <- 0.0
            }
            cen[1, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in (pos[1]+N_z):J){
              cen[1, j] <- 0.0
            }
            
            #primary type 2
            
            for(j in 1:(pos[1]-1)){
              cen[2, j] <- 0.0
            }
            cen[2, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in pos[2]:(pos[2]+N_z-1)){
              cen[2, j] <- 0.0
            }
            cen[2, (pos[2]+N_z):(pos[2]+2*N_z-1)] <- mux
            for(j in (pos[2]+2*N_z):J){
              cen[2, j] <- 0.0
            }
            
            #secondary 1 type 1
            for(j in 1:(pos[2]-1)){
              cen[3, j] <- 0.0
            }
            cen[3, pos[2]:(pos[2]+N_z-1)] <- mux
            for(j in (pos[2]+N_z):J){
              cen[3, j] <- 0.0
            }
            
            #secondary 1 type 2
            for(j in 1:(pos[2]-1)){
              cen[4, j] <- 0.0
            }
            cen[4, pos[2]:(pos[2]+N_z-1)] <- mux
            for(j in (pos[2]+N_z):(pos[2]+2*N_z-1)){
              cen[4, j] <- 0.0
            }
            cen[4, (pos[2]+2*N_z):(pos[2]+3*N_z-1)] <- mux
            for(j in (pos[2]+3*N_z):J){
              cen[4, j] <- 0.0
            }
            
            #centering values
            for (j in 1:n_1){
              X[j, pos[1]:(N_z+pos[1]-1)] ~ dmnorm(mux[], E[,])
            }
            
            # Prior for mux
            for(j in 1:N_z){
              mux[j] ~ dnorm(0.0, 1.0E-8)
            }
            
            #prior for E
            E[1:N_z, 1:N_z] ~ dwish(Omega[,], N_z)
            #######
            
            # Prior for beta
            for(j in 1:J){
              beta[j] ~ dnorm(mub[j],1/(var_b))
            }
            
            # Prior for the inverse variance
            for(j in 1:nsource){
              inv.var[j]   ~ dgamma(a[j], b[j])
              sigma[j] <- 1/inv.var[j]
            }
            for(i in 1:n_tot){
              prec[i] <- inv.var[df[i]]
            }
          }
          pos=c(2+(n_sources), 2+(n_sources)+length(beta_m))
        } else if(n_sources==3){
          
          #JAGS model:
          M <- function() {
            # Likelihood
            for(i in 1:n_tot){
              Y[i]   ~ dnorm(mu[i],prec[i])
              #mu[i] <- X[i, ] %*% beta
              mu[i] <- (X[i, ]-cen[type[i], ]) %*% beta
            }
            
            
            ###### 4
            
            #primary type 1
            for(j in 1:(pos[1]-1)){
              cen[1, j] <- 0.0
            }
            cen[1, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in (pos[1]+N_z):J){
              cen[1, j] <- 0.0
            }
            
            #primary type 2
            
            for(j in 1:(pos[1]-1)){
              cen[2, j] <- 0.0
            }
            cen[2, pos[1]:(pos[1]+N_z-1)] <- mux
            for(j in (pos[1]+N_z):(pos[2]-1)){
              cen[2, j] <- 0.0
            }
            cen[2, (pos[2]):(pos[2]+N_z-1)] <- mux
            for(j in (pos[2]+N_z):J){
              cen[2, j] <- 0.0
            }
            
            #secondary 1 type 1
            for(j in 1:(pos[3]-1)){
              cen[3, j] <- 0.0
            }
            cen[3, pos[3]:(pos[3]+N_z-1)] <- mux
            for(j in (pos[3]+N_z):J){
              cen[3, j] <- 0.0
            }
            
            #secondary 1 type 2
            for(j in 1:(pos[3]-1)){
              cen[4, j] <- 0.0
            }
            cen[4, pos[3]:(pos[3]+N_z-1)] <- mux
            for(j in (pos[3]+N_z):(pos[4]-1)){
              cen[4, j] <- 0.0
            }
            cen[4, pos[4]:(pos[4]+N_z-1)] <- mux
            for(j in (pos[4]+N_z):J){
              cen[4, j] <- 0.0
            }
            
            #secondary 2 type 1
            for(j in 1:(pos[5]-1)){
              cen[5, j] <- 0.0
            }
            cen[5, pos[5]:(pos[5]+N_z-1)] <- mux
            for(j in (pos[5]+N_z):J){
              cen[5, j] <- 0.0
            }
            
            #secondary 2 type 2
            for(j in 1:(pos[5]-1)){
              cen[6, j] <- 0.0
            }
            cen[6, pos[5]:(pos[5]+N_z-1)] <- mux
            for(j in (pos[5]+N_z):(pos[6]-1)){
              cen[6, j] <- 0.0
            }
            cen[6, pos[6]:(pos[6]+N_z-1)] <- mux
            for(j in (pos[6]+N_z):J){
              cen[6, j] <- 0.0
            }
            
            
            #centering values
            
            for (j in 1:n_1){
              X[j, pos[1]:(N_z+pos[1]-1)] ~ dmnorm(mux[], E[,])
            }
            
            # Prior for mux
            for(j in 1:N_z){
              mux[j] ~ dnorm(0.0, 1.0E-8)
            }
            
            #prior for E
            E[1:N_z, 1:N_z] ~ dwish(Omega[,], N_z)
            #######
            
            # Prior for beta
            for(j in 1:J){
              beta[j] ~ dnorm(mub[j],1/(var_b))
            }
            
            # Prior for the inverse variance
            for(j in 1:nsource){
              inv.var[j]   ~ dgamma(a[j], b[j])
              sigma[j] <- 1/inv.var[j]
            }
            for(i in 1:n_tot){
              prec[i] <- inv.var[df[i]]
            }
          }
          
          pos=c(2+(n_sources), 2+(n_sources)+(n_sources)*length(beta_m), 2+(n_sources)+length(beta_m), 2+(n_sources)+(n_sources)*length(beta_m) + length(beta_m), 2+(n_sources)+2*length(beta_m), 2+(n_sources)+2*length(beta_m) + (n_sources-1)*length(beta_m)+1)
          
        } 
      }
      
      out <- jags(list(Y=data$Y, X=X, n_tot=sum(n), df=data$df, type=data$type, nsource=n_sources, mub = mub, var_b= var_b, a = a, b =b, pos= pos,
                       J=ncol(X), N_z= length(beta_m), n_1= n[1], Omega= diag(length(beta_m))),  inits = jags.inits, n.chains=1, n.thin=2, n.iter=ndraw, n.burnin = burnin, DIC=F, parameters.to.save =  c('beta'), M)
      
      out_mcmc <- as.mcmc(out)
      temp <- data.frame(unlist(out_mcmc[[1]]))
      names(temp) <- gsub("\\.", "", names(temp))
      names(temp) <- gsub("beta", "", names(temp))
      temp <- temp[ , order(as.numeric(names(temp)))]
      betas <- cbind(betas, temp[,2])
      
      if(marginal == "BIC" | marginal == "AIC") {
        #BIC weights:
        #iterated weighted least squares
        fit2 <- allModelsResults2[[i]]
        X_cen <- as.matrix(model.matrix(fit2))
        
        cf <- coef(fit2)
        s <- rep(NA, n_sources)
        for(k in 1:BICiterations) {
          old <- cf
          for(j in 1:n_sources){
            s[j] <- sum((data$Y[data$df==j] - X_cen[data$df==j,] %*% cf)^2)/(length(data$Y[data$df==j]))
          }
         
          cf <- coef(lm(data$Y ~ X_cen + 0, weights = sqrt(rep(s, times=n))))
          if (sum(abs(old - cf)) < 0.0000000000001)
          {
            break;
          }
        }
        
        #BIC:
        
        bic <- c(bic, - log(nrow(data))*(ncol(X_cen)+n_sources) + 2*dmvnorm(data$Y, mean= X_cen %*% cf, sigma = diag(rep(s, times=n)), log=TRUE))
        aic <- c(aic, - 2*(ncol(X_cen)+n_sources) + 2*dmvnorm(data$Y, mean= X_cen %*% cf, sigma = diag(rep(s, times=n)), log=TRUE))
        
        
      }
      
    }
    
    if(marginal == "BIC") {
      wts <- exp(0.5*(bic-max(bic)))/sum(exp(0.5*(bic-max(bic))))
    }else if(marginal == "AIC") {
      wts <- exp(0.5*(aic-max(aic)))/sum(exp(0.5*(aic-max(aic))))
    }
    
    #posterior: draw from individual posteriors according to weights
    samples <- rmixture(wts=wts, betas=betas)
    CI <- c(quantile(samples, 0.025), quantile(samples, 0.975))
    pointest <- mean(samples)
    vars <- diag(var(betas))
  } 
  else if(method == "PRIMARY") {
    #fit model only based on primary data
    fmla <- as.formula(
      paste(c("Y ~ 1 + trt",paste("Z", rep(1:1, each=n_m), rep(1:n_m, 1), sep="")),
            collapse=" + "))
    fit <- lm(fmla, data=data[data$secondary==0,])
    
    CI <- confint(fit)[2,]
    pointest <- fit$coef[2]
    wts=rep(0, 2^(n_sources-1))
    vars <- rep(0, 2^(n_sources-1))
    
  }
  
  return(list(CI=CI, mean= pointest, weights= wts, vars= vars))
  
}


#scenario 1: no trt effect modifiers

muz=list(c(18,-22), c(15,-18))
beta_m=c(0)
nsim=1000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(13-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(11-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(9-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(8-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(7-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(6-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(5-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(4-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(2-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(1-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(0-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-1-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-2-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-5-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-7-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),

  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(13-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(11-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(9-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(8-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(7-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(6-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(5-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(4-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(2-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(1-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(0-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-1-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-2-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-5-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0, 0), p=c(0.5, 0.5), beta_t=c(-7-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[1]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18,-22), c(15,-18)), varz=list(diag(c(2^2,3^2)), diag(c(4^2,2^2))), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2")

)


#summary matrix: first column is coverage, then bias, then MSE, then power to detect nonzero trt eff.
summ <- matrix(nrow=length(pars),ncol=4)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}


wt <- NULL
sdratio <- NULL
pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####

  #### Do the simulation with these settings
  sim.results <- t(mclapply(1:nsim,doSim,n=par$n, mu0=par$mu0, beta_z=par$beta_z, p=par$p, beta_t=par$beta_t, beta_m=par$beta_m, muz=par$muz, varz=par$varz, sd=par$sd, var_b=par$var_b, a=par$a, b=par$b, burnin=par$burnin, ndraw=par$ndraw, BICiterations=par$BICiterations, method = par$method, marginal=par$marginal, mc.cores=cores, mc.silent=T))
  result <- unlist(sim.results)
  result <- matrix(result, ncol=3 +2^(length(par$n)-1)+ 2^(length(par$n)-1), byrow=T)
  #target of estimation: trt effect at mean values of modifiers in primary study, not par$muz[[1]][1:length(par$beta_m)] but result[,4:(4+length(par$beta_m))]
  cvrg <- is.between(as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)), result[,1], result[,2])
  pwr <- is.between(0, result[,1], result[,2])
  ## Summarize the simulation results
  summ[i,1] <- sum(cvrg==T)/nsim
  summ[i,2] <- mean(result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))
  summ[i,3] <- sum((result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))**2)/nsim
  summ[i,4] <- sum(pwr==F)/nsim
  #wt <- cbind(wt, result[,4])
  #sdratio <- c(sdratio, mean(result[,8]/result[,9]))

  sink('output2not.txt')
  print(summ[1:i,])
  setTxtProgressBar(pb, i)
  sink()
})
close(pb)

summsc1 <- summ
save(summsc1, file = "summsc1.RData")

save.image(file='scenario1.RData')


#scenario 4: no effect modifiers, one confounder we wrongly think is a modifier , vary muz in primary study: now marginally and conditionally exchangeable at all points

muz=list(c(99), c(12))
beta_m=c(0)
nsim=1000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(32), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(28), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(24), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(20), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(16), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(12), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(10), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(6), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(0), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(-4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(-8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),

  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(32), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(28), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(24), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(20), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(18), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(16), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(12), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(10), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(6), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(0), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(-4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0), muz=list(c(-8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2")
)

#summary matrix: first column is coverage, then bias, then MSE, then power to detect nonzero trt eff.
summ <- matrix(nrow=length(pars),ncol=4)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}


wt <- NULL
sdratio <- NULL
pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####

  #### Do the simulation with these settings
  sim.results <- t(mclapply(1:nsim,doSim,n=par$n, mu0=par$mu0, beta_z=par$beta_z, p=par$p, beta_t=par$beta_t, beta_m=par$beta_m, muz=par$muz, varz=par$varz, sd=par$sd, var_b=par$var_b, a=par$a, b=par$b, burnin=par$burnin, ndraw=par$ndraw, BICiterations=par$BICiterations, method = par$method, marginal=par$marginal, mc.cores=cores, mc.silent=T))
  result <- unlist(sim.results)
  result <- matrix(result, ncol=3 +2^(length(par$n)-1)+ 2^(length(par$n)-1), byrow=T)
  #target of estimation: trt effect at mean values of modifiers in primary study, not par$muz[[1]][1:length(par$beta_m)] but result[,4:(4+length(par$beta_m))]
  cvrg <- is.between(as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)), result[,1], result[,2])
  pwr <- is.between(0, result[,1], result[,2])
  ## Summarize the simulation results
  summ[i,1] <- sum(cvrg==T)/nsim
  summ[i,2] <- mean(result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))
  summ[i,3] <- sum((result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))**2)/nsim
  summ[i,4] <- sum(pwr==F)/nsim
  #wt <- cbind(wt, result[,4])
  #sdratio <- c(sdratio, mean(result[,8]/result[,9]))

  sink('output2not.txt')
  print(summ[1:i,])
  setTxtProgressBar(pb, i)
  sink()
})
close(pb)

summsc2 <- summ
save(summsc2, file = "summsc2.RData")

save.image(file='scenario2.RData')


#scenario 2: one strong modifier, vary muz in primary study: now conditionally exchangeable at all points, marginally at only one

muz=list(c(99), c(12))
beta_m=c(0.5)
nsim=1000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(32), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(28), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(24), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(20), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(18), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(16), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(12), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(10), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(6), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(0), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),

  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(32), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(28), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(24), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(20), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(18), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(16), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(14), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(12), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(10), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(6), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(0), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-4), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 400), mu0=c(5, -2),  beta_z=c(0.3), p=c(0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-8), c(12)), varz=list(as.matrix(2^2), as.matrix(3^2)), sd=c(5,4), var_b=100^2, a=c(0.001, 0.001), b=c(0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2")
)


#summary matrix: first column is coverage, then bias, then MSE, then power to detect nonzero trt eff.
summ <- matrix(nrow=length(pars),ncol=4)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}


wt <- NULL
sdratio <- NULL
pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####

  #### Do the simulation with these settings
  sim.results <- t(mclapply(1:nsim,doSim,n=par$n, mu0=par$mu0, beta_z=par$beta_z, p=par$p, beta_t=par$beta_t, beta_m=par$beta_m, muz=par$muz, varz=par$varz, sd=par$sd, var_b=par$var_b, a=par$a, b=par$b, burnin=par$burnin, ndraw=par$ndraw, BICiterations=par$BICiterations, method = par$method, marginal=par$marginal, mc.cores=cores, mc.silent=T))
  result <- unlist(sim.results)
  result <- matrix(result, ncol=3 +2^(length(par$n)-1)+ 2^(length(par$n)-1), byrow=T)
  #target of estimation: trt effect at mean values of modifiers in primary study, not par$muz[[1]][1:length(par$beta_m)] but result[,4:(4+length(par$beta_m))]
  cvrg <- is.between(as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)), result[,1], result[,2])
  pwr <- is.between(0, result[,1], result[,2])
  ## Summarize the simulation results
  summ[i,1] <- sum(cvrg==T)/nsim
  summ[i,2] <- mean(result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))
  summ[i,3] <- sum((result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))**2)/nsim
  summ[i,4] <- sum(pwr==F)/nsim
  #wt <- cbind(wt, result[,4])
  #sdratio <- c(sdratio, mean(result[,8]/result[,9]))

  sink('output2not.txt')
  print(summ[1:i,])
  setTxtProgressBar(pb, i)
  sink()
})
close(pb)

summsc3 <- summ
save(summsc3, file = "summsc3.RData")

save.image(file='scenario3.RData')

#scenario 3: one strong modifier, vary muz in primary study, two supplemental sources

muz=list(c(99), c(2), c(16))
beta_m=c(0.5)
nsim=1000
RNGkind("L'Ecuyer-CMRG")
pars=list(
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(36), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(33), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(30), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(26), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(24), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(22), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(20), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(18), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(16), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(14), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(12), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(10), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(8), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(6), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(4), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(2), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(0), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-2), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-4), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-6), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-10), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-13), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-16), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM1"),
  
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(36), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(33), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(30), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(26), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(24), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(22), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(20), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(18), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(16), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(14), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(12), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(10), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(8), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(6), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(4), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(2), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(0), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-2), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-4), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-6), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-10), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-13), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2"),
  list(n=c(100, 200, 200), mu0=c(5, -2, -2),  beta_z=c(0.3), p=c(0.5, 0.5, 0.5), beta_t=c(3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m), 3-muz[[2]][1:length(beta_m)]%*%as.matrix(beta_m)), beta_m=c(0.5), muz=list(c(-16), c(2), c(16)), varz=list(as.matrix(2^2), as.matrix(3^2), as.matrix(3^2)), sd=c(5,4,5), var_b=100^2, a=c(0.001, 0.001, 0.001), b=c(0.001, 0.001, 0.001), burnin=1000, ndraw=4000, BICiterations=1000, marginal="BIC", method="MEM2")
)


#summary matrix: first column is coverage, then bias, then MSE, then power to detect nonzero trt eff.
summ <- matrix(nrow=length(pars),ncol=4)

is.between <- function(x, a, b) {
  (x - a)  *  (b - x) >= 0
}


wt <- NULL
sdratio <- NULL
pb= txtProgressBar(min = 0, max = length(pars), style = 3, char=":)")
system.time(for(i in 1:length(pars)) {
  ## Define the parameters
  par <- pars[[i]] ####
  
  #### Do the simulation with these settings
  sim.results <- t(mclapply(1:nsim,doSim,n=par$n, mu0=par$mu0, beta_z=par$beta_z, p=par$p, beta_t=par$beta_t, beta_m=par$beta_m, muz=par$muz, varz=par$varz, sd=par$sd, var_b=par$var_b, a=par$a, b=par$b, burnin=par$burnin, ndraw=par$ndraw, BICiterations=par$BICiterations, method = par$method, marginal=par$marginal, mc.cores=cores, mc.silent=T))
  result <- unlist(sim.results)
  result <- matrix(result, ncol=3 +2^(length(par$n)-1)+ 2^(length(par$n)-1), byrow=T)
  #target of estimation: trt effect at mean values of modifiers in primary study, not par$muz[[1]][1:length(par$beta_m)] but result[,4:(4+length(par$beta_m))]
  cvrg <- is.between(as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)), result[,1], result[,2])
  pwr <- is.between(0, result[,1], result[,2])
  ## Summarize the simulation results
  summ[i,1] <- sum(cvrg==T)/nsim
  summ[i,2] <- mean(result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))
  summ[i,3] <- sum((result[,3] - as.numeric(par$beta_t[1] + par$muz[[1]][1:length(par$beta_m)]%*%as.matrix(par$beta_m)))**2)/nsim
  summ[i,4] <- sum(pwr==F)/nsim
  #wt <- cbind(wt, result[,4])
  #sdratio <- c(sdratio, mean(result[,8]/result[,9]))
  
  sink('output2not.txt')
  print(summ[1:i,])
  setTxtProgressBar(pb, i)
  sink()
})
close(pb)

summsc4 <- summ
save(summsc4, file = "summsc4.RData")

save.image(file='scenario4.RData')
