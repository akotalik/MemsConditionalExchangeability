#Illustration og IPCW bagging code
#Application to KMsurv datasets:

#load necessary libraries
library(data.table)
library(dplyr)
library(boot)
library(survival) # for Kaplan-Meier estimators
library(KMsurv)
library(gam) # for generalized additive models

# number of cores for multicore boot
cores=2


data(kidtran)
SURVTIME=5

# Define the outcome of interest (time to Hard CVD, Soft CVD, etc.) 
kidtran <- mutate(kidtran, T.use = time, C.use = delta, gender_f = as.factor(gender), race_f = as.factor(race), agegender = (kidtran$gender-1)*kidtran$age,
                  eventsurv.ind = as.factor(C.use == 1 & T.use <= (SURVTIME*365)))


#  Features to use in analysis
varnms_kidtran <- c("gender_f", "age", "race_f", "agegender")

#Get IPC Weights- Cox:
fmla= as.formula(paste0("Surv(T.use, 1-C.use)~",paste0(varnms_kidtran,collapse="+")))
cox.C <- coxph(fmla, data=kidtran)

newdata_zero=kidtran[1,]
newdata_zero$gender=newdata_zero$race=as.factor(1)
newdata_zero$age=0
res_zero=survfit(cox.C, newdata=newdata_zero)
beta = cox.C$coef

kidtran2 <- mutate(kidtran, 
                   gender2 = ifelse(gender == 2, 1, 0),
                   race2 = ifelse(race == 2, 1, 0),
                   nonzero=(1-(T.use < (SURVTIME*365) & C.use==0))
                   
)

kidtran2=as.data.frame(kidtran2)

wts=NULL
for (i in 1:nrow(kidtran2)){
  newdata=kidtran2[i,]
  newdata_cov = kidtran2[i, c("gender2", "age", "race2")]
  res=res_zero$surv ^(exp(sum(newdata_cov * beta )))
  wts=c(wts, newdata$nonzero / (res[min(which(res_zero$time > min(newdata$T.use, SURVTIME*365)))]) )
}

kidtran$wts=wts
kidtran2=NULL

set.seed(1281985)
frac.train <- 0.75
train.set <- sample(1:nrow(kidtran), floor(frac.train*nrow(kidtran)), replace=FALSE)
test.set <- setdiff(1:nrow(kidtran), train.set)

kidtran.train <- data.frame(kidtran[train.set,])
kidtran.test <- data.frame(kidtran[test.set,])



######################################################################
##########Functions###################################################
######################################################################

#example function:

#logistic regression backward selection where you can specify no. of steps (for crossvalidation):
backselection = function(sampleddata,steps,fmla, cens_method) {
  if(cens_method!="NATIVE") {
    fit <- glm(fmla, data = sampleddata,family = "binomial")
  }
  else{
    wts <- sampleddata$wts
    fit <- glm(fmla, data = sampleddata,family = "binomial", weights=wts)
  }
  if(!is.null(steps)){
    if(steps%%1==0 & steps>0) {
      for (i in 1:steps) {
        out<-which.max(summary(fit)$coef[-1,grep("Pr(>|z|)", colnames(summary(fit)$coef))])
        if(grepl("+_f[[:digit:]]", x=names(out))==TRUE){names(out)<-substr(names(out), 1, nchar(names(out))-1)}
        currentvars<-attr(fit$terms, "term.labels")
        newvars<-currentvars[!is.element(currentvars, names(out))]
        newfmla<-as.formula(paste0("eventsurv.ind~",paste0(newvars,collapse="+")))
        if(cens_method!="NATIVE") {
          fit <- glm(newfmla, data = sampleddata,family = "binomial")
        }
        else{
          wts <- sampleddata$wts
          fit <- glm(newfmla, data = sampleddata,family = "binomial", weights=wts)
        }
      }
    }
  }
  return(fit)
}

#in general, only need a function that fits a model and produces predictions:

#function that does logistic regression of y on x and returns predictions
logfun = function(data,fmla,testdata,pars,i,cens_method) {
  #pars are desired parameters
  #data is the training set
  #testdata is the test set
  #i is a vector that boot will pass (sampled rows)
  chosen <- as.data.frame(data[i, ])
  fit<-backselection(sampleddata=chosen, steps=pars$step, fmla=fmla, cens_method=cens_method)
  if(is.null(testdata)) {
    predicted <- predict(fit, newdata=data, type = "response", na.action=na.omit)    
  } else {
    predicted <- predict(fit, newdata=testdata, type = "response", na.action=na.omit)    
  }
  #returns survival probability, i.e. 1-risk
  return(1-predicted)
}


#function that does GAM and returns predictions
GAMfun = function(data,fmla,testdata,pars,i,cens_method) {
  #pars are parameters
  #data is the training set
  #testdata is the test set
  #i is a vector that boot will pass (sampled rows)
  
  #update formula with df:
  #grep( "s\\(.*\\)", attr(terms(fmla), "term.labels") )
  newterms <- gsub( "\\)", paste(" ,df=", pars$df, ")"), attr(terms(fmla), "term.labels"))
  newfmla <- reformulate(newterms, fmla[[2]])
  
  chosen <- data[i, ]
  if(cens_method!="NATIVE") {
    fit <- gam(formula=newfmla, data = chosen, family = binomial)
  }
  else{
    wts <- chosen$wts
    fit <- gam(formula=newfmla, data = chosen, family = binomial, weights=wts)
  }
  if(is.null(testdata)) {
    predicted <- predict(fit, newdata=data, type = "response", na.action=na.omit) 
  } else {
    predicted <- predict(fit, newdata=testdata, type = "response", na.action=na.omit) 
  }
  return(1-predicted)
}


##########################################################################################
########################k-fold crossvalidation fn  ################################################
##########################################################################################
crossval_function <- function(fun, parts, calibrationpars, formula, B=NULL, data, cens_method) {
  data$index<-as.integer(runif(nrow(data),1,parts+1))
  data=data.table(data, key="index")
  result <- vector("list", length(calibrationpars))
  # create progress bar
  pb= txtProgressBar(min = 0, max = parts, style = 3, char=":)")
  for (k in 1:parts) {
    test.set=data[J(k)]
    train.set=data[!J(k)]
    for (j in 1:length(calibrationpars)) {
      if(cens_method == "IPCW") {
        b <-boot(data=train.set, statistic=fun, cens_method=cens_method, R=B, fmla=formula, pars=calibrationpars[[j]], testdata=test.set, weights = train.set$wts, parallel="multicore", ncpus=cores)
        result[[j]]<-c(result[[j]], colMeans(b$t, na.rm=T))
      }
      
      else if(cens_method == "ZERO") {
        probs <- fun(data=train.set, fmla=formula, cens_method=cens_method, testdata=test.set, pars=calibrationpars[[j]], i=seq.int(1,nrow(train.set)))
        result[[j]]<-c(result[[j]], probs)
      }
      
      else if(cens_method == "DISCARD") {
        train.set.disc <- filter(train.set, !(T.use <= (SURVTIME*365) & C.use == 0))
        probs <- fun(data=train.set.disc, fmla=formula, cens_method=cens_method, testdata=test.set, pars=calibrationpars[[j]], i=seq.int(1,nrow(train.set.disc)))
        result[[j]]<-c(result[[j]], probs)
      }
      
      else if(cens_method == "NATIVE") {
        probs <- fun(data=train.set, fmla=formula, cens_method=cens_method, testdata=test.set, pars=calibrationpars[[j]], i=seq.int(1,nrow(train.set)))
        result[[j]]<-c(result[[j]], probs)
      }
      
      else{
        stop("Need to use censoring method IPCW, ZERO, NATIVE or DISCARD.")
      }
      
    }
    # update progress bar
    setTxtProgressBar(pb, k)
  }
  matrix = do.call(cbind, result)
  matrix = cbind(as.numeric(data$eventsurv.ind)-1, data$wts, matrix)
  matrix = as.data.frame(matrix)
  names(matrix) = c("eventsurv.ind", "wts", paste("Set", 1:(ncol(matrix)-2), sep=""))
  close(pb)
  return(matrix)
}

########################################################################
########################Prediction Function################################
########################################################################
pred_function <- function(parameters, formula, fun, B=NULL, cens_method = "IPCW", train.dat, test.dat) {
  result<-NULL
  if(cens_method == "IPCW") {
    # create progress bar
    pb= txtProgressBar(min = 0, max = length(parameters), style = 3, char=":)")
    for (j in 1:length(parameters)) {
      b <- boot(data=train.dat, statistic=fun, cens_method=cens_method, fmla=formula, R=B, pars=parameters[[j]], testdata=test.dat, weights = train.dat$wts, parallel="multicore", ncpus=cores)
      result <- cbind(result, colMeans(b$t, na.rm=T))
      # update progress bar
      setTxtProgressBar(pb, j)
    }
    close(pb)
    return(result)
  }
  else if(cens_method == "ZERO") {
    # create progress bar
    pb= txtProgressBar(min = 0, max = length(parameters), style = 3, char=":)")
    for (j in 1:length(parameters)) {
      probs <- fun(data=train.dat, fmla=formula, cens_method=cens_method, testdata=test.dat, pars=parameters[[j]], i=seq.int(1,nrow(train.dat)))
      result <- cbind(result, probs)
      # update progress bar
      setTxtProgressBar(pb, j)
    }
    close(pb)
    return(result)
  }
  else if(cens_method == "DISCARD") {
    # create progress bar
    pb= txtProgressBar(min = 0, max = length(parameters), style = 3, char=":)")
    train.dat.disc <- filter(train.dat, !(T.use <= (SURVTIME*365) & C.use == 0))
    for (j in 1:length(parameters)) {
      probs <- fun(data=train.dat.disc, fmla=formula, cens_method=cens_method, testdata=test.dat, pars=parameters[[j]], i=seq.int(1,nrow(train.dat.disc)))
      result <- cbind(result, probs)
      # update progress bar
      setTxtProgressBar(pb, j)
    }
    close(pb)
    return(result)
  }
  else if(cens_method == "NATIVE") {
    # create progress bar
    pb= txtProgressBar(min = 0, max = length(parameters), style = 3, char=":)")
    for (j in 1:length(parameters)) {
      probs <- fun(data=train.dat, fmla=formula, cens_method=cens_method, testdata=test.dat, pars=parameters[[j]], i=seq.int(1,nrow(train.dat)))
      result <- cbind(result, probs)
      # update progress bar
      setTxtProgressBar(pb, j)
    }
    close(pb)
    return(result)
  }
  else{
    stop("Need to use censoring method IPCW, ZERO, NATIVE or DISCARD.")
  }
}

##Function to get coefficients
getCoefficients <- function(res, cens_method) {
  for (i in 3:ncol(res)) {
    res[,i][res[,i] == 1] <- 0.999
    res[,i][res[,i] == 0] <- 0.001
    res[,i]=log(res[,i]/(1-res[,i]))
  }
  
  if(cens_method=="NATIVE" | cens_method=="IPCW") {
    fit=glm(eventsurv.ind~.-1-wts, data=res, weights=res$wts, family="quasibinomial")
  }
  else if (cens_method=="DISCARD") {
    res.disc=res[res$wts!=0,]
    fit=glm(eventsurv.ind~.-1-wts, data=res.disc, family="quasibinomial")
  }
  else {
    fit=glm(eventsurv.ind~.-1-wts, data=res, family="quasibinomial")
  }
  coefficients <- coef(fit)
  coefficients[is.na(coefficients)] <- 0
  return(coef(fit))
}

#combining predictions

combine <- function(predictions, coefficients) {
  for (i in 1:ncol(predictions)) {
    predictions[,i][predictions[,i] == 1] <- 0.999
    predictions[,i][predictions[,i] == 0] <- 0.001
    predictions[,i]=log(predictions[,i]/(1-predictions[,i]))
  }
  a=as.matrix(predictions) %*% coefficients
  result=exp(a)/(1+exp(a))
  return(1-result)
}


########################################################################
####################### Wrapper Function ###############################
########################################################################
wrapperfn <- function(fmla, parameters, fun, cens_method, parts, traindata, testdata, B=NULL) {
  res<-crossval_function(fun=fun, cens_method="DISCARD", formula=fmla, parts=parts, data=traindata, calibrationpars=parameters)
  sse<-NULL
  for (i in 3:ncol(res)) {
    sse<-cbind(sse, sum(((1-res[,i])-as.logical(res$eventsurv.ind))**2))
  }
  threebest<-pars[which(sse %in% sort(sse)[1:3])]
  parameters<-list(threebest[[1]],threebest[[2]],threebest[[3]])
  
  res<-crossval_function(fun=fun, cens_method=cens_method, formula=fmla, parts=parts, data=traindata, calibrationpars=parameters, B=B)
  coef<-getCoefficients(res, cens_method=cens_method)
  pred <- pred_function(parameters=parameters, formula=fmla, fun=fun, B=B, cens_method = cens_method, train.dat=traindata, test.dat=testdata)
  predicted <- combine(pred, coef)
  return(predicted)
}

set.seed(1101985)
fmla= as.formula(paste0("eventsurv.ind~",paste0(varnms_kidtran,collapse="+")))
pars=list(list(step=0), list(step=1), list(step=2), list(step=3))
logipc=wrapperfn(fmla=fmla, parameters=pars, fun=logfun, cens_method = "IPCW", parts=4, traindata=kidtran.train, testdata=kidtran.test, B=10)
logdisc=wrapperfn(fmla=fmla, parameters=pars, fun=logfun, cens_method = "DISCARD", parts=4, traindata=kidtran.train, testdata=kidtran.test)
logzero=wrapperfn(fmla=fmla, parameters=pars, fun=logfun, cens_method = "ZERO", parts=4, traindata=kidtran.train, testdata=kidtran.test)
lognative=wrapperfn(fmla=fmla, parameters=pars, fun=logfun, cens_method = "NATIVE", parts=4, traindata=kidtran.train, testdata=kidtran.test)

fmla= as.formula("eventsurv.ind~gender_f + race_f + s(age) + s(agegender)")
pars=list(list(df=1),list(df=2),list(df=3),list(df=4),list(df=5))
gamipc=wrapperfn(fmla=fmla, parameters=pars, fun=GAMfun, cens_method = "IPCW", parts=4, traindata=kidtran.train, testdata=kidtran.test, B=10)
gamdisc=wrapperfn(fmla=fmla, parameters=pars, fun=GAMfun, cens_method = "DISCARD", parts=4, traindata=kidtran.train, testdata=kidtran.test)
gamzero=wrapperfn(fmla=fmla, parameters=pars, fun=GAMfun, cens_method = "ZERO", parts=4, traindata=kidtran.train, testdata=kidtran.test)
gamnative=wrapperfn(fmla=fmla, parameters=pars, fun=GAMfun, cens_method = "NATIVE", parts=4, traindata=kidtran.train, testdata=kidtran.test)



#results are predicted survival probabilities (1-risk) for test set
results<-cbind(logipc, logdisc, logzero, lognative, gamipc, gamdisc, gamzero, gamnative)
results<-as.data.frame(results)
names(results)<-c("logipc", "logdisc", "logzero", "lognative", "gamipc", "gamdisc", "gamzero", "gamnative")



####################################################################################################
###################################################Performance metrics#########################
####################################################################################################


#Code illustration only: results do not necessarily have any meaning

#SSE

brierscore<-apply(results, 2, function(predicted) sum(kidtran.test$wts*((1-predicted)-as.logical(kidtran.test$eventsurv.ind))**2)/length(kidtran.test$wts))


#calibration function:
calib.stat <- function(p,
                       T.test,
                       C.test,
                       cutpts,
                       t) {
  ## p = vector of predicted survival probabilities (1-risk) for test set
  ## T.test = Observation times for test set
  ## C.test = Event indicators for test set
  ## cutpts = Cut points for risk categories
  ## t = Time at which calibration is to be evaluated
  risk.class <- cut(p,cutpts,labels=FALSE)
  lev.stats <- sapply(1:(length(cutpts)-1),function(f) {
    ind <- which(risk.class==f)
    sf <- summary( survfit(Surv(T.test[ind],C.test[ind])~1) )
    ind.surv <- max(which(sf$time<=t))
    p.KM <- sf$surv[ind.surv]
    ##print(c(cutpts[f],p.KM,mean(p[ind],na.rm=TRUE),sqrt(S.KM$SKMV[ind.surv])))
    (mean(p[ind],na.rm=TRUE) - p.KM)^2/(sf$std.err[ind.surv])^2
  })
  
  return( sum(lev.stats,na.rm=TRUE) )
}


calib<-NULL
for (i in 1:ncol(results)){
  if(all(is.na(results[,i]))==FALSE)
  {
    c<-NULL
    ct<-NULL
    p<-NULL
    p<-results[,i]
    ct<-unique(round(quantile(p, probs=seq(0,1,0.1), names=F), digits=8))
    #ct<-quantile(p, probs=seq(0,1,0.1), names=F)
    c<-calib.stat(p=p, T.test=kidtran.test$T.use, C.test=kidtran.test$C.use, cutpts=ct, t=365*SURVTIME)
  }
  else {
    c<-NA
  }
  calib<-c(calib, c)
}
calib


#C-index:
## C-index for censored data due to Harrell
compute.cIndex <- function(pred,T,C,t) {
  Tt <- T
  Ct <- C
  Tt[T>t] <- t
  Ct[T>t&C==1] <- 0
  
  survConcordance(Surv(Tt,Ct)~pred)
}

cindex<-NULL
for (i in 1:ncol(results)){
  p<-1-results[,i]
  c<-compute.cIndex(p, kidtran.test$T.use, kidtran.test$C.use, 365*SURVTIME)$concordance
  cindex<-c(cindex, unname(c))
}
cindex

