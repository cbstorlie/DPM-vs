
library("survival")
library("mclust")
library("doMC")
library("gbm")
library("glmnet")
library("missForest")
library("MASS")
library(rjags)


###############################
##### Utility Functions #######
###############################

Rsq <- function(y,yhat)
  return(1-sum((y-yhat)^2)/sum((y-mean(y))^2))



run.jags.par <- function(model.str, data, n.chains, n.adapt, n.burn, n.iter, thin, inits, variable.names){

  samples.list <- foreach(l=1:n.chains)%dopar%{
    init0 <- inits[[l]]
    init0$.RNG.name <- "base::Wichmann-Hill"
    init0$.RNG.seed <- 100+l
    bc.model <- jags.model(textConnection(model.str), data=data, n.chains=1, n.adapt=n.adapt, inits=init0)
    update(bc.model, n.burn)
    bc.samples <- coda.samples(bc.model, variable.names=variable.names, n.iter=n.iter, thin=thin)
  }
  all.samples <- mcmc.list()
  for(l in 1:n.chains)
     all.samples[[l]] <- samples.list[[l]][[1]]  
  return(all.samples)
}



my.summary <- function(ans.mcmc){

  foo <- summary(ans.mcmc)
  Rhat <- gelman.diag(ans.mcmc, multivariate=FALSE)$psrf[,1]
  n.eff <- effectiveSize(ans.mcmc)

  foo$statistics <- cbind(foo$statistics, Rhat, n.eff)
  return(foo)
}

plot.mcmc <- function(ans.mcmc, vnames=NULL, ask=TRUE){

  if(is.null(vnames))
    vnames <- colnames(ans.mcmc[[1]])
  n.chains <- length(ans.mcmc)
  sub.samples <- mcmc.list()
  for(l in 1:n.chains)
     sub.samples[[l]] <- ans.mcmc[[l]][,grep(vnames, colnames(ans.mcmc[[1]]))]
  plot(sub.samples, ask=ask)
}




logit <- function(p)
  return(log(p)-log(1-p))



inv.logit <- function(x)
  return(exp(x)/(1+exp(x)))


plot.ROC <- function(ans.ROC, xlim=NULL, ylim=NULL, add=FALSE, col=1, lty=1, lwd=1.5, smooth=FALSE){
  if(is.null(xlim))
    xlim <- c(0,1)
  if(is.null(ylim))
    ylim <- c(0,1)
  acc1 <- ans.ROC$acc1
  if(smooth)
    acc1 <- predict(smooth.spline(1-ans.ROC$acc0, acc1), 1-ans.ROC$acc0)$y
  if(!add)
    plot(c(0,1),c(0,1), ylab="sensitivity", xlab="1-specificity", col=0, xlim=xlim, ylim=ylim)
  lines(1-ans.ROC$acc0, acc1, col=col, lwd=lwd, lty=lty)
  invisible()
}


plot.PPV <- function(ans.ROC, xlim=NULL, ylim=NULL, add=FALSE, col=1, lty=1, lwd=1.5, smooth=FALSE){
  if(is.null(xlim))
    xlim <- c(0,1)
  if(is.null(ylim))
    ylim <- c(0,1)
  ppv <- ans.ROC$ppv
  if(smooth)
    ppv <- predict(smooth.spline(1-ans.ROC$acc0, ppv), 1-ans.ROC$acc0)$y
  ng <- length(ppv)
  ppv <- ppv[ng:1]
  ppv[is.na(ppv)] <- 0
  for(i in 1:(ng-1)){
    if(ppv[i]<ppv[i+1])
      ppv[1:i][ppv[1:i]<ppv[i+1]] <- ppv[i+1]
  }
  ppv <- ppv[ng:1]
  if(!add)
    plot(c(0,1),c(0,1), ylab="PPV", xlab="1-specificity", col=0, xlim=xlim, ylim=ylim)
  lines(1-ans.ROC$acc0, ppv, col=col, lwd=lwd, lty=lty)
  invisible()
}



get.ROC <- function(y, score){

  n0 <- sum(y==0)
  n1 <- sum(y==1)
  n <- n0+n1
  score0 <- score[y==0]
  score1 <- score[y==1]
  thresh <- sort(c(unique(score)-1E-6, unique(score), unique(score)+1E-6))
  nt <- length(thresh)
  acc0 <- acc1 <- acc <- ppv <- npv <- rep(0,nt)
  for(t in 1:nt){
    acc0[t] <- mean(score0<=thresh[t])
    acc1[t] <- mean(score1>thresh[t])
    acc[t] <- (n0*acc0[t] + n1*acc1[t])/n
    ind0 <- score<thresh[t]
    ppv[t] <- mean(y[!ind0])
    npv[t] <- 1-mean(y[ind0])
  }
  ppv[nt] <- ppv[nt-1]
  npv[1] <- npv[2]
  maxacc <- max(acc)
  auc <- sum((acc1[-1]+acc1[-nt])/2*(acc0[-1]-acc0[-nt]))
  return(list(acc=maxacc, auc=auc, acc0=acc0, acc1=acc1, ppv=ppv, npv=npv))
}



get.ROC.boot <- function(y, score.mat, ns=1000, B=getDoParWorkers()){

  S <- ncol(score.mat)
  n <- length(y)
  blocks <- list()
  ind.now <- 1
  inc <- ceiling(ns/B)

  acc <- auc <- rep(0,S)
  for(s in 1:S){
    ansROC.s <- get.ROC(y, score.mat[,s])
    acc[s] <- ansROC.s$acc
    auc[s] <- ansROC.s$auc
  }


  for(b in 1:B){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, ns)
    ind.now <- ind.now+inc
    if(ind.now>ns){
      B <- b
      break
    }
  }
  accauc <- foreach(b=1:B, .combine=rbind)%dopar%{
    n.b <- length(blocks[[b]])
    acc.b <- matrix(0, n.b, S)
    auc.b <- matrix(0, n.b, S)
    for(k in 1:n.b){
      i <- blocks[[b]][k]
      ind.i <- sample(1:n, n, replace=TRUE)
      y.i <- y[ind.i]
      score.i <- score.mat[ind.i,]
      for(s in 1:S){
        ansROC.s <- get.ROC(y.i, score.i[,s])
	acc.b[k,s] <- ansROC.s$acc
	auc.b[k,s] <- ansROC.s$auc
      }
    }
    cbind(acc.b, auc.b)
  }
  acc.mat <- accauc[,1:S]
  auc.mat <- accauc[,(S+1):(2*S)]
  colnames(acc.mat) <- colnames(auc.mat) <- colnames(score.mat)

#  acc.comp <- auc.comp <- array(0,c(ns,S,S))
#  for(s in 1:(S-1)){
#    for(t in (s+1):S){
#      acc.comp[,s,t] <- acc[,s] - acc[,t]
#      acc.comp[,t,s] <- -acc.comp[,s,t]
#      auc.comp[,s,t] <- auc[,s] - auc[,t]
#      auc.comp[,t,s] <- -auc.comp[,s,t]
#    }
#  }
  return(list(acc=acc, auc=auc, acc.bs=acc.mat, auc.bs=auc.mat))
}



get.CIpval <- function(boot.mat, alpha=.05){

  S <- dim(boot.mat)[2]
  results <- matrix(0,S,2)
  pw.results <- array(1, c(S,S,3))
  for(s in 1:S)
    results[s,1:2] <- quantile(boot.mat[,s],prob=c(alpha/2, 1-alpha/2))
  
  for(s in 1:(S-1)){
    for(t in (s+1):S){
      diff.st <- boot.mat[,s]-boot.mat[,t]
      pw.results[s,t,1:2] <- quantile(diff.st,prob=c(alpha/2, 1-alpha/2))
      pw.results[t,s,1:2] <- quantile(-diff.st,prob=c(alpha/2, 1-alpha/2))
      pw.results[s,t,3] <- pw.results[t,s,3] <- 2*min(mean(diff.st < 0), mean(-diff.st < 0))
    }
  }
  dimnames(pw.results)[[1]] <- dimnames(pw.results)[[2]] <- score.names
  dimnames(results)[[1]] <- score.names
  dimnames(results)[[2]] <- paste(100*(1-alpha), "%", c(" LB", " UB"), sep="")
  return(list(CI=results, pw.results=pw.results))
}



get.foldid <- function(n, nfolds=10, seed=220){

  replace.seed <- T
  if(missing(seed))
    replace.seed <- F

  if(replace.seed){
   ## set seed to specified value
    if(!any(ls(name='.GlobalEnv', all.names=T)=='.Random.seed')){
      set.seed(1)
    }
    save.seed <- .Random.seed
    set.seed(seed)  
  }

  perm <- sample(1:n, n)
  n.cv <- rep(floor(n/nfolds),nfolds)
  rem <- n - n.cv[1]*nfolds
  n.cv[index(1,rem)] <- n.cv[index(1,rem)]+1
  foldid <- rep(0,n)
  ind2 <- 0
  
  for(i in 1:nfolds){
    ind1 <- ind2+1
    ind2 <- ind2+n.cv[i]
    foldid[perm[ind1:ind2]] <- i
  }
  if(replace.seed){
   ## restore random seed to previous value
    .Random.seed <<- save.seed
  }
  return(foldid)
}






makeprobs <- function(v){ 
   #compute the stick-breaking weights
   M <- length(v)
   if(M==1)
     return(1)
   probs <- 1-v
   probs[2:M] <- probs[2:M]*cumprod(v[2:M-1])
return(probs)}



probs2v <- function(probs){ 
   #compute the stick-breaking weights
   M <- length(probs)
   v <- probs
   v[2:M] <- probs[2:M]/(1-cumsum(probs[2:M-1]))
   v[v>=1] <- 1-1E-16
   v[probs==0] <- .5
   return(1-v)
}




## Get ‘nsample' samples from the prior distribution of the number of clusters (M) in a
## data set of size ’n' that would be generated by a DP(delta), with delta ~ gamma(‘A','B')

sampleM <- function(A,B,n, nsamp=1000){
  Mmax <- min(n,500)
  Msamp <- rep(0,nsamp)
  for(i in 1:nsamp){
    delta <- rgamma(1,A,B)
    v <- c(rbeta(Mmax-1,delta,1),0)
    omega <- makeprobs(v)
    Msamp[i] <- length(unique(sample(1:Mmax, prob=omega, replace=TRUE)))
  }
  return(Msamp)
}




myrtruncnorm <- function(n, mu, sd, T, greater=TRUE){

  ans <- rep(0,n)
  if(length(greater)==1)
    greater <- rep(greater,n)
  less <- !greater
  if(length(mu)==1)
    mu <- rep(mu,n)
  if(length(sd)==1)
    sd <- rep(sd,n)
  if(length(T)==1)
    T <- rep(T,n)
  if(any(greater)){
    ng <- sum(greater)
    lp <- log(runif(ng))+pnorm(T[greater],mu[greater],sd[greater],lower.tail=FALSE,log.p=TRUE)
    ans[greater] <- qnorm(lp,mu[greater],sd[greater],log.p=TRUE,lower.tail=FALSE)
  }
  if(any(less)){
    nl <- sum(less)
    lp <- log(runif(nl))+pnorm(T[less],mu[less],sd[less],lower.tail=TRUE,log.p=TRUE)
    ans[less] <- qnorm(lp,mu[less],sd[less],lower.tail=TRUE,log.p=TRUE)
  }
  return(ans)
}


myrtruncnorm.old <- function(n, mu, sd, T, greater=TRUE){

  if(greater)
    x <- 1-runif(n)*(1-pnorm(T,mu,sd))
  else
    x <- runif(n)*pnorm(T,mu,sd)
  return(qnorm(x,mu,sd))
}



index <- function(m,n){
  if(m<=n) return(m:n)
  else return(numeric(0))
}

which.equal <- function(x, y){

  n <- length(x)
  ans <- rep(0,n)
  for(i in 1:n){
    ans[i] <- any(approx.equal(y,x[i]))
  }
  return(as.logical(ans))
}

approx.equal <- function(x, y, tol=1E-9){

  return(abs(x - y) < tol)
}


rwish <- function(v, S){
  S <- as.matrix(S)
  if (v < nrow(S)){
      stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if(p > 1){
    pseq <- 1:(p - 1)
    Z[rep(p*pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p*(p-1)/2)
  }
  return(crossprod(Z%*%CC))
}



riwish <- function(v, S){
  Sinv <- ginv.gp(S)$inv
  W <- rwish(v, Sinv)
  Winv <- ginv.gp(W)$inv
  return(Winv)
}



dwish <- function (W, v, S) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (!is.matrix(W)) 
        W <- matrix(W)
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in  dwish()\n\n")
    }
    k <- nrow(S)
    gammapart <- 0
    for (i in 1:k) {
        gammapart <- gammapart + lgamma((v + 1 - i)/2)
    }
    logdenom <- gammapart + (v * k/2)*log(2) + (k * (k - 1)/4)*log(pi)
    ans.S <- ginv.gp(S)
    ans.W <- ginv.gp(W)
    logdetS <- ans.S$log.det
    logdetW <- ans.W$log.det
    hold <- ans.S$inv%*%W
    tracehold <- sum(hold[row(hold) == col(hold)])
    lognum <- (-v/2)*logdetS + ((v - k - 1)/2)*logdetW - 1/2*tracehold
    return(lognum-logdenom)
}



diwish <- function (W, v, S, Winv=NULL, logdetW=NULL, logdetS=NULL) 
{
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (!is.matrix(W)) 
        W <- matrix(W)
    if (v < nrow(S)) {
        stop("v is less than the dimension of S in  diwish().\n")
    }
    k <- nrow(S)
    gammapart <- 0
    for (i in 1:k) {
        gammapart <- gammapart + lgamma((v + 1 - i)/2)
    }
    logdenom <- gammapart + (v * k/2)*log(2) + (k * (k - 1)/4)*log(pi)
    if(is.null(logdetS)){
      ans.S <- ginv.gp(S)
      logdetS <- ans.S$log.det
    }    
    if(is.null(logdetW) || is.null(Winv)){
      ans.W <- ginv.gp(W)
      logdetW <- ans.W$log.det
      Winv <- ans.W$inv
    }
    hold <- S%*%Winv
    tracehold <- sum(hold[row(hold) == col(hold)])
    lognum <- (v/2)*logdetS - (v + k + 1)/2*logdetW -1/2*tracehold
    return(lognum - logdenom)
}




gamma2 <- function(x)
  return(digamma(x)*gamma(x))

gamma3 <- function(x)
  return(gamma(x)*(trigamma(x)+1/(gamma(x)^2)*gamma2(x)^2))


ginv.gp.old <- function (X, eps = 1e-12){
    if(any(X==Inf))
      return(list(inv = diag(nrow(X)), log.det = Inf))
    eig.X <- eigen(X, symmetric = TRUE)
    P <- eig.X[[2]]
    lambda <- eig.X[[1]]
    ind <- lambda > eps
    lambda.inv <- lambda.sqrt <- lambda
    lambda.inv[ind] <- 1/lambda[ind]
    lambda.sqrt[ind] <- sqrt(lambda[ind])
    lambda.inv[!ind] <- lambda.sqrt[!ind] <- 0
    Q <- P%*%sqrt(diag(lambda.inv, nrow=length(lambda)))
    inv <- Q%*%t(Q)
    sqrt <- P %*% diag(lambda.sqrt, nrow=length(lambda)) %*% t(P)
    sqrt.inv <- P %*% diag(sqrt(lambda.inv), nrow=length(lambda)) %*% t(P)
    if(all(lambda>0))
      log.det <- sum(log(lambda))
    else
      log.det <- -Inf
    return(list(inv = inv, log.det = log.det, sqrt=sqrt, sqrt.inv=sqrt.inv))
}


ginv.gp <- function (X, eps = 1e-12){
    if(any(X==Inf))
      return(list(inv = diag(nrow(X)), log.det = Inf))
    L <- chol(X)
    inv <- chol2inv(L)
    sqrt <- t(L)
    sqrt.inv <- backsolve(L,diag(1,nrow(L)))
    log.det <- as.numeric(determinant(L)$modulus)*2
    return(list(inv = inv, log.det = log.det, sqrt=sqrt, sqrt.inv=sqrt.inv))
}



rmvnorm <- function(n=1, mu=0, Sigma, S.sqrt=NA){

  if(is.na(S.sqrt[1])){
    ans <- ginv.gp(Sigma)
    S.sqrt <- ans$sqrt

#S.sqrt..<<-S.sqrt

  }
  p <- nrow(S.sqrt)
  if(length(mu)==1) 
    mu <- rep(mu,p)
  X <- matrix(0,n,p)
  for(i in 1:n)
    X[i,] <- as.numeric(S.sqrt%*%rnorm(length(mu),0,1) + mu)
  return(X[1:n,])
}


dmvnorm <- function(x, mu, Sigma, S.inv=NULL, log.det=NULL){

  if(is.null(S.inv[1]) || is.null(log.det)){
    ans <- ginv.gp(Sigma)
    S.inv <- ans$inv
    log.det <- ans$log.det
  }
  if(is.null(dim(x))){
    n <- 1
    p <- length(x)
    x <- matrix(x,1,p)
  }
  else{
    n <- nrow(x)
    p <- ncol(x)
  }
  dens <- -n*p/2*log(2*pi)-n/2*log.det
  for(i in 1:n)
    dens <- dens - 0.5*as.numeric(t(x[i,]-mu)%*%S.inv%*%(x[i,]-mu))
  return(dens)
}



rniw <- function(nu,lambda,P,eta){

  Q <- rwish(eta, ginv.gp(P)$inv)
  ans.Q <- ginv.gp(Q)
  Sigma <- ans.Q$inv
  S.sqrt <-  1/sqrt(lambda)*ans.Q$sqrt.inv
  mu <- rmvnorm(1,nu,S.sqrt=1/sqrt(lambda)*ans.Q$sqrt.inv)
  return(list(mu=mu,Sigma=Sigma,Q=Q))
}



rmatnorm <- function(M,U,V,U.sqrt=NULL,V.sqrt=NULL){

  n <- nrow(M)
  p <- ncol(M)
  if(is.null(U.sqrt))
    U.sqrt <- ginv.gp(U)$sqrt
  if(is.null(V.sqrt))
    V.sqrt <- ginv.gp(V)$sqrt
  Z <- matrix(rmvnorm(1,as.numeric(M),S.sqrt=V.sqrt%x%U.sqrt),n,p)
  return(Z)
}




dmvmixnorm <- function(z, omega, mu, S.inv, log.det, y, outcome, x, beta, kappa){
  p <- length(z)
  M <- nrow(mu)
  dens <- prob <- rep(0,M)
  for(m in 1:M)
    dens[m] <- log(omega[m])-p/2*log(2*pi)-1/2*log.det[m] - 0.5*as.numeric(t(z-mu[m,])%*%S.inv[m,,]%*%(z-mu[m,])) + get.1obs.like(y, outcome, x, beta[m,], kappa)

  for(m in 1:M){
    if(dens[m]== -Inf)
      prob[m] <- 0
    else
      prob[m] <- 1/(sum(exp(dens - dens[m])))
  }
  return(prob)
}



dmvmixnorm.probit <- function(z, omega, mu, S.inv, log.det, yc, x, beta){
  p <- length(z)
  M <- nrow(mu)
  dens <- prob <- rep(0,M)
  for(m in 1:M)
    dens[m] <- log(omega[m])-p/2*log(2*pi)-1/2*log.det[m] - 0.5*as.numeric(t(z-mu[m,])%*%S.inv[m,,]%*%(z-mu[m,])) + dnorm(yc, sum(x*beta), 1, log=TRUE)

  for(m in 1:M){
    if(dens[m]== -Inf)
      prob[m] <- 0
    else
      prob[m] <- 1/(sum(exp(dens - dens[m])))
  }
  return(prob)
}




dmvmixnorm_noy <- function(z, omega, mu, S.inv, log.det){
  p <- length(z)
  M <- nrow(mu)
  dens <- prob <- rep(0,M)
  for(m in 1:M)
    dens[m] <- log(omega[m])-p/2*log(2*pi)-1/2*log.det[m] - 0.5*as.numeric(t(z-mu[m,])%*%S.inv[m,,]%*%(z-mu[m,]))

  for(m in 1:M){
    if(dens[m]== -Inf)
      prob[m] <- 0
    else
      prob[m] <- 1/(sum(exp(dens - dens[m])))
  }
  return(prob)
}



get.best.phi <- function(x, bin.X, omega, muX, SigmaX){

  ind.obs <- bin.X==0 & !is.na(x[-1])
  if(sum(ind.obs)==0){
    prob <- omega
  }
  else{
    p <- sum(ind.obs)
    M <- nrow(muX)
    dens <- prob <- rep(0,M)
    for(m in 1:M){
      ans.inv.m <- ginv.gp(SigmaX[m,ind.obs,ind.obs])
      dens[m] <- log(omega[m])-p/2*log(2*pi)-1/2*ans.inv.m$log.det - 0.5*as.numeric(t(x[-1][ind.obs]-muX[m,ind.obs])%*%ans.inv.m$inv%*%(x[-1][ind.obs]-muX[m,ind.obs]))
    }
    for(m in 1:M){
      if(dens[m]== -Inf)
        prob[m] <- 0
      else
        prob[m] <- 1/(sum(exp(dens - dens[m])))
    }
  }
  return(which(prob==max(prob))[1])
}



get.kmeans.BIC <- function(ans.kmeans, Z){

  p <- ncol(Z)
  phi <- ans.kmeans$cluster
  n <- length(phi)
  M <- max(phi)
  tab.phi <- table(phi.now)
  ord.tab <- order(as.numeric(names(tab.phi)))
  tab.phi <- tab.phi[ord.tab]
  omega <- as.numeric(tab.phi/n)
  mu <- matrix(0,M,p)
  Sigma <- array(0,c(M,p,p))
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    mu[m,] <- colMeans(Z[ind.m,])
    Sigma[m,,] <- cov(Z[ind.m,])
  }
  dens <- 0
  for(m in 1:M){
    dens.m <- try(dmvnorm(Z, mu[m,], Sigma[m,,]), silent=TRUE)
    if(is.character(dens.m[1])){
      dens <- -Inf
      break
    }
    else{
      dens <- dens + log(omega[m])+dens.m
    }
  }
  BIC <- dens - (p/2)*log(n)
  return(BIC)
}



kmeansBIC <- function(fit){

  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(D + log(n)*m*k)
}


mykmeans <- function(Z, Mvec){

  BIC <- foreach(m=Mvec, .combine=c)%dopar%{
    foo.m <- kmeans(Z,m)
    kmeansBIC(foo.m)
  }
  Mopt <- Mvec[which(BIC==min(BIC))[1]]
  fopt <- kmeans(Z,Mopt)
  return(fopt)
}


######################################################################
############# Factorize missingness and use glmnet ###################
######################################################################

convert2factors <- function(X, bin.X){

  Xf <- X
  Xd <- matrix(0,nrow(X),0)
  p <- ncol(X)
  for(j in 1:p){
    if(bin.X[j]==1){
      Xf[is.na(Xf[,j]),j] <- -1
      Xf[,j] <- Xf[,j]+1
    }
    else{
      cuts <- quantile(X[,j],prob=c(.33,.67), na.rm=TRUE)
      Xf[X[,j]<=cuts[1],j] <- 1
      Xf[X[,j]>cuts[1] & X[,j]<=cuts[2],j] <- 2
      Xf[X[,j]>cuts[2],j] <- 3
      Xf[is.na(X[,j]),j] <- 0
    }
  }
  Xf <- as.data.frame(Xf)
#  for(j in 1:p){
#    Xf[,j] <- as.factor(Xf[,j])
#    Fj <- outer(Xf[,j], levels(Xf[,j]), `==`)*1
#    Xd <- cbind(Xd,Fj)
#  }
  return(list(Xf=Xf, Xd=Xd))
}





convertframe2factors <- function(X){

  Xf <- X
  p <- ncol(X)
  for(j in 1:p){
    if(any(is.na(Xf[,j]))){
      if(is.factor(Xf[,j])){
        levels(Xf[,j]) <- c(levels(Xf[,j]), "missing")
        Xf[is.na(Xf[,j]),j] <- "missing"
      }
      else{
        cuts <- quantile(X[,j],prob=c(.33,.67), na.rm=TRUE)
        blah <- rep("missing", n)
        blah <- as.factor(blah)
        levels(blah) <- c("1", "2", "3", "missing")
        blah[X[,j]<=cuts[1]] <- "1"
        blah[X[,j]>cuts[1] & X[,j]<=cuts[2]] <- "2"
        blah[X[,j]>cuts[2]] <- "3"
        blah[is.na(X[,j])] <- "missing"
        Xf[,j] <- blah
      }
    }
  }
  Xf <- as.data.frame(Xf)
  return(Xf)
}



add.missing.ind.cols <- function(X){

  Xf <- X
  xnames <- names(Xf)
  p <- ncol(X)
  for(j in 1:p){
    if(any(is.na(Xf[,j]))){
      if(is.factor(Xf[,j])){
        levels(Xf[,j]) <- c(levels(Xf[,j]), "missing")
        Xf[is.na(Xf[,j]),j] <- "missing"
      }
      else{
        Xf[is.na(Xf[,j]),j] <- 0
	Xf <- cbind(Xf, is.na(X[,j])*1)
        xnames <- c(xnames, paste(xnames[j],".miss",sep=""))
      }
    }
  }
  Xf <- as.data.frame(Xf)
  names(Xf) <- xnames
  return(Xf)
}






#############################################################
############# RF interpolation and GLMNET ###################
#############################################################

RFEN <- function(X, y, outcome=NULL, bin.X=NULL, maxiter=5, parallelize="no", family="cox", alpha=.5){

 # Fill in X with RF
  if(is.null(bin.X)){
    bin.X <- rep(0,ncol(X))
    for(j in 1:ncol(X))
      bin.X[j] <- is.factor(X[,j])
  }
  Xd <- as.data.frame(X)
  for(j in 1:ncol(X)){
    if(bin.X[j]==1)
      Xd[,j] <- as.factor(Xd[,j])
  }
  if(!is.null(outcome))
    Xd <- cbind(Xd,y,outcome)
  else
    Xd <- cbind(Xd,y)
  Xh <- missForest(Xd, strata=as.list(c(1:ncol(Xd))), maxiter=maxiter, parallelize=parallelize)$ximp
  Xh <- matrix(as.numeric(as.matrix(Xh)), nrow(Xh), ncol(Xh))[,1:ncol(X)]

 # Fit glmnet
  if(family=="cox"){
    yS <- cbind(y,outcome)
    colnames(yS) <- c("time", "status")
  }
  else
    yS <- y
  ans.fit <- cv.glmnet(Xh, yS, family=family, alpha=alpha)
  return(list(fit=ans.fit, X=X, Xh=Xh, y=y, outcome=outcome, bin.X=bin.X))
}



predict.RFEN <- function(obj, Xp, maxiter=5, K=1, parallelize="no"){

  np <- nrow(Xp)
  if(K > np)
    K <- np
  p <- ncol(Xp)
  X <- obj$X
  n <- nrow(X)
  bin.X <- obj$bin.X
  
 # Fill in Xp using both X and Xp
  Xd <- as.data.frame(rbind(Xp,X))
  for(j in 1:ncol(X)){
    if(bin.X[j]==1)
      Xd[,j] <- as.factor(Xd[,j])
  }
  ind.list <- list()
  ind.now <- 1
  inc <- ceiling(np/K)
  for(k in 1:K){
    ind.list[[k]] <- ind.now:(ind.now+inc-1)
    ind.now <- ind.now + inc
  }
  Xhp <- matrix(0,np,p)
  for(k in 1:K){
    np.k <- length(ind.list[[k]])
    ind.k <- c(ind.list[[k]],(np+1):(np+n))
    Xd.k <- Xd[ind.k,]
    Xh.k <- missForest(Xd.k, strata=as.list(c(1:ncol(X))), maxiter=maxiter, parallelize=parallelize)$ximp
    Xhp[ind.list[[k]],] <- matrix(as.numeric(as.matrix(Xh.k[1:np.k,])), np.k, p)
  }
  lamhat.RFEN <- as.numeric(predict(obj$fit, Xhp))
  return(lamhat.RFEN)
}







##################################################################
############# Tuning by CV concordance for GBM ###################
##################################################################


cv.gbm <- function(X, y, outcome=NULL, n.trees=seq(150,500,by=10), interaction.depth, n.minobsinnode = 5, shrinkage=10^(seq(-3,-1,length=10)), bag.fraction = .5, distribution="bernoulli", foldid=NULL, nfolds=10, seed=220, verbose=FALSE){

  n <- nrow(X)
  if(is.null(foldid))
    foldid <- get.foldid(n,10, seed=seed)
  nfolds <- max(foldid)
  ns <- length(shrinkage)
  nnt <- length(n.trees)
  concord.mat <- foreach(s=1:ns, .combine=rbind)%dopar%{
    concord.s <- rep(0,nnt)
    yhat <- matrix(0,n,nnt)
    for(k in 1:nfolds){
      X.k <- X[foldid==k,,drop=F]
      y.k <- y[foldid==k]
      X.mk <- X[foldid!=k,,drop=F]
      y.mk <- y[foldid!=k]
      if(distribution=="poisson"){
        outcome.mk <- outcome[foldid!=k]
        ans.sk <- gbm.fit(X.mk, outcome.mk, offset=y.mk, n.trees=n.trees[nnt], interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, shrinkage=shrinkage[s], bag.fraction=bag.fraction, distribution=distribution, verbose=verbose)
      }
      else{
        ans.sk <- gbm.fit(X.mk, y.mk, n.trees=n.trees[nnt], interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, shrinkage=shrinkage[s], bag.fraction=bag.fraction, distribution=distribution, verbose=verbose)
      }
      if(distribution=="multinomial"){
        for(t in 1:nnt){
	  foo <- predict(ans.sk, X.k, n.trees=n.trees[t],type="response", verbose=verbose)
	  fuh <- foo[cbind(1:length(y.k),as.numeric(y.k),1)]
	  fuh[is.na(fuh)] <- 0
          yhat[foldid==k,t] <- fuh
	}
      }
      else{
        for(t in 1:nnt)
          yhat[foldid==k,t] <- predict(ans.sk, X.k, n.trees=n.trees[t])
      }
    }
    for(t in 1:nnt){
      if(distribution=="coxph")
        concord.s[t] <- survConcordance(y~yhat[,t])$concord
      else if(distribution=="poisson")
        concord.s[t] <- survConcordance(Surv(y,outcome)~yhat[,t])$concord
      else if(distribution=="multinomial")
        concord.s[t] <- sum(log(yhat[,t]))
	#        concord.s[t] <- mean(y==yhat[,t])
      else if(distribution=="gaussian")
        concord.s[t] <- 1-sum((y-yhat[,t])^2)/sum((y-mean(y))^2)
      else
        concord.s[t] <- get.ROC(y, yhat[,t])$auc
    }
    concord.s
  }

  inds <- which(concord.mat==max(concord.mat), arr.ind=TRUE)
  opt.shrink <- shrinkage[inds[1]]
  opt.n.trees <- n.trees[inds[2]]
  opt.concord <- concord.mat[inds[1], inds[2]]
  return(list(shrinkage=opt.shrink, n.trees=opt.n.trees, concord=opt.concord))
}



cv.gbm.serial <- function(X, y, n.trees=seq(150,500,by=10), interaction.depth, n.minobsinnode = 5, shrinkage=10^(seq(-3,-1,length=10)), bag.fraction = .5, distribution="bernoulli", foldid=NULL, nfolds=10, seed=220){

  n <- nrow(X)
  if(is.null(foldid))
    foldid <- get.foldid(n,10, seed=seed)
  nfolds <- max(foldid)
  ns <- length(shrinkage)
  nnt <- length(n.trees)
  concord.mat <- matrix(0, ns, nnt)
  for(s in 1:ns){
    yhat <- matrix(0,n,nnt)
    for(k in 1:nfolds){

print(s)
print(k)

      X.k <- X[foldid==k,]
      y.k <- y[foldid==k]
      X.mk <- X[foldid!=k,]
      y.mk <- y[foldid!=k]
      ans.sk <- gbm.fit(X.mk, y.mk, n.trees=n.trees[nnt], interaction.depth=interaction.depth, n.minobsinnode=n.minobsinnode, shrinkage=shrinkage[s], bag.fraction=bag.fraction, distribution=distribution)
      for(t in 1:nnt)
        yhat[foldid==k,t] <- exp(predict(ans.sk, X.k, n.trees=n.trees[t]))
    }
    for(t in 1:nnt){
      if(distribution=="coxph")
        concord.mat[s,t] <- survConcordance(y~yhat[,t])$concord
      else
        concord.mat[s,t] <- get.ROC(y, yhat[,t])$auc
    }
  }
  inds <- which(concord.mat==max(concord.mat), arr.ind=TRUE)
  opt.shrink <- shrinkage[inds[1]]
  opt.n.trees <- n.trees[inds[2]]
  opt.concord <- concord.mat[inds[1], inds[2]]
  return(list(shrinkage=opt.shrink, n.trees=opt.n.trees, concord=opt.concord))
}





##################################################################
################### MCMC Updating Functions ######################
##################################################################


#### Positive proposal functions for rho and kappa ####

dlogt <- function(x, df=1, sigma=1, mu=0){
  num <- gamma((df+1)/2)
  denom <- x*sigma*sqrt(pi*df)*gamma(df/2)*(1+1/df*((log(x)-mu)/sigma)^2)^((df+1)/2)
  return(num/denom)
}

ldlogt <- function(x, df=1, sigma=1, mu=0){
  num <- lgamma((df+1)/2)
  denom <- log(x)+log(sigma)+.5*log(pi*df)+lgamma(df/2)+((df+1)/2)*log(1+1/df*((log(x)-mu)/sigma)^2)
  return(num-denom)
}


rlogt <- function(n, df=1, sigma=1, mu=0){
  return(exp(rt(n,df)*sigma+mu))
}


r.pos.proposal <- function(prev, df, sigma){
  return(rlogt(1, df=df, sigma=sigma, mu=log(prev)))
}

d.pos.proposal <- function(prop, prev, df, sigma){
  return(ldlogt(prop, df=df, sigma=sigma, mu=log(prev)))
}




get.like <- function(y, outcome, X, beta.now, kappa.now, phi.now, lambda.now=NULL){

  n <- length(y)
  if(is.null(lambda.now)){
    lambda.now <- rep(0,n)
    M <- nrow(beta.now)
    for(m in 1:M){
      ind.m <- which(phi.now==m)
      lambda.now[ind.m] <- exp(as.numeric(X[ind.m,]%*%beta.now[m,]))
    }
  }
  ind.cens <- (outcome==0)
  ans <- 0
  if(any(ind.cens))
    ans <- ans + sum(pweibull(y[ind.cens], kappa.now, 1/lambda.now[ind.cens], lower.tail=FALSE, log.p=TRUE))
  if(any(!ind.cens))
    ans <- ans + sum(dweibull(y[!ind.cens], kappa.now, 1/lambda.now[!ind.cens], log=TRUE))
  return(ans)
}



get.like.probit <- function(y, X, beta.now, phi.now){

  n <- length(y)
  p.now <- rep(0,n)
  M <- nrow(beta.now)
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    p.now[ind.m] <- pnorm(as.numeric(X[ind.m,]%*%beta.now[m,]))
  }
  ans <- sum(dbinom(y,1,p.now, log=TRUE))
  return(ans)
}



get.1obs.like <- function(y, outcome, x, beta.now, kappa.now){

  lambda.now <- exp(sum(beta.now*x))
  if(outcome==0)
    ans <- pweibull(y, kappa.now, 1/lambda.now, lower.tail=FALSE, log.p=TRUE)
  else
    ans <- dweibull(y, kappa.now, 1/lambda.now, log=TRUE)
  return(ans)
}





update.y <- function(y, outcome, Xc.now, beta.now, phi.now, kappa.now){

  M <- nrow(beta.now)
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    outcome.m <- outcome[ind.m]
    if(length(ind.m)>0){
      ind.cens.m <- ind.m[outcome.m==0]
      if(length(ind.cens.m)>0){
        y.cens.m <- y[ind.cens.m]
        lambda.cens.m <- exp(as.numeric(Xc.now[ind.cens.m,]%*%beta.now[m,]))
        lprob <- log(runif(length(ind.cens.m)))+pweibull(y.cens.m, kappa.now, 1/lambda.cens.m, lower.tail=FALSE,log.p=TRUE)
        y.cens.m <- qweibull(lprob, kappa.now, 1/lambda.cens.m, lower.tail=FALSE,log.p=TRUE)
        y[ind.cens.m] <- y.cens.m
      }
    }
  }
  return(y)
}




update.y.probit <- function(y, Xc.now, beta.now, phi.now){

  M <- nrow(beta.now)
  yc <- rep(0,length(y))
  for(m in 1:M){
    ind.m <- which(phi.now==m & !is.na(y))
    if(length(ind.m)>0){
      n.m <- length(ind.m)
      mu.m <- as.numeric(Xc.now[ind.m,]%*%beta.now[m,])
      y.m <- y[ind.m]
      yc[ind.m] <- myrtruncnorm(n.m, mu.m, sd=1, T=0, greater=as.logical(y.m))
    }
    ind.m.na <- which(phi.now==m & is.na(y))
    if(length(ind.m.na)>0){
      n.m.na <- length(ind.m.na)
      mu.m.na <- as.numeric(Xc.now[ind.m.na,]%*%beta.now[m,])
      y.m.na <- y[ind.m.na]
      yc[ind.m.na] <- rnorm(n.m.na, mu.m.na, sd=1)
    }
  }
  return(yc)
}





update.kappa <- function(y, outcome, Xc.now, beta.now, kappa.now, phi.now, A.kappa, B.kappa, rho.kappa, prop.sd.kappa, kappa.prop.01, kappa.prop.10, dfp){

  accept <- 0
  B.k.now <- as.numeric(kappa.now!=1)
  prob.p.g.now <- ifelse(B.k.now==0, kappa.prop.01, 1-kappa.prop.10)
  B.k.p <- rbinom(1,1,prob.p.g.now)
  prob.now.g.p <- ifelse(B.k.p==0, kappa.prop.01, 1-kappa.prop.10)

  if(B.k.p==1)
    kappa.p <- r.pos.proposal(kappa.now, df=dfp, sigma=prop.sd.kappa)
  else
    kappa.p <- 1
  like.p <- get.like(y, outcome, Xc.now, beta.now, kappa.p, phi.now)
  like.now <- get.like(y, outcome, Xc.now, beta.now, kappa.now, phi.now)
  pi.now <- B.k.now*(log(rho.kappa)+ dgamma(kappa.now, A.kappa, B.kappa, log=TRUE)) + (1-B.k.now)*log(1-rho.kappa)
  pi.p <- B.k.p*(log(rho.kappa)+ dgamma(kappa.p, A.kappa, B.kappa, log=TRUE)) + (1-B.k.p)*log(1-rho.kappa)
  dprop.now <- B.k.p*(log(prob.p.g.now)+d.pos.proposal(kappa.now, kappa.p, df=dfp, sigma=prop.sd.kappa)) + (1-B.k.p)*log(1-prob.p.g.now)
  dprop.p <- B.k.now*(log(prob.now.g.p)+d.pos.proposal(kappa.p, kappa.now, df=dfp, sigma=prop.sd.kappa)) + (1-B.k.now)*log(1-prob.now.g.p)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    kappa.now <- kappa.p
    accept <- 1
  }
  return(list(kappa=kappa.now, accept=accept))
}




update.Beta <- function(beta.now, groups, tau2.now, one.tau, M.Beta, S2.Beta, one.beta){

  nG <- length(tau2.now)
  M <- nrow(beta.now)
  p <- ncol(beta.now)-1
  Beta.now <- rep(0,p+1)
  if(!one.beta){
    for(k in 1:nG){
      for(l in 1:length(groups[[k]])){
        j <- groups[[k]][l]
        s2.s <- 1/(1/S2.Beta + M/tau2.now[k])
        mu.s <- s2.s*(M.Beta/S2.Beta + sum(beta.now[,j]/tau2.now[k]))
        Beta.now[j] <- rnorm(1,mu.s,sqrt(s2.s))
      }
    }
  }
  return(Beta.now)
}



update.tau2 <- function(beta.now, groups, Beta.now, A.tau2, B.tau2, one.tau, one.beta){

  nG <- length(groups)
  p <- length(Beta.now)-1
  M <- nrow(beta.now)
  if(one.tau){
    indneq0 <- which(beta.now[1,]!=0)
    nneq0 <- length(indneq0)
    if(one.beta){
      A.s <- A.tau2 + nneq0/2
      B.s <- B.tau2 + sum(beta.now[1,indneq0]^2)/2      
    }
    else{
      A.s <- A.tau2 + M*nneq0/2
      B.s <- B.tau2 + sum((beta.now[,indneq0]-matrix(Beta.now[indneq0],M,nneq0,byrow=TRUE))^2)/2
    }
    tau2.now <- rep(1/rgamma(1,A.s,B.s), nG)
  }
  else{
    M <- nrow(beta.now)
    tau2.now <- rep(0,nG)
    for(k in 1:nG){
      indneq0 <- intersect(groups[[k]],which(beta.now[1,]!=0))
      nneq0 <- length(indneq0)
      if(one.beta){
        A.s <- A.tau2 + nneq0/2
        B.s <- B.tau2 + sum(beta.now[1,indneq0]^2)/2
      }
      else{
        A.s <- A.tau2 + nneq0*M/2
        B.s <- B.tau2 + sum((beta.now[,indneq0] - matrix(Beta.now[indneq0],M,nneq0,byrow=TRUE))^2)/2
      }
      tau2.now[k] <- 1/rgamma(1,A.s,B.s)
    }
  }
  return(tau2.now)
}





update.beta <- function(y, yc.now, outcome, Xc.now, beta.now, kappa.now, phi.now, Beta.now, tau2.now, groups, mult1, rho, rho.prop.01, rho.prop.10, one.beta){

  n <- length(y)
  ys <- -log(yc.now)
  nG <- length(groups)
  M <- Mf <- nrow(beta.now)
  if(one.beta){
    Mf <- 1
    phi.now <- rep(1,n)
  }
  accept <- rep(0,nG)
  lambda.now <- rep(0,n)
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    lambda.now[ind.m] <- exp(as.numeric(Xc.now[ind.m,]%*%beta.now[m,]))
  }
  for(k in 1:nG){
    sig2 <- 1/kappa.now^2*(gammapp.1 - gammap.1^2)*mult1[k]
    J.k <- length(groups[[k]])
    B.k.now <- as.numeric(beta.now[1,groups[[k]]][1]!=0)
    prob.p.g.now <- ifelse(B.k.now==0, rho.prop.01[k], 1-rho.prop.10[k])
    B.k.p <- rbinom(1,1,prob.p.g.now)
    prob.now.g.p <- ifelse(B.k.p==0, rho.prop.01[k], 1-rho.prop.10[k])
    Q0 <- diag(1/tau2.now[k], J.k)
    M0 <- Beta.now[groups[[k]]]
    Q0M0 <- Q0%*%M0
    dp.g.now <- dnow.g.p <- pi.p <- pi.now <- 0
    beta.p <- beta.now
    lambda.p <- lambda.now

    for(m in 1:Mf){
      ind.m <- which(phi.now==m)
      beta.k.now.m <- beta.now[m,groups[[k]]]
      if(length(ind.m) > 0){
        X.k.m <- Xc.now[ind.m,groups[[k]],drop=FALSE]
        log.lambda.mk.m <- as.numeric(log(lambda.now[ind.m]) - X.k.m%*%beta.k.now.m)
        yt.m <- as.numeric(ys[ind.m] + 1/kappa.now*gammap.1 - log.lambda.mk.m)
        Q.s.m <- as.matrix(1/sig2*(crossprod(X.k.m) + sig2*Q0))
        ans.inv <- ginv.gp(Q.s.m)
        S.s.m <- ans.inv$inv
        mu.s.m <- as.numeric(1/sig2*S.s.m%*%( sig2*Q0M0 + crossprod(X.k.m,yt.m) ))
      }
      else{
        Q.s.m <- Q0
        ans.inv <- ginv.gp(Q.s.m)
	mu.s.m <- M0
      }
      
      if(B.k.p==1)
        beta.k.p.m <- rmvnorm(1,mu=mu.s.m, S.sqrt=ans.inv$sqrt.inv)
      else
        beta.k.p.m <- rep(0,J.k)

      dp.g.now <- dp.g.now + B.k.p*(1/M*log(prob.p.g.now)+dmvnorm(beta.k.p.m, mu.s.m, S.inv=Q.s.m, log.det=-ans.inv$log.det)) + (1-B.k.p)*1/M*log(1-prob.p.g.now)
      dnow.g.p <- dnow.g.p + B.k.now*(1/M*log(prob.now.g.p)+dmvnorm(beta.k.now.m, mu.s.m, S.inv=Q.s.m, log.det=-ans.inv$log.det)) + (1-B.k.now)*1/M*log(1-prob.now.g.p)

      pi.p <- pi.p + B.k.p*(1/M*log(rho[k])+sum(dnorm(beta.k.p.m,M0,sqrt(tau2.now[k]),log=TRUE))) + (1-B.k.p)*1/M*log(1-rho[k])
      pi.now <- pi.now + B.k.now*(1/M*log(rho[k])+sum(dnorm(beta.k.now.m,M0,sqrt(tau2.now[k]),log=TRUE))) + (1-B.k.now)*1/M*log(1-rho[k])
      beta.p[m,groups[[k]]] <- beta.k.p.m
      lambda.p[ind.m] <- as.numeric(exp(log.lambda.mk.m + X.k.m%*%beta.k.p.m))
    }

    if(one.beta && M>1)
      beta.p[2:M,] <- matrix(beta.p[1,],M-2+1,ncol(beta.now),byrow=TRUE)

    like.p <- get.like(y, outcome, kappa.now=kappa.now, lambda.now=lambda.p)
    like.now <- get.like(y, outcome, kappa.now=kappa.now, lambda.now=lambda.now)
    MH.ratio <- exp(like.p+pi.p+dnow.g.p-like.now-pi.now-dp.g.now)
    if(runif(1) < MH.ratio){
      beta.now <- beta.p
      lambda.now <- lambda.p
      accept[k] <- 1
    }
  }
  return(list(beta=beta.now, accept=accept))
}




update.beta.probit <- function(yc.now, Xc.now, beta.now, phi.now, Beta.now, tau2.now, groups, mult1, rho, rho.prop.01, rho.prop.10, one.beta){

  n <- length(yc.now)
  nG <- length(groups)
  M <- Mf <- nrow(beta.now)
  if(one.beta){
    Mf <- 1
    phi.now <- rep(1,n)
  }
  accept <- rep(0,nG)
  mu.now <- rep(0,n)
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    mu.now[ind.m] <- as.numeric(Xc.now[ind.m,]%*%beta.now[m,])
  }
  for(k in 1:nG){
    sig2 <- 1
    J.k <- length(groups[[k]])
    B.k.now <- any(beta.now[1,groups[[k]]]!=0)*1
    prob.p.g.now <- ifelse(B.k.now==0, rho.prop.01[k], 1-rho.prop.10[k])
    B.k.p <- rbinom(1,1,prob.p.g.now)
    prob.now.g.p <- ifelse(B.k.p==0, rho.prop.01[k], 1-rho.prop.10[k])
    Q0 <- diag(1/tau2.now[k], J.k)
    M0 <- Beta.now[groups[[k]]]
    Q0M0 <- Q0%*%M0
    dp.g.now <- dnow.g.p <- pi.p <- pi.now <- 0
    beta.p <- beta.now
    mu.p <- mu.now

    for(m in 1:Mf){
      ind.m <- which(phi.now==m)
      beta.k.now.m <- beta.now[m,groups[[k]]]
      if(length(ind.m) > 0){
        X.k.m <- Xc.now[ind.m,groups[[k]],drop=FALSE]
        mu.mk.m <- as.numeric(mu.now[ind.m] - X.k.m%*%beta.k.now.m)
        yt.m <- yc.now[ind.m] - mu.mk.m
        Q.s.m <- as.matrix(1/sig2*(crossprod(X.k.m) + sig2*Q0))
        ans.inv <- ginv.gp(Q.s.m)
        S.s.m <- ans.inv$inv
        mu.s.m <- as.numeric(1/sig2*S.s.m%*%( sig2*Q0M0 + crossprod(X.k.m,yt.m) ))
      }
      else{
        Q.s.m <- Q0
        ans.inv <- ginv.gp(Q.s.m)
	mu.s.m <- M0
      }
      
      if(B.k.p==1)
        beta.k.p.m <- rmvnorm(1,mu=mu.s.m, S.sqrt=ans.inv$sqrt.inv)
      else
        beta.k.p.m <- rep(0,J.k)

      dp.g.now <- dp.g.now + B.k.p*(1/M*log(prob.p.g.now)+dmvnorm(beta.k.p.m, mu.s.m, S.inv=Q.s.m, log.det=-ans.inv$log.det)) + (1-B.k.p)*1/M*log(1-prob.p.g.now)
      dnow.g.p <- dnow.g.p + B.k.now*(1/M*log(prob.now.g.p)+dmvnorm(beta.k.now.m, mu.s.m, S.inv=Q.s.m, log.det=-ans.inv$log.det)) + (1-B.k.now)*1/M*log(1-prob.now.g.p)

      pi.p <- pi.p + B.k.p*(1/M*log(rho[k])+sum(dnorm(beta.k.p.m,M0,sqrt(tau2.now[k]),log=TRUE))) + (1-B.k.p)*1/M*log(1-rho[k])
      pi.now <- pi.now + B.k.now*(1/M*log(rho[k])+sum(dnorm(beta.k.now.m,M0,sqrt(tau2.now[k]),log=TRUE))) + (1-B.k.now)*1/M*log(1-rho[k])
      beta.p[m,groups[[k]]] <- beta.k.p.m
      mu.p[ind.m] <- as.numeric(mu.mk.m + X.k.m%*%beta.k.p.m)
    }

    if(one.beta && M>1)
      beta.p[2:M,] <- matrix(beta.p[1,],M-2+1,ncol(beta.now),byrow=TRUE)

    like.p <- sum(dnorm(yc.now, mu.p, 1, log=TRUE))
    like.now <- sum(dnorm(yc.now, mu.now, 1, log=TRUE))
    MH.ratio <- exp(like.p+pi.p+dnow.g.p-like.now-pi.now-dp.g.now)
    if(runif(1) < MH.ratio){
      beta.now <- beta.p
      mu.now <- mu.p
      accept[k] <- 1
    }
  }
  return(list(beta=beta.now, accept=accept))
}






initialize.XZ.kmeans <- function(X, bin.X, M, Mr=4){
  
  Xc.now <- X
  Z.now <- X[,-1]
  p <- ncol(X)-1
  n <- nrow(X)
  muX.now <- rep(0,p)
  SigmaX.now <- diag(1,p)
  for(j in 1:p){
    ind.miss.j <- which(is.na(X[,j+1]))
    if(bin.X[j]==1)
      Z.now[Xc.now[,j+1]==0,j] <- -1
    muX.now[j] <- mean(Z.now[,j], na.rm=TRUE)
    SigmaX.now[j,j] <- ifelse(bin.X[j]==1, 1, var(Z.now[,j], na.rm=TRUE))
    Z.now[ind.miss.j,j] <- rnorm(length(ind.miss.j), muX.now[j], sqrt(SigmaX.now[j,j]))
    if(bin.X[j]==1)
      Xc.now[,j+1] <- as.numeric(Z.now[,j]>0)
    else
      Xc.now[,j+1] <- Z.now[,j]
  }
 ## initial clustering
  ans.kmeans <- kmeans(Z.now,Mr)
  phi.now <- ans.kmeans$cluster
  ord.m <- order(-table(phi.now))
  phi.now <- match(phi.now, ord.m)
  muX.now <- matrix(0,M,p)
  SigmaX.now <- array(0,c(M,p,p))
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    if(length(ind.m)>0){
      muX.now[m,] <- colMeans(Z.now[ind.m,])
      vars.m <- apply(Z.now[ind.m,],2,var)
      vars.m[vars.m==0] <- min(vars.m[vars.m!=0])
      SigmaX.now[m,,] <- diag(vars.m)
    }
    else
      SigmaX.now[m,,] <- diag(rep(1,p))
  }
  return(list(Xc.now=Xc.now, Z.now=Z.now, muX.now=muX.now, SigmaX.now=SigmaX.now, phi.now=phi.now))
}






initialize.XZ <- function(X, y, outcome, bin.X, M, G=1:M, phi.init, beta.init, kappa.now, groups, rho, use.mclust=FALSE){

  Xc.now <- X
  Z.now <- X[,-1]
  p <- ncol(X)-1
  n <- nrow(X)
  muX.now <- rep(0,p)
  SigmaX.now <- diag(1,p)
  for(j in 1:p){
    ind.miss.j <- which(is.na(X[,j+1]))
    if(bin.X[j]==1){
      pj <- mean(X[,j+1]==1, na.rm=TRUE)
      muX.now[j] <- max(-3,min(qnorm(pj),3))
      ind0.j <- which(Xc.now[,j+1]==0)
      n0.j <- length(ind0.j)
      ind1.j <- which(Xc.now[,j+1]==1)
      n1.j <- length(ind1.j)
      Z.now[ind0.j,j] <- myrtruncnorm(n0.j, muX.now[j], 1, T=0, greater=FALSE)
      Z.now[ind1.j,j] <- myrtruncnorm(n1.j, muX.now[j], 1, T=0, greater=TRUE)
    }
    muX.now[j] <- mean(Z.now[,j], na.rm=TRUE)
    SigmaX.now[j,j] <- ifelse(bin.X[j]==1, 1, var(Z.now[,j], na.rm=TRUE))
    Z.now[ind.miss.j,j] <- rnorm(length(ind.miss.j), muX.now[j], sqrt(SigmaX.now[j,j]))
    if(bin.X[j]==1)
      Xc.now[,j+1] <- as.numeric(Z.now[,j]>0)
    else
      Xc.now[,j+1] <- Z.now[,j]
  }

 ## initial clustering
  if(is.null(phi.init)){
    if(use.mclust){
      foo2 <- Mclust(Z.now,G,modelNames="VVV")
      phi.now <- foo2$classification
    }
    else{
      ans.kmeans <- mykmeans(Z.now,G)
      phi.now <- ans.kmeans$cluster
    }
    ord.m <- order(-table(phi.now))
    phi.now <- match(phi.now, ord.m)
  }
  else{
    phi.now <- phi.init
  }
  K.foo <- max(phi.now)
  muX.now <- matrix(0,M,p)
  SigmaX.now <- array(0,c(M,p,p))
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    if(length(ind.m)>0)
      muX.now[m,] <- colMeans(Z.now[ind.m,,drop=FALSE])
    if(length(ind.m) >= 2*p){
      SigmaX.now[m,,] <- cov(Z.now[ind.m,])
      vars.m <- diag(SigmaX.now[m,,])
      vars.m[vars.m==0] <- min(vars.m[vars.m!=0])
      SigmaX.now[m,,] <- diag(vars.m)
    }
    else if(length(ind.m)>1){
      vars.m <- apply(Z.now[ind.m,],2,var)
      vars.m[vars.m==0] <- min(vars.m[vars.m!=0])
      SigmaX.now[m,,] <- diag(vars.m)
    }
    else
      SigmaX.now[m,,] <- diag(rep(1,p))
  }
  if(is.null(beta.init)){
    rhofoo <- rep(rho, sapply(groups, length))
    indkeep <- which(rhofoo > .0001)[-1] - 1
    yS <- cbind(y,outcome)
    yS[yS[,1]==0,] <- 1E-6
    colnames(yS) <- c("time", "status")
    beta.1 <- rep(0,p+1)
    beta.1[1] <- -mean(log(y[y>0])) + 1/kappa.now*gammap.1
    ans.fit <- cv.glmnet(Xc.now[,indkeep], yS, family="cox", alpha=.5, nfolds=5)
    ind <- which(ans.fit$glmnet.fit$lambda==ans.fit$lambda.1se)
    beta.1[indkeep+1] <- ans.fit$glmnet.fit$beta[,ind]
    beta.now <- matrix(beta.1, M, p+1, byrow=TRUE)
  }
  else{
    if(is.null(dim(beta.init)))
      beta.now <- matrix(beta.init, M, p+1, byrow=TRUE)
    else
      beta.now <- beta.init
  }  
  return(list(Xc.now=Xc.now, Z.now=Z.now, muX.now=muX.now, SigmaX.now=SigmaX.now, phi.now=phi.now, beta.now=beta.now))
}








initialize.XZ.probit <- function(X, y, bin.X, M, G=1:M, phi.init, beta.init, groups, rho, use.mclust=FALSE){

  Xc.now <- X
  Z.now <- X[,-1]
  p <- ncol(X)-1
  n <- nrow(X)
  muX.now <- rep(0,p)
  SigmaX.now <- diag(1,p)
  for(j in 1:p){
    ind.miss.j <- which(is.na(X[,j+1]))
    if(bin.X[j]==1){
      pj <- mean(X[,j+1]==1, na.rm=TRUE)
      muX.now[j] <- max(-3,min(qnorm(pj),3))
      ind0.j <- which(Xc.now[,j+1]==0)
      n0.j <- length(ind0.j)
      ind1.j <- which(Xc.now[,j+1]==1)
      n1.j <- length(ind1.j)
      Z.now[ind0.j,j] <- myrtruncnorm(n0.j, muX.now[j], 1, T=0, greater=FALSE)
      Z.now[ind1.j,j] <- myrtruncnorm(n1.j, muX.now[j], 1, T=0, greater=TRUE)
    }
    muX.now[j] <- mean(Z.now[,j], na.rm=TRUE)
    SigmaX.now[j,j] <- ifelse(bin.X[j]==1, 1, var(Z.now[,j], na.rm=TRUE))
    Z.now[ind.miss.j,j] <- rnorm(length(ind.miss.j), muX.now[j], sqrt(SigmaX.now[j,j]))
    if(bin.X[j]==1)
      Xc.now[,j+1] <- as.numeric(Z.now[,j]>0)
    else
      Xc.now[,j+1] <- Z.now[,j]
  }

 ## initial clustering
  if(is.null(phi.init)){
    if(use.mclust){
      foo2 <- Mclust(Z.now,G,modelNames="VVV")
      phi.now <- foo2$classification
    }
    else{
      ans.kmeans <- mykmeans(Z.now,G)
      phi.now <- ans.kmeans$cluster
    }
    ord.m <- order(-table(phi.now))
    phi.now <- match(phi.now, ord.m)
  }
  else{
    phi.now <- phi.init
  }
  K.foo <- max(phi.now)
  muX.now <- matrix(0,M,p)
  SigmaX.now <- array(0,c(M,p,p))
  for(m in 1:M){
    ind.m <- which(phi.now==m)
    if(length(ind.m)>0)
      muX.now[m,] <- colMeans(Z.now[ind.m,,drop=FALSE])
    if(length(ind.m) >= 2*p){
      SigmaX.now[m,,] <- cov(Z.now[ind.m,])
      vars.m <- diag(SigmaX.now[m,,])
      vars.m[vars.m==0] <- min(vars.m[vars.m!=0])
      SigmaX.now[m,,] <- diag(vars.m)
    }
    else if(length(ind.m)>1){
      vars.m <- apply(Z.now[ind.m,],2,var)
      vars.m[vars.m==0] <- min(vars.m[vars.m!=0])
      SigmaX.now[m,,] <- diag(vars.m)
    }
    else
      SigmaX.now[m,,] <- diag(rep(1,p))
  }
  if(is.null(beta.init)){
    ind.com <- which(!is.na(y))
    rhofoo <- rep(rho, sapply(groups, length))
    indkeep <- which(rhofoo > .0001)[-1] - 1
    ans.fit <- cv.glmnet(Xc.now[ind.com,indkeep], y[ind.com], family="binomial", alpha=.5, nfolds=5)
    ind <- which(ans.fit$glmnet.fit$lambda==ans.fit$lambda.1se)
    beta.1 <- rep(0,p+1)
    beta.1[1] <- ans.fit$glmnet.fit$a0[ind]
    beta.1[indkeep] <- ans.fit$glmnet.fit$beta[,ind]
    beta.now <- matrix(beta.1, M, p+1, byrow=TRUE)
  }
  else{
    if(is.null(dim(beta.init)))
      beta.now <- matrix(beta.init, M, p+1, byrow=TRUE)
    else
      beta.now <- beta.init
  }  
  return(list(Xc.now=Xc.now, Z.now=Z.now, muX.now=muX.now, SigmaX.now=SigmaX.now, phi.now=phi.now, beta.now=beta.now))
}







update.XZ <- function(X, Xc.now, Z.now, bin.X, muX.now, SigmaX.now, y, outcome, beta.now, kappa.now, mult2, omega.now, nb=getDoParWorkers()){

  n <- nrow(X)
  p <- ncol(X)-1
  M <- length(omega.now)
  ans.QX <- foreach(m=1:M)%dopar%{
    ans.inv <- ginv.gp(SigmaX.now[m,,])
    list(ans.inv$inv, ans.inv$log.det)
  }
  QX <- SigmaX.now
  ldet.vec <- rep(0,M)
  for(m in 1:M){
    QX[m,,] <- ans.QX[[m]][[1]]
    ldet.vec[m] <- ans.QX[[m]][[2]]
  }
  
## for each missing j, draw Zij (and effectively Xij) | rest, and accept/reject
  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/nb)
  for(b in 1:nb){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
  }

  XZ <- foreach(b=1:nb, .combine=rbind)%dopar%{
    Xc.now.b <- matrix(0, length(blocks[[b]]), p+1)
    Z.now.b <- matrix(0, length(blocks[[b]]), p)
    accept.b <- matrix(NA, length(blocks[[b]]), p)
    phi.now.b <- rep(0, length(blocks[[b]]))
    for(k in 1:length(blocks[[b]])){
      i <- blocks[[b]][k]
      
     ## Gibbs sample a new phi.now[i]
      xc.now.i <- Xc.now[i,]
      z.now.i <- Z.now[i,]
      probs.i <- dmvmixnorm(z.now.i, omega.now, muX.now, QX, ldet.vec, y[i], outcome[i], xc.now.i, beta.now, kappa.now)
      phi.i <- sample(1:M, 1, prob=probs.i)

     ## for each missing j, draw Xij, Zij | rest 
      ind.miss.i <- which(is.na(X[i,-1]))
      accept.i <- rep(0,p)
      if(length(ind.miss.i)>0){
        for(j in ind.miss.i){
          mu.miss.ij <- muX.now[phi.i,j] - 1/QX[phi.i,j,j]*sum(QX[phi.i,j,-j]*(z.now.i[-j] - muX.now[phi.i,-j]))
          z.miss.p <- rnorm(1, mu.miss.ij, sqrt(mult2/QX[phi.i,j,j]))
	  if(bin.X[j]==1)
	    xc.miss.p <- as.numeric(z.miss.p > 0)
	  else
	    xc.miss.p <- z.miss.p
          xc.p.i <- xc.now.i
          z.p.i <- z.now.i
          xc.p.i[j+1] <- xc.miss.p
          z.p.i[j] <- z.miss.p
	  if(mult2==1){
	    dZ.p <- dZ.now <- pi.p <- pi.now <- 0
	  }
	  else{
            dZ.p <- dnorm(z.miss.p, mu.miss.ij, sqrt(mult2/QX[phi.i,j,j]), log=TRUE)
            dZ.now <- dnorm(z.now.i[j], mu.miss.ij, sqrt(mult2/QX[phi.i,j,j]), log=TRUE)
            pi.p <- dnorm(z.miss.p, mu.miss.ij, sqrt(1/QX[phi.i,j,j]), log=TRUE)
            pi.now <- dnorm(z.now.i[j], mu.miss.ij, sqrt(1/QX[phi.i,j,j]), log=TRUE)
	  }
          likey.p <- get.1obs.like(y[i], outcome[i], xc.p.i, beta.now[phi.i,], kappa.now)
          likey.now <- get.1obs.like(y[i], outcome[i], xc.now.i, beta.now[phi.i,], kappa.now)
          MH.ratio <- exp(likey.p + pi.p - dZ.p - likey.now - pi.now + dZ.now)
          if(runif(1) < MH.ratio){
            xc.now.i <- xc.p.i
	    z.now.i <- z.p.i
            accept.i[j] <- 1
          }
	  else
	    accept.i[j] <- 0
        }
      }
    ## for each binary, non-missing j, draw Zij | rest 
      ind.bin.nomiss.i <- which(!is.na(X[i,-1]) & bin.X==1)
      if(length(ind.bin.nomiss.i)>0){
        for(j in ind.bin.nomiss.i){
          mu.ij <- muX.now[phi.i,j] - 1/QX[phi.i,j,j]*sum(QX[phi.i,j,-j]*(z.now.i[-j]-muX.now[phi.i,-j]))
          s2.ij <- 1/QX[phi.i,j,j]
          z.now.i[j] <- myrtruncnorm(1, mu.ij, sqrt(s2.ij), T=0, greater=as.logical(xc.now.i[j+1]))
        }
      }
      Xc.now.b[k,] <- xc.now.i
      Z.now.b[k,] <- z.now.i
      accept.b[k,] <- accept.i
      phi.now.b[k] <- phi.i
    }
    cbind(Xc.now.b, Z.now.b, accept.b, phi.now.b)
  }
  Xc.now <- XZ[,1:(p+1)]
  Z.now <- XZ[,(p+2):(2*p+1)]
  accept <- XZ[,(2*p+2):(3*p+1)]
  phi.now <- XZ[,3*p+2]
  accept <- accept[is.na(X[,-1])]

  return(list(Xc.now=Xc.now, Z.now=Z.now, phi.now=phi.now, accept=accept))
}





update.XZ.probit <- function(X, Xc.now, Z.now, bin.X, muX.now, SigmaX.now, yc.now, beta.now, mult2, omega.now, groups, inv.groups, nb=getDoParWorkers()){

  n <- nrow(X)
  p <- ncol(X)-1
  M <- length(omega.now)
  ans.QX <- foreach(m=1:M)%dopar%{
    ans.inv <- ginv.gp(SigmaX.now[m,,])
    list(ans.inv$inv, ans.inv$log.det)
  }
  QX <- SigmaX.now
  ldet.vec <- rep(0,M)
  for(m in 1:M){
    QX[m,,] <- ans.QX[[m]][[1]]
    ldet.vec[m] <- ans.QX[[m]][[2]]
  }
  
## for each missing j, draw Zij (and effectively Xij) | rest, and accept/reject
  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/nb)
  for(b in 1:nb){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
  }

  XZ <- foreach(b=1:nb, .combine=rbind)%dopar%{
    Xc.now.b <- matrix(0, length(blocks[[b]]), p+1)
    Z.now.b <- matrix(0, length(blocks[[b]]), p)
    accept.b <- matrix(NA, length(blocks[[b]]), p)
    phi.now.b <- rep(0, length(blocks[[b]]))
    for(k in 1:length(blocks[[b]])){
      i <- blocks[[b]][k]
      
     ## Gibbs sample a new phi.now[i]
      xc.now.i <- Xc.now[i,]
      z.now.i <- Z.now[i,]
      probs.i <- dmvmixnorm.probit(z.now.i, omega.now, muX.now, QX, ldet.vec, yc.now[i], xc.now.i, beta.now)
      phi.i <- sample(1:M, 1, prob=probs.i)

     ## for each missing j, draw Xij, Zij | rest 
      ind.miss.i <- which(is.na(X[i,-1]))
      accept.i <- rep(0,p)
      if(length(ind.miss.i)>0){
        for(j in ind.miss.i){
          mu.miss.ij <- muX.now[phi.i,j] - 1/QX[phi.i,j,j]*sum(QX[phi.i,j,-j]*(z.now.i[-j] - muX.now[phi.i,-j]))
          z.miss.p <- rnorm(1, mu.miss.ij, sqrt(mult2/QX[phi.i,j,j]))
          z.p.i <- z.now.i
          z.p.i[j] <- z.miss.p
          gind <- groups[[inv.groups[j+1]]]-1
	  if(bin.X[j]==1){
	    ngj <- length(gind) 
            if(ngj>1){
	      foo <- which(z.p.i[gind]==max(z.p.i[gind]))
	      xc.miss.p <- rep(0,ngj)
	      if(z.p.i[gind[foo]] > 0)
  	        xc.miss.p[foo] <- 1
	    }
	    else
              xc.miss.p <- as.numeric(z.miss.p > 0)
          }
	  else
	    xc.miss.p <- z.miss.p
          xc.p.i <- xc.now.i
          xc.p.i[gind+1] <- xc.miss.p
	  if(mult2==1){
	    dZ.p <- dZ.now <- pi.p <- pi.now <- 0
	  }
	  else{
            dZ.p <- dnorm(z.miss.p, mu.miss.ij, sqrt(mult2/QX[phi.i,j,j]), log=TRUE)
            dZ.now <- dnorm(z.now.i[j], mu.miss.ij, sqrt(mult2/QX[phi.i,j,j]), log=TRUE)
            pi.p <- dnorm(z.miss.p, mu.miss.ij, sqrt(1/QX[phi.i,j,j]), log=TRUE)
            pi.now <- dnorm(z.now.i[j], mu.miss.ij, sqrt(1/QX[phi.i,j,j]), log=TRUE)
	  }
          likey.p <- dnorm(yc.now[i], sum(xc.p.i*beta.now[phi.i,]), 1, log=TRUE)
          likey.now <- dnorm(yc.now[i], sum(xc.now.i*beta.now[phi.i,]), 1, log=TRUE)
          MH.ratio <- exp(likey.p + pi.p - dZ.p - likey.now - pi.now + dZ.now)
          if(runif(1) < MH.ratio){
            xc.now.i <- xc.p.i
	    z.now.i <- z.p.i
            accept.i[j] <- 1
          }
	  else
	    accept.i[j] <- 0
        }
      }
    ## for each binary, non-missing j, draw Zij | rest 
      ind.bin.nomiss.i <- which(!is.na(X[i,-1]) & bin.X==1)
      if(length(ind.bin.nomiss.i)>0){
        for(j in ind.bin.nomiss.i){
          gind <- groups[[inv.groups[j+1]]]-1
          mu.ij <- muX.now[phi.i,j] - 1/QX[phi.i,j,j]*sum(QX[phi.i,j,-j]*(z.now.i[-j]-muX.now[phi.i,-j]))
          s2.ij <- 1/QX[phi.i,j,j]
	  foo <- which(xc.now.i[gind+1]==1)
	  if(length(foo)==0){
            greater.j <- FALSE
	    T.j <- 0
	  }
	  else if(gind[foo]==j){
	    greater.j <- TRUE
	    T.j <- max(0,z.now.i[gind[-foo]])
	  }
	  else{
	    greater.j <- FALSE
	    T.j <- max(z.now.i[gind[foo]])
          }	    
          z.now.i[j] <- myrtruncnorm(1, mu.ij, sqrt(s2.ij), T=T.j, greater=greater.j)
        }
      }
      Xc.now.b[k,] <- xc.now.i
      Z.now.b[k,] <- z.now.i
      accept.b[k,] <- accept.i
      phi.now.b[k] <- phi.i
    }
    cbind(Xc.now.b, Z.now.b, accept.b, phi.now.b)
  }
  Xc.now <- XZ[,1:(p+1)]
  Z.now <- XZ[,(p+2):(2*p+1)]
  accept <- XZ[,(2*p+2):(3*p+1)]
  phi.now <- XZ[,3*p+2]
  accept <- accept[is.na(X[,-1])]

  return(list(Xc.now=Xc.now, Z.now=Z.now, phi.now=phi.now, accept=accept))
}







update.muX <- function(Z.now, SigmaX.now, phi.now, M.mu, S.mu, ind.scale){

  M <- dim(SigmaX.now)[1]
  p <- ncol(Z.now)
  muX.now <- foreach(m=1:M, .combine=rbind)%dopar%{
    scale.m <- diag(sqrt(diag(SigmaX.now[m,,])))
    diag(scale.m)[!ind.scale] <- 1
    M.mu.m <- scale.m%*%M.mu
    S.mu.m <- scale.m%*%S.mu%*%scale.m
    ind.m <- which(phi.now==m)
    n.m <- length(ind.m)
    if(n.m > 0){
      zbar.m <- colMeans(Z.now[ind.m,,drop=FALSE])
      Q0.m <- ginv.gp(S.mu.m)$inv
      QX.m <- ginv.gp(SigmaX.now[m,,])$inv
      Qs.m <- Q0.m + n.m*QX.m
      ans.inv <- ginv.gp(Qs.m)
      Ss.m <- ans.inv$inv
      Ss.sqrt.m <- ans.inv$sqrt.inv
      mus.m <- Ss.m%*%(Q0.m%*%M.mu.m + n.m*QX.m%*%zbar.m)
    }
    else{
      mus.m <- M.mu.m
      Ss.sqrt.m <- ginv.gp(S.mu.m)$sqrt
    }
    rmvnorm(1,mu=mus.m,S.sqrt=Ss.sqrt.m)
  }
  muX.now <- matrix(muX.now,M,p)
  return(muX.now)
}



update.SigmaX <- function(SigmaX.now, Z.now, muX.now, phi.now, P.Sigma, N.Sigma, bin.X){

  M <- dim(muX.now)[1]
  p <- ncol(Z.now)
  SigmaX.list <- foreach(m=1:M)%dopar%{
    ind.m <- which(phi.now==m)
    n.m <- length(ind.m)
    Ns.m <- N.Sigma + n.m
    if(n.m > 0){
      mu.mat.m <- matrix(muX.now[m,], n.m, p, byrow=TRUE)
      CP.nnd <- crossprod(Z.now[ind.m,,drop=FALSE] - mu.mat.m)
#      CP <- crossprod(Z.now[ind.m,,drop=FALSE] - mu.mat.m)
#      blah <- eigen(CP)
#      lam <- blah$val
#      lam[lam<0] <- 1E-6
#      CP.nnd <- blah$vec%*%diag(lam)%*%t(blah$vec)
      Ps.m <- P.Sigma + CP.nnd
    }
    else{
      Ps.m <- P.Sigma      
    }
    ans <- try(riwish(Ns.m, Ps.m), silent=TRUE)
    if(is.character(ans[1])){
      ans <- SigmaX.now[m,,]
      warning("non PD update matrix in update.SigmaX")
    }
    ans
  }  
  SigmaX.now <- array(0,c(M,p,p))
  for(m in 1:M)
    SigmaX.now[m,,] <- SigmaX.list[[m]]
  return(SigmaX.now)
}
    


update.vee <- function(phi.now, delta.now, M){

  if(M==1){
    vee <- 0
    omega <- 1
  }
  else{
    vee <- rep(0,M)
    for(m in 1:M){
      a.star <- sum(phi.now==m)+1
      b.star <- sum(phi.now>m)+delta.now
      vee[m] <- rbeta(1,b.star, a.star)
    }
    vee[M] <- 0
    omega <- makeprobs(vee)
  }
  return(list(vee=vee,omega=omega))
}



update.delta <- function(vee.now, A.delta, B.delta){

  M <- length(vee.now)
  a.star <- M + A.delta
  b.star <- -sum(log(vee.now[-M])) + B.delta
  delta <- rgamma(1, a.star, b.star)
  return(delta)
}






update.Omega <- function(SigmaX.now, eta.now, P.Omega, N.Omega){

#SigmaX.now<<-SigmaX.now
#eta.now<<-eta.now
#P.Omega<<-P.Omega
#N.Omega<<-N.Omega

  M <- dim(SigmaX.now)[1]
  p <- dim(SigmaX.now)[2]
  QX.list <- foreach(m=1:M)%dopar%{
    ans.inv <- ginv.gp(SigmaX.now[m,,])
    ans.inv$inv
  }
  sumQX <- matrix(0,p,p)
  for(m in 1:M)
    sumQX <- sumQX + QX.list[[m]]
  Ps.inv <- sumQX + ginv.gp(P.Omega)$inv
  Ps <- ginv.gp(Ps.inv)$inv
  Ns <- M*eta.now + N.Omega
  Omega.now <- rwish(Ns, Ps)
}




update.eta <- function(SigmaX.now, eta.now, Omega.now, A.eta, B.eta, prop.sd.eta, dfp, etagp){

  accept <- 0
  M <- dim(SigmaX.now)[1]
  p <- dim(SigmaX.now)[2]
  eta.p <- r.pos.proposal(eta.now-p-etagp, df=dfp, sigma=prop.sd.eta) + p + etagp
  like.p <- foreach(m=1:M, .combine=sum)%dopar%{
    diwish(SigmaX.now[m,,], eta.p, Omega.now)
  }
  like.now <- foreach(m=1:M, .combine=sum)%dopar%{
    diwish(SigmaX.now[m,,], eta.now, Omega.now)
  }
  pi.p <- dgamma(eta.p-p-etagp, A.eta, B.eta, log=TRUE)
  pi.now <- dgamma(eta.now-p-etagp, A.eta, B.eta, log=TRUE)
  dprop.p <- d.pos.proposal(eta.p-p-etagp, eta.now-p-etagp, df=dfp, sigma=prop.sd.eta)
  dprop.now <- d.pos.proposal(eta.now-p-etagp, eta.p-p-etagp, df=dfp, sigma=prop.sd.eta)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    eta.now <- eta.p
    accept <- 1
  }
  return(list(eta=eta.now, accept=accept))
}





update.nu <- function(muX.now, Psi.now, M.nu, S.nu){
  M <- nrow(muX.now)
  p <- ncol(muX.now)
  mubar <- colMeans(muX.now)
  Q0 <- ginv.gp(S.nu)$inv
  QX <- ginv.gp(Psi.now)$inv
  Qs <- Q0 + M*QX
  ans.inv <- ginv.gp(Qs)
  Ss <- ans.inv$inv
  Ss.sqrt <- ans.inv$sqrt.inv
  mus <- Ss%*%(Q0%*%M.nu + M*QX%*%mubar)
  nu.now <- rmvnorm(1,mu=mus,S.sqrt=Ss.sqrt)
  return(nu.now)
}



update.Psi <- function(muX.now, nu.now, P.Psi, N.Psi){

#muX.now<<-muX.now
#nu.now<<-nu.now
#P.Psi<<-P.Psi
#N.Psi<<-N.Psi

  M <- nrow(muX.now)
  p <- ncol(muX.now)
  Ns <- N.Psi + M
  nu.mat <- matrix(nu.now, M, p, byrow=TRUE)
  Ps <- P.Psi + crossprod(muX.now - nu.mat)
  Psi.now <- riwish(Ns, Ps)
  return(Psi.now)
}


##################################################################
################# Weibull Regression MCMC  #######################
##################################################################


weibull.MCMC <- function(X, y, outcome=NULL, rho=.5, rho.prop.01=.3, rho.prop.10=.3, A.tau2=1, B.tau2=1, rho.kappa=.5, A.kappa=1, B.kappa=1, prop.sd.kappa=.05, kappa.prop.01=.3, kappa.prop.10=.3, M.nu=0, S.nu=100, P.Psi=NULL, N.Psi=NULL, M.Beta=0, S2.Beta=.001, P.Omega=NULL, N.Omega=NULL, A.eta=1, B.eta=.1, prop.sd.eta=.5, A.delta=1, B.delta=.1, N.mcmc=5000, every=1, nplot=10, nback=1000, dfp=20, groups=NULL, mult1=2, mult2=2, bin.X=NULL, everyXZ=1, maxplot=82, one.tau=TRUE, one.beta=TRUE, beta.init=NULL, phi.init=NULL, omega.init=NULL, M=20, G=NULL, begin=0, nb=getDoParWorkers(), etagp=2){

 ## Create needed variables ##
  n <- length(y)
  p <- ncol(X)
  if(is.null(groups))
    groups <- as.list(1:(p+1))
  else
    groups <- c(list(1),sapply(groups,"+",1))
  nG <- length(groups)
  if(is.null(outcome))
    outcome <- rep(1,n)
  y.Surv <- Surv(y, event=outcome)

  if(length(rho.prop.01)==1)
    rho.prop.01 <- rep(rho.prop.01, nG)
  if(length(rho.prop.10)==1)
    rho.prop.10 <- rep(rho.prop.10, nG)
  if(length(rho)==1)
    rho <- rep(rho, nG)
  rho.prop.01[1] <- .9
  rho.prop.10[1] <- .1
  rho[1] <- .999999

  if(is.null(G))
    G <- 1:M
  if(length(mult1)==1)
    mult1 <- rep(mult1, nG)

  if(length(M.nu)==1)
    M.nu <- rep(M.nu,p)
  if(length(S.nu)==1)
    S.nu <- diag(S.nu,p)
  if(is.null(N.Omega))
    N.Omega <- (p+1)*1.1
  if(is.null(P.Omega))
    P.Omega <- diag(1,p)*(etagp+A.eta/B.eta-1)/N.Omega
  if(is.null(N.Psi))
    N.Psi <- (p+1)*1.1
  if(is.null(P.Psi))
    P.Psi <- diag(1,p)*(N.Psi-p-1)

  if(is.null(bin.X)){
    bin.X <- rep(0,p)
    for(j in 1:p){
      bin.X[j] <- as.numeric(length(unique(X[!is.na(X[,j]),j]))<3)
    }
  }

  ind.scale <- rep(FALSE,p)
  for(j in 2:nG){
    bin.j <- bin.X[groups[[j]][1]-1]
    if(bin.j)
      ind.scale[groups[[j]][1]-1] <- TRUE
  }

  gammap.1 <<- gamma2(1)
  gammapp.1 <<- gamma3(1)

 ## From here on out X is a p+1 column matrix (1st col is intercept) ##
  X <- cbind(1,X)

## initialize kappa
  kappa.now <- 1

#X<<-X
#y<<-y
#outcome<<-outcome
#bin.X<<-bin.X
#M<<-M
#G<<-G
#phi.init<<-phi.init
#beta.init<<-beta.init
#kappa.now<<-kappa.now
#groups<<-groups
#rho<<-rho
#take.dump

 ## initialize Xc, Z, cluster probs, means, Sigmas, and beta
cat("\nInitializing Clusters \n")
  ans.XZ <- initialize.XZ(X, y, outcome, bin.X, M, G, phi.init, beta.init, kappa.now, groups, rho)
  Xc.now <- ans.XZ$Xc
  Z.now <- ans.XZ$Z
  muX.now <- ans.XZ$muX
  SigmaX.now <- ans.XZ$SigmaX
  phi.now <- ans.XZ$phi.now
  beta.now <- ans.XZ$beta.now
  Beta.now <- colMeans(beta.now)
  nu.now <- apply(muX.now,2,mean)
  if(M==1){
    Psi.now <- N.Psi*P.Psi
  }
  else{
    Psi.now <- cov(muX.now)
    while(min(eigen(Psi.now)$val) < 1E-8)
      diag(Psi.now) <- 1.5*diag(Psi.now)+.01
  }
  Omega.now <- apply(SigmaX.now,c(2,3),mean)*(N.Omega-p-1)
  eta.now <- A.eta/B.eta + p + etagp

  delta.now <- A.delta/B.delta
  ans.vee <- update.vee(phi.now, delta.now, M)
  vee.now <- ans.vee$vee
  omega.now <- ans.vee$omega
  
  tau2.now <- rep(B.tau2/A.tau2, nG)
  risk.now <- rep(0,n)

cat("\nAllocating Memory for Posterior Objects \n")

 ## Allocate posterior objects for which to store MCMC samples
  Beta <- matrix(0, (N.mcmc-begin)%/%every, p+1)
  beta <- array(0, c((N.mcmc-begin)%/%every, M, p+1))
  kappa <- rep(0, (N.mcmc-begin)%/%every)
  tau2 <- matrix(0, (N.mcmc-begin)%/%every, nG)
  muX <- array(0, c((N.mcmc-begin)%/%every, M, p))
  SigmaX <- array(0, c((N.mcmc-begin)%/%every,M,p,p))
  vee <- matrix(0, (N.mcmc-begin)%/%every, M)
  omega <- matrix(0, (N.mcmc-begin)%/%every, M)
  lambda <- matrix(0, (N.mcmc-begin)%/%every, n)
  
  accept.beta <- matrix(0, N.mcmc, nG)
  accept.kappa <- rep(0, N.mcmc)
  accept.eta <- rep(0, N.mcmc)
  accept.XZ <- rep(0, sum(is.na(X[,-1])))

  concord1.avg <- .8
  concord2.avg <- .8
  h.avg <- .05

 
 ################
 ## Begin MCMC ##
 ################
 cat("\n")
  for(it in 1:N.mcmc){

#if(it%%nplot==0){
#
#if(it%%100==0)
  cat("\nIteration", it, "out of", N.mcmc)
#}

#print("y update")

#y<<-y
#outcome<<-outcome
#Xc.now<<-Xc.now
#beta.now<<-beta.now
#kappa.now<<-kappa.now
#phi.now<<-phi.now

    yc.now <- update.y(y, outcome, Xc.now, beta.now, phi.now, kappa.now)


#print("beta update")

#y<<-y
#yc.now<<-yc.now
#outcome<<-outcome
#Xc.now<<-Xc.now
#beta.now<<-beta.now
#kappa.now<<-kappa.now
#phi.now<<-phi.now
#Beta.now<<-Beta.now
#tau2.now<<-tau2.now
#groups<<-groups
#mult1<<-mult1
#rho<<-rho
#rho.prop.01<<-rho.prop.01
#rho.prop.10<<-rho.prop.10


    ans.beta <- update.beta(y, yc.now, outcome, Xc.now, beta.now, kappa.now, phi.now, Beta.now, tau2.now, groups, mult1, rho, rho.prop.01, rho.prop.10, one.beta)
    beta.now <- ans.beta$beta
    accept.beta[it,] <- ans.beta$accept

#print("B update")

    Beta.now <- update.Beta(beta.now, groups, tau2.now, one.tau, M.Beta, S2.Beta, one.beta)

#print("tau^2 update")

#beta.now<<-beta.now
#groups<<-groups
#Beta.now<<-Beta.now
#A.tau2<<-A.tau2
#B.tau2<<-B.tau2
#one.tau<<-one.tau


    tau2.now <- update.tau2(beta.now, groups, Beta.now, A.tau2, B.tau2, one.tau, one.beta)


#print("kappa update")

#y<<-y
#outcome<<-outcome
#Xc.now<<-Xc.now
#beta.now<<-beta.now
#kappa.now<<-kappa.now
#phi.now<<-phi.now
#A.kappa<<-A.kappa
#B.kappa<<-B.kappa
#rho.kappa<<-rho.kappa
#prop.sd.kappa<<-prop.sd.kappa
#kappa.prop.01<<-kappa.prop.01
#kappa.prop.10<<-kappa.prop.10
#dfp<<-dfp

    ans.kappa <- update.kappa(y, outcome, Xc.now, beta.now, kappa.now, phi.now, A.kappa, B.kappa, rho.kappa, prop.sd.kappa, kappa.prop.01, kappa.prop.10, dfp)
    kappa.now <- ans.kappa$kappa
    accept.kappa[it] <- ans.kappa$accept

#print("X and Z updates")

#X<<-X
#Xc.now<<-Xc.now
#Z.now<<-Z.now
#bin.X<<-bin.X
#muX.now<<-muX.now
#SigmaX.new<<-SigmaX.now
#y<<-y
#beta.now<<-beta.now
#kappa.now<<-kappa.now
#mult2<<-mult2
#omega.now<<-omega.now
#nb<<-nb

  ## Update X and Z only every everyXZ iterations
    if(it%%everyXZ==0){
      ans.XZ <- update.XZ(X, Xc.now, Z.now, bin.X, muX.now, SigmaX.now, y, outcome, beta.now, kappa.now, mult2, omega.now, nb)
      Xc.now <- ans.XZ$Xc
      Z.now <- ans.XZ$Z
      phi.now <- ans.XZ$phi.now
      accept.XZ <- accept.XZ + ans.XZ$accept
    }

#tab.phi <- table(phi.now)
#foo <- as.numeric(names(tab.phi[tab.phi==max(tab.phi)]))
#DS <- diag(SigmaX.now[foo,,])
#indp <- which(DS==max(DS[bin.X==0]))
#Z.foo <- Z.now[phi.now==foo,indp]
#cat("\nm =",foo)
#cat("\nj =",indp)
#cat("\n summary(Z.foo) =\n")
#print(summary(Z.foo))
#print(summary(Z.foo[!is.na(X[phi.now==foo,indp+1])]))
#print(summary(Z.foo[is.na(X[phi.now==foo,indp+1])]))
#cat("\n var(Z.foo) =",var(Z.foo), "\n")
#cat("\n SigmaX.now[foo,indp,indp] =",SigmaX.now[foo,indp,indp], "\n")


#phi.now<<-phi.now
#delta.now<<-delta.now
#M<<-M

## Update vee/omega
    ans.vee <- update.vee(phi.now, delta.now, M)
    vee.now <- ans.vee$vee
    omega.now <- ans.vee$omega


  ## Update delta
    delta.now <- update.delta(vee.now, A.delta, B.delta)



#print("Omega and eta update")


#SigmaX.now<<-SigmaX.now
#eta.now<<-eta.now
#P.Omega<<-P.Omega
#N.Omega<<-N.Omega

    Omega.now <- update.Omega(SigmaX.now, eta.now, P.Omega, N.Omega)
    ans.eta <- update.eta(SigmaX.now, eta.now, Omega.now, A.eta, B.eta, prop.sd.eta, dfp, etagp)
    eta.now <- ans.eta$eta
    accept.eta[it] <- ans.eta$accept


#print("SigmaX update")

#SigmaX.now<<-SigmaX.now
#Z.now<<-Z.now
#muX.now<<-muX.now
#phi.now<<-phi.now
#Omega.now<<-Omega.now
#eta.now<<-eta.now
#bin.X<<-bin.X

    SigmaX.now <- update.SigmaX(SigmaX.now, Z.now, muX.now, phi.now, Omega.now, eta.now, bin.X)


#print("muX update")

#Z.now<<-Z.now
#SigmaX.new<<-SigmaX.now
#phi.now<<-phi.now
#nu.now<<-nu.now
#Psi.now<<-Psi.now
#bin.X<<-bin.X

    muX.now <- update.muX(Z.now, SigmaX.now, phi.now, nu.now, Psi.now, ind.scale)


   ## scale muX and SigmaX
    OmegaX.now <- SigmaX.now
    nuX.now <- muX.now
    for(m in 1:M){
      scale.m <- diag(1/sqrt(diag(SigmaX.now[m,,])))
      diag(scale.m)[!ind.scale] <- 1
      OmegaX.now[m,,] <- scale.m%*%SigmaX.now[m,,]%*%t(scale.m)
      nuX.now[m,] <- muX.now[m,]%*%scale.m
    }

#print("nu and Psi update")

    nu.now <- update.nu(nuX.now, Psi.now, M.nu, S.nu)
    Psi.now <- update.Psi(nuX.now, nu.now, P.Psi, N.Psi)



#print("End of Updates")

    for(m in 1:M){
      ind.m <- which(phi.now==m)
      risk.now[ind.m] <- as.numeric(Xc.now[ind.m,]%*%beta.now[m,])
    }

   ## record params.now in posterior sample
    if(it>begin && (it-begin)%%every==0){

      Beta[(it-begin)/every,] <- Beta.now*as.numeric(beta.now[1,]!=0)
      beta[(it-begin)/every,,] <- beta.now
      tau2[(it-begin)/every,] <- tau2.now
      kappa[(it-begin)/every] <- kappa.now
      omega[(it-begin)/every,] <- omega.now
      vee[(it-begin)/every,] <- vee.now
      muX[(it-begin)/every,,] <- nuX.now
      SigmaX[(it-begin)/every,,,] <- OmegaX.now
      lambda[(it-begin)/every,] <- risk.now
    }
    

   ## Summarize and Plot posterior
    if((it-begin)%%nplot==0){
    
      concord1 <- survConcordance(y.Surv~risk.now)$concord
      concord2 <- survConcordance(Surv(yc.now, event=rep(1,n))~risk.now)$concord
      concord1.avg <- (1-h.avg)*concord1.avg + h.avg*concord1
      concord2.avg <- (1-h.avg)*concord2.avg + h.avg*concord2

      pp <- min(maxplot,p)
      if(one.beta)
        N.plots <- (pp+1)+nG*(1-one.tau)+one.tau+1+2
      else
        N.plots <- 2*(pp+1)+nG*(1-one.tau)+one.tau+1+2
      cols <- min(13, ceiling(sqrt(N.plots)))
      rows <- min(13, ceiling(N.plots/cols))
      it.e <- floor((it-begin)/every)
      ind.now <- max(1,floor(it.e/2),it.e-nback+1):max(1,it.e)
      par(mfrow=c(rows,cols), mar=c(2,2,2,1))

     ## Print Rsq and recent average
      cat("\n\nConcordance  =",concord1)
      cat("\nConcordance (recent avg) =",concord1.avg)
      cat("\nConcordance2  =",concord2)
      cat("\nConcordance2 (recent avg) =",concord2.avg,"\n")
      like.y <- get.like(y, outcome, Xc.now, beta.now, kappa.now, phi.now)
      cat("\nlog(like) =",like.y,"\n")

     ## Print acceptance %
      cat("\nkappa acceptance = ", mean(accept.kappa[1:it]))
      cat("\neta acceptance = ", mean(accept.eta[1:it]))
      cat("\nbeta acceptance = \n")
      print(summary(colMeans(accept.beta[1:it,,drop=FALSE])))
      cat("\nnX/Z acceptance = \n")
      print(summary(accept.XZ/floor(it/everyXZ)))
      cat("\n\n\n")
    }
    if((it-begin)%%nplot==0 && it>begin){
     ## Plot Betas
      if(one.beta){
        Beta.ord <- c(1,order(-abs(colMeans(beta[ind.now,1,-1])))+1)
        for(k in 1:(pp+1))
          plot(beta[ind.now,1,Beta.ord[k]],ylab="",main=paste("Beta_",Beta.ord[k]-1,sep=""),cex=.5)
      }
      else{
        Beta.ord <- c(1,order(-colMeans(abs(beta.now[,-1])*omega.now))+1)
        for(k in 1:(pp+1)){
          plot(Beta[ind.now,Beta.ord[k]],ylab="",main=paste("Beta_",Beta.ord[k]-1,sep=""),cex=.5)
          beta.k.hist <- rep(beta.now[,Beta.ord[k]], round(omega.now*100))
	  hist(beta.k.hist, main=paste("beta_",Beta.ord[k]-1,sep=""),cex=.5)
	}
      }

#     ## Plot tau2's
#      if(one.tau){
#        plot(tau2[ind.now,1],ylab="",main="tau2",cex=.5)
#      }
#      else{
#        for(k in 1:nG)
#          plot(tau2[ind.now,k],ylab="",main=paste("tau2_",k-1,sep=""),cex=.5)
#      }
     
     ## Plot kappa
      plot(kappa[ind.now],ylab="",main="kappa",cex=.5)

    ## Bar plot of omega.now
      barplot(omega.now, main="omega")
      abline(h=exp(-6), col=2)
      abline(h=exp(-4), col=4)

      lomega <- -sort(-log(omega.now+1E-300))
      barplot(lomega, main="omega", yaxt='n')
      at <- ceiling(1/10*round(seq(min(lomega),0,length=4),0))*10
      labels <- format(exp(at),nsmall=1,digits=1)
      axis(2,at=at, labels=labels)
      abline(h=-6, col=2)
      abline(h=-4, col=4)
           
    }
  }
  return(list(beta=beta, Beta=Beta, tau2=tau2, kappa=kappa, muX=muX, SigmaX=SigmaX, omega=omega, vee=vee, X=X, y=y, outcome=outcome, bin.X=bin.X, lambda=lambda))
}












############################################
####### Weibull Prediction Functions #######
############################################




sample.xz.predict <- function(x, bin.X, muX, SigmaX, QX, omega, ldet, N.mcmc, post.ind){

  p <- length(x)-1
  M <- length(omega)
  nreal <- length(post.ind)
  xc <- matrix(0,N.mcmc,p+1)
  phi <- rep(0,N.mcmc)
  ind.miss <- which(is.na(x[-1]))
#  if(length(ind.miss)==0)
#    return(list(xc=matrix(x,nrow=nreal,p+1,byrow=TRUE), phi=rep(1,nreal)))
  nm <- length(ind.miss)
  nb <- sum(bin.X==1)
  S.miss <- S.miss.sqrt <- array(0, c(M, nm, nm))
  if(nm>0){
    for(m in 1:M){
      Q.miss.m <- QX[m,ind.miss,ind.miss]
      ans.inv <- ginv.gp(Q.miss.m)
      S.miss[m,,] <- ans.inv$inv
      S.miss.sqrt[m,,] <- ans.inv$sqrt.inv
    }
  }
  
 # initialize xc.now, z.now, phi
    phi[1] <- get.best.phi(x, bin.X, omega, muX, SigmaX)
    xc.now <- x[-1]
    z.now <- xc.now
    z.now[bin.X==1] <- 2*xc.now[bin.X==1]-1
    z.now[is.na(x[-1])] <- muX[phi[1],is.na(x[-1])]
    xc.now[bin.X==0] <- z.now[bin.X==0]
    xc.now[bin.X==1] <- as.numeric(z.now[bin.X==1]>0)
    xc.now <- c(1,xc.now)

  for(it in 1:N.mcmc){

   ## Gibbs sample a new phi.now[i]
    if(it>1){
      probs <- dmvmixnorm_noy(z.now, omega, muX, QX, ldet)
      phi[it] <- sample(1:M, 1, prob=probs)
    }
    QX.it <- QX[phi[it],,]
    muX.it <- muX[phi[it],]

   ## for each missing j, draw Xij, Zij | rest
    if(nm>0){
      S.miss.it <- matrix(S.miss[phi[it],,],nm,nm)
      S.miss.sqrt.it <- matrix(S.miss.sqrt[phi[it],,],nm,nm)
      mu.miss <- muX.it[ind.miss] - S.miss.it%*%QX.it[ind.miss,-ind.miss,drop=FALSE]%*%(z.now[-ind.miss] - muX.it[-ind.miss])
      z.miss.now <- rmvnorm(1,mu=mu.miss, S.sqrt=S.miss.sqrt.it)
      xc.miss.now <- z.miss.now
      xc.miss.now[bin.X[ind.miss]==1] <- as.numeric(z.miss.now[bin.X[ind.miss]==1] > 0)
      xc.now[ind.miss+1] <- xc.miss.now
      z.now[ind.miss] <- z.miss.now
    }
   ## for each binary, non-missing j, draw Zij | rest
    if(nb>0){
      ind.bin.nomiss <- which(!is.na(x[-1]) & bin.X==1)
      if(length(ind.bin.nomiss)>0){
        for(j in ind.bin.nomiss){
          mu.j <- muX.it[j] - 1/QX.it[j,j]*sum(QX.it[j,-j]*(z.now[-j]-muX.it[-j]))
          s2.j <- 1/QX.it[j,j]
          z.now[j] <- myrtruncnorm(1, mu.j, sqrt(s2.j), T=0, greater=as.logical(xc.now[j+1]))
        }
      }
    }
    xc[it,] <- xc.now
  }
  xc <- xc[post.ind,]
  phi <- phi[post.ind]

  return(list(xc=xc, phi=phi))
}







### Predict at theta hat function ###



predict.at.theta.hat <- function(obj, Xp, alpha=.05, nrealX=1000, N.mcmcX=NULL, post.indX=NULL, post.ind=NULL, B=getDoParWorkers()){

  Xp <- cbind(1,Xp)
  beta <- obj$beta
  muX <- obj$muX
  SigmaX <- obj$SigmaX
  omega <- obj$omega
  bin.X <- obj$bin.X
  N.mcmc <- nrow(beta)
  if(is.null(N.mcmcX))
    N.mcmcX <- max(2*nrealX,200)
  if(is.null(post.indX))
    post.indX <- sort(sample(floor(N.mcmcX/2+1):N.mcmcX, nrealX, replace=nrealX>N.mcmcX/2))
  if(is.null(post.ind))
    post.ind <- floor(N.mcmc/2):N.mcmc

#  beta.hat <- apply(beta[post.ind,,,drop=FALSE],c(2,3),median)
  beta.hat <- apply(beta[post.ind,,,drop=FALSE],c(2,3),mean)
  muX.hat <- apply(muX[post.ind,,,drop=FALSE], c(2,3), mean)
  SigmaX.hat <- apply(SigmaX[post.ind,,,,drop=FALSE], c(2,3,4), mean)
  omega.hat <- apply(omega[post.ind,,drop=FALSE],2,mean)

  n <- nrow(Xp)
  p <- ncol(Xp)-1
  M <- length(omega.hat)

  QX.hat <- SigmaX.hat
  ldet.hat <- rep(0,M)
  for(m in 1:M){
    ans.inv <- ginv.gp(SigmaX.hat[m,,])
    QX.hat[m,,] <- ans.inv$inv
    ldet.hat[m] <- ans.inv$log.det
  }

  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/B)
  for(b in 1:B){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
    if(ind.now>n){
      B <- b
      break
    }
  }

  lambda <- foreach(b=1:B, .combine=rbind)%dopar%{
    n.b <- length(blocks[[b]])
    Xc.b <- matrix(0, nrealX*n.b, p+1)
    phi.b <- matrix(0,nrealX, n.b)
    lambda.hat.b <- rep(0, n.b)
    lambda.b <- matrix(0, n.b, nrealX)
    for(k in 1:length(blocks[[b]])){

if(b==1 && k%%10==0)
cat("\nk = ",k," out of ",length(blocks[[b]]), sep="")

      i <- blocks[[b]][k]
      ans.i <- sample.xz.predict(Xp[i,], bin.X, muX.hat, SigmaX.hat, QX.hat, omega.hat, ldet.hat, N.mcmcX, post.indX)
      Xc.i <- ans.i$xc
      phi.b[,k] <- ans.i$phi
      lambda.b[k,] <- rowSums(Xc.i*beta.hat[ans.i$phi,])
      tab.k <- table(phi.b[,k])
      phi.hat.k <- as.numeric(names(sort(-tab.k)))[1]
      X.hat.k <- colMeans(Xc.i[phi.b[,k]==phi.hat.k,])
      lambda.hat.b[k] <- sum(X.hat.k*beta.hat[phi.hat.k,])
    }
    cbind(lambda.b,lambda.hat.b)
  }

  lambda.hat <- rowMeans(lambda[,-(nrealX+1)])
  lambda.hat2 <- lambda[,nrealX+1]
  lambda.CI <- t(apply(lambda[,-(nrealX+1)],1,quantile,prob=c(alpha/2,1-alpha/2)))
  return(list(pred=lambda.hat, pred2=lambda.hat2, CI=lambda.CI, lambda=lambda))
}

  




predict.pm <- function(obj, Xp, alpha=.05, nrealX=1000, N.mcmcX=NULL, post.indX=NULL, post.ind=NULL, B=getDoParWorkers()){

  Xp <- cbind(1,Xp)
  beta <- obj$beta
  muX <- obj$muX
  SigmaX <- obj$SigmaX
  omega <- obj$omega
  bin.X <- obj$bin.X
  N.mcmc <- nrow(beta)
  if(is.null(N.mcmcX))
    N.mcmcX <- max(2*nrealX,200)
  if(is.null(post.indX))
    post.indX <- sort(sample(floor(N.mcmcX/2+1):N.mcmcX, nrealX, replace=nrealX>N.mcmcX/2))
  if(is.null(post.ind))
    post.ind <- floor(N.mcmc/2):N.mcmc

  beta.hat <- apply(beta[post.ind,,,drop=FALSE],c(2,3),mean)
  ind.beta <- sample(post.ind, nrealX, replace=nrealX>length(post.ind))
  beta <- beta[ind.beta,,,drop=FALSE]
  muX.hat <- apply(muX[post.ind,,,drop=FALSE], c(2,3), mean)
  SigmaX.hat <- apply(SigmaX[post.ind,,,,drop=FALSE], c(2,3,4), mean)
  omega.hat <- apply(omega[post.ind,,drop=FALSE],2,mean)

  n <- nrow(Xp)
  p <- ncol(Xp)-1
  M <- length(omega.hat)

  QX.hat <- SigmaX.hat
  ldet.hat <- rep(0,M)
  for(m in 1:M){
    ans.inv <- ginv.gp(SigmaX.hat[m,,])
    QX.hat[m,,] <- ans.inv$inv
    ldet.hat[m] <- ans.inv$log.det
  }

  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/B)
  for(b in 1:B){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
    if(ind.now>n){
      B <- b
      break
    }
  }

  lambda <- foreach(b=1:B, .combine=rbind)%dopar%{
    n.b <- length(blocks[[b]])
    Xc.b <- matrix(0, nrealX*n.b, p+1)
    phi.b <- matrix(0,nrealX, n.b)
    lambda.hat.b <- rep(0, n.b)
    lambda.b <- matrix(0, n.b, nrealX)
    for(k in 1:length(blocks[[b]])){

if(b==1 && k%%10==0)
cat("\nk = ",k," out of ",length(blocks[[b]]), sep="")

      i <- blocks[[b]][k]
      ans.i <- sample.xz.predict(Xp[i,], bin.X, muX.hat, SigmaX.hat, QX.hat, omega.hat, ldet.hat, N.mcmcX, post.indX)
      Xc.i <- ans.i$xc
      phi.b[,k] <- ans.i$phi
      beta.phi <- 
      lambda.b[k,] <- rowSums(Xc.i*beta[,1,])
      tab.k <- table(phi.b[,k])
      phi.hat.k <- as.numeric(names(sort(-tab.k)))[1]
      X.hat.k <- colMeans(Xc.i[phi.b[,k]==phi.hat.k,])
      lambda.hat.b[k] <- sum(X.hat.k*beta.hat[phi.hat.k,])
    }
    cbind(lambda.b,lambda.hat.b)
  }

  lambda.hat <- rowMeans(lambda[,-(nrealX+1)])
  lambda.hat2 <- lambda[,nrealX+1]
  lambda.CI <- t(apply(lambda[,-(nrealX+1)],1,quantile,prob=c(alpha/2,1-alpha/2)))
  return(list(pred=lambda.hat, pred2=lambda.hat2, CI=lambda.CI, lambda=lambda))
}

  






predict.at.true <- function(Xp, omega, muX, SigmaX, beta, bin.X, alpha=.05, nrealX=1000, N.mcmcX=NULL, post.indX=NULL, B=getDoParWorkers()){

  Xp <- cbind(1,Xp)
  if(is.null(N.mcmcX))
    N.mcmcX <- max(2*nrealX,200)
  if(is.null(post.indX))
    post.indX <- sort(sample(floor(N.mcmcX/2+1):N.mcmcX, nrealX, replace=nrealX>N.mcmcX/2))

  beta.hat <- beta
  muX.hat <- muX
  SigmaX.hat <- SigmaX
  omega.hat <- omega

  n <- nrow(Xp)
  p <- ncol(Xp)-1
  M <- length(omega.hat)

  QX.hat <- SigmaX.hat
  ldet.hat <- rep(0,M)
  for(m in 1:M){
    ans.inv <- ginv.gp(SigmaX.hat[m,,])
    QX.hat[m,,] <- ans.inv$inv
    ldet.hat[m] <- ans.inv$log.det
  }

  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/B)
  for(b in 1:B){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
    if(ind.now>n){
      B <- b
      break
    }
  }

  lambda <- foreach(b=1:B, .combine=rbind)%dopar%{
    n.b <- length(blocks[[b]])
    Xc.b <- matrix(0, nrealX*n.b, p+1)
    phi.b <- matrix(0,nrealX, n.b)
    lambda.hat.b <- rep(0, n.b)
    lambda.b <- matrix(0, n.b, nrealX)
    for(k in 1:length(blocks[[b]])){
      i <- blocks[[b]][k]
      ans.i <- sample.xz.predict(Xp[i,], bin.X, muX.hat, SigmaX.hat, QX.hat, omega.hat, ldet.hat, N.mcmcX, post.indX)
      Xc.i <- ans.i$xc
      phi.b[,k] <- ans.i$phi
      lambda.b[k,] <- rowSums(Xc.i*beta.hat[ans.i$phi,])
      tab.k <- table(phi.b[,k])
#if(is.na(Xp[i,2]))
#cat("\nk =",k, "  ", tab.k)
      phi.hat.k <- as.numeric(names(sort(-tab.k)))[1]
      X.hat.k <- colMeans(Xc.i[phi.b[,k]==phi.hat.k,])
      lambda.hat.b[k] <- sum(X.hat.k*beta.hat[phi.hat.k,])
    }
    cbind(lambda.b,lambda.hat.b)
  }
  lambda.hat <- rowMeans(lambda[,-(nrealX+1)])
  lambda.hat2 <- lambda[,nrealX+1]
  lambda.CI <- t(apply(lambda[,-(nrealX+1)],1,quantile,prob=c(alpha/2,1-alpha/2)))
  return(list(pred=lambda.hat, pred2=lambda.hat2, CI=lambda.CI, lambda=lambda))
}

  










sample.xz.predict.LOS <- function(x, bin.X, LOS.ind, dayt, muX, SigmaX, QX, omega, ldet, N.mcmc, post.ind){

  p <- length(x)-1
  M <- length(omega)
  nreal <- length(post.ind)
  xc <- matrix(0,N.mcmc,p+1)
  phi <- rep(0,N.mcmc)
  ind.miss <- which(is.na(x[-1]))
  LOS <- x[LOS.ind+1]
  nm <- length(ind.miss)
  nb <- sum(bin.X==1)
  S.miss <- S.miss.sqrt <- array(0, c(M, nm, nm))
  if(nm>0){
    for(m in 1:M){
      Q.miss.m <- QX[m,ind.miss,ind.miss]
      ans.inv <- ginv.gp(Q.miss.m)
      S.miss[m,,] <- ans.inv$inv
      S.miss.sqrt[m,,] <- ans.inv$sqrt.inv
    }
  }
  
 # initialize xc.now, z.now, phi
    phi[1] <- get.best.phi(x, bin.X, omega, muX, SigmaX)
    xc.now <- x[-1]
    z.now <- xc.now
    z.now[bin.X==1] <- 2*xc.now[bin.X==1]-1
    z.now[is.na(x[-1])] <- muX[phi[1],is.na(x[-1])]
    xc.now[bin.X==0] <- z.now[bin.X==0]
    xc.now[bin.X==1] <- as.numeric(z.now[bin.X==1]>0)
    xc.now <- c(1,xc.now)

  for(it in 1:N.mcmc){

   ## Gibbs sample a new phi.now[i]
    if(it>1){
      probs <- dmvmixnorm_noy(z.now, omega, muX, QX, ldet)
      phi[it] <- sample(1:M, 1, prob=probs)
    }
    QX.it <- QX[phi[it],,]
    muX.it <- muX[phi[it],]

   ## for each missing j, draw Xij, Zij | rest
    if(nm>0){
      S.miss.it <- matrix(S.miss[phi[it],,],nm,nm)
      S.miss.sqrt.it <- matrix(S.miss.sqrt[phi[it],,],nm,nm)
      mu.miss <- muX.it[ind.miss] - S.miss.it%*%QX.it[ind.miss,-ind.miss,drop=FALSE]%*%(z.now[-ind.miss] - muX.it[-ind.miss])
      z.miss.now <- rmvnorm(1,mu=mu.miss, S.sqrt=S.miss.sqrt.it)
      xc.miss.now <- z.miss.now
      xc.miss.now[bin.X[ind.miss]==1] <- as.numeric(z.miss.now[bin.X[ind.miss]==1] > 0)
      xc.now[ind.miss+1] <- xc.miss.now
      z.now[ind.miss] <- z.miss.now
    }
   ## for each binary, non-missing j, draw Zij | rest
    if(nb>0){
      ind.bin.nomiss <- which(!is.na(x[-1]) & bin.X==1)
      if(length(ind.bin.nomiss)>0){
        for(j in ind.bin.nomiss){
          mu.j <- muX.it[j] - 1/QX.it[j,j]*sum(QX.it[j,-j]*(z.now[-j]-muX.it[-j]))
          s2.j <- 1/QX.it[j,j]
          z.now[j] <- myrtruncnorm(1, mu.j, sqrt(s2.j), T=0, greater=as.logical(xc.now[j+1]))
        }
      }
    }
   ## Sample LOS variable, truncated s.t. LOS > dayt
    if(!is.na(LOS) && LOS > dayt){
      j <- LOS.ind
      mu.j <- muX.it[j] - 1/QX.it[j,j]*sum(QX.it[j,-j]*(z.now[-j]-muX.it[-j]))
      s2.j <- 1/QX.it[j,j]
      xc.now[j+1] <- z.now[j] <- myrtruncnorm(1, mu.j, sqrt(s2.j), T=dayt, greater=TRUE)
    }
    xc[it,] <- xc.now
  }
  xc <- xc[post.ind,]
  phi <- phi[post.ind]

#print(summary(exp(xc[,96])-.1))
#print(phi)

  return(list(xc=xc, phi=phi))
}






predict.at.theta.hat.LOS <- function(obj, Xp, LOS.ind=NULL, dayt=-5, alpha=.05, nrealX=1000, N.mcmcX=NULL, post.indX=NULL, post.ind=NULL, B=getDoParWorkers(), SigmaX.hat=NULL){

  Xp <- cbind(1,Xp)
#  beta <- obj$beta
#  muX <- obj$muX
#  SigmaX <- obj$SigmaX
#  omega <- obj$omega
  bin.X <- obj$bin.X
  N.mcmc <- nrow(obj$beta)
  if(is.null(N.mcmcX))
    N.mcmcX <- max(2*nrealX,200)
  if(is.null(post.indX))
    post.indX <- sort(sample(floor(N.mcmcX/2+1):N.mcmcX, nrealX, replace=nrealX>N.mcmcX/2))
  if(is.null(post.ind))
    post.ind <- floor(N.mcmc/2):N.mcmc

#  beta.hat <- apply(obj$beta[post.ind,,,drop=FALSE],c(2,3),median)
  beta.hat <- apply(obj$beta[post.ind,,,drop=FALSE],c(2,3),mean)
  muX.hat <- apply(obj$muX[post.ind,,,drop=FALSE], c(2,3), mean)
  if(is.null(SigmaX.hat))
    SigmaX.hat <- apply(obj$SigmaX[post.ind,,,,drop=FALSE], c(2,3,4), mean)
  omega.hat <- apply(obj$omega[post.ind,,drop=FALSE],2,mean)

  n <- nrow(Xp)
  p <- ncol(Xp)-1
  M <- length(omega.hat)

  QX.hat <- SigmaX.hat
  ldet.hat <- rep(0,M)
#  for(m in 1:M){
#    ans.inv <- ginv.gp(SigmaX.hat[m,,])
#    QX.hat[m,,] <- ans.inv$inv
#    ldet.hat[m] <- ans.inv$log.det
#  }

  ans.Q <- foreach(m=1:M)%dopar%{
    ginv.gp(SigmaX.hat[m,,])
  }
  for(m in 1:M){
    QX.hat[m,,] <- ans.Q[[m]]$inv
    ldet.hat[m] <- ans.Q[[m]]$log.det
  }


  blocks <- list()
  ind.now <- 1
  inc <- ceiling(n/B)
  for(b in 1:B){
    blocks[[b]] <- ind.now:min(ind.now+inc-1, n)
    ind.now <- ind.now+inc
    if(ind.now>n){
      B <- b
      break
    }
  }

  lambda <- foreach(b=1:B, .combine=rbind)%dopar%{
    n.b <- length(blocks[[b]])
    Xc.b <- matrix(0, nrealX*n.b, p+1)
    phi.b <- matrix(0,nrealX, n.b)
    lambda.hat.b <- rep(0, n.b)
    lambda.b <- matrix(0, n.b, nrealX)
    for(k in 1:length(blocks[[b]])){

if(b==1 && k%%10==0)
cat("\nk = ",k," out of ",length(blocks[[b]]), sep="")

      i <- blocks[[b]][k]
      ans.i <- sample.xz.predict.LOS(Xp[i,], bin.X, LOS.ind, dayt, muX.hat, SigmaX.hat, QX.hat, omega.hat, ldet.hat, N.mcmcX, post.indX)
      Xc.i <- ans.i$xc
      phi.b[,k] <- ans.i$phi
      lambda.b[k,] <- rowSums(Xc.i*beta.hat[ans.i$phi,])
      tab.k <- table(phi.b[,k])
      phi.hat.k <- as.numeric(names(sort(-tab.k)))[1]
      X.hat.k <- colMeans(Xc.i[phi.b[,k]==phi.hat.k,])
      lambda.hat.b[k] <- sum(X.hat.k*beta.hat[phi.hat.k,])
    }
    cbind(lambda.b,lambda.hat.b)
  }

  lambda.hat <- rowMeans(lambda[,-(nrealX+1),drop=FALSE])
  lambda.hat2 <- lambda[,nrealX+1]
  lambda.CI <- t(apply(lambda[,-(nrealX+1),drop=FALSE],1,quantile,prob=c(alpha/2,1-alpha/2)))
  return(list(pred=lambda.hat, pred2=lambda.hat2, CI=lambda.CI, lambda=lambda[,-(nrealX+1),drop=FALSE]))
}

  







##################################################################
################# Probit Regression MCMC  #######################
##################################################################

## Needed additional functions
#
#  initialize.XZ.probit
#  update.y.probit
#  update.beta.probit
#  update.XZ.probit
#  get.like.probit
#  dmvmixnorm.probit
#


probit.MCMC <- function(X, y, rho=.5, rho.prop.01=.3, rho.prop.10=.3, A.tau2=1, B.tau2=1, M.nu=0, S.nu=100, P.Psi=NULL, N.Psi=NULL, M.Beta=0, S2.Beta=.001, P.Omega=NULL, N.Omega=NULL, A.eta=1, B.eta=.1, prop.sd.eta=.5, A.delta=1, B.delta=1, N.mcmc=5000, every=1, nplot=10, nback=1000, dfp=20, groups=NULL, mult1=2, mult2=2, bin.X=NULL, everyXZ=1, maxplot=82, one.tau=TRUE, one.beta=TRUE, beta.init=NULL, phi.init=NULL, omega.init=NULL, M=20, G=NULL, begin=0, nb=getDoParWorkers(), etagp=10, yNA=NULL, fold=NULL){

 ## Create needed variables ##
  n <- length(y)
  p <- ncol(X)
  if(is.null(groups))
    groups <- as.list(1:(p+1))
  else
    groups <- c(list(1),sapply(groups,"+",1))
  nG <- length(groups)
  inv.groups <- rep(0,p)
  for(j in 1:nG)
    inv.groups[groups[[j]]] <- j

  if(length(rho.prop.01)==1)
    rho.prop.01 <- rep(rho.prop.01, nG)
  if(length(rho.prop.10)==1)
    rho.prop.10 <- rep(rho.prop.10, nG)
  if(length(rho)==1)
    rho <- rep(rho, nG)
  rho.prop.01[1] <- .9
  rho.prop.10[1] <- .1
  rho[1] <- .999999

  if(is.null(G))
    G <- 1:M
  if(length(mult1)==1)
    mult1 <- rep(mult1, nG)

  if(length(M.nu)==1)
    M.nu <- rep(M.nu,p)
  if(length(S.nu)==1)
    S.nu <- diag(S.nu,p)
  if(is.null(P.Omega))
    P.Omega <- diag(1,p)
  if(is.null(N.Omega))
    N.Omega <- p
  if(is.null(P.Psi))
    P.Psi <- diag(1,p)
  if(is.null(N.Psi))
    N.Psi <- p

  if(is.null(bin.X)){
    bin.X <- rep(0,p)
    for(j in 1:p){
      bin.X[j] <- as.numeric(length(unique(X[!is.na(X[,j]),j]))<3)
    }
  }

  ind.scale <- rep(FALSE,p)
  for(j in 2:nG){
    bin.j <- bin.X[groups[[j]][1]-1]
    if(bin.j)
      ind.scale[groups[[j]][1]-1] <- TRUE
  }

  gammap.1 <<- gamma2(1)
  gammapp.1 <<- gamma3(1)

 ## From here on out X is a p+1 column matrix (1st col is intercept) ##
  X <- cbind(1,X)


#X<<-X
#y<<-y
#bin.X<<-bin.X
#M<<-M
#G<<-G
#phi.init<<-phi.init
#beta.init<<-beta.init
#groups<<-groups
#rho<<-rho
#take.dump

 ## initialize Xc, Z, cluster probs, means, Sigmas, and beta
cat("\nInitializing Clusters \n")
  ans.XZ <- initialize.XZ.probit(X, y, bin.X, M, G, phi.init, beta.init, groups, rho)
  Xc.now <- ans.XZ$Xc
  Z.now <- ans.XZ$Z
  muX.now <- ans.XZ$muX
  SigmaX.now <- ans.XZ$SigmaX
  phi.now <- ans.XZ$phi.now
  beta.now <- ans.XZ$beta.now
  Beta.now <- colMeans(beta.now)
  nu.now <- apply(muX.now,2,mean)
  if(M==1){
    Psi.now <- N.Psi*P.Psi
  }
  else{
    Psi.now <- cov(muX.now)
    while(min(eigen(Psi.now)$val) < 1E-8)
      diag(Psi.now) <- 1.5*diag(Psi.now)+.01
  }
  Omega.now <- apply(SigmaX.now,c(2,3),mean)*(N.Omega-p-1)
  eta.now <- A.eta/B.eta + p + etagp

  delta.now <- A.delta/B.delta
  ans.vee <- update.vee(phi.now, delta.now, M)
  vee.now <- ans.vee$vee
  omega.now <- ans.vee$omega
  
  tau2.now <- rep(B.tau2/A.tau2, nG)
  risk.now <- rep(0,n)

cat("\nAllocating Memory for Posterior Objects \n")

 ## Allocate posterior objects for which to store MCMC samples
  Beta <- matrix(0, (N.mcmc-begin)%/%every, p+1)
  beta <- array(0, c((N.mcmc-begin)%/%every, M, p+1))
  tau2 <- matrix(0, (N.mcmc-begin)%/%every, nG)
  muX <- array(0, c((N.mcmc-begin)%/%every, M, p))
  SigmaX <- array(0, c((N.mcmc-begin)%/%every,M,p,p))
  vee <- matrix(0, (N.mcmc-begin)%/%every, M)
  omega <- matrix(0, (N.mcmc-begin)%/%every, M)
  lambda <- matrix(0, (N.mcmc-begin)%/%every, n)

  accept.beta <- matrix(0, N.mcmc, nG)
  accept.eta <- rep(0, N.mcmc)
  accept.XZ <- rep(0, sum(is.na(X[,-1])))

  AUC.avg <- .8
  ACC.avg <- .8
  h.avg <- .05

 
 ################
 ## Begin MCMC ##
 ################
 cat("\n")
  for(it in 1:N.mcmc){

#if(it%%nplot==0){
#
#if(it%%10==0){
  cat("\nIteration", it, "out of", N.mcmc)
#}

#print("y update")

#y<<-y
#Xc.now<<-Xc.now
#beta.now<<-beta.now
#phi.now<<-phi.now

    yc.now <- update.y.probit(y, Xc.now, beta.now, phi.now)


#print("beta update")

#y<<-y
#yc.now<<-yc.now
#Xc.now<<-Xc.now
#beta.now<<-beta.now
#phi.now<<-phi.now
#Beta.now<<-Beta.now
#tau2.now<<-tau2.now
#groups<<-groups
#mult1<<-mult1
#rho<<-rho
#rho.prop.01<<-rho.prop.01
#rho.prop.10<<-rho.prop.10
#one.beta<<-one.beta

    ans.beta <- update.beta.probit(yc.now, Xc.now, beta.now, phi.now, Beta.now, tau2.now, groups, mult1, rho, rho.prop.01, rho.prop.10, one.beta)
    beta.now <- ans.beta$beta
    accept.beta[it,] <- ans.beta$accept

#print("B update")

    Beta.now <- update.Beta(beta.now, groups, tau2.now, one.tau, M.Beta, S2.Beta, one.beta)

#print("tau^2 update")

#beta.now<<-beta.now
#groups<<-groups
#Beta.now<<-Beta.now
#A.tau2<<-A.tau2
#B.tau2<<-B.tau2
#one.tau<<-one.tau


    tau2.now <- update.tau2(beta.now, groups, Beta.now, A.tau2, B.tau2, one.tau, one.beta)


#print("X and Z updates")

#X<<-X
#Xc.now<<-Xc.now
#Z.now<<-Z.now
#bin.X<<-bin.X
#muX.now<<-muX.now
#SigmaX.now<<-SigmaX.now
#yc.now<<-yc.now
#beta.now<<-beta.now
#mult2<<-mult2
#omega.now<<-omega.now
#groups<<-groups
#inv.groups<<-inv.groups
#nb<<-nb

  ## Update X and Z only every everyXZ iterations
    if(it%%everyXZ==0){
      ans.XZ <- update.XZ.probit(X, Xc.now, Z.now, bin.X, muX.now, SigmaX.now, yc.now, beta.now, mult2, omega.now, groups, inv.groups, nb)
      Xc.now <- ans.XZ$Xc
      Z.now <- ans.XZ$Z
      phi.now <- ans.XZ$phi.now
      accept.XZ <- accept.XZ + ans.XZ$accept
    }


#phi.now<<-phi.now
#delta.now<<-delta.now
#M<<-M

## Update vee/omega
    ans.vee <- update.vee(phi.now, delta.now, M)
    vee.now <- ans.vee$vee
    omega.now <- ans.vee$omega


  ## Update delta
    delta.now <- update.delta(vee.now, A.delta, B.delta)



#print("Omega and eta update")

    Omega.now <- update.Omega(SigmaX.now, eta.now, P.Omega, N.Omega)

#SigmaX.now<<-SigmaX.now
#eta.now<<-eta.now
#Omega.now<<-Omega.now
#A.eta<<-A.eta
#B.eta<<-B.eta
#prop.sd.eta<<-prop.sd.eta
#dfp<<-dfp

    ans.eta <- update.eta(SigmaX.now, eta.now, Omega.now, A.eta, B.eta, prop.sd.eta, dfp, etagp)
    eta.now <- ans.eta$eta
    accept.eta[it] <- ans.eta$accept

#print("SigmaX update")


#SigmaX.now<<-SigmaX.now
#Z.now<<-Z.now
#muX.now<<-muX.now
#phi.now<<-phi.now
#Omega.now<<-Omega.now
#eta.now<<-eta.now
#bin.X<<-bin.X

    SigmaX.now <- update.SigmaX(SigmaX.now, Z.now, muX.now, phi.now, Omega.now, eta.now, bin.X)



#print("muX update")

#Z.now<<-Z.now
#SigmaX.now<<-SigmaX.now
#phi.now<<-phi.now
#nu.now<<-nu.now
#Psi.now<<-Psi.now

    muX.now <- update.muX(Z.now, SigmaX.now, phi.now, nu.now, Psi.now, ind.scale)


   ## scale muX and SigmaX
    OmegaX.now <- SigmaX.now
    nuX.now <- muX.now
    for(m in 1:M){
      scale.m <- diag(1/sqrt(diag(SigmaX.now[m,,])))
      diag(scale.m)[!ind.scale] <- 1
      OmegaX.now[m,,] <- scale.m%*%SigmaX.now[m,,]%*%t(scale.m)
      nuX.now[m,] <- muX.now[m,]%*%scale.m
    }

#print("nu and Psi update")

    nu.now <- update.nu(nuX.now, Psi.now, M.nu, S.nu)
    Psi.now <- update.Psi(nuX.now, nu.now, P.Psi, N.Psi)


#print("End of Updates")

    for(m in 1:M){
      ind.m <- which(phi.now==m)
      risk.now[ind.m] <- as.numeric(Xc.now[ind.m,]%*%beta.now[m,])
    }

beta..<<-beta
beta.now..<<-beta.now


   ## record params.now in posterior sample
    if(it>begin && (it-begin)%%every==0){

      Beta[(it-begin)/every,] <- Beta.now*as.numeric(beta.now[1,]!=0)
      beta[(it-begin)/every,,] <- beta.now
      tau2[(it-begin)/every,] <- tau2.now
      omega[(it-begin)/every,] <- omega.now
      vee[(it-begin)/every,] <- vee.now
      muX[(it-begin)/every,,] <- nuX.now
      SigmaX[(it-begin)/every,,,] <- OmegaX.now
      lambda[(it-begin)/every,] <- risk.now
    }
    

   ## Summarize and Plot posterior
    if(!is.null(yNA) && it%%5==0){
      ind.na <- which(is.na(y))
      ans.ROC <- get.ROC(yNA, risk.now[ind.na])
      cat("\n\nk =",fold,", AUC (CV) =",ans.ROC$auc)
    }
    if(it>=nplot && (it-begin)%%nplot==0){
      ind.com <- which(!is.na(y))
      ans.ROC <- get.ROC(y[ind.com], risk.now[ind.com])
      AUC <- ans.ROC$auc
      ACC <- ans.ROC$acc
      AUC.avg <- (1-h.avg)*AUC.avg + h.avg*AUC
      ACC.avg <- (1-h.avg)*ACC.avg + h.avg*ACC

      pp <- min(maxplot,p)
      if(one.beta)
        N.plots <- (pp+1)+nG*(1-one.tau)+one.tau+1+2
      else
        N.plots <- 2*(pp+1)+nG*(1-one.tau)+one.tau+1+2
      cols <- min(13, ceiling(sqrt(N.plots)))
      rows <- min(13, ceiling(N.plots/cols))
      it.e <- floor((it-begin)/every)
      ind.now <- max(1,floor(it.e/2),it.e-nback+1):max(1,it.e)
      par(mfrow=c(rows,cols), mar=c(2,2,2,1))

     ## Print Rsq and recent average
      cat("\n\nAUC =",AUC)
      cat("\nAUC (recent avg) =",AUC.avg)
      cat("\nACC  =",ACC)
      cat("\nACC (recent avg) =",ACC.avg,"\n")
      like.y <- get.like.probit(y[ind.com], Xc.now[ind.com,], beta.now, phi.now[ind.com])
      cat("\nlog(like) =",like.y,"\n")

     ## Print acceptance %
      cat("\neta acceptance = ", mean(accept.eta[1:it]))
      cat("\nbeta acceptance = \n")
      print(summary(colMeans(accept.beta[1:it,,drop=FALSE])))
      cat("\nnX/Z acceptance = \n")
      print(summary(accept.XZ/floor(it/everyXZ)))
      cat("\n\n\n")
    }
    if((it-begin)%%nplot==0 && it>begin){
     ## Plot Betas
      if(one.beta){
        Beta.ord <- c(1,order(-abs(colMeans(beta[ind.now,1,-1])))+1)
        for(k in 1:(pp+1))
          plot(beta[ind.now,1,Beta.ord[k]],ylab="",main=paste("Beta_",Beta.ord[k]-1,sep=""),cex=.5)
      }
      else{
        Beta.ord <- c(1,order(-colMeans(abs(beta.now[,-1])*omega.now))+1)
        for(k in 1:(pp+1)){
          plot(Beta[ind.now,Beta.ord[k]],ylab="",main=paste("Beta_",Beta.ord[k]-1,sep=""),cex=.5)
          beta.k.hist <- rep(beta.now[,Beta.ord[k]], round(omega.now*100))
	  hist(beta.k.hist, main=paste("beta_",Beta.ord[k]-1,sep=""),cex=.5)
	}
      }
     

    ## Bar plot of omega.now
      barplot(omega.now, main="omega")
      abline(h=exp(-6), col=2)
      abline(h=exp(-4), col=4)

      lomega <- -sort(-log(omega.now+1E-300))
      barplot(lomega, main="omega", yaxt='n')
      at <- ceiling(1/10*round(seq(min(lomega),0,length=4),0))*10
      labels <- format(exp(at),nsmall=1,digits=1)
      axis(2,at=at, labels=labels)
      abline(h=-6, col=2)
      abline(h=-4, col=4)
           
    }
  }
  return(list(beta=beta, Beta=Beta, tau2=tau2, muX=muX, SigmaX=SigmaX, omega=omega, vee=vee, X=X, y=y, bin.X=bin.X, lambda=lambda))
}
















#############################################################################
############# missForest altered to allow for variable labels ###############
#############################################################################



mymissForest <- function (xmis, maxiter = 10, ntree = 100, variablewise = FALSE, 
    decreasing = FALSE, verbose = FALSE, mtry = floor(sqrt(ncol(xmis))), 
    replace = TRUE, classwt = NULL, cutoff = NULL, strata = NULL, 
    sampsize = NULL, nodesize = NULL, maxnodes = NULL, xtrue = NA, 
    parallelize = c("no", "variables", "forests")) 
{
    n <- nrow(xmis)
    p <- ncol(xmis)
    if (!is.null(classwt)) 
        stopifnot(length(classwt) == p, typeof(classwt) == "list")
    if (!is.null(cutoff)) 
        stopifnot(length(cutoff) == p, typeof(cutoff) == "list")
    if (!is.null(strata)) 
        stopifnot(length(strata) == p, typeof(strata) == "list")
    if (!is.null(nodesize)) 
        stopifnot(length(nodesize) == 2)
    if (any(apply(is.na(xmis), 2, sum) == n)) {
        indCmis <- which(apply(is.na(xmis), 2, sum) == n)
        xmis <- xmis[, -indCmis]
        p <- ncol(xmis)
        cat("  removed variable(s)", indCmis, "due to the missingness of all entries\n")
    }
    parallelize <- match.arg(parallelize)
    if (parallelize %in% c("variables", "forests")) {
        if (getDoParWorkers() == 1) {
            stop("You must register a 'foreach' parallel backend to run 'missForest' in parallel. Set 'parallelize' to 'no' to compute serially.")
        }
        else if (verbose) {
            if (parallelize == "variables") {
                cat("  parallelizing over the variables of the input data matrix 'xmis'\n")
            }
            else {
                cat("  parallelizing computation of the random forest model objects\n")
            }
        }
        if (getDoParWorkers() > p) {
            stop("The number of parallel cores should not exceed the number of variables (p=", 
                p, ")")
        }
    }
    ximp <- xmis
    xAttrib <- lapply(xmis, attributes)
    varType <- character(p)
    for (t.co in 1:p) {
        if (is.null(xAttrib[[t.co]]$levels) && is.null(xAttrib[[t.co]]$levels) ) {
            varType[t.co] <- "numeric"
            ximp[is.na(xmis[, t.co]), t.co] <- mean(xmis[, t.co], 
                na.rm = TRUE)
        }
        else {
            varType[t.co] <- "factor"
            max.level <- max(table(ximp[, t.co]))
            class.assign <- sample(names(which(max.level == summary(ximp[, 
                t.co]))), 1)
            if (class.assign != "NA's") {
                ximp[is.na(xmis[, t.co]), t.co] <- class.assign
            }
            else {
                while (class.assign == "NA's") {
                  class.assign <- sample(names(which(max.level == 
                    summary(ximp[, t.co]))), 1)
                }
                ximp[is.na(xmis[, t.co]), t.co] <- class.assign
            }
        }
    }
    NAloc <- is.na(xmis)
    noNAvar <- apply(NAloc, 2, sum)
    sort.j <- order(noNAvar)
    if (decreasing) 
        sort.j <- rev(sort.j)
    sort.noNAvar <- noNAvar[sort.j]
    nzsort.j <- sort.j[sort.noNAvar > 0]
    if (parallelize == "variables") {
        "%cols%" <- get("%dopar%")
        idxList <- as.list(isplitVector(nzsort.j, chunkSize = getDoParWorkers()))
    }
    Ximp <- vector("list", maxiter)
    iter <- 0
    k <- length(unique(varType))
    convNew <- rep(0, k)
    convOld <- rep(Inf, k)
    OOBerror <- numeric(p)
    names(OOBerror) <- varType
    if (k == 1) {
        if (unique(varType) == "numeric") {
            names(convNew) <- c("numeric")
        }
        else {
            names(convNew) <- c("factor")
        }
        convergence <- c()
        OOBerr <- numeric(1)
    }
    else {
        names(convNew) <- c("numeric", "factor")
        convergence <- matrix(NA, ncol = 2)
        OOBerr <- numeric(2)
    }
    stopCriterion <- function(varType, convNew, convOld, iter, 
        maxiter) {
        k <- length(unique(varType))
        if (k == 1) {
            (convNew < convOld) & (iter < maxiter)
        }
        else {
            ((convNew[1] < convOld[1]) | (convNew[2] < convOld[2])) & 
                (iter < maxiter)
        }
    }
    while (stopCriterion(varType, convNew, convOld, iter, maxiter)) {
        if (iter != 0) {
            convOld <- convNew
            OOBerrOld <- OOBerr
        }
        cat("  missForest iteration", iter + 1, "in progress...")
        t.start <- proc.time()
        ximp.old <- ximp
        if (parallelize == "variables") {
            for (idx in idxList) {
                results <- foreach(varInd = idx, .packages = "randomForest") %cols% 
                  {
                    obsi <- !NAloc[, varInd]
                    misi <- NAloc[, varInd]
                    obsY <- ximp[obsi, varInd]
                    obsX <- ximp[obsi, seq(1, p)[-varInd]]
                    misX <- ximp[misi, seq(1, p)[-varInd]]
                    typeY <- varType[varInd]
                    if (typeY == "numeric") {
                      RF <- randomForest(x = obsX, y = obsY, 
                        ntree = ntree, mtry = mtry, replace = replace, 
                        sampsize = if (!is.null(sampsize)) 
                          sampsize[[varInd]]
                        else if (replace) 
                          nrow(obsX)
                        else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                          nodesize[1]
                        else 1, maxnodes = if (!is.null(maxnodes)) 
                          maxnodes
                        else NULL)
                      oerr <- RF$mse[ntree]
                      misY <- predict(RF, misX)
                    }
                    else {
                      obsY <- factor(obsY)
                      summarY <- summary(obsY)
                      if (length(summarY) == 1) {
                        oerr <- 0
                        misY <- factor(rep(names(summarY), length(misi)))
                      }
                      else {
                        RF <- randomForest(x = obsX, y = obsY, 
                          ntree = ntree, mtry = mtry, replace = replace, 
                          classwt = if (!is.null(classwt)) 
                            classwt[[varInd]]
                          else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                            cutoff[[varInd]]
                          else rep(1/nlevels(obsY), nlevels(obsY)), 
                          strata = if (!is.null(strata)) 
                            strata[[varInd]]
                          else obsY, sampsize = if (!is.null(sampsize)) 
                            sampsize[[varInd]]
                          else if (replace) 
                            nrow(obsX)
                          else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                            nodesize[2]
                          else 5, maxnodes = if (!is.null(maxnodes)) 
                            maxnodes
                          else NULL)
                        oerr <- RF$err.rate[[ntree, 1]]
                        misY <- predict(RF, misX)
                      }
                    }
                    list(varInd = varInd, misY = misY, oerr = oerr)
                  }
                for (res in results) {
                  misi <- NAloc[, res$varInd]
                  ximp[misi, res$varInd] <- res$misY
                  OOBerror[res$varInd] <- res$oerr
                }
            }
        }
        else {
            for (s in 1:p) {
                varInd <- sort.j[s]
                if (noNAvar[[varInd]] != 0) {
                  obsi <- !NAloc[, varInd]
                  misi <- NAloc[, varInd]
                  obsY <- ximp[obsi, varInd]
                  obsX <- ximp[obsi, seq(1, p)[-varInd]]
                  misX <- ximp[misi, seq(1, p)[-varInd]]
                  typeY <- varType[varInd]
                  if (typeY == "numeric") {
                    if (parallelize == "forests") {
                      xntree <- NULL
                      RF <- foreach(xntree = idiv(ntree, chunks = getDoParWorkers()), 
                        .combine = "combine", .multicombine = TRUE, 
                        .packages = "randomForest") %dopar% {
                        randomForest(x = obsX, y = obsY, ntree = xntree, 
                          mtry = mtry, replace = replace, sampsize = if (!is.null(sampsize)) 
                            sampsize[[varInd]]
                          else if (replace) 
                            nrow(obsX)
                          else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                            nodesize[1]
                          else 1, maxnodes = if (!is.null(maxnodes)) 
                            maxnodes
                          else NULL)
                      }
                      OOBerror[varInd] <- mean((predict(RF) - 
                        RF$y)^2, na.rm = TRUE)
                    }
                    else {
                      RF <- randomForest(x = obsX, y = obsY, 
                        ntree = ntree, mtry = mtry, replace = replace, 
                        sampsize = if (!is.null(sampsize)) 
                          sampsize[[varInd]]
                        else if (replace) 
                          nrow(obsX)
                        else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                          nodesize[1]
                        else 1, maxnodes = if (!is.null(maxnodes)) 
                          maxnodes
                        else NULL)
                      OOBerror[varInd] <- RF$mse[ntree]
                    }
                    misY <- predict(RF, misX)
                  }
                  else {
                    obsY <- factor(obsY)
                    summarY <- summary(obsY)
                    if (length(summarY) == 1) {
                      misY <- factor(rep(names(summarY), sum(misi)))
                    }
                    else {
                      if (parallelize == "forests") {
                        RF <- foreach(xntree = idiv(ntree, chunks = getDoParWorkers()), 
                          .combine = "combine", .multicombine = TRUE, 
                          .packages = "randomForest") %dopar% 
                          {
                            randomForest(x = obsX, y = obsY, 
                              ntree = xntree, mtry = mtry, replace = replace, 
                              classwt = if (!is.null(classwt)) 
                                classwt[[varInd]]
                              else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                                cutoff[[varInd]]
                              else rep(1/nlevels(obsY), nlevels(obsY)), 
                              strata = if (!is.null(strata)) 
                                strata[[varInd]]
                              else obsY, sampsize = if (!is.null(sampsize)) 
                                sampsize[[varInd]]
                              else if (replace) 
                                nrow(obsX)
                              else ceiling(0.632 * nrow(obsX)), 
                              nodesize = if (!is.null(nodesize)) 
                                nodesize[2]
                              else 5, maxnodes = if (!is.null(maxnodes)) 
                                maxnodes
                              else NULL)
                          }
                        ne <- as.integer(predict(RF)) != as.integer(RF$y)
                        ne <- ne[!is.na(ne)]
                        OOBerror[varInd] <- sum(ne)/length(ne)
                      }
                      else {
                        RF <- randomForest(x = obsX, y = obsY, 
                          ntree = ntree, mtry = mtry, replace = replace, 
                          classwt = if (!is.null(classwt)) 
                            classwt[[varInd]]
                          else rep(1, nlevels(obsY)), cutoff = if (!is.null(cutoff)) 
                            cutoff[[varInd]]
                          else rep(1/nlevels(obsY), nlevels(obsY)), 
                          strata = if (!is.null(strata)) 
                            strata[[varInd]]
                          else obsY, sampsize = if (!is.null(sampsize)) 
                            sampsize[[varInd]]
                          else if (replace) 
                            nrow(obsX)
                          else ceiling(0.632 * nrow(obsX)), nodesize = if (!is.null(nodesize)) 
                            nodesize[2]
                          else 5, maxnodes = if (!is.null(maxnodes)) 
                            maxnodes
                          else NULL)
                        OOBerror[varInd] <- RF$err.rate[[ntree, 
                          1]]
                      }
                      misY <- predict(RF, misX)
                    }
                  }
                  ximp[misi, varInd] <- misY
                }
            }
        }
        cat("done!\n")
        iter <- iter + 1
        Ximp[[iter]] <- ximp
        t.co2 <- 1
        for (t.type in names(convNew)) {
            t.ind <- which(varType == t.type)
            if (t.type == "numeric") {
                convNew[t.co2] <- sum((ximp[, t.ind] - ximp.old[, 
                  t.ind])^2)/sum(ximp[, t.ind]^2)
            }
            else {
                dist <- sum(as.character(as.matrix(ximp[, t.ind])) != 
                  as.character(as.matrix(ximp.old[, t.ind])))
                convNew[t.co2] <- dist/(n * sum(varType == "factor"))
            }
            t.co2 <- t.co2 + 1
        }
        if (!variablewise) {
            NRMSE <- sqrt(mean(OOBerror[varType == "numeric"])/var(as.vector(as.matrix(xmis[, 
                varType == "numeric"])), na.rm = TRUE))
            PFC <- mean(OOBerror[varType == "factor"])
            if (k == 1) {
                if (unique(varType) == "numeric") {
                  OOBerr <- NRMSE
                  names(OOBerr) <- "NRMSE"
                }
                else {
                  OOBerr <- PFC
                  names(OOBerr) <- "PFC"
                }
            }
            else {
                OOBerr <- c(NRMSE, PFC)
                names(OOBerr) <- c("NRMSE", "PFC")
            }
        }
        else {
            OOBerr <- OOBerror
            names(OOBerr)[varType == "numeric"] <- "MSE"
            names(OOBerr)[varType == "factor"] <- "PFC"
        }
        if (any(!is.na(xtrue))) {
            err <- suppressWarnings(mixError(ximp, xmis, xtrue))
        }
        if (verbose) {
            delta.start <- proc.time() - t.start
            if (any(!is.na(xtrue))) {
                cat("    error(s):", err, "\n")
            }
            cat("    estimated error(s):", OOBerr, "\n")
            cat("    difference(s):", convNew, "\n")
            cat("    time:", delta.start[3], "seconds\n\n")
        }
    }
    if (iter == maxiter) {
        if (any(is.na(xtrue))) {
            out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr)
        }
        else {
            out <- list(ximp = Ximp[[iter]], OOBerror = OOBerr, 
                error = err)
        }
    }
    else {
        if (any(is.na(xtrue))) {
            out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld)
        }
        else {
            out <- list(ximp = Ximp[[iter - 1]], OOBerror = OOBerrOld, 
                error = suppressWarnings(mixError(Ximp[[iter - 
                  1]], xmis, xtrue)))
        }
    }
    class(out) <- "missForest"
    return(out)
}

















