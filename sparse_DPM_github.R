
source("weibull_probit_DPM_4.R")
library("sas7bdat")
library("mclust")
library("clustvarsel")
library("rgl")
library("truncnorm")
library("doMC")
library("flexclust")
library("gtools")
registerDoMC() 



##############################################################
######### DPM CLustering w/ Variable Selection MCMC ##########
##############################################################


DPMvs.MCMC <- function(y, rho=.5, p.ad=.45, p.swap=.45, A.lambda=1, B.lambda=1, prop.sd.lambda=.5, P.Psi=NULL, N.Psi=NULL, A.eta=1, B.eta=.1, prop.sd.eta=.5, A.delta=1, B.delta=.1, prop.sd.delta=NULL, N.mcmc=5000, every=1, nplot=10, nback=1000, dfp=20, maxplot=100, gamma.init=NULL, phi.init=NULL, G=2:5, begin=0, bsy=3, prop.sd.y.m=.25, prop.sd.y.d=.25, prop.sd.y.c=.25, etagp=1, TT=2, Lg=200, Lgp=20, Lgy=0, Lpy=0, Lgpy=0, discrete=NULL, limits=NULL, prop.sd.y1=NULL, prop.sd.y2=NULL, yc.init=NULL, Ly.init=100, groups=NULL){

#############
## Inputs: ##
#############
#
# See "Clustering and Variable Selection in the Presence of Mixed Variable Types and Missing Data" by Storlie et.al. for more detailed definitions of variable names.
#
# y - a matrix or data frame n (observations) rows by p (variables) columns.  Recommend that y is standardized to have unit variances for each column.
# rho - the prior probability that a variable is informative.  Can be a scalar or a vector of values length p.
# limits - a p x 2 matrix of limits for each variable.  Defaults to -Inf, Inf.
# discrete - a vector of length p of 0s and 1s, to indicate continuous or discrete, respectively, for each variable.
# p.ad - probability of proposing the adding/deleting of a variable in the gamma update for informative variables
# p.swap - probability of proposing a swap move in the gamma update for informative variables
# A.lambda, B.lambda - Prior on lambda is Gamma(A.lambda, B.lambda)
# prop.sd.lambda - scale parameter for lambda proposals
# P.Psi, N.Psi - Prior for Psi is Wishart(P.Psi, N.Psi), defaults to a diagonal matrix for P.Psi of the variances for y (or identity if using standardized y as recommended).  N.Psi defaults to (p+1)*1.1
# 
# TT - number of restricted Gibbs sample updates to perform before proposing a split in the split/merge update
# A.eta, B.eta - Prior on eta is: eta = p + etagp + Gamma(A.lambda, B.lambda)
# prop.sd.eta - scale parameter for eta proposals
# A.delta, B.delta - Prior on delta is Gamma(A.lambda, B.lambda)
# prop.sd.delta - scale parameter for delta proposals
# N.mcmc - number of MCMC iterations to run
# begin - begin saving MCMC iterations at the 'begin'-th iteration.
# every - record/save the parameters every 'every' iterations.
# nplot - plot progress every 'nplot' iterations.
# nback - use samples from current iteration back to 'nback' iterations for progress plots
# dfp - degrees of freedom for scaled t distribution proposals.
# maxplot - maximum number of frames to plot in progress plots.
# gamma.init - starting values for gamma, defaults to updating Z several (Ly.init) times assuming one cluster, then running clustvarsel on the resulting Z.
# phi.init - starting values for phi; defaults to the resulting phi from Mclust after following the process above to obtain gamma.init
# G - the possible number of clusters to try for inititializing gamma and phi.  Ignored if phi.init is specified.
# bsy - number of observations to block update at one time when updating Z.
# prop.sd.y.m - scale parameter for proposals of Z for missing data elements.
# prop.sd.y.d - scale parameter for proposals for latent Z for observed, but discrete variables
# prop.sd.y.c - scale parameter for proposals for latent Z for observed, continuous, but censored variables
# etagp - see definition for A.eta, B.eta and the prior on eta above.
# Lg - the number of times to perform an MH update of gamma each sweep of the MCMC.
# Lgp - the number of times to perform an MH update of (gamma, phi), jointly, each sweep of the MCMC.
# Lgy - the number of times to perform an MH update of (gamma, Z), jointly, each sweep of the MCMC.
# Lpy - the number of times to perform an MH update of (phi, Z), jointly, each sweep of the MCMC.
# Lgpy - the number of times to perform an MH update of (gamma, phi, Z), jointly, each sweep of the MCMC.
# yc.init - initial values for Z (yc in the code is the latent Z in the paper).  Defaults to updating Z several (Ly.init) times assuming one cluster, with each missing value initialied as N(0, sd(y[,j])), each discrete variable set to its discrete value in y, and each censored value set to its cenored value in y.
# Ly.init - see yc.init and gamma.init
# groups - not supported, intended to specify distinct groups of variables, and at least one variable from each group must be chosen as informative.

##############
## Outputs: ##
##############

# gamma - N x p matrix, posterior sample of gamma, where N is the number of recorded MCMC iterations.
# phi - N x n matrix, posterior sample of phi.
# yc - N x M matrix, posterior sample for the latent values of y (i.e., \tilde{Z}).  M is the number of missing, discrete, or censored values in y.
# lambda - N-vector, posterior sample for lambda
# Psi - N x p x p array, posterior sample for Psi
# eta - N-vector, posterior sample for eta
# delta - N-vector, posterior sample for delta
# y, rho, p.ad, p.swap, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, prop.sd.y2, bsy, - the values that were input for these variables, respectively.



 ## Create needed variables ##
  n <- nrow(y)
  p <- ncol(y)
  if(is.null(prop.sd.y1))
    prop.sd.y1 <- prop.sd.y.m
  if(is.null(prop.sd.y2))
    prop.sd.y2 <- prop.sd.y.d
  if(is.null(discrete))
    discrete <- rep(FALSE,p)
  if(is.null(limits))
    limits <- cbind(rep(-Inf,p), rep(Inf,p))
  if(is.null(colnames(y)))
    colnames(y) <- paste("y_",1:p,sep="")
  if(length(rho)==1)
    rho <- rep(rho, p)
  if(is.null(groups))
    groups <- rep(1,p)

  if(is.null(N.Psi))
    N.Psi <- (p+1)*1.1
  if(is.null(P.Psi))
    P.Psi <- diag(apply(y,2,var,na.rm=TRUE))
  lambda.now <- A.lambda/B.lambda
  eta.now <- A.eta/B.eta + p + etagp
  delta.now <- A.delta/B.delta
  Psi.now <- P.Psi*N.Psi

 ## initialize yc, phi, gamma, mu, Sigma, Beta, Gamma
  cat("\nInitializing Clusters \n")
  ans.init <- initialize.clusters(y, G, lambda.now, eta.now, Psi.now, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, gamma.init, phi.init, discrete, limits, yc.init, Ly.init, groups)
  yc.now <- ans.init$yc.now
  phi.now <- ans.init$phi.now
  gamma.now <- ans.init$gamma.now

cat("\nAllocating Memory for Posterior Objects \n")

 ## Allocate posterior objects for which to store MCMC samples
  gamma <- matrix(0, (N.mcmc-begin)%/%every, p)
  phi <- matrix(0, (N.mcmc-begin)%/%every, n)
  lambda <- rep(0, (N.mcmc-begin)%/%every)
  Psi <- array(0, c((N.mcmc-begin)%/%every, p, p))
  eta <- rep(0, (N.mcmc-begin)%/%every)
  delta <- rep(0, (N.mcmc-begin)%/%every)
  NN <- length(y)
  yc <- matrix(0, (N.mcmc-begin)%/%every, NN)
  
  accept.gamma <- rep(0, N.mcmc)
  accept.phi <- rep(0, N.mcmc)
  accept.gp <- rep(0, N.mcmc)
  accept.gy <- rep(0, N.mcmc)
  accept.py <- rep(0, N.mcmc)
  accept.gpy <- rep(0, N.mcmc)
  accept.lambda <- rep(0, N.mcmc)
  accept.Psi <- rep(0, N.mcmc)
  accept.eta <- rep(0, N.mcmc)
  accept.delta <- rep(0, N.mcmc)

  nb <- ceiling(n/bsy)
  accept.yc.m <- rep(0, n)

 ################
 ## Begin MCMC ##
 ################
 cat("\n")
  for(it in 1:N.mcmc){

  cat("\nIteration", it, "out of", N.mcmc)


discrete..<<-discrete

#print("y update")
    ans.yc <- update.yc(yc.now, y, phi.now, gamma.now, lambda.now, eta.now, Psi.now, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, discrete, limits)
    yc.now <- ans.yc$yc.now
    accept.yc.m <- accept.yc.m + ans.yc$accept.m

#print("gamma update")
  ## Update gamma
  ans.gamma <- update.gamma(gamma.now, p.ad, p.swap, yc.now, phi.now, lambda.now, eta.now, Psi.now, rho, Lg, groups)
  gamma.now <- ans.gamma$gamma
  accept.gamma[it] <- ans.gamma$accept



#print("phi split/merge update")
#phi.now<<-phi.now
#yc.now<<-yc.now
#gamma.now<<-gamma.now
#lambda.now<<-lambda.now
#eta.now<<-eta.now
#Psi.now<<-Psi.now
#delta.now<<-delta.now

    ans.phi <- update.phi.sm(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, TT)
    phi.now <- ans.phi$phi
    accept.phi[it] <- ans.phi$accept

#print("gamma & y joint update")
    ans.gy <- update.gamma.y(gamma.now, p.ad, p.swap, yc.now, y, phi.now, lambda.now, eta.now, Psi.now, rho, Lgy, discrete, limits, prop.sd.y1, prop.sd.y2, groups)
    gamma.now <- ans.gy$gamma
    yc.now <- ans.gy$yc
    accept.gy[it] <- ans.gy$accept


#print("phi & y joint update")
    ans.py <- update.phi.y(gamma.now, phi.now, p.ad, p.swap, yc.now, y, lambda.now, eta.now, Psi.now, rho, Lpy, discrete, limits, prop.sd.y1, prop.sd.y2)
    phi.now <- ans.py$phi
    yc.now <- ans.py$yc
    accept.py[it] <- ans.py$accept


#print("gamma phi & y joint update")
#  if(runif(1)<1)
#    ans.gpy <- update.gamma.phi.y(gamma.now, phi.now, p.ad, p.swap, yc.now, y, lambda.now, eta.now, Psi.now, rho, delta.now, Lgpy, discrete, limits, prop.sd.y1, prop.sd.y2, groups)
#  else
#    ans.gpy <- update.gamma.phi.sm.y(gamma.now, phi.now, p.ad, p.swap, yc.now, y, lambda.now, eta.now, Psi.now, rho, delta.now, TT, discrete, limits, prop.sd.y1, prop.sd.y2, groups)
#    gamma.now <- ans.gpy$gamma
#    phi.now <- ans.gpy$phi
#    yc.now <- ans.gpy$yc
#    accept.gpy[it] <- ans.gpy$accept


#print("gamma & phi joint update")
    if(runif(1)<.5)
      ans.gp <- update.gamma.phi(gamma.now, phi.now, p.ad, p.swap, yc.now, lambda.now, eta.now, Psi.now, rho, delta.now, Lgp, groups)
    else
      ans.gp <- update.gamma.phi.sm(gamma.now, phi.now, p.ad, p.swap, yc.now, lambda.now, eta.now, Psi.now, rho, delta.now, TT, groups)
    gamma.now <- ans.gp$gamma
    phi.now <- ans.gp$phi
    accept.gp[it] <- ans.gp$accept



#print("phi Gibbs update")
    phi.now <- update.phi(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now)


#print("delta update")
    ans.delta <- update.delta(delta.now, phi.now, A.delta, B.delta, prop.sd.delta, dfp)
    delta.now <- ans.delta$delta.now
    accept.delta[it] <- ans.delta$accept


#print("lambda update")
    ans.lambda <- update.lambda(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, A.lambda, B.lambda, prop.sd.lambda, dfp)
    lambda.now <- ans.lambda$lambda.now
    accept.lambda[it] <- ans.lambda$accept


#print("Psi update")
    ans.Psi <- update.Psi(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, P.Psi, N.Psi)
    Psi.now <- ans.Psi$Psi.now
    accept.Psi[it] <- ans.Psi$accept


#print("eta update")
    ans.eta <- update.eta(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, A.eta, B.eta, prop.sd.eta, dfp, etagp)
    eta.now <- ans.eta$eta
    accept.eta[it] <- ans.eta$accept
    

#print("End of Updates")

   ## record params.now in posterior sample
    if(it>begin && (it-begin)%%every==0){

      yc[(it-begin)/every,] <- yc.now
      gamma[(it-begin)/every,] <- gamma.now
      phi[(it-begin)/every,] <- phi.now
      lambda[(it-begin)/every] <- lambda.now
      Psi[(it-begin)/every,,] <- Psi.now
      eta[(it-begin)/every] <- eta.now
      delta[(it-begin)/every] <- delta.now
     }

yc.now..<<-yc.now

   ## Summarize and Plot posterior
    if((it-begin)%%nplot==0){
      ncv <- sum(gamma.now)
      N.plots <- min(maxplot, max(1,choose(ncv,2)) + 2)
      cols <- min(13, ceiling(sqrt(N.plots)))
      rows <- min(13, ceiling(N.plots/cols))
      it.e <- floor((it-begin)/every)
      ind.now <- max(1,floor(it.e/2),it.e-nback+1):max(1,it.e)
      par(mfrow=c(rows,cols), mar=c(2,2,2,1))

     ## Print likelihood
      like.y <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
      cat("\nlog(like) =",like.y,"\n")

     ## Print model
      cat("\ngamma =",gamma.now,"\n")

     ## Print acceptance %
      cat("\ngamma acceptance = ", mean(accept.gamma[1:it]))
      cat("\nphi acceptance = ", mean(accept.phi[1:it]))
      cat("\ngamma & y acceptance = ", mean(accept.gy[1:it]))
      cat("\nphi & y acceptance = ", mean(accept.py[1:it]))
      cat("\ngamma phi & y acceptance = ", mean(accept.gpy[1:it]))
      cat("\ngamma & phi acceptance = ", mean(accept.gp[1:it]))
      cat("\nlambda acceptance = ", mean(accept.lambda[1:it]))
      cat("\nPsi acceptance = ", mean(accept.Psi[1:it]))
      cat("\neta acceptance = ", mean(accept.eta[1:it]))
      cat("\ndelta acceptance = ", mean(accept.delta[1:it]))
      cat("\nyc acceptance = ", summary(accept.yc.m/it))
      cat("\n\n\n")
    }
    if((it-begin)%%nplot==0 && it>begin){
      par(mar=c(1,1.5,1.5,.5))
     ## Plot clusters in pairwise scatter of gamma=1 vars
      ind.cv <- which(gamma.now==1)
      if(length(ind.cv) <= 6)
        names.cv <- colnames(y)[ind.cv]
      else
        names.cv <- paste("y",ind.cv,sep="")
      n.cv <- length(ind.cv)
      indp <- 1
      if(n.cv==1){
        ind2 <- ifelse(ind.cv==1,2,1)
        plot(yc.now[,ind.cv], yc.now[,ind2], col=phi.now, main=paste(names.cv[1],"by",colnames(y)[ind2]),pch=16)
      }
      else{
        for(j in 1:(n.cv-1)){
          for(k in (j+1):n.cv){
            if(indp <= N.plots-2)
	      plot(yc.now[,ind.cv[j]], yc.now[,ind.cv[k]], col=phi.now, pch=16, cex=.7, cex.axis=.8, mgp=c(2,.35,0))
	      title(main=paste(names.cv[j],"by",names.cv[k]),line=.5, cex.main=.8)
            indp <- indp + 1	  
	  }
        }
      }
     ## Bar plot of omega.now
      omega.now <- update.vee(phi.now, delta.now, M=max(phi.now)+1)$omega
      barplot(omega.now, main="omega")
      abline(h=exp(-6), col=2)
      abline(h=exp(-4), col=4)           
    }
  }
  return(list(gamma=gamma, phi=phi, yc=yc, lambda=lambda, Psi=Psi, eta=eta, delta=delta, y=y, rho=rho, p.ad=p.ad, p.swap=p.swap, prop.sd.y.m=prop.sd.y.m, prop.sd.y.d=prop.sd.y.d, prop.sd.y.c=prop.sd.y.c, prop.sd.y2=prop.sd.y2, bsy=bsy))
}




####################################################################
########### Get gamma hat and phi.hat for an MCMC output ###########
####################################################################


get.gamma.probs <- function(ans.DPMvs, post.ind=NULL){

  p <- dim(ans.DPMvs$gamma)[2]
  N.mcmc <- dim(ans.DPMvs$gamma)[1]
  if(is.null(post.ind))
    post.ind <- 1:N.mcmc
  gamma.probs <- apply(ans.DPMvs$gamma[post.ind,],2,mean)
  return(gamma.probs)
}



get.phi.hat <- function(ans.DPMvs, post.ind=NULL, M=NULL, pw.prob=NULL){

  n <- dim(ans.DPMvs$phi)[2]
  N.mcmc <- dim(ans.DPMvs$phi)[1]
  if(is.null(post.ind))
    post.ind <- 1:N.mcmc
  N.mcmc <- length(post.ind)
  phi <- ans.DPMvs$phi[post.ind,]
  if(!is.null(M)){
    foo <- apply(allDPM$phi,1,max)
    ind.M <- which(foo==M)
  }
  else{
    ind.M <- 1:dim(allDPM$phi)[1]
  }
  if(is.null(pw.prob)){
    pw.prob <- foreach(i=1:n,.combine=rbind)%dopar%{
      if(i%%10==0)
        print(i)
      ans.i <- rep(0,n)
      for(j in 1:n)
        ans.i[j] <- mean(phi[ind.M,i]==phi[ind.M,j])
      ans.i
    }
  }
  dist <- foreach(k=1:N.mcmc,.combine=rbind)%dopar%{
    if(k%%100==0)
      print(k)
    if(!is.null(M) && max(phi[k,])!=M){
      A.it <- 1E6
    }
    else{
      A.it <- matrix(1,n,n)
      for(i in 1:(n-1))
        for(j in (i+1):n)
          A.it[i,j] <- A.it[j,i] <- (phi[k,i]==phi[k,j])*1
    }
    c(sum(abs(A.it-pw.prob)), sum((A.it-pw.prob)^2))
  }
  ind1.min <- which(dist[,1]==min(dist[,1]))[1]
  ind2.min <- which(dist[,2]==min(dist[,2]))[1]
  phi.hat1 <- phi[ind1.min,]
  phi.hat2 <- phi[ind2.min,]
  
  phi.hat2r <- phi.hat2
  phi.hat1r <- try(refine.phi.hat(phi.hat1, pw.prob, ans.DPMvs, maxit=100, Nmax=5, pow=1), silent=TRUE)
  if(is.character(phi.hat1r[1]))
    phi.hat1r <- phi.hat1
  phi.hat2r <- try(refine.phi.hat(phi.hat2, pw.prob, ans.DPMvs, maxit=100, Nmax=5, pow=2), silent=TRUE)
  if(is.character(phi.hat2r[1]))
    phi.hat2r <- phi.hat2
    
  return(list(phi.hat1=phi.hat1, phi.hat2=phi.hat2, phi.hat1r=phi.hat1r, phi.hat2r=phi.hat2r, ind1.min=post.ind[ind1.min], ind2.min=post.ind[ind2.min], pw.prob=pw.prob))
}




refine.phi.hat <- function(phi.hat, pw.prob, ans.DPMvs, maxit=20, Nmax=5, pow=2){

## loop through each obs and consider it's best possible switch, according to sum((A.it-pw.prob)^2)
## keep the best switch , iterate to convergence

  n <- dim(ans.DPMvs$phi)[2]
  M <- max(phi.hat)
  phi.now <- phi.best <- phi.hat
  Nneg <- 0
  avail <- 1:n
  acc.dist <- 0
  
  for(it in 1:maxit){

print(it)

    A.now <- matrix(1,n,n)
      for(i in 1:(n-1))
        for(j in (i+1):n)
          A.now[i,j] <- A.now[j,i] <- (phi.now[i]==phi.now[j])*1
    dist.now <- sum(abs(A.now-pw.prob)^pow)
    m.dist <- foreach(i=1:n, .combine=rbind)%dopar%{
      if(!any(avail==i))
        ans.i <- c(phi.now[i], -Inf)
      else{
        dist.i.prev <- sum(abs(A.now[i,]-pw.prob[i,])^pow)
        dist.i <- Inf
        for(m in (1:M)[-phi.now[i]]){
          a.im <- (phi.now==m)*1
          a.im[i] <- 1
          dist.im <- sum(abs(a.im-pw.prob[i,])^pow)
          if(dist.im < dist.i){
            mi.now <- m
	    dist.i <- dist.im
          }
        }
        ans.i <- c(mi.now, dist.i.prev-dist.i)
      }
      ans.i
    }
    ind.i <- which(m.dist[,2]==max(m.dist[,2]))[1]
    phi.now[ind.i] <- m.dist[ind.i,1]
    acc.dist <- acc.dist + m.dist[ind.i,2]
    if(acc.dist > 0){
      phi.best <- phi.now
      Nneg <- 0
      avail <- 1:n
      acc.dist <- 0
    }
    else{
      Nneg <- Nneg+1
      avail <- avail[-ind.i]
    }
    if(Nneg==Nmax)
      break
  }
  return(phi.best)
}




get.phi.hat2 <- function(phi.probs){

  n <- nrow(phi.probs)
  phi.hat2 <- rep(0,n)
  for(i in 1:n)
    phi.hat2[i] <- which(phi.probs[i,]==max(phi.probs[i,]))
  return(phi.hat2)
}



get.phi.probs.hat <- function(ans.DPMvs, y.init, phi.init=NULL, post.ind=NULL, probs.now=NULL, maxit=100, Mmin=1, Mmax=8, tol=1E-8){

  n <- dim(ans.DPMvs$phi)[2]
  N.mcmc <- dim(ans.DPMvs$phi)[1]
  if(is.null(post.ind))
    post.ind <- 1:N.mcmc
  foo <- apply(ans.DPMvs$phi[post.ind,],1,max)
  ind.keep <- which(foo<=Mmax & foo>=Mmin)
  post.ind <- post.ind[ind.keep]
  phi.now <- ans.DPMvs$phi[post.ind,]
  if(is.null(phi.init))
    phi.init <- phi.now[1,]
  m.rank <- order(-table(phi.init))
  phi.init <- m.rank[phi.init]

  M.init <- max(phi.init)
  Mmax <- max(phi.now)

  if(is.null(probs.now)){
   ## initialize probs.now
   ## initialize permutations of phi.now
    mu.init <- matrix(0,M.init,dim(y.init)[2])
    pi.init <- table(phi.init)/n
    for(m in 1:M.init)
      mu.init[m,] <- apply(y.init[phi.init==m,,drop=FALSE],2,mean)
    phi.now <- foreach(k=1:length(post.ind), .combine=rbind)%dopar%{
if(k%%1000==0)
  print(k)
      phi.k <- phi.now[k,]
      M.k <- max(phi.k)
      mu.k <- matrix(0,M.k,dim(y.init)[2])
      pi.k <- table(phi.k)/n
      for(m in 1:M.k)
        mu.k[m,] <- apply(y.init[phi.k==m,,drop=FALSE],2,mean)
      D1.k <- as.matrix(dist(rbind(mu.init,mu.k)))[1:M.init,(M.init+1):(M.init+M.k)]
      D2.k <- as.matrix(dist(cbind(c(pi.init,pi.k))))[1:M.init,(M.init+1):(M.init+M.k)]
      D1.k <- D1.k/max(D1.k)
      D2.k <- D2.k/max(D2.k)
      D.k <- D1.k+D2.k
      foo.k <- rep(0,M.k)
      avail <- rep(TRUE,M.k)
      for(m in 1:min(M.k,M.init)){
        foo.k[m] <- which(D.k[m,]==min(D.k[m,avail]))
        avail[foo.k[m]] <- FALSE
      }
      if(any(avail))
        foo.k <- c(foo.k,which(avail))
      perm.k <- 1:M.k
      perm.k[foo.k] <- 1:M.k
      perm.k[phi.k]
    }    
    probs.now <- matrix(0,n,Mmax)
    for(m in 1:Mmax)
      probs.now[,m] <- apply(phi.now==m, 2, mean)
  }

  for(it in index(1,maxit)){

print(it)

    phi.now <- foreach(k=1:length(post.ind), .combine=rbind)%dopar%{
if(k%%1000==0)
  print(k)
      phi.k <- phi.k.new <- phi.now[k,]
      M.k <- max(phi.k)
      perms.k <- permutations(n=Mmax, r=M.k)
      np <- nrow(perms.k)
      like.best <- -Inf
      for(l in 1:np){
        phi.kl <- perms.k[l,][phi.k]
	like.l <- sum(log(probs.now[cbind(1:n,phi.kl)]))
        if(like.l>like.best){
	  phi.k.new <- phi.kl
	  like.best <- like.l
	}
      }
      phi.k.new
    }
    probs.new <- matrix(0,n,Mmax)
    for(m in 1:Mmax)
      probs.new[,m] <- apply(phi.now==m, 2, mean)
    dist <- sum((probs.new-probs.now)^2)
    probs.now <- probs.new
    if(dist<tol)
      break
    if(it==maxit)
      warning("reached maxit")
  }
  return(list(phi.now=phi.now, probs.now=probs.now))
}












####################################
######### Stepwise Mclust ##########
####################################


get.acc <- function(phi, phi.hat, maxp=1E6){

  phi1 <- match(phi, unique(phi))
  phi2 <- match(phi.hat, unique(phi.hat))
  n <- length(phi1)
  M1 <- max(phi)
  M2 <- max(phi.hat)
  if(M1 > M2){
    phi1t <- phi1
    phi1 <- phi2
    phi2 <- phi1t
    M1t <- M1
    M1 <- M2
    M2 <- M1t
  }
  perms <- permutations(n=M2, r=M1)
  np <- nrow(perms)
  if(np > maxp)
    stop("np > maxp")
  acc.best <- 0
  for(i in 1:np){
    phi1.i <- perms[i,][phi1]
    acc.i <- sum(phi1.i==phi2)
    if(acc.i>acc.best)
      acc.best <- acc.i
  }
  return(acc.best/n)
}


mydtruncnorm <- function(x,a,b,mu,sd){
  ans0 <- dnorm(x,mu,sd,log=TRUE)
  ans1 <- pnorm(a,mu,sd,log=TRUE) - ans0
  ans2 <- pnorm(b,mu,sd,log=TRUE) - ans0
  ans <- 1/(exp(ans1)-exp(ans2))
  return(ans)
}




myjitter.old <- function(x, lcens=FALSE, rcens=FALSE){
  ux <- sort(unique(x))
  nu <- length(ux)
  n <- length(x)
  xt <- x
  for(i in 1:n){
    if(is.na(x[i]))
      next
    indu <- which(ux==x[i])
    a <- ifelse(indu==1, ux[1], (ux[indu-1]+ux[indu])/2) 
    b <- ifelse(indu==nu, ux[nu], (ux[indu]+ux[indu+1])/2)
    a <- ifelse(indu==nu && rcens, ux[nu], a)
    b <- ifelse(indu==1 && lcens, ux[1], b)    
    a <- ifelse(indu==2 && lcens, ux[1], a)
    b <- ifelse(indu==(nu-1) && rcens, ux[nu], b)
    xt[i] <- runif(1,a,b)
  }
  return(xt)
}



myjitter <- function(x, lcens=FALSE, rcens=FALSE){
  ux <- sort(unique(x))
  nu <- length(ux)
  n <- length(x)
  xt <- x
  for(i in 1:n){
    if(is.na(x[i]))
      next
    indu <- which(ux==x[i])
    a <- ifelse(indu==1, 2*ux[1]-ux[2], ux[indu-1]) 
    b <- ifelse(indu==nu, 2*ux[nu]-ux[nu-1], ux[indu+1])
    a <- ifelse(indu==nu && rcens, ux[nu], a)
    b <- ifelse(indu==1 && lcens, ux[1], b)
    a <- ifelse(indu==1 && lcens, ux[1], a)
    b <- ifelse(indu==nu && rcens, ux[nu], b)
    if(a==b)
      xt[i] <- a
    else
      xt[i] <- rtruncnorm(1,a,b,(a+b)/2,(b-a)/2)
  }
  return(xt)
}



lmgamma <- function(p,x){
  ans <- p*(p-1)/4*log(pi) + sum(lgamma(x + (1-1:p)/2))
  return(ans)
}





get.bic <- function(ans.lm){
  p <- length(ans.lm$coefficients)+1
  n <- length(ans.lm$fitted)
  RSS <- sum(ans.lm$residuals^2)
  bic <- -n*log(2*pi) - n*log(RSS/n) - n - p*log(n)
  return(bic)
}





sw.mclust.par <- function(X, G=1:9, modelNames=NULL, maxp=20, maxiter=30){

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  Xr.now <- matrix(0,n,0)
  ind.keep <- c()
  ind.avail <- 1:p
  next.var <- 0

## Select the first clustering vars
  bic.inc <- -Inf
  for(j in 1:p){
    ans.j <- Mclust(X[,j],G=G,modelNames="V")
    lm.j <- lm.fit(cbind(1,Xr.now), X[,j])
    inc.j <- ans.j$bic - get.bic(lm.j)
    if(inc.j > bic.inc){
      next.var <- j
      bic.inc <- inc.j
      bic.prev <- ans.j$bic
    }
  }
  Xr.now <- X[,next.var,drop=FALSE]
  ind.avail <- ind.avail[ind.avail!=next.var]
  ind.keep <- c(ind.keep, next.var)
  print(ind.keep)
  
## Stepwise addition/removal from now on
  flag1 <- flag2 <- 1
  iter <- 1
  while(flag1 || flag2){
    flag1 <- 0
    flag2 <- 0
    if(iter > maxiter || length(ind.keep)==maxp)
      break
    if(iter==1)
      bic.inc <- -Inf
    else
      bic.inc <- 0
   ## possible addition
    incs <- foreach(j=ind.avail, .combine=rbind)%dopar%{
      X.j <- cbind(Xr.now, X[,j])
      ans.j <- Mclust(X.j,G=G,modelNames=modelNames)
      lm.j <- lm.fit(cbind(1,Xr.now), X[,j])
      inc.j <- ans.j$bic - (bic.prev + get.bic(lm.j))
      cbind(inc.j, ans.j$bic)
    }
    if(length(incs)>0 && max(incs[,1]) > bic.inc){
      next.var <- ind.avail[which(incs[,1]==max(incs[,1]))[1]]
      bic.next <- incs[which(incs[,1]==max(incs[,1]))[1],2]
      flag1 <- 1
    }
    if(flag1){
      ind.avail <- ind.avail[ind.avail!=next.var]
      ind.keep <- sort(c(ind.keep, next.var))
      Xr.now <- X[,ind.keep,drop=FALSE]
      bic.prev <- bic.next      
    }
    print(ind.keep)
    
   ## possible deletion
    bic.inc <- 0
    for(j in 1:length(ind.keep)){
      X.mj <- as.matrix(Xr.now[,-j])
      if(ncol(X.mj)>1)
        ans.j <- Mclust(X.mj,G=G,modelNames=modelNames)
      else
        ans.j <- Mclust(X.mj,G=G,modelNames=NULL)
      lm.j <- lm.fit(cbind(1,X.mj), Xr.now[,j])
      inc.j <- ans.j$bic + get.bic(lm.j) - bic.prev
      if(inc.j > bic.inc){
        rm.var <- j
        bic.inc <- inc.j
        bic.next <- ans.j$bic
	flag2 <- 1
      }
    }
    if(flag2){
      ind.avail <- sort(c(ind.avail,ind.keep[rm.var]))
      ind.keep <- ind.keep[-rm.var]
      Xr.now <- X[,ind.keep,drop=FALSE]
      bic.prev <- bic.next
    }
    print(ind.keep)
    iter <- iter + 1
  }
  return(ind.keep)
}








sw.mclust <- function(X, G=1:9, modelNames=NULL, maxp=20, maxiter=30, back=TRUE){

  n <- nrow(X)
  p <- ncol(X)
  X <- as.matrix(X)

  Xr.now <- matrix(0,n,0)
  ind.keep <- c()
  ind.avail <- 1:p
  next.var <- 0

## Select the first clustering vars
  bic.inc <- -Inf
  for(j in 1:p){
    ans.j <- Mclust(X[,j],G=G,modelNames=NULL)
    lm.j <- lm.fit(cbind(1,Xr.now), X[,j])
    inc.j <- ans.j$bic - get.bic(lm.j)
    if(inc.j > bic.inc){
      next.var <- j
      bic.inc <- inc.j
      bic.prev <- ans.j$bic
    }
  }
  Xr.now <- X[,next.var,drop=FALSE]
  ind.avail <- ind.avail[ind.avail!=next.var]
  ind.keep <- next.var
  print(ind.keep)
 
## Stepwise addition/removal from now on
  flag1 <- flag2 <- 1
  iter <- 1
  while(flag1 || flag2){
    flag1 <- 0
    flag2 <- 0
    if(iter > maxiter || length(ind.keep)==maxp)
      break
    if(iter==1)
      bic.inc <- -Inf
    else
      bic.inc <- 0

   ## possible addition
    for(j in ind.avail){
#print(j)
      X.j <- cbind(X[,j], Xr.now)
      ans.j <- Mclust(X.j,G=G,modelNames=modelNames)
      lm.j <- lm.fit(cbind(1,Xr.now), X[,j])
      inc.j <- ans.j$bic - (bic.prev + get.bic(lm.j))
      if(inc.j > bic.inc){
        next.var <- j
        bic.inc <- inc.j
	bic.next <- ans.j$bic
	flag1 <- 1
      }
    }
    if(flag1){
      ind.avail <- ind.avail[ind.avail!=next.var]
      ind.keep <- sort(c(ind.keep, next.var))
      Xr.now <- X[,ind.keep,drop=FALSE]
      bic.prev <- bic.next      
    }
    print(ind.keep)
    
   ## possible deletion
    if(length(ind.keep)>2 && back==TRUE){
      bic.inc <- 0
      for(j in 1:length(ind.keep)){
        X.mj <- as.matrix(Xr.now[,-j])
        if(ncol(X.mj)>1)
          ans.j <- Mclust(X.mj,G=G,modelNames=modelNames)
        else
          ans.j <- Mclust(X.mj,G=G,modelNames=NULL)
        lm.j <- lm.fit(cbind(1,X.mj), Xr.now[,j])
        inc.j <- ans.j$bic + get.bic(lm.j) - bic.prev
        if(inc.j > bic.inc){
          rm.var <- j
          bic.inc <- inc.j
          bic.next <- ans.j$bic
	  flag2 <- 1
        }
      }
    }
    if(flag2){
      ind.avail <- sort(c(ind.avail,ind.keep[rm.var]))
      ind.keep <- ind.keep[-rm.var]
      Xr.now <- X[,ind.keep,drop=FALSE]
      bic.prev <- bic.next
    }
    print(ind.keep)
    iter <- iter + 1
  }
  return(ind.keep)
}





############################################
######### MCMC Updating Functions ##########
############################################


get.like <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now){

  n <- nrow(yc.now)
  p <- ncol(yc.now)
  phi.now <- match(phi.now, sort(unique(phi.now)))
  M <- max(phi.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  y1 <- yc.now[,c.vars,drop=FALSE]
  y2 <- yc.now[,nc.vars,drop=FALSE]
  p1 <- ncol(y1)
  p2 <- ncol(y2)
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  P12 <- Psi.now[c.vars,nc.vars,drop=FALSE]
  P22 <- Psi.now[nc.vars,nc.vars,drop=FALSE]
  ans.P11 <- ginv.gp(P11)
  P22.1 <- P22 - t(P12)%*%ans.P11$inv%*%P12
  ans.P22.1 <- ginv.gp(P22.1)
  y1bar <- apply(y1,2,mean)
  y1barM <- matrix(y1bar,n,p1,byrow=TRUE)
  y2bar <- apply(y2,2,mean)
  y2barM <- matrix(y2bar,n,p2,byrow=TRUE)
  V11 <- crossprod(y1-y1barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y1bar) + P11
  V22 <- crossprod(y2-y2barM) + n*lambda.now/(n+lambda.now)*y2bar%*%t(y2bar) + P22
  V12 <- crossprod(y1-y1barM,y2-y2barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y2bar) + P12
  ans.V11 <- ginv.gp(V11)
  V22.1 <- V22 - t(V12)%*%ans.V11$inv%*%V12
  ans.V22.1 <- ginv.gp(V22.1)
  ldet.P11 <- ans.P11$log.det
  ldet.V11 <- ans.V11$log.det
  ldet.P22.1 <- ans.P22.1$log.det
  ldet.V22.1 <- ans.V22.1$log.det
  
  like <- -n*p/2*log(2*pi) + p2/2*log(lambda.now/(n+lambda.now)) + p2/2*ldet.P11 - p2/2*ldet.V11 + eta.now/2*ldet.P22.1 - (n+eta.now)/2*ldet.V22.1 + lmgamma(p2, (n+eta.now)/2) - lmgamma(p2, eta.now/2)

  for(m in 1:M){
    ym1 <- y1[phi.now==m,,drop=FALSE]
    nm <- nrow(ym1)
    if(nm > 0){
      ym1bar <- apply(ym1,2,mean)
      ym1barM <- matrix(ym1bar,nm,p1,byrow=TRUE)
      Vm11 <- crossprod(ym1-ym1barM) + nm*lambda.now/(nm+lambda.now)*ym1bar%*%t(ym1bar) + P11
      ldet.Vm11 <- ginv.gp(Vm11)$log.det
      like <- like + p1/2*log(lambda.now/(nm+lambda.now)) + (eta.now-p2)/2*ldet.P11 - (nm+eta.now-p2)/2*ldet.Vm11 + lmgamma(p1,(nm+eta.now-p2)/2) - lmgamma(p1,(eta.now-p2)/2)
    }
  }
  return(like)
}







update.yc <- function(yc.now, y, phi.now, gamma.now, lambda.now, eta.now, Psi.now, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, discrete, limits, use.mean=FALSE, force=FALSE, obs=NULL){

##  Draw mu_m1, Sigma_m11, b_2, Q21, Q22 | y_{-i}, phi, gamma
##  ... then draw missing to complete y_i
##  Also block y_i


  n <- nrow(y)
  p <- ncol(y)
  if(is.null(obs))
    obs <- sample(1:n)
  M <- max(phi.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  na.mat <- is.na(y)
  na.obs <- which(apply(na.mat,1,sum)>0)
  n.na <- length(na.obs)
  nd <- sum(discrete)
  udis <- list()
  for(j in 1:p){
    if(discrete[j]){
      udis[[j]] <- sort(unique(y[,j]))
      udis[[j]] <- c(udis[[j]], Inf)
    }
    else
      udis[[j]] <- c()
  }
#  if(n.na==0 && nd==0 &&)
#    return(list(yc.now=y, accept=rep(1,0)))

  blocks <- list()
  ind.now <- 1
  nb <- ceiling(n/bsy)
  for(b in 1:nb){
    blocks[[b]] <- ind.now:min(ind.now+bsy-1, n)
    ind.now <- ind.now+bsy
  }
  accept.m <- accept.d <- accept.c <- rep(0,n)
  
  for(b in 1:nb){
    yc.p <- yc.now
    ind.b <- obs[blocks[[b]]]
    yc.mb <- yc.now[-ind.b,]
    ans.muQ <- rmuQ(yc.mb, phi.now[-ind.b], gamma.now, lambda.now, eta.now, Psi.now, M)
    mu.b <- ans.muQ$mu
    Q.b <- ans.muQ$Q
    
    d.p <- d.now <- 0
    for(i in ind.b){
      m <- phi.now[i]
      ind.na.i <- is.na(y[i,])
      ind.c.i <- !ind.na.i
      
     ## sample yc for missing y's    
      if(any(ind.na.i)){
        Q.na.i <- Q.b[m,ind.na.i,ind.na.i]
        ans.Q.na.i <- ginv.gp(Q.na.i*1/prop.sd.y.m^2)
        S.na.i.sqrt <- ans.Q.na.i$sqrt.inv
        if(use.mean){
          mu.na.i <- mu.b[m,ind.na.i] - ginv.gp(Q.na.i)$inv%*%Q.b[m,ind.na.i,ind.c.i]%*%(yc.now[i,ind.c.i] - mu.b[m,ind.c.i])
          yc.p[i,ind.na.i] <- rmvnorm(1,mu.na.i,S.sqrt=S.na.i.sqrt)
          d.p <- d.p + dmvnorm(yc.p[i,ind.na.i], mu.na.i, S.inv=1/prop.sd.y.m^2*Q.na.i,log.det=-ans.Q.na.i$log.det)
          d.now <- d.now + dmvnorm(yc.now[i,ind.na.i], mu.na.i, S.inv=1/prop.sd.y.m^2*Q.na.i, log.det=-ans.Q.na.i$log.det)
        }
        else{
          mu.na.i <- yc.now[i,ind.na.i] 
          yc.p[i,ind.na.i] <- rmvnorm(1,mu.na.i,S.sqrt=S.na.i.sqrt)
          d.p <- d.now <- 0
        }
      }

     ## sample latent yc for observed, but discrete y's
      ind.dobs.i <- which(!ind.na.i & discrete)
      ndobs.i <- length(ind.dobs.i)
      for(l in index(1,ndobs.i)){
        j <- ind.dobs.i[l]
        ind.mj <- (1:p)[-j]
#        mu.j <- mu.b[m,j] - 1/Q.b[m,j,j]*Q.b[m,j,ind.mj]%*%(yc.now[i,ind.mj] - mu.b[m,ind.mj])
        mu.j <- yc.now[i,j]
	indu <- which(udis[[j]]==y[i,j])
	aa <- ifelse(indu==1, -Inf, y[i,j])
	bb <- udis[[j]][indu+1]
        yc.p[i,j] <- rtruncnorm(1,a=aa, b=bb, mu.j, prop.sd.y.d*sqrt(1/Q.b[m,j,j]))
#        mu.jp <- mu.j
	mu.jp <- yc.p[i,j]
	d.p <- d.p + log(dtruncnorm(yc.p[i,j], a=aa, b=bb, mu.j, prop.sd.y.d*sqrt(1/Q.b[m,j,j])))
	d.now <- d.now + log(dtruncnorm(yc.now[i,j], a=aa, b=bb, mu.jp, prop.sd.y.d*sqrt(1/Q.b[m,j,j])))
      }

     ## sample latent yc for observed, but left or right censored y's
      ind.lcens <- y[i,]<=limits[,1]
      ind.rcens <- y[i,]>=limits[,2]
      ind.lrobs.i <- which(!discrete & (ind.lcens | ind.rcens))
      nlrobs.i <- length(ind.lrobs.i)
      for(l in index(1,nlrobs.i)){
        j <- ind.lrobs.i[l]
        ind.mj <- (1:p)[-j]
#        mu.j <- mu.b[m,j] - 1/Q.b[m,j,j]*Q.b[m,j,ind.mj]%*%(yc.now[i,ind.mj] - mu.b[m,ind.mj])
        mu.j <- yc.now[i,j]
	aa <- ifelse(ind.rcens[j], limits[j,2], -Inf)
	bb <- ifelse(ind.lcens[j], limits[j,1], Inf)
        yc.p[i,j] <- rtruncnorm(1,a=aa, b=bb, mu.j, prop.sd.y.c*sqrt(1/Q.b[m,j,j]))
#        mu.jp <- mu.j
	mu.jp <- yc.p[i,j]
	d.p <- d.p + log(dtruncnorm(yc.p[i,j], a=aa, b=bb, mu.j, prop.sd.y.c*sqrt(1/Q.b[m,j,j])))
	d.now <- d.now + log(dtruncnorm(yc.now[i,j], a=aa, b=bb, mu.jp, prop.sd.y.c*sqrt(1/Q.b[m,j,j])))
      }
    }
    if(force)
      MH.ratio <- 1
    else{
      like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
      like.p <- get.like(yc.p, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
      MH.ratio <- exp(like.p - d.p - like.now + d.now)
    }
    if(runif(1) < MH.ratio){
      yc.now <- yc.p
      accept.m[ind.b] <- 1
    }
  }
  return(list(yc.now=yc.now, accept.m=accept.m))
}





rmuQ <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, M, use.mean=FALSE){

  n <- nrow(yc.now)
  p <- ncol(yc.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  y1 <- yc.now[,c.vars,drop=FALSE]
  y2 <- yc.now[,nc.vars,drop=FALSE]
  p1 <- ncol(y1)
  p2 <- ncol(y2)
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  P12 <- Psi.now[c.vars,nc.vars,drop=FALSE]
  P22 <- Psi.now[nc.vars,nc.vars,drop=FALSE]
  ans.P11 <- ginv.gp(P11)
  P22.1 <- P22 - t(P12)%*%ans.P11$inv%*%P12
  ans.P22.1 <- ginv.gp(P22.1)
  y1bar <- apply(y1,2,mean)
  y1barM <- matrix(y1bar,n,p1,byrow=TRUE)
  y2bar <- apply(y2,2,mean)
  y2barM <- matrix(y2bar,n,p2,byrow=TRUE)
  V11 <- crossprod(y1-y1barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y1bar) + P11
  V22 <- crossprod(y2-y2barM) + n*lambda.now/(n+lambda.now)*y2bar%*%t(y2bar) + P22
  V12 <- crossprod(y1-y1barM,y2-y2barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y2bar) + P12
  ans.V11 <- ginv.gp(V11)
  V22.1 <- V22 - t(V12)%*%ans.V11$inv%*%V12
  ans.V22.1 <- ginv.gp(V22.1)
  ldet.P11 <- ans.P11$log.det
  ldet.V11 <- ans.V11$log.det
  ldet.P22.1 <- ans.P22.1$log.det
  ldet.V22.1 <- ans.V22.1$log.det

 ## Draw Q22
  if(use.mean)
    Q22 <- ans.V22.1$inv*(n+eta.now)
  else
    Q22 <- rwish(n+eta.now, ans.V22.1$inv)
  ans.Q22 <- ginv.gp(Q22)
  if(use.mean)
    Q21 <- -Q22%*%t(V12)%*%ans.V11$inv
  else
    Q21 <- rmatnorm(-Q22%*%t(V12)%*%ans.V11$inv, U.sqrt=ans.Q22$sqrt, V.sqrt=ans.V11$sqrt.inv)
  mu.b2.s <- n/(n+lambda.now)*(Q22%*%y2bar+Q21%*%y1bar)
  S.b2.sqrt <- 1/sqrt(n+lambda.now)*ans.Q22$sqrt
  if(use.mean)
    b2 <- mu.b2.s
  else
    b2 <- rmvnorm(1,mu.b2.s,S.sqrt=S.b2.sqrt)

  Q <- array(0,c(M,p,p))
  mu <- matrix(0,M,p)
  for(m in 1:M){
    ym1 <- y1[phi.now==m,,drop=FALSE]
    nm <- nrow(ym1)
    if(nm > 0){
      ym1bar <- apply(ym1,2,mean)
      ym1barM <- matrix(ym1bar,nm,p1,byrow=TRUE)
      P.s <- crossprod(ym1-ym1barM) + nm*lambda.now/(nm+lambda.now)*ym1bar%*%t(ym1bar) + P11
      mu.s <- nm/(nm+lambda.now)*ym1bar
    }
    else{
      mu.s <- 0
      P.s <- P11
    }
    eta.s <- nm + eta.now - p2
    lambda.s <- nm + lambda.now
    if(use.mean)
      ans.m <- list(mu=mu.s, Q=ginv.gp(P.s)$inv*eta.s)
    else
      ans.m <- rniw(mu.s,lambda.s,P.s,eta.s)
    mu.m1 <- ans.m$mu
    Q.m11 <- ans.m$Q + t(Q21)%*%ans.Q22$inv%*%Q21
    mu.m2 <- ans.Q22$inv%*%(b2-Q21%*%mu.m1)
    mu[m,c.vars] <- mu.m1
    mu[m,nc.vars] <- mu.m2
    Q[m,c.vars,c.vars] <- Q.m11
    Q[m,c.vars,nc.vars] <- t(Q21)
    Q[m,nc.vars,c.vars] <- Q21
    Q[m,nc.vars,nc.vars] <- Q22
  }
  return(list(mu=mu, Q=Q))
}





rmu1S11 <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, M, draw=FALSE){

  n <- nrow(yc.now)
  p <- ncol(yc.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  y1 <- yc.now[,c.vars,drop=FALSE]
  p1 <- ncol(y1)
  p2 <- p-p1
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  ans.P11 <- ginv.gp(P11)
  y1bar <- apply(y1,2,mean)
  y1barM <- matrix(y1bar,n,p1,byrow=TRUE)
  V11 <- crossprod(y1-y1barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y1bar) + P11
  ans.V11 <- ginv.gp(V11)
  ldet.P11 <- ans.P11$log.det
  ldet.V11 <- ans.V11$log.det

  S11 <- array(0,c(M,p1,p1))
  mu1 <- matrix(0,M,p1)
  for(m in 1:M){
    ym1 <- y1[phi.now==m,,drop=FALSE]
    nm <- nrow(ym1)
    if(nm > 0){
      ym1bar <- apply(ym1,2,mean)
      ym1barM <- matrix(ym1bar,nm,p1,byrow=TRUE)
      P.s <- crossprod(ym1-ym1barM) + nm*lambda.now/(nm+lambda.now)*ym1bar%*%t(ym1bar) + P11
      mu.s <- nm/(nm+lambda.now)*ym1bar
    }
    else{
      mu.s <- 0
      P.s <- P11
    }
    eta.s <- nm + eta.now - p2
    lambda.s <- nm + lambda.now
    if(draw){
      ans.m <- rniw(mu.s,lambda.s,P.s,eta.s)
      mu1[m,] <- ans.m$mu
      S11[m,,] <- ans.m$Sigma
    }
    else{
      mu1[m,] <- mu.s
      S11[m,,] <- P.s/(eta.s-p1-1)
    }
  }
  return(list(mu1=mu1, S11=S11))
}







update.gamma <- function(gamma.now, p.ad, p.swap, yc.now, phi.now, lambda.now, eta.now, Psi.now, rho, Lg, groups){

  accept <- 0
  p <- length(gamma.now)
  sum.as <- p.ad+p.swap
  p.ad <- p.ad/sum.as
  p.swap <- p.swap/sum.as
  for(l in 1:Lg){
    gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)
    d.g.pgnow <- dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    d.g.nowgp <- dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)

    like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
    like.p <- get.like(yc.now, phi.now, gamma.p, lambda.now, eta.now, Psi.now)
    pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
    pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
    MH.ratio <- exp(like.p + pi.g.p - d.g.pgnow - like.now - pi.g.now + d.g.nowgp)
    if(runif(1) < MH.ratio){
      gamma.now <- gamma.p
      accept <- 1
    }
  }
  return(list(gamma=gamma.now, accept=accept))
}





rgamma.prop <- function(gamma.now, p.ad, p.swap, groups){

  p <- length(gamma.now)
  gamma.p <- gamma.now
  probs <- c(p.ad, p.swap, 1-p.ad-p.swap)
  probs[probs<0] <- 0
  move <- sample(1:3, 1, prob=probs)
  ng.vec <- as.numeric(table(groups))
  G <- length(ng.vec)
  g <- sample(1:G, 1, prob=ng.vec)
  ind.g <- which(groups==g)
  ncv <- sum(gamma.now==1)
  ind0.g <- which(gamma.now==0 & groups==g)
  ind1.g <- which(gamma.now==1 & groups==g)
  if(move==1){
    if(length(ind1.g)==1)
      ind <- sample(ind0.g, 1)
    else if(ncv==(p-1))
      ind <- sample(ind1.g, 1)
    else
      ind <- sample(ind.g, 1)
    gamma.p[ind] <- 1-gamma.now[ind]
  }
  if(move==2){
    ind120 <- ifelse(length(ind1.g)==1, ind1.g, sample(ind1.g, 1))
    ind021 <- ifelse(length(ind0.g)==1, ind0.g, sample(ind0.g, 1))
    gamma.p[ind120] <- 0
    gamma.p[ind021] <- 1
  }
  return(gamma.p)
}




dgamma.prop <- function(gamma.p, gamma.now, p.ad, p.swap, groups){

  ndiff <- sum(abs(gamma.p-gamma.now))
  if(ndiff==0)
    ans <- log(1-p.ad-p.swap)
  else{
    ng.vec <- as.numeric(table(groups))
    G <- length(ng.vec)
    g <- groups[which(gamma.now!=gamma.p)][1]
    p.g <- ng.vec[g]/sum(ng.vec)
    ind.g <- which(groups==g)
    ncv <- sum(gamma.now==1)
    ind0.g <- which(gamma.now==0 & groups==g)
    ind1.g <- which(gamma.now==1 & groups==g)
    if(ndiff==2)
      ans <- log(p.swap*p.g/length(ind1.g)/length(ind0.g))
    if(ndiff==1){
      if(length(ind1.g)==1 || length(ind0.g)==1)
        ans <- p.ad*p.g/(ng.vec[g]-1)
      else
        ans <- p.ad*p.g/ng.vec[g]
    }
  }
  return(ans)
}






update.phi.sm <- function(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, TT){

  accept <- 0
  n <- nrow(yc.now)
  M <- max(phi.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  i <- sample(1:n,1)
  j <- sample((1:n)[-i],1)
  y1 <- yc.now[,c.vars,drop=FALSE]

  split <- phi.now[i]==phi.now[j]
  obs <- which(phi.now==phi.now[i] | phi.now==phi.now[j])
  obs.s <- obs[obs!=i & obs!=j]

#i<<-i
#j<<-j

  phi.launch <- get.phi.launch(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs.s, i, j, TT)

#phi.launch<<-phi.launch

  if(split){
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, m.ind=c(phi.now[i],M+1))
    phi.p <- ans.s$phi
    d.p <- ans.s$dens
    d.now <- 0
    like.p <- get.like.sm(yc.now, phi.p, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=M+1)
    like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[i])
    n1 <- sum(phi.p==phi.now[i])
    n2 <- sum(phi.p==M+1)
    pi.pdivnow <- log(delta.now) + lfactorial(n1-1) + lfactorial(n2-1) - lfactorial(n1+n2-1)
  }
  else{
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs=obs, m.ind=c(phi.now[i],phi.now[j]), phi.p=phi.now)
    phi.p <- phi.now
    phi.p[obs] <- phi.now[i]
    u.phi.p <- sort(unique(phi.p))
    phi.p <- match(phi.p, u.phi.p)
    phi.p.i <- which(u.phi.p==phi.now[i])
    d.p <- 0
    d.now <- ans.s$dens
    like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[j])
    like.p <- get.like.sm(yc.now, phi.p, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.p.i, phi.j=phi.p.i)
    n1 <- sum(phi.now==phi.now[i])
    n2 <- sum(phi.now==phi.now[j])
    pi.pdivnow <- -log(delta.now) - lfactorial(n1-1) - lfactorial(n2-1) + lfactorial(n1+n2-1)
  }
  MH.ratio <- exp(like.p - like.now + pi.pdivnow + d.now - d.p)
  if(runif(1) < MH.ratio){
    phi.now <- match(phi.p, sort(unique(phi.p)))
    accept <- 1
  }

  return(list(phi.now=phi.now, accept=accept))
}





get.phi.launch <- function(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs.s, i, j, TT){

  M <- max(phi.now)
  split <- phi.now[i]==phi.now[j]
  phi.i <- phi.now[i]
  phi.j <- ifelse(split, M+1, phi.now[j])
  phi.s0 <- phi.now
  phi.s0[j] <- phi.j 
  if(length(obs.s)==0)
    return(phi.s0)
    
  c.vars <- (gamma.now==1)
  y1.s <- yc.now[obs.s,c.vars,drop=FALSE]
  n.s <- nrow(y1.s)
  p1 <- ncol(y1.s)
  y1.i <- yc.now[i,c.vars]
  y1.j <- yc.now[j,c.vars]
  dist.i <- apply((y1.s - matrix(y1.i,n.s,p1, byrow=TRUE))^2,1,sum)
  dist.j <- apply((y1.s - matrix(y1.j,n.s,p1, byrow=TRUE))^2,1,sum)
  kmc <- apply(cbind(dist.i,dist.j),1,order)[1,]
  phi.s0[obs.s] <- phi.i*(kmc==1) + phi.j*(kmc==2) 

  for(t in index(1,TT)){
    ans.s <- update.phi.2G(phi.s0, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, m.ind=c(phi.i,phi.j), dens=FALSE)
    phi.s0 <- ans.s$phi
  }
  return(phi.s0)
}




get.like.sm <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i, phi.j){

  n <- nrow(yc.now)
  p <- ncol(yc.now)
  M <- max(phi.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  ans.P11 <- ginv.gp(P11)
  ldet.P11 <- ans.P11$log.det  
  y1 <- yc.now[,c.vars,drop=FALSE]
  p1 <- ncol(y1)
  p2 <- p-p1
  m.ind <- unique(c(phi.i,phi.j))

  like <- 0
  for(m in m.ind){
    ym1 <- y1[phi.now==m,,drop=FALSE]
    nm <- nrow(ym1)
    if(nm > 0){
      ym1bar <- apply(ym1,2,mean)
      ym1barM <- matrix(ym1bar,nm,p1,byrow=TRUE)
      Vm11 <- crossprod(ym1-ym1barM) + nm*lambda.now/(nm+lambda.now)*ym1bar%*%t(ym1bar) + P11
      ldet.Vm11 <- ginv.gp(Vm11)$log.det
      like <- like + p1/2*log(lambda.now/(nm+lambda.now)) + (eta.now-p2)/2*ldet.P11 - (nm+eta.now-p2)/2*ldet.Vm11 + lmgamma(p1,(nm+eta.now-p2)/2) - lmgamma(p1,(eta.now-p2)/2)
    }
  }
  return(like)
}




update.phi <- function(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now){

  n <- nrow(yc.now)
  obs <- sample(1:n,n,replace=FALSE)
  p <- ncol(yc.now)
  c.vars <- (gamma.now==1)
  y1 <- yc.now[,c.vars,drop=FALSE]
  dens <- 0
  for(i in obs){
    phi.prev.i <- phi.now[i]
    probs.i <- get.phi.probs(i, phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, m.ind=NULL)
    phi.now[i] <- sample(1:length(probs.i),1,prob=probs.i)
    if(!any(phi.now==phi.prev.i))
      phi.now <- match(phi.now, sort(unique(phi.now)))
  }
  return(phi.now)
}




update.phi.2G <- function(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs, m.ind, phi.p=NULL, dens=TRUE){

  draw <- is.null(phi.p)
  n <- nrow(yc.now)
  p <- ncol(yc.now)
  c.vars <- (gamma.now==1)
  y1 <- yc.now[,c.vars,drop=FALSE]
  dens <- 0
  phi.now2 <- phi.now
  phi.now2[phi.now==m.ind[1]] <- m.ind[2]
  phi.now2[phi.now==m.ind[2]] <- m.ind[1]  
  for(i in obs){
    probs.i <- get.phi.probs(i, phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, m.ind)
    if(draw)
      phi.now[i] <- sample(1:length(probs.i),1,prob=probs.i)
    else
      phi.now[i] <- phi.p[i]
    dens <- dens + log(probs.i[phi.now[i]])
  }
  phi.p <- phi.now
  return(list(phi.now=phi.now, dens=dens))
}



get.phi.probs <- function(i, phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, m.ind=NULL){

  n <- nrow(yc.now)
  M <- max(phi.now[-i])
  if(is.null(m.ind))
    m.ind <- 1:(M+1)
  c.vars <- (gamma.now==1)
  p <- ncol(yc.now)
  p1 <- sum(c.vars)
  p2 <- p - p1
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  y1.mi <- yc.now[-i,c.vars,drop=FALSE]
  y1.i <- yc.now[i,c.vars]
  phi.mi <- phi.now[-i]
  like <- rep(-Inf,M+1)
  prob <- rep(0,M+1)
  for(m in m.ind){
    ym1.mi <- y1.mi[phi.mi==m,,drop=FALSE]
    ym1 <- rbind(ym1.mi,y1.i)
    nm <- nrow(ym1)
    if(nm > 1){
      ym1bar <- apply(ym1,2,mean)
      ym1barM <- matrix(ym1bar,nm,p1,byrow=TRUE)
      ym1bar.mi <- apply(ym1.mi,2,mean)
      ym1barM.mi <- matrix(ym1bar.mi,nm-1,p1,byrow=TRUE)
      Vm11 <- crossprod(ym1-ym1barM) + nm*lambda.now/(nm+lambda.now)*ym1bar%*%t(ym1bar) + P11
      Vm11.mi <- crossprod(ym1.mi-ym1barM.mi) + (nm-1)*lambda.now/(nm-1+lambda.now)*ym1bar.mi%*%t(ym1bar.mi) + P11
      ldet.Vm11 <- ginv.gp(Vm11)$log.det
      ldet.Vm11.mi <- ginv.gp(Vm11.mi)$log.det
      like[m] <-  log(nm-1) + p1/2*log((nm-1+lambda.now)/(nm+lambda.now)) + (nm-1+eta.now-p2)/2*ldet.Vm11.mi - (nm+eta.now-p2)/2*ldet.Vm11 + lmgamma(p1,(nm+eta.now-p2)/2) - lmgamma(p1,(nm-1+eta.now-p2)/2)
    }
    else{
      Vm11 <- lambda.now/(1+lambda.now)*t(ym1)%*%ym1 + P11
      Vm11.mi <- P11
      ldet.Vm11 <- ginv.gp(Vm11)$log.det
      ldet.Vm11.mi <- ginv.gp(Vm11.mi)$log.det

p1..<<-p1
p2..<<-p2
eta.now..<<-eta.now


      like[m] <-  log(delta.now) + p1/2*log(lambda.now/(1+lambda.now)) + (eta.now-p2)/2*ldet.Vm11.mi - (1+eta.now-p2)/2*ldet.Vm11 + lmgamma(p1,(1+eta.now-p2)/2) - lmgamma(p1,(eta.now-p2)/2)
    }
  }

#phi.now..<<-phi.now
#m.ind<<-m.ind
#like<<-like

  for(m in m.ind){
    prob[m] <- 1/(sum(exp(like - like[m])))
  }
  return(prob)
}







update.phi.AG <- function(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=NULL, mu1=NULL, Q11=NULL, ldet.vec=NULL, obs){

  n <- nrow(yc.now)
  if(is.null(obs))
    obs <- 1:n
  p <- ncol(yc.now)
  c.vars <- (gamma.now==1)
  y1 <- yc.now[,c.vars,drop=FALSE]
  p1 <- ncol(y1)
  M <- max(phi.now)

  if(is.null(Q11)){
    ans.mu1S11 <- rmu1S11(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, M)
    mu1 <- ans.mu1S11$mu1
    S11 <- ans.mu1S11$S11
    ldet.vec <- rep(0,M)
    Q11 <- S11
    for(m in 1:M){
      ans.m <- ginv.gp(S11[m,,])
      Q11[m,,] <- ans.m$inv
      ldet.vec[m] <- ans.m$log.det
    }
  }

  omega <- update.vee.hat(phi.now, delta.now, M=max(phi.now)+1)$omega[1:M]
  dens <- 0
  draw <- is.null(phi.p)
  prob <- rep(0,M)

  for(i in obs){
    like.vec <- rep(0,M)
    if(sum(phi.now==phi.now[i])==1){
      prob <- rep(0,M)
      prob[phi.now[i]] <- 1
    }
    else{
      for(m in 1:M)
        like.vec[m] <- log(omega[m]) + dmvnorm(y1[i,],mu1[m,],S.inv=Q11[m,,],log.det=ldet.vec[m])
      for(m in 1:M)
        prob[m] <- 1/(sum(exp(like.vec - like.vec[m])))
    }
    if(draw)
      phi.now[i] <- sample(1:M,1,prob=prob)
    else
      phi.now[i] <- phi.p[i]
    dens <- dens + log(prob[phi.now[i]])

#print(i)
#print(prob[phi.now[i]])

  }    
  return(list(phi.now=phi.now, dens=dens, mu1=mu1, Q11=Q11, ldet.vec=ldet.vec))
}




update.vee.hat <- function(phi.now, delta.now, M){

  if(M==1){
    vee <- 0
    omega <- 1
  }
  else{
    vee <- rep(0,M)
    for(m in 1:M){
      a.star <- sum(phi.now==m)+1
      b.star <- sum(phi.now>m)+delta.now
      vee[m] <- b.star/(a.star+b.star)
    }
    vee[M] <- 0
    omega <- makeprobs(vee)
  }
  return(list(vee=vee,omega=omega))
}




update.yg <- function(yc.now, y, phi.p, gamma.p, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2, mw=1000){

  n <- nrow(y)
  M <- max(phi.p)
  c.vars <- (gamma.p==1)
  nc.vars <- !c.vars
  gamma.all <- (gamma.p+gamma.now)>0
  ycv <- yc.now
  ycv[,gamma.all] <- y[,gamma.all]
  na.mat <- is.na(ycv)
  na.obs <- which(apply(na.mat,1,sum)>0)
  n.na <- length(na.obs)
  nd <- sum(discrete)
  udis <- list()
  for(j in 1:p){
    if(discrete[j]){
      udis[[j]] <- sort(unique(y[,j]))
      udis[[j]] <- c(udis[[j]], Inf)
    }
    else
      udis[[j]] <- c()
  }
#  if(n.na==0 && nd==0)
#    return(list(yc=yc.now, dens=0))
  
  ans.muQ <- rmuQ(yc.now, phi.p, gamma.all, lambda.now, eta.now, Psi.now, M)  #, use.mean=TRUE)
  mu <- ans.muQ$mu
  Q <- ans.muQ$Q
  dens <- 0
  draw <- is.null(yc.p)
  yc.p2 <- yc.now

if(0){
  for(i in na.obs){
    m <- phi.p[i]
    ind.na.i <- is.na(ycv[i,])
    ind.c.i <- !ind.na.i
    Q.na.i <- Q[m,ind.na.i,ind.na.i]
    ans.Q.na.i <- ginv.gp(Q.na.i*1/prop.sd.y1^2)
    S.na.i.sqrt <- ans.Q.na.i$sqrt.inv
    mu.na.i <- mu[m,ind.na.i] - ginv.gp(Q.na.i)$inv%*%Q[m,ind.na.i,ind.c.i]%*%(yc.now[i,ind.c.i] - mu[m,ind.c.i])
    mu.na.i <- (mw*mu.na.i+yc.now[i,ind.na.i])/(1+mw)
    if(draw)
      yc.p2[i,ind.na.i] <- rmvnorm(1,mu.na.i,S.sqrt=S.na.i.sqrt)
    else
      yc.p2[i,ind.na.i] <- yc.p[i,ind.na.i]
    dens <- dens + dmvnorm(yc.p2[i,ind.na.i], mu.na.i, S.inv=1/prop.sd.y1^2*Q.na.i, log.det=-ans.Q.na.i$log.det)
  }
}

  dis.cv <- discrete
  dis.cv[!gamma.all] <- FALSE
  for(i in 1:n){
    m <- phi.p[i]
    ind.na.i <- is.na(ycv[i,])
    ind.dobs.i <- which(!ind.na.i & dis.cv)
    ndobs.i <- length(ind.dobs.i)
   ## sample latent yc for observed, but discrete y's
    for(l in index(1,ndobs.i)){
      j <- ind.dobs.i[l]
      ind.mj <- (1:p)[-j]
      mu.j <- mu[m,j] - 1/Q[m,j,j]*Q[m,j,ind.mj]%*%(yc.p2[i,ind.mj] - mu[m,ind.mj])
      mu.j <- (mw*mu.j+yc.now[i,j])/(1+mw)
      indu <- which(udis[[j]]==y[i,j])
      aa <- ifelse(indu==1, -Inf, y[i,j])
      bb <- udis[[j]][indu+1]
      if(draw)
        yc.p2[i,j] <- rtruncnorm(1,a=aa, b=bb, mu.j, sqrt(prop.sd.y2/Q[m,j,j]))
      else
        yc.p2[i,j] <- yc.p[i,j]
      dens <- dens + log(dtruncnorm(yc.p2[i,j], a=aa, b=bb, mu.j, sqrt(prop.sd.y2/Q[m,j,j])))
    }
   ## sample latent yc for observed, but left or right censored y's
    ind.lcens <- y[i,]<=limits[,1]
    ind.rcens <- y[i,]>=limits[,2]
    ind.lrobs.i <- which(gamma.all & (ind.lcens | ind.rcens))
    nlrobs.i <- length(ind.lrobs.i)
    for(l in index(1,nlrobs.i)){
      j <- ind.lrobs.i[l]
      ind.mj <- (1:p)[-j]
      mu.j <- mu[m,j] - 1/Q[m,j,j]*Q[m,j,ind.mj]%*%(yc.p2[i,ind.mj] - mu[m,ind.mj])
      mu.j <- (mw*mu.j+yc.now[i,j])/(1+mw)
      aa <- ifelse(ind.rcens[j], limits[j,2], -Inf)
      bb <- ifelse(ind.lcens[j], limits[j,1], Inf)
      if(draw)
        yc.p2[i,j] <- rtruncnorm(1,a=aa, b=bb, mu.j, sqrt(prop.sd.y2/Q[m,j,j]))
      else
        yc.p2[i,j] <- yc.p[i,j]
      dens <- dens + log(dtruncnorm(yc.p2[i,j], a=aa, b=bb, mu.j, sqrt(prop.sd.y2/Q[m,j,j])))
    }
  }

  return(list(yc.now=yc.p2, dens=dens))
}





update.gamma.y <- function(gamma.now, p.ad, p.swap, yc.now, y, phi.now, lambda.now, eta.now, Psi.now, rho, Lgy, discrete, limits, prop.sd.y1, prop.sd.y2, groups){

  accept <- 0
  p <- length(gamma.now)
  sum.ads <- p.ad+p.swap
  p.stay <- 1-(sum.ads)
  sum.new <- (1-p.stay)^(1/Lgy)
  p.ad <- p.ad/sum.ads*sum.new
  p.swap <- p.swap/sum.ads*sum.new
  for(l in index(1,Lgy)){
    gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)
    d.now <- dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    d.p <- dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)

    ans.y.p <- update.yg(yc.now, y, phi.now, gamma.p, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2)
    yc.p <- ans.y.p$yc
    d.p <- d.p + ans.y.p$dens
    ans.y.now <- update.yg(yc.p, y, phi.now, gamma.now, gamma.p, lambda.now, eta.now, Psi.now, yc.p=yc.now, discrete, limits, prop.sd.y1, prop.sd.y2)
    d.now <- d.now + ans.y.now$dens
    
    like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
    like.p <- get.like(yc.p, phi.now, gamma.p, lambda.now, eta.now, Psi.now)
    pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
    pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
    MH.ratio <- exp(like.p + pi.g.p - d.p - like.now - pi.g.now + d.now)
    if(runif(1) < MH.ratio){
      gamma.now <- gamma.p
      yc.now <- yc.p
      accept <- 1
    }
  }
  return(list(gamma=gamma.now, yc.now=yc.now, accept=accept))
}





update.phi.y <- function(gamma.now, phi.now, p.ad, p.swap, yc.now, y, lambda.now, eta.now, Psi.now, rho, Lpy, discrete, limits, prop.sd.y1, prop.sd.y2){

  accept <- 0
  n <- nrow(yc.now)
  obs <- sample(1:n)
  for(l in index(1,Lpy)){
    ans.p <- update.phi.AG(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=NULL, obs=obs)
    phi.p <- ans.p$phi
    d.p <- ans.p$dens
    ans.now <- update.phi.AG(phi.p, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=phi.now, obs=obs[n:1])
    d.now <- ans.now$dens

    ans.y.p <- update.yg(yc.now, y, phi.p, gamma.now, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2)
    yc.p <- ans.y.p$yc
    d.p <- d.p + ans.y.p$dens
    ans.y.now <- update.yg(yc.p, y, phi.now, gamma.now, gamma.now, lambda.now, eta.now, Psi.now, yc.p=yc.now, discrete, limits, prop.sd.y1, prop.sd.y2)
    d.now <- d.now + ans.y.now$dens
    
    like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
    like.p <- get.like(yc.p, phi.p, gamma.now, lambda.now, eta.now, Psi.now)

    n.p <- as.numeric(table(phi.p))
    n.now <- as.numeric(table(phi.now))
    pi.pdivnow <- (length(n.p)-length(n.now))*log(delta.now) + sum(lfactorial(n.p-1)) - sum(lfactorial(n.now-1))
    
    MH.ratio <- exp(like.p - like.now + d.now - d.p + pi.pdivnow)

    if(runif(1) < MH.ratio){
      phi.now <- phi.p
      yc.now <- yc.p
      accept <- 1
    }
  }
  return(list(phi.now=phi.now, yc.now=yc.now, accept=accept))
}







update.gamma.phi.y <- function(gamma.now, phi.now, p.ad, p.swap, yc.now, y, lambda.now, eta.now, Psi.now, rho, delta.now, Lgp, discrete, limits, prop.sd.y1, prop.sd.y2, groups){

  accept <- 0
  n <- nrow(yc.now)
  sum.ads <- p.ad+p.swap
  p.stay <- 1-(sum.ads)
  sum.new <- (1-p.stay)^(1/Lgp)
  p.ad <- p.ad/sum.ads*sum.new
  p.swap <- p.swap/sum.ads*sum.new
  obs <- sample(1:n)
  
  for(l in index(1,Lgp)){
    gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)  

    ans.p <- update.phi.AG(phi.now, yc.now, gamma.p, lambda.now, eta.now, Psi.now, delta.now, phi.p=NULL, obs=obs)
    phi.p <- ans.p$phi
    d.p <- ans.p$dens + dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    ans.now <- update.phi.AG(phi.p, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=phi.now, obs=obs[n:1])
    d.now <- ans.now$dens + dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)

    ans.y.p <- update.yg(yc.now, y, phi.p, gamma.p, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2)
    yc.p <- ans.y.p$yc
    d.p <- d.p + ans.y.p$dens
    ans.y.now <- update.yg(yc.p, y, phi.now, gamma.now, gamma.p, lambda.now, eta.now, Psi.now, yc.p=yc.now, discrete, limits, prop.sd.y1, prop.sd.y2)
    d.now <- d.now + ans.y.now$dens
    
    like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
    like.p <- get.like(yc.p, phi.p, gamma.p, lambda.now, eta.now, Psi.now)

    pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
    pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
    n.p <- as.numeric(table(phi.p))
    n.now <- as.numeric(table(phi.now))
    pi.pdivnow <- pi.g.p - pi.g.now + (length(n.p)-length(n.now))*log(delta.now) + sum(lfactorial(n.p-1)) - sum(lfactorial(n.now-1))
    
    MH.ratio <- exp(like.p - like.now + d.now - d.p + pi.pdivnow)

    if(runif(1) < MH.ratio){
      gamma.now <- gamma.p
      phi.now <- phi.p
      yc.now <- yc.p
      accept <- 1
    }
  }
  return(list(gamma.now=gamma.now, phi.now=phi.now, yc.now=yc.now, accept=accept))
}






update.delta <- function(delta.now, phi.now, A.delta, B.delta, prop.sd.delta, dfp){

  accept <- 0
  M <- max(phi.now)
  n <- length(phi.now)
  delta.p <- r.pos.proposal(delta.now, df=dfp, sigma=prop.sd.delta)
  like.p <- M*log(delta.p) + lgamma(delta.p) - lgamma(delta.p+n)
  like.now <- M*log(delta.now) + lgamma(delta.now) - lgamma(delta.now+n)
  pi.p <- dgamma(delta.p, A.delta, B.delta, log=TRUE)
  pi.now <- dgamma(delta.now, A.delta, B.delta, log=TRUE)
  dprop.p <- d.pos.proposal(delta.p, delta.now, df=dfp, sigma=prop.sd.delta)
  dprop.now <- d.pos.proposal(delta.now, delta.p, df=dfp, sigma=prop.sd.delta)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    delta.now <- delta.p
    accept <- 1
  }
  return(list(delta.now=delta.now, accept=accept))
}




update.lambda <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, A.lambda, B.lambda, prop.sd.lambda, dfp){
 
  accept <- 0
  M <- max(phi.now)
  n <- length(phi.now)
  lambda.p <- r.pos.proposal(lambda.now, df=dfp, sigma=prop.sd.lambda)
  like.p <- get.like(yc.now, phi.now, gamma.now, lambda.p, eta.now, Psi.now)
  like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
  pi.p <- dgamma(lambda.p, A.lambda, B.lambda, log=TRUE)
  pi.now <- dgamma(lambda.now, A.lambda, B.lambda, log=TRUE)
  dprop.p <- d.pos.proposal(lambda.p, lambda.now, df=dfp, sigma=prop.sd.lambda)
  dprop.now <- d.pos.proposal(lambda.now, lambda.p, df=dfp, sigma=prop.sd.lambda)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    lambda.now <- lambda.p
    accept <- 1
  }
  return(list(lambda.now=lambda.now, accept=accept))
}




update.eta <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, A.eta, B.eta, prop.sd.eta, dfp, etagp){
 
  accept <- 0
  p <- ncol(yc.now)
  M <- max(phi.now)
  n <- length(phi.now)
  eta.p <- eta.p <- r.pos.proposal(eta.now-p-etagp, df=dfp, sigma=prop.sd.eta) + p + etagp
  like.p <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.p, Psi.now)
  like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
  pi.p <- dgamma(eta.p-p-etagp, A.eta, B.eta, log=TRUE)
  pi.now <- dgamma(eta.now-p-etagp, A.eta, B.eta, log=TRUE)
  dprop.p <- d.pos.proposal(eta.p-p-etagp, eta.now-p-etagp, df=dfp, sigma=prop.sd.eta)
  dprop.now <- d.pos.proposal(eta.now-p-etagp, eta.p-p-etagp, df=dfp, sigma=prop.sd.eta)
  MH.ratio <- exp((like.p + pi.p - dprop.p) - (like.now + pi.now - dprop.now))
  if(runif(1) < MH.ratio){
    eta.now <- eta.p
    accept <- 1
  }
  return(list(eta.now=eta.now, accept=accept))
}




update.Psi <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, P.Psi, N.Psi){
 
  accept <- 0
  M <- max(phi.now)
  n <- length(phi.now)
  
  Psi.now <- P.Psi*N.Psi
  accept <- 1
  return(list(Psi.now=Psi.now, accept=accept))
}




initialize.clusters <- function(y, G, lambda.now, eta.now, Psi.now, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, gamma.init=NULL, phi.init=NULL, discrete, limits, yc.now=NULL, Ly=200, groups){

  n <- nrow(y)
  p <- ncol(y)
  if(is.null(phi.init))
    phi.now <- rep(1,n)
  else
    phi.now <- phi.init
  if(is.null(gamma.init))
    gamma.now <- c(1,rep(0,p-1))
  else
    gamma.now <- gamma.init
  if(is.null(yc.now)){
    accept.m <- accept.d <- accept.c <- rep(0,n)
    yc.now <- y
    yc.now[is.na(y)] <- rnorm(sum(is.na(y)),0,.5)
    yc.now <- update.yc(yc.now, y, phi.now, gamma.now, lambda.now, eta.now, Psi.now, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, discrete, limits,force=TRUE)$yc.now
    for(i in 1:Ly){
      if(i%%1==0)
        print(i)
      ans.i <- update.yc(yc.now, y, phi.now, gamma.now, lambda.now, eta.now, Psi.now, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, discrete, limits)
      accept.m <- accept.m + ans.i$accept.m
      print(summary(accept.m/i))
      yc.now <- ans.i$yc.now
    }
  }
#yc.now..<<-yc.now
#save(yc.now..,file="yc.RData")

  if(is.null(gamma.init)){
#    ind.cv <- clustvarsel(yc.now, G=G, emModels2="VVV", parallel=getDoParWorkers())$subset
    ind.cv <- clustvarsel(yc.now, G=G, parallel=FALSE)$subset
    if(is.null(ind.cv))
      ind.cv <- 1
#    ind.cv <- try(sw.mclust(yc.now, G=G, modelNames="VVV", back=FALSE), silent=TRUE)
#    if(is.character(ind.cv[1]))
#      ind.cv <- sw.mclust(yc.now, G=G, modelNames="VVI", back=FALSE)

    gamma.now <- rep(1,p)
    gamma.now[-ind.cv] <- 0
    if(all(gamma.now==1))
      gamma.now[p] <- 0
    G <- length(unique(groups))
    for(g in 1:G){
      if(sum(gamma.now[groups==g])==0)
        gamma.now[which(groups==g)[1]] <- 1
    }
  }
  else
    ind.cv <- which(gamma.now==1)

gamma.now..<<-gamma.now
print(gamma.now)

  if(is.null(phi.init)){
    if(length(ind.cv)>1)
      ans.mclust <- Mclust(yc.now[,ind.cv,drop=FALSE],G=G)
    else
      ans.mclust <- Mclust(yc.now[,ind.cv,drop=FALSE],G=G)    

    phi.now <- ans.mclust$class 
  }

  return(list(yc.now=yc.now, phi.now=phi.now, gamma.now=gamma.now))
}






get.phi.probs.hat.old <- function(ans.DPMvs, post.ind=NULL, maxit=100, post.ind2=NULL, Mmax=8){

  n <- dim(ans.DPMvs$phi)[2]
  N.mcmc <- dim(ans.DPMvs$phi)[1]
  if(is.null(post.ind))
    post.ind <- 1:N.mcmc
  foo <- apply(ans.DPMvs$phi[post.ind,],1,max)
  ind.keep <- which(foo<=Mmax)
  post.ind <- post.ind[ind.keep]
  if(is.null(post.ind2))
    post.ind2 <- post.ind
  phi.now <- ans.DPMvs$phi[post.ind,]
  Mmax <- max(phi.now)
  probs.now <- matrix(0,n,Mmax)
  for(m in 1:Mmax)
    probs.now[,m] <- apply(ans.DPMvs$phi[post.ind2,]==m, 2, mean)
  for(it in 1:maxit){

print(it)

    phi.now <- foreach(k=1:length(post.ind), .combine=rbind)%dopar%{
      phi.k <- phi.k.new <- phi.now[k,]
      M.k <- max(phi.k)
      perms.k <- permutations(n=Mmax, r=M.k)
      np <- nrow(perms.k)
      like.best <- -Inf
      for(l in 1:np){
        phi.kl <- perms.k[l,][phi.k]
	like.l <- sum(log(probs.now[cbind(1:n,phi.kl)]))
        if(like.l>like.best)
	  phi.k.new <- phi.kl
      }
      phi.k.new
    }
    probs.new <- matrix(0,n,Mmax)
    for(m in 1:Mmax)
      probs.new[,m] <- apply(phi.now==m, 2, mean)
    dist <- sum((probs.new-probs.now)^2)
    probs.now <- probs.new
    if(dist<1E-6)
      break
  }
  if(it==maxit)
    warning("reached maxit")
  return(list(phi.now=phi.now, probs.now=probs.now))
}





get.phi.hat.probs.old <- function(ans.DPMvs, phi.hat, post.ind=NULL, N.post=100){

 ## for each i, sample a gamma, ymiss, and hyperparams, with phi[-i] fixed at phi.hat
 ## then sample phi[i]; record Gibbs prob for phi[i].  Report "posterior" average of Gibbs prob.
  N.mcmc <- dim(ans.DPMvs$phi)[1]
  if(is.null(post.ind))
    post.ind <- sample(1:N.mcmc, N.post, replace=TRUE)
  N.post <- length(post.ind)

  M <- max(phi.hat)
  y <- yc.now <- ans.DPMvs$y
  ind.miss <- is.na(y)
  n <- dim(y)[1]
  
  phi.probs.list <- foreach(kk=1:N.post)%dopar%{
    if(kk%%10==0)
      print(kk)
    it <- post.ind[kk]
    gamma <- ans.DPMvs$gamma[it,]
    delta <- ans.DPMvs$delta[it]
    lambda <- ans.DPMvs$lambda[it]
    eta <- ans.DPMvs$eta[it]
    Psi <- ans.DPMvs$Psi[it,,]
    yc.now[ind.miss] <- ans.DPMvs$ymiss[it,]
    phi.probs.it <- matrix(0,n,M+1)
    
    for(i in 1:n){
      phi.probs.it[i,] <- get.phi.probs(i, phi.hat, yc.now, gamma, lambda, eta, Psi, delta, m.ind=1:(M+1))
    }
    phi.probs.it
  }
  phi.probs.ar <- array(0,c(N.post,n,M+1))
  for(kk in 1:N.post)
    phi.probs.ar[kk,,] <- phi.probs.list[[kk]]
    
  phi.probs <- apply(phi.probs.ar, c(2,3), mean)
  return(phi.probs)
}












#############################################################
################ DPM-ind (Kim06 approach) ###################
#############################################################




DPMind.MCMC <- function(y, rho=.5, p.ad=.3, p.swap=.2, lambda=.1, lambda0=.1, Psi=NULL, eta=NULL, aa=1.5, bb=0.5, delta=.2, N.mcmc=5000, every=1, nplot=10, nback=1000, dfp=20, maxplot=100, gamma.init=NULL, phi.init=NULL, G=2:5, begin=0, bsy=1, prop.sd.y.m=.4, prop.sd.y.d=.4, prop.sd.y.c=.4, TT=2, Lg=200, discrete=NULL, limits=NULL, prop.sd.y1=NULL, prop.sd.y2=NULL, yc.init=NULL, Ly.init=100, groups=NULL){

 ## Create needed variables ##
  n <- nrow(y)
  p <- ncol(y)
  if(is.null(prop.sd.y1))
    prop.sd.y1 <- prop.sd.y.m
  if(is.null(prop.sd.y2))
    prop.sd.y2 <- prop.sd.y.d
  if(is.null(discrete))
    discrete <- rep(FALSE,p)
  if(is.null(limits))
    limits <- cbind(rep(-Inf,p), rep(Inf,p))
  if(is.null(colnames(y)))
    colnames(y) <- paste("y_",1:p,sep="")
  if(length(rho)==1)
    rho <- rep(rho, p)
  if(is.null(groups))
    groups <- rep(1,p)

  if(is.null(Psi))
    Psi <- diag(p)
  if(is.null(eta))
    eta <- p+2      ##(i.e., delta=3 in Kim06)

 ## initialize yc, phi, gamma
  cat("\nInitializing Clusters \n")
  ans.init <- initialize.clusters(y, G, lambda, eta, Psi, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, gamma.init, phi.init, discrete, limits, yc.init, Ly.init, groups)
  yc.now <- ans.init$yc.now
  phi.now <- ans.init$phi.now
  gamma.now <- ans.init$gamma.now

cat("\nAllocating Memory for Posterior Objects \n")

 ## Allocate posterior objects for which to store MCMC samples
  gamma <- matrix(0, (N.mcmc-begin)%/%every, p)
  phi <- matrix(0, (N.mcmc-begin)%/%every, n)
  nmiss <- sum(is.na(y))
  ymiss <- matrix(0, (N.mcmc-begin)%/%every, nmiss)
  
  accept.gamma <- rep(0, N.mcmc)
  accept.phi <- rep(0, N.mcmc)

  nb <- ceiling(n/bsy)
  accept.yc.m <- rep(0, n)

 ################
 ## Begin MCMC ##
 ################
 cat("\n")
  for(it in 1:N.mcmc){

  cat("\nIteration", it, "out of", N.mcmc)

#print("y update")
    ans.yc <- update.yc.ind(yc.now, y, phi.now, gamma.now, lambda, eta, Psi, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, discrete, limits, aa, bb, lambda0)
    yc.now <- ans.yc$yc.now
    accept.yc.m <- accept.yc.m + ans.yc$accept.m

#print("gamma update")
  ## Update gamma
  ans.gamma <- update.gamma.ind(gamma.now, p.ad, p.swap, yc.now, phi.now, lambda, eta, Psi, rho, Lg, groups, aa, bb, lambda0)
  gamma.now <- ans.gamma$gamma
  accept.gamma[it] <- ans.gamma$accept



#print("phi split/merge update")
#phi.now<<-phi.now
#yc.now<<-yc.now
#gamma.now<<-gamma.now
#lambda.now<<-lambda.now
#eta.now<<-eta.now
#Psi.now<<-Psi.now
#delta.now<<-delta.now

    ans.phi <- update.phi.sm(phi.now, yc.now, gamma.now, lambda, eta, Psi, delta, TT)
    phi.now <- ans.phi$phi
    accept.phi[it] <- ans.phi$accept


#print("phi Gibbs update")
    phi.now <- update.phi(phi.now, yc.now, gamma.now, lambda, eta, Psi, delta)


#print("End of Updates")

   ## record params.now in posterior sample
    if(it>begin && (it-begin)%%every==0){

      ymiss[(it-begin)/every,] <- yc.now[is.na(y)]
      gamma[(it-begin)/every,] <- gamma.now
      phi[(it-begin)/every,] <- phi.now
     }
    

   ## Summarize and Plot posterior
    if((it-begin)%%nplot==0){
      ncv <- sum(gamma.now)
      N.plots <- min(maxplot, max(1,choose(ncv,2)) + 2)
      cols <- min(13, ceiling(sqrt(N.plots)))
      rows <- min(13, ceiling(N.plots/cols))
      it.e <- floor((it-begin)/every)
      ind.now <- max(1,floor(it.e/2),it.e-nback+1):max(1,it.e)
      par(mfrow=c(rows,cols), mar=c(2,2,2,1))

     ## Print likelihood
      like.y <- get.like(yc.now, phi.now, gamma.now, lambda, eta, Psi)
      cat("\nlog(like) =",like.y,"\n")

     ## Print model
      cat("\ngamma =",gamma.now,"\n")

     ## Print acceptance %
      cat("\ngamma acceptance = ", mean(accept.gamma[1:it]))
      cat("\nphi acceptance = ", mean(accept.phi[1:it]))
      cat("\nyc acceptance = ", summary(accept.yc.m/it))
      cat("\n\n\n")
    }
    if((it-begin)%%nplot==0 && it>begin){
      par(mar=c(1,1.5,1.5,.5))
     ## Plot clusters in pairwise scatter of gamma=1 vars
      ind.cv <- which(gamma.now==1)
      if(length(ind.cv) <= 6)
        names.cv <- colnames(y)[ind.cv]
      else
        names.cv <- paste("y",ind.cv,sep="")
      n.cv <- length(ind.cv)
      indp <- 1
      if(n.cv==1){
        ind2 <- ifelse(ind.cv==1,2,1)
        plot(yc.now[,ind.cv], yc.now[,ind2], col=phi.now, main=paste(names.cv[1],"by",colnames(y)[ind2]),pch=16)
      }
      else{
        for(j in 1:(n.cv-1)){
          for(k in (j+1):n.cv){
            if(indp <= N.plots-2)
	      plot(yc.now[,ind.cv[j]], yc.now[,ind.cv[k]], col=phi.now, pch=16, cex=.7, cex.axis=.8, mgp=c(2,.35,0))
	      title(main=paste(names.cv[j],"by",names.cv[k]),line=.5, cex.main=.8)
            indp <- indp + 1	  
	  }
        }
      }
     ## Bar plot of omega.now
      omega.now <- update.vee(phi.now, delta, M=max(phi.now)+1)$omega
      barplot(omega.now, main="omega")
      abline(h=exp(-6), col=2)
      abline(h=exp(-4), col=4)           
    }
  }
  return(list(gamma=gamma, phi=phi, ymiss=ymiss, lambda=lambda, Psi=Psi, eta=eta, delta=delta, y=y, rho=rho, p.ad=p.ad, p.swap=p.swap, prop.sd.y.m=prop.sd.y.m, prop.sd.y.d=prop.sd.y.d, prop.sd.y.c=prop.sd.y.c, prop.sd.y2=prop.sd.y2, bsy=bsy))
}






update.yc.ind <- function(yc.now, y, phi.now, gamma.now, lambda, eta, Psi, bsy, prop.sd.y.m, prop.sd.y.d, prop.sd.y.c, discrete, limits, Aa, Bb, lambda0, use.mean=FALSE, force=FALSE, obs=NULL){

##  Draw mu_m1, Sigma_m11, b_2, Q21, Q22 | y_{-i}, phi, gamma
##  ... then draw missing to complete y_i
##  Also block y_i


  n <- nrow(y)
  p <- ncol(y)
  if(is.null(obs))
    obs <- sample(1:n)
  M <- max(phi.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  na.mat <- is.na(y)
  na.obs <- which(apply(na.mat,1,sum)>0)
  n.na <- length(na.obs)
  nd <- sum(discrete)
  udis <- list()
  for(j in 1:p){
    if(discrete[j]){
      udis[[j]] <- sort(unique(y[,j]))
      udis[[j]] <- c(udis[[j]], Inf)
    }
    else
      udis[[j]] <- c()
  }
#  if(n.na==0 && nd==0 &&)
#    return(list(yc.now=y, accept=rep(1,0)))

  blocks <- list()
  ind.now <- 1
  nb <- ceiling(n/bsy)
  for(b in 1:nb){
    blocks[[b]] <- ind.now:min(ind.now+bsy-1, n)
    ind.now <- ind.now+bsy
  }
  accept.m <- accept.d <- accept.c <- rep(0,n)
  
  for(b in 1:nb){
    yc.p <- yc.now
    ind.b <- obs[blocks[[b]]]
    yc.mb <- yc.now[-ind.b,]
    ans.muQ <- rmuQ(yc.mb, phi.now[-ind.b], gamma.now, lambda, eta, Psi, M)
    mu.b <- ans.muQ$mu
    Q.b <- ans.muQ$Q
    
    d.p <- d.now <- 0
    for(i in ind.b){
      m <- phi.now[i]
      ind.na.i <- is.na(y[i,])
      ind.c.i <- !ind.na.i
      
     ## sample yc for missing y's    
      if(any(ind.na.i)){
        Q.na.i <- Q.b[m,ind.na.i,ind.na.i]
        ans.Q.na.i <- ginv.gp(Q.na.i*1/prop.sd.y.m^2)
        S.na.i.sqrt <- ans.Q.na.i$sqrt.inv
        if(use.mean){
          mu.na.i <- mu.b[m,ind.na.i] - ginv.gp(Q.na.i)$inv%*%Q.b[m,ind.na.i,ind.c.i]%*%(yc.now[i,ind.c.i] - mu.b[m,ind.c.i])
          yc.p[i,ind.na.i] <- rmvnorm(1,mu.na.i,S.sqrt=S.na.i.sqrt)
          d.p <- d.p + dmvnorm(yc.p[i,ind.na.i], mu.na.i, S.inv=1/prop.sd.y.m^2*Q.na.i,log.det=-ans.Q.na.i$log.det)
          d.now <- d.now + dmvnorm(yc.now[i,ind.na.i], mu.na.i, S.inv=1/prop.sd.y.m^2*Q.na.i, log.det=-ans.Q.na.i$log.det)
        }
        else{
          mu.na.i <- yc.now[i,ind.na.i] 
          yc.p[i,ind.na.i] <- rmvnorm(1,mu.na.i,S.sqrt=S.na.i.sqrt)
          d.p <- d.now <- 0
        }
      }

     ## sample latent yc for observed, but discrete y's
      ind.dobs.i <- which(!ind.na.i & discrete)
      ndobs.i <- length(ind.dobs.i)
      for(l in index(1,ndobs.i)){
        j <- ind.dobs.i[l]
        ind.mj <- (1:p)[-j]
#        mu.j <- mu.b[m,j] - 1/Q.b[m,j,j]*Q.b[m,j,ind.mj]%*%(yc.now[i,ind.mj] - mu.b[m,ind.mj])
        mu.j <- yc.now[i,j]
	indu <- which(udis[[j]]==y[i,j])
	aa <- ifelse(indu==1, -Inf, y[i,j])
	bb <- udis[[j]][indu+1]
        yc.p[i,j] <- rtruncnorm(1,a=aa, b=bb, mu.j, prop.sd.y.d*sqrt(1/Q.b[m,j,j]))
#        mu.jp <- mu.j
	mu.jp <- yc.p[i,j]
	d.p <- d.p + log(dtruncnorm(yc.p[i,j], a=aa, b=bb, mu.j, prop.sd.y.d*sqrt(1/Q.b[m,j,j])))
	d.now <- d.now + log(dtruncnorm(yc.now[i,j], a=aa, b=bb, mu.jp, prop.sd.y.d*sqrt(1/Q.b[m,j,j])))
      }

     ## sample latent yc for observed, but left or right censored y's
      ind.lcens <- y[i,]<=limits[,1]
      ind.rcens <- y[i,]>=limits[,2]
      ind.lrobs.i <- which(!discrete & (ind.lcens | ind.rcens))
      nlrobs.i <- length(ind.lrobs.i)
      for(l in index(1,nlrobs.i)){
        j <- ind.lrobs.i[l]
        ind.mj <- (1:p)[-j]
#        mu.j <- mu.b[m,j] - 1/Q.b[m,j,j]*Q.b[m,j,ind.mj]%*%(yc.now[i,ind.mj] - mu.b[m,ind.mj])
        mu.j <- yc.now[i,j]
	aa <- ifelse(ind.rcens[j], limits[j,2], -Inf)
	bb <- ifelse(ind.lcens[j], limits[j,1], Inf)
        yc.p[i,j] <- rtruncnorm(1,a=aa, b=bb, mu.j, prop.sd.y.c*sqrt(1/Q.b[m,j,j]))
#        mu.jp <- mu.j
	mu.jp <- yc.p[i,j]
	d.p <- d.p + log(dtruncnorm(yc.p[i,j], a=aa, b=bb, mu.j, prop.sd.y.c*sqrt(1/Q.b[m,j,j])))
	d.now <- d.now + log(dtruncnorm(yc.now[i,j], a=aa, b=bb, mu.jp, prop.sd.y.c*sqrt(1/Q.b[m,j,j])))
      }
    }
    if(force)
      MH.ratio <- 1
    else{
      like.now <- get.like.ind(yc.now, phi.now, gamma.now, lambda, eta, Psi, Aa, Bb, lambda)
      like.p <- get.like.ind(yc.p, phi.now, gamma.now, lambda, eta, Psi, Aa, Bb, lambda)
      MH.ratio <- exp(like.p - d.p - like.now + d.now)
    }
    if(runif(1) < MH.ratio){
      yc.now <- yc.p
      accept.m[ind.b] <- 1
    }
  }
  return(list(yc.now=yc.now, accept.m=accept.m))
}





update.gamma.ind <- function(gamma.now, p.ad, p.swap, yc.now, phi.now, lambda, eta, Psi, rho, Lg, groups, aa, bb, lambda0){

  accept <- 0
  p <- length(gamma.now)
  sum.as <- p.ad+p.swap
  p.ad <- p.ad/sum.as
  p.swap <- p.swap/sum.as
  for(l in 1:Lg){
    gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)
    d.g.pgnow <- dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    d.g.nowgp <- dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)

    like.now <- get.like.ind(yc.now, phi.now, gamma.now, lambda, eta, Psi, aa, bb, lambda0)
    like.p <- get.like.ind(yc.now, phi.now, gamma.p, lambda, eta, Psi, aa, bb, lambda0)
    pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
    pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
    MH.ratio <- exp(like.p + pi.g.p - d.g.pgnow - like.now - pi.g.now + d.g.nowgp)
    if(runif(1) < MH.ratio){
      gamma.now <- gamma.p
      accept <- 1
    }
  }
  return(list(gamma=gamma.now, accept=accept))
}


get.like.ind <- function(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, aa, bb, lambda0){

  n <- nrow(yc.now)
  p <- ncol(yc.now)
  phi.now <- match(phi.now, sort(unique(phi.now)))
  M <- max(phi.now)
  c.vars <- (gamma.now==1)
  nc.vars <- !c.vars
  y1 <- yc.now[,c.vars,drop=FALSE]
  y2 <- yc.now[,nc.vars,drop=FALSE]
  p1 <- ncol(y1)
  p2 <- ncol(y2)
  P11 <- Psi.now[c.vars,c.vars,drop=FALSE]
  ans.P11 <- ginv.gp(P11)
  y1bar <- apply(y1,2,mean)
  y1barM <- matrix(y1bar,n,p1,byrow=TRUE)
  y2bar <- apply(y2,2,mean)
  y2barM <- matrix(y2bar,n,p2,byrow=TRUE)
  V11 <- crossprod(y1-y1barM) + n*lambda.now/(n+lambda.now)*y1bar%*%t(y1bar) + P11
  ans.V11 <- ginv.gp(V11)
  ldet.P11 <- ans.P11$log.det
  ldet.V11 <- ans.V11$log.det
  sum2.2 <- apply(y2,2,var)*(n-1)
  y2bar <- apply(y2,2,mean)
  D2 <- prod(bb+1/2*(sum2.2+n*lambda0/(n+lambda0)*y2bar^2))
  
  like <- -n*p/2*log(2*pi) + p2/2*log(lambda.now/(n+lambda.now)) + aa*p2*log(bb) - (aa+n/2)*log(D2) + p2*(lgamma(aa+n/2) - lgamma(aa))

  for(m in 1:M){
    ym1 <- y1[phi.now==m,,drop=FALSE]
    nm <- nrow(ym1)
    if(nm > 0){
      ym1bar <- apply(ym1,2,mean)
      ym1barM <- matrix(ym1bar,nm,p1,byrow=TRUE)
      Vm11 <- crossprod(ym1-ym1barM) + nm*lambda.now/(nm+lambda.now)*ym1bar%*%t(ym1bar) + P11
      ldet.Vm11 <- ginv.gp(Vm11)$log.det
      like <- like + p1/2*log(lambda.now/(nm+lambda.now)) + (eta.now-p2)/2*ldet.P11 - (nm+eta.now-p2)/2*ldet.Vm11 + lmgamma(p1,(nm+eta.now-p2)/2) - lmgamma(p1,(eta.now-p2)/2)
    }
  }
  return(like)
}




####################################################################
########################### OLD CODE ###############################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################















####################################################################
####################################################################
####################################################################
















####################################################################
########################### OLD CODE ###############################
####################################################################






get.phi.hat.probs.old <- function(ans.DPMvs, ans.hat=NULL, N.mcmc.yg=1000, post.ind=NULL, Lg=5, bsy=NULL, each.y=10, ynew=NULL){

 ## for each i, sample ymiss, gamma, and hyperparams, with phi[-i] fixed at phi.hat
 ## then sample phi[i]; record Gibbs prob for phi[i].  Report "posterior" average of Gibbs prob.
  N.mcmc <- dim(ans.DPMvs$phi)[1]
  if(is.null(post.ind))
    post.ind <- 1:N.mcmc
  if(is.null(ans.hat))
    ans.hat <- get.phi.hat(ans.DPMvs, post.ind=post.ind)
  phi.hat <- ans.hat$phi.hat
  ind.min <- ans.hat$ind.min

  p.ad <- ans.DPMvs$p.ad
  p.swap <- ans.DPMvs$p.swap
  discrete <- ans.DPMvs$discrete
  limits <- ans.DPMvs$limits
  if(is.null(bsy))
    bsy <- ans.DPMvs$bsy
  prop.sd.y <- ans.DPMvs$prop.sd.y
  y <- ans.DPMvs$y
  yc.now <- y
  ind.na <- is.na(y)
  yc.now[ind.na] <- ans.DPMvs$ymiss[ind.min,]
  phi.now <- phi.hat
  n <- nrow(y)
  if(!is.null(ynew)){
    y <- rbind(y,ynew)
    nnew <- nrow(ynew)
    obs <- (n+1):(n+nnew)
    yc.new <- ynew
    yc.new[is.na(yc.new)] <- 0
    yc.now <- rbind(yc.now,yc.new)
    phi.now <- c(phi.now,rep(1,nnew))
  }
  else
    obs <- 1:n
  gamma.now <- ans.DPMvs$gamma[ind.min,]

  M <- max(phi.hat)
  probs.now <- matrix(0,N.mcmc.yg,n)
  rho <- ans.DPMvs$rho
  delta.hat <- mean(ans.DPMvs$delta[post.ind])
  lambda.hat <- mean(ans.DPMvs$lambda[post.ind])
  eta.hat <- mean(ans.DPMvs$eta[post.ind])
  Psi.hat <- apply(ans.DPMvs$Psi[post.ind,,],c(2,3),mean)
  phi.probs <- foreach(i=obs, .combine=rbind)%dopar%{
    if(i%%10==0)
      print(i)
    accept.i <- rep(0,N.mcmc.yg)
    M.i <- max(phi.hat[-i])+1
    phi.probs.i <- matrix(0,N.mcmc.yg,M+1)
    for(it in 1:N.mcmc.yg){
      for(l in index(1,Lg))
        gamma.now <- update.gamma(gamma.now, p.ad, p.swap,yc.now,phi.now,lambda.hat,eta.hat,Psi.hat,rho)$gamma

        ans.yp <- update.y.phi.i(i, yc.now, y, phi.now, gamma.now, lambda.hat, eta.hat, Psi.hat, delta.hat, bsy, prop.sd.y=1)
        yc.now <- ans.yp$yc.now
	phi.now <- ans.yp$phi.now
	accept.i[it] <- ans.yp$accept

print(phi.now[i])
print(yc.now[i,1:2])

      if(it%%each.y==0)
        yc.now <- update.yc(yc.now, y, phi.now, gamma.now, lambda.hat, eta.hat, Psi.hat, bsy, prop.sd.y, discrete, limits)$yc.now
      phi.probs.i[it,1:M.i] <- get.phi.probs(i, phi.now, yc.now, gamma.now, lambda.hat, eta.hat, Psi.hat, delta.hat, m.ind=1:M.i)
      phi.now[i] <- sample(1:M.i,1,prob=phi.probs.i[it,1:M.i])
    }
    c(mean(accept.i), apply(phi.probs.i,2,mean))
  }
  phi.probs <- matrix(phi.probs,nnew,M+2)
  return(list(accept=phi.probs[,1], phi.probs=phi.probs[,-1]))
}



update.y.phi.i <- function(i, yc.now, y, phi.now, gamma.now, lambda.hat, eta.hat, Psi.hat, delta.hat, bsy, prop.sd.y){

  M <- max(phi.now[-i])+1
  prob.p <- rep(1/M,M)
  prob.p[phi.now[i]] <- .5
  phi.p.i <- sample(1:M,1,prob=prob.p)
  phi.p <- phi.now
  phi.p[i] <- phi.p.i
  yc.p <- yc.now
  accept <- 0
  ind.na.i <- is.na(y[i,])
  if(sum(ind.na.i)>0){
    ans.muQ <- rmuQ(yc.now[-i,], phi.now[-i], gamma.now, lambda.hat, eta.hat, Psi.hat, M)
    mu.i <- ans.muQ$mu
    Q.i <- ans.muQ$Q
    
    ind.c.i <- !ind.na.i
    Q.na.i.p <- Q.i[phi.p[i],ind.na.i,ind.na.i]
    Q.na.i.now <- Q.i[phi.now[i],ind.na.i,ind.na.i]
    ans.Q.na.i.p <- ginv.gp(Q.na.i.p/prop.sd.y^2)
    ans.Q.na.i.now <- ginv.gp(Q.na.i.now/prop.sd.y^2)
    S.na.i.sqrt.p <- ans.Q.na.i.p$sqrt.inv
    S.na.i.sqrt.now <- ans.Q.na.i.now$sqrt.inv
    mu.na.i.p <- mu.i[phi.p[i],ind.na.i] - ginv.gp(Q.na.i.p)$inv%*%Q.i[phi.p[i],ind.na.i,ind.c.i]%*%(y[i,ind.c.i] - mu.i[phi.p[i],ind.c.i])
    mu.na.i.now <- mu.i[phi.now[i],ind.na.i] - ginv.gp(Q.na.i.now)$inv%*%Q.i[phi.now[i],ind.na.i,ind.c.i]%*%(y[i,ind.c.i] - mu.i[phi.now[i],ind.c.i])
    yc.p[i,ind.na.i] <- rmvnorm(1,mu.na.i.p,S.sqrt=S.na.i.sqrt.p)
    d.p <- dmvnorm(yc.p[i,ind.na.i], mu.na.i.p, S.inv=1/prop.sd.y^2*Q.na.i.p, log.det=-ans.Q.na.i.p$log.det)
    d.now <- dmvnorm(yc.now[i,ind.na.i], mu.na.i.now, S.inv=1/prop.sd.y^2*Q.na.i.now, log.det=-ans.Q.na.i.now$log.det)
  }
  else
    d.p <- d.now <- 0
    
  like.now <- get.like(yc.now, phi.now, gamma.now, lambda.hat, eta.hat, Psi.hat)
  like.p <- get.like(yc.p, phi.p, gamma.now, lambda.hat, eta.hat, Psi.hat)
  n1 <- sum(phi.p==phi.now[i])
  n2 <- sum(phi.now==phi.p[i])
  pi.pdivnow <- log(max(n2,delta.hat)) - log(max(n1,delta.hat))
   
  MH.ratio <- exp(like.p - d.p - like.now + d.now + pi.pdivnow)
  if(runif(1) < MH.ratio){
    yc.now <- yc.p
    phi.now <- phi.p
    accept <- 1
  }
  return(list(yc.now=yc.now, phi.now=phi.now, accept=accept))
}



update.phi.AG.old <- function(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=NULL, obs=NULL){

  n <- nrow(yc.now)
  M <- max(phi.now)
  if(M==1)
    return(list(phi.now=rep(1,n), dens=0))
  if(is.null(obs))
    obs <- 1:n
  draw <- is.null(phi.p)  
  dens <- 0
  for(i in obs){
    phi.prev.i <- phi.now[i]
    if(sum(phi.now==phi.prev.i)>1){
      probs.i <- get.phi.probs(i, phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now,m.ind=1:M)
      if(draw)
        phi.now[i] <- sample(1:length(probs.i),1,prob=probs.i)
      else
        phi.now[i] <- phi.p[i]
      dens <- dens + log(probs.i[phi.now[i]])
    }
    else if(!draw){
      if(phi.now[i]!=phi.p[i])
        return(list(phi.now=phi.p, dens=-Inf))
    }
  }
  return(list(phi.now=phi.now, dens=dens))
}



update.gamma.phi.old <- function(gamma.now, phi.now, p.ad, p.swap, yc.now, lambda.now, eta.now, Psi.now, rho, groups){

  accept <- 0
  n <- nrow(yc.now)
  obs <- sample(1:n,n,replace=FALSE)
  gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)
  ans.p <- update.phi.AG(phi.now, yc.now, gamma.p, lambda.now, eta.now, Psi.now, delta.now, phi.p=NULL, obs=obs)
  phi.p <- ans.p$phi
  d.p <- ans.p$dens + dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
  ans.now <- update.phi.AG(phi.p, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=phi.now, obs=obs)
  d.now <- ans.now$dens + dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)

  like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
  like.p <- get.like(yc.now, phi.p, gamma.p, lambda.now, eta.now, Psi.now)

  pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
  pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
  n.p <- as.numeric(table(phi.p))
  n.now <- as.numeric(table(phi.now))
  pi.pdivnow <- pi.g.p - pi.g.now + sum(lfactorial(n.p-1)) - sum(lfactorial(n.now-1))

  MH.ratio <- exp(like.p - like.now + d.now - d.p + pi.pdivnow)

  if(runif(1) < MH.ratio){
    gamma.now <- gamma.p
    phi.now <- phi.p
    accept <- 1
  }
  return(list(gamma.now=gamma.now, phi.now=phi.now, accept=accept))
}









update.phi.AG.old <- function(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=NULL, obs=NULL, thresh=.025){

  n <- nrow(yc.now)
  n.now <- as.numeric(table(phi.now))
  M <- sum(n.now/n > thresh)
#  M <- max(phi.now)
  p <- ncol(yc.now)
  if(is.null(obs))
    obs <- 1:n
  c.vars <- (gamma.now==1)
  y1 <- yc.now[,c.vars,drop=FALSE]
  p1 <- ncol(y1)
  if(p1==1)
    phi.0 <- Mclust(y1, G=M, model="V")$class
  else
    phi.0 <- Mclust(y1, G=M, model="VVV")$class


#phi.p <<- phi.p

 ## match cluster id in phi.0 to closet match in phi.p
  if(!is.null(phi.p)){
    n.vec <- as.numeric(table(phi.p))
    ord.m <- order(-n.vec)
    phi.00 <- phi.0
    m.avail <- rep(TRUE,M)
    for(m1 in ord.m){
      if(all(!m.avail))
        break
      best <- 0
      for(m2 in (1:M)[m.avail]){
        n.match.12 <- sum(phi.0==m2 & phi.p==m1)
        if(n.match.12 > best){
	  best <- n.match.12
	  best.m2 <- m2
	}
      }
      phi.00[phi.0==best.m2] <- m1
      m.avail[best.m2] <- FALSE
    }
    phi.0 <- phi.00
  }

 ## Gibbs update and density
  draw <- is.null(phi.p)  
  dens <- 0
  phi.now <- phi.0
  if(draw){
    for(i in obs){
      phi.prev.i <- phi.now[i]
      probs.i <- get.phi.probs(i, phi.now, yc.now, gamma.now, lambda.now,eta.now,Psi.now,delta.now,m.ind=NULL)
      m.max <- which(probs.i==max(probs.i))
      p.max <- probs.i[m.max]
      probs.i[m.max] <- p.max*.97
      probs.i[-m.max] <- probs.i[-m.max]*(1-probs.i[m.max])/sum(probs.i[-m.max])
      phi.now[i] <- sample(1:length(probs.i),1,prob=probs.i)
      dens <- dens + log(probs.i[phi.now[i]])
      if(!any(phi.now==phi.prev.i))
        phi.now <- match(phi.now, sort(unique(phi.now)))
    }
  }
  else{
    for(i in obs){
      phi.prev.i <- phi.now[i]
      m.ind <- unique(phi.now)
      if(any(m.ind==phi.p[i]))
        m.ind <- c(m.ind, max(m.ind)+1)
      else
        m.ind <- c(m.ind, phi.p[i])
      probs.i <- get.phi.probs(i, phi.now, yc.now, gamma.now, lambda.now,eta.now,Psi.now,delta.now,m.ind=m.ind)
      m.max <- which(probs.i==max(probs.i))
      p.max <- probs.i[m.max]
      probs.i[m.max] <- p.max*.97
      probs.i[-m.max] <- probs.i[-m.max]*(1-probs.i[m.max])/sum(probs.i[-m.max])
      phi.now[i] <- phi.p[i]
      dens <- dens + log(probs.i[phi.now[i]])
    }
  }
  return(list(phi.now=phi.now, dens=dens))
}






update.gamma.old <- function(gamma.now, p.ad, p.swap, yc.now, phi.now, lambda.now, eta.now, Psi.now, rho, Lg, groups){

  accept <- 0
  p <- length(gamma.now)
  sum.as <- p.ad+p.swap
  p.ad <- p.ad/sum.as
  p.swap <- p.swap/sum.as
  while(Lg>0){
    ans <- foreach(l=1:Lg, .combine=rbind)%dopar%{
      gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)
      d.g.pgnow <- dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
      d.g.nowgp <- dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)

      like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
      like.p <- get.like(yc.now, phi.now, gamma.p, lambda.now, eta.now, Psi.now)
      pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
      pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
      MH.ratio <- exp(like.p + pi.g.p - d.g.pgnow - like.now - pi.g.now + d.g.nowgp)
      if(runif(1) < MH.ratio)
        foo <- 1
      else
        foo <- 0
      c(gamma.p, foo)
    }
    if(any(ans[,p+1]==1)){
      ind <- which(ans[,p+1]==1)[1]
      gamma.now <- ans[ind,1:p]
      Lg <- Lg - ind
      accept <- 1
    }
    else
      Lg <- 0

    Lg <- 0
  }
  return(list(gamma=gamma.now, accept=accept))
}







update.gamma.phi <- function(gamma.now, phi.now, p.ad, p.swap, yc.now, lambda.now, eta.now, Psi.now, rho, delta.now, Lgp, groups){

  accept <- 0
  n <- nrow(yc.now)
  sum.ads <- p.ad+p.swap
  p.stay <- 1-(sum.ads)
  sum.new <- (1-p.stay)^(1/Lgp)
  p.ad <- p.ad/sum.ads*sum.new
  p.swap <- p.swap/sum.ads*sum.new
  obs <- sample(1:n)
  
  for(l in index(1,Lgp)){
    gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)  

    ans.p <- update.phi.AG(phi.now, yc.now, gamma.p, lambda.now, eta.now, Psi.now, delta.now, phi.p=NULL, obs=obs)
    phi.p <- ans.p$phi
    d.p <- ans.p$dens + dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    ans.now <- update.phi.AG(phi.p, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, phi.p=phi.now, obs=obs[n:1])
    d.now <- ans.now$dens + dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)

    like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
    like.p <- get.like(yc.now, phi.p, gamma.p, lambda.now, eta.now, Psi.now)

    pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
    pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
    n.p <- as.numeric(table(phi.p))
    n.now <- as.numeric(table(phi.now))
    pi.pdivnow <- pi.g.p - pi.g.now + (length(n.p)-length(n.now))*log(delta.now) + sum(lfactorial(n.p-1)) - sum(lfactorial(n.now-1))
    
    MH.ratio <- exp(like.p - like.now + d.now - d.p + pi.pdivnow)

    if(runif(1) < MH.ratio){
      gamma.now <- gamma.p
      phi.now <- phi.p
      accept <- 1
    }
  }
  return(list(gamma.now=gamma.now, phi.now=phi.now, accept=accept))
}








update.gamma.phi.sm <- function(gamma.now, phi.now, p.ad, p.swap, yc.now, lambda.now, eta.now, Psi.now, rho, delta.now, TT, groups){

  accept <- 0
  n <- nrow(yc.now)
  M <- max(phi.now)
  i <- sample(1:n,1)
  j <- sample((1:n)[-i],1)

  gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)
  pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
  pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
  change <- any(gamma.p!=gamma.now)

  split <- phi.now[i]==phi.now[j]
  obs <- which(phi.now==phi.now[i] | phi.now==phi.now[j])
  obs.s <- sample(obs[obs!=i & obs!=j])

#  gamma.launch <- (gamma.p==1 | gamma.now==1)*1
#  phi.launch <- get.phi.launch(phi.now, yc.now, gamma.launch, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=TT)

  
  if(split){
    phi.launch <- get.phi.launch(phi.now, yc.now, gamma.p, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=0)
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.p, lambda.now, eta.now, Psi.now, delta.now, obs.s, m.ind=c(phi.now[i],M+1))
    phi.p <- ans.s$phi
    d.p <- ans.s$dens + dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    d.now <- dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)
    if(change){
      like.p <- get.like(yc.now, phi.p, gamma.p, lambda.now, eta.now, Psi.now)
      like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
    }
    else{
      like.p <- get.like.sm(yc.now, phi.p, gamma.p, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=M+1)
      like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[i])
    }
    n1 <- sum(phi.p==phi.now[i])
    n2 <- sum(phi.p==M+1)
    pi.pdivnow <- pi.g.p - pi.g.now + log(delta.now) + lfactorial(n1-1) + lfactorial(n2-1) - lfactorial(n1+n2-1)
  }
  else{
    phi.launch <- get.phi.launch(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=0)
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs.s, m.ind=c(phi.now[i],phi.now[j]), phi.p=phi.now)
    phi.p <- phi.now
    phi.p[obs] <- phi.now[i]
    u.phi.p <- sort(unique(phi.p))
    phi.p <- match(phi.p, u.phi.p)
    phi.p.i <- which(u.phi.p==phi.now[i])
    d.p <- dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    d.now <- ans.s$dens + dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)
    if(change){
      like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
      like.p <- get.like(yc.now, phi.p, gamma.p, lambda.now, eta.now, Psi.now)
    }
    else{
      like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[j])
      like.p <- get.like.sm(yc.now, phi.p, gamma.p, lambda.now, eta.now, Psi.now, phi.i=phi.p.i, phi.j=phi.p.i)
    }
    n1 <- sum(phi.now==phi.now[i])
    n2 <- sum(phi.now==phi.now[j])
    pi.pdivnow <- pi.g.p - pi.g.now - log(delta.now) - lfactorial(n1-1) - lfactorial(n2-1) + lfactorial(n1+n2-1)
  }
  MH.ratio <- exp(like.p - like.now + pi.pdivnow + d.now - d.p)
  if(runif(1) < MH.ratio){

#i<<-i
#j<<-j
#yc.now<<-yc.now
#lambda.now<<-lambda.now
#delta.now<<-delta.now
#eta.now<<-eta.now
#Psi.now<<-Psi.now
#phi.now<<-phi.now
#phi.p<<-phi.p
#phi.p.i<<-phi.p.i
#gamma.now<<-gamma.now
#gamma.p<<-gamma.p
#phi.launch<<-phi.launch
#pi.pdivnow<<-pi.pdivnow
#like.p<<-like.p
#like.now<<-like.now
#d.now<<-d.now
#d.p<<-d.p
#MH.ratio<<-MH.ratio

    phi.now <- match(phi.p, sort(unique(phi.p)))
    gamma.now <- gamma.p
    accept <- 1
  }
  return(list(phi.now=phi.now, gamma.now=gamma.now, accept=accept))
}







update.gamma.phi.y.sm <- function(gamma.now, phi.now, p.ad, p.swap, yc.now, lambda.now, eta.now, Psi.now, rho, delta.now, TT, discrete, limits, prop.sd.y1, prop.sd.y2, groups){

  accept <- 0
  n <- nrow(yc.now)
  M <- max(phi.now)
  i <- sample(1:n,1)
  j <- sample((1:n)[-i],1)

  gamma.p <- rgamma.prop(gamma.now, p.ad, p.swap, groups)
  pi.g.p <- sum(gamma.p*log(rho))+sum((1-gamma.p)*log(1-rho))
  pi.g.now <- sum(gamma.now*log(rho))+sum((1-gamma.now)*log(1-rho))
  change <- any(gamma.p!=gamma.now)

  split <- phi.now[i]==phi.now[j]
  obs <- which(phi.now==phi.now[i] | phi.now==phi.now[j])
  obs.s <- sample(obs[obs!=i & obs!=j])

#  gamma.launch <- (gamma.p==1 | gamma.now==1)*1
#  phi.launch <- get.phi.launch(phi.now, yc.now, gamma.launch, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=TT)

  
  if(split){
    phi.launch <- get.phi.launch(phi.now, yc.now, gamma.p, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=0)
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.p, lambda.now, eta.now, Psi.now, delta.now, obs.s, m.ind=c(phi.now[i],M+1))
    phi.p <- ans.s$phi
    d.p <- ans.s$dens + dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    d.now <- dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)
    ans.y.p <- update.yg(yc.now, y, phi.p, gamma.p, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2)
    yc.p <- ans.y.p$yc
    d.p <- d.p + ans.y.p$dens
    ans.y.now <- update.yg(yc.p, y, phi.now, gamma.now, gamma.p, lambda.now, eta.now, Psi.now, yc.p=yc.now, discrete, limits, prop.sd.y1, prop.sd.y2)
    d.now <- d.now + ans.y.now$dens    
    if(change){
      like.p <- get.like(yc.p, phi.p, gamma.p, lambda.now, eta.now, Psi.now)
      like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
    }
    else{
      like.p <- get.like.sm(yc.p, phi.p, gamma.p, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=M+1)
      like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[i])
    }
    n1 <- sum(phi.p==phi.now[i])
    n2 <- sum(phi.p==M+1)
    pi.pdivnow <- pi.g.p - pi.g.now + log(delta.now) + lfactorial(n1-1) + lfactorial(n2-1) - lfactorial(n1+n2-1)
  }
  else{
    phi.launch <- get.phi.launch(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=0)
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs.s, m.ind=c(phi.now[i],phi.now[j]), phi.p=phi.now)
    phi.p <- phi.now
    phi.p[obs] <- phi.now[i]
    u.phi.p <- sort(unique(phi.p))
    phi.p <- match(phi.p, u.phi.p)
    phi.p.i <- which(u.phi.p==phi.now[i])
    d.p <- dgamma.prop(gamma.p, gamma.now, p.ad, p.swap, groups)
    d.now <- ans.s$dens + dgamma.prop(gamma.now, gamma.p, p.ad, p.swap, groups)
    ans.y.p <- update.yg(yc.now, y, phi.p, gamma.p, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2)
    yc.p <- ans.y.p$yc
    d.p <- d.p + ans.y.p$dens
    ans.y.now <- update.yg(yc.p, y, phi.now, gamma.now, gamma.p, lambda.now, eta.now, Psi.now, yc.p=yc.now, discrete, limits, prop.sd.y1, prop.sd.y2)
    d.now <- d.now + ans.y.now$dens
    if(change){
      like.now <- get.like(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now)
      like.p <- get.like(yc.p, phi.p, gamma.p, lambda.now, eta.now, Psi.now)
    }
    else{
      like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[j])
      like.p <- get.like.sm(yc.p, phi.p, gamma.p, lambda.now, eta.now, Psi.now, phi.i=phi.p.i, phi.j=phi.p.i)
    }
    n1 <- sum(phi.now==phi.now[i])
    n2 <- sum(phi.now==phi.now[j])
    pi.pdivnow <- pi.g.p - pi.g.now - log(delta.now) - lfactorial(n1-1) - lfactorial(n2-1) + lfactorial(n1+n2-1)
  }

  MH.ratio <- exp(like.p - like.now + pi.pdivnow + d.now - d.p)
  if(runif(1) < MH.ratio){

#i<<-i
#j<<-j
#yc.now<<-yc.now
#lambda.now<<-lambda.now
#delta.now<<-delta.now
#eta.now<<-eta.now
#Psi.now<<-Psi.now
#phi.now<<-phi.now
#phi.p<<-phi.p
#phi.p.i<<-phi.p.i
#gamma.now<<-gamma.now
#gamma.p<<-gamma.p
#phi.launch<<-phi.launch
#pi.pdivnow<<-pi.pdivnow
#like.p<<-like.p
#like.now<<-like.now
#d.now<<-d.now
#d.p<<-d.p
#MH.ratio<<-MH.ratio

    phi.now <- match(phi.p, sort(unique(phi.p)))
    gamma.now <- gamma.p
    yc.now <- yc.p
    accept <- 1
  }
  return(list(phi.now=phi.now, gamma.now=gamma.now, yc.now=yc.now, accept=accept))
}







update.phi.y.sm <- function(gamma.now, phi.now, yc.now, lambda.now, eta.now, Psi.now, delta.now, TT, discrete, limits, prop.sd.y1, prop.sd.y2){

  accept <- 0
  n <- nrow(yc.now)
  M <- max(phi.now)
  i <- sample(1:n,1)
  j <- sample((1:n)[-i],1)

  split <- phi.now[i]==phi.now[j]
  obs <- which(phi.now==phi.now[i] | phi.now==phi.now[j])
  obs.s <- sample(obs[obs!=i & obs!=j])
  
  if(split){
    phi.launch <- get.phi.launch(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=TT)
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs.s, m.ind=c(phi.now[i],M+1))
    phi.p <- ans.s$phi
    d.p <- ans.s$dens
    d.now <- 0
    ans.y.p <- update.yg(yc.now, y, phi.p, gamma.now, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2)
    yc.p <- ans.y.p$yc
    d.p <- d.p + ans.y.p$dens
    ans.y.now <- update.yg(yc.p, y, phi.now, gamma.now, gamma.now, lambda.now, eta.now, Psi.now, yc.p=yc.now, discrete, limits, prop.sd.y1, prop.sd.y2)
    d.now <- d.now + ans.y.now$dens    
    like.p <- get.like.sm(yc.p, phi.p, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=M+1)
    like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[i])
    n1 <- sum(phi.p==phi.now[i])
    n2 <- sum(phi.p==M+1)
    pi.pdivnow <- log(delta.now) + lfactorial(n1-1) + lfactorial(n2-1) - lfactorial(n1+n2-1)
  }
  else{
    phi.launch <- get.phi.launch(phi.now, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs=obs.s, i, j, TT=TT)
    ans.s <- update.phi.2G(phi.launch, yc.now, gamma.now, lambda.now, eta.now, Psi.now, delta.now, obs.s, m.ind=c(phi.now[i],phi.now[j]), phi.p=phi.now)
    phi.p <- phi.now
    phi.p[obs] <- phi.now[i]
    u.phi.p <- sort(unique(phi.p))
    phi.p <- match(phi.p, u.phi.p)
    phi.p.i <- which(u.phi.p==phi.now[i])
    d.p <- 0
    d.now <- ans.s$dens
    ans.y.p <- update.yg(yc.now, y, phi.p, gamma.now, gamma.now, lambda.now, eta.now, Psi.now, yc.p=NULL, discrete, limits, prop.sd.y1, prop.sd.y2)
    yc.p <- ans.y.p$yc
    d.p <- d.p + ans.y.p$dens
    ans.y.now <- update.yg(yc.p, y, phi.now, gamma.now, gamma.now, lambda.now, eta.now, Psi.now, yc.p=yc.now, discrete, limits, prop.sd.y1, prop.sd.y2)
    d.now <- d.now + ans.y.now$dens
    like.now <- get.like.sm(yc.now, phi.now, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.now[i], phi.j=phi.now[j])
    like.p <- get.like.sm(yc.p, phi.p, gamma.now, lambda.now, eta.now, Psi.now, phi.i=phi.p.i, phi.j=phi.p.i)
    n1 <- sum(phi.now==phi.now[i])
    n2 <- sum(phi.now==phi.now[j])
    pi.pdivnow <-  -log(delta.now) - lfactorial(n1-1) - lfactorial(n2-1) + lfactorial(n1+n2-1)
  }
  MH.ratio <- exp(like.p - like.now + pi.pdivnow + d.now - d.p)
print(MH.ratio)
  if(runif(1) < MH.ratio){
    phi.now <- match(phi.p, sort(unique(phi.p)))
    yc.now <- yc.p
    accept <- 1
  }
  return(list(phi.now=phi.now, yc.now=yc.now, accept=accept))
}

