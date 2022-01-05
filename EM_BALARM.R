


## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Misc functions
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## Function for computing some statistical moments / parameters for an ALARM(1)
## (theoretical transition matrix as a Markov chain, marginal probability,
## marginal variance, lag-1 autocorrelation)
## ----------------------------------------------------------------------------------------------------------------

## For the calculations leading to these expressions, see handwritten notes ALARM + Teugels.
mcchar.fun <- function(b, c)
   {
   	t.tmp <- matrix(nrow = 2, ncol = 2)
   	t.tmp[1,1] <- exp(b+c) / (1 + exp(b+c))
   	t.tmp[1,2] <- 1 / (1 + exp(b+c))
   	t.tmp[2,1] <- exp(c) / (1 + exp(c))
   	t.tmp[2,2] <- 1 / (1 + exp(c))
   	p.tmp <- (1 + exp(c)) / (1 + 2*exp(c) + exp(b+2*c)) 
   	return(list(transition.mat = t.tmp, stationary.distr = c(1 - p.tmp, p.tmp),
   	   variance = p.tmp*(1-p.tmp), lag1ac = (t.tmp[1,1] - (1-p.tmp)) / p.tmp))
   }


## ----------------------------------------------------------------------------------------------------------------
## BIC
## ----------------------------------------------------------------------------------------------------------------

bic.fun <- function(lst, K, G, D, Tlength)
   {
   	i.tmp <- min(which(lst$ll == max(lst$ll, na.rm = T)))
   	ll <- lst$ll[i.tmp]
   	nobs <- ncol(lst$membership.prob)*(Tlength - K)
   	p <- (K+D+2)*G
    -2*ll + p*log(nobs)
   }

# bic.fun <- function(lst)
   # {
   	# ll <- unique(lst$ll[which(lst$ll == max(lst$ll, na.rm = T))])
   	# K <- ncol(lst$est) - 1
   	# G <- nrow(lst$est)
   	# nobs <- ncol(lst$membership.prob) - K
   	# p <- (K+2)*G
    # -2*ll + p*log(nobs)
   # }


## ----------------------------------------------------------------------------------------------------------------
## Initial values for the AR parameters
## ----------------------------------------------------------------------------------------------------------------

## Get a reasonable initial value (all ar coeffs 0, constant is the logit of some
## reasonable marginal probs):
initval_ar.fun <- function(x, G, K)
   {
	mp.tmp <- apply(x, 1, function(vec) sum(vec)/length(vec))
	if(min(mp.tmp) > 0)
	   {
		q.tmp <- min(mp.tmp) + (0:(G-1)) * (max(mp.tmp) - min(mp.tmp)) / G
	   } else {
	   	q.tmp <- (1:G) * max(mp.tmp) / (G+1)
	   }
	init.tmp <- matrix(0, nrow = G, ncol = K+1)
	init.tmp[,K+1] <- log(q.tmp / (1-q.tmp))
	init.tmp
   }



## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Fitting functions for BALARM with only AR terms
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Qmultinom
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## Log-likelihood contribution to Q
## ----------------------------------------------------------------------------------------------------------------

## Function for computing the contibution from the multinomial part to the expected complete log-likelihood Q 
## Arguments:
##   p.cl: cluster probabilities (weight of mixture components); length: (nb of clusters) 
##   w: weight matrix; dimensions: (nb of clusters) x (nb of potential edges) 
multinomll.fun <- function(p.cl, w)
   {
	sumw.tmp <- apply(w, MAR = 1, sum)
	return(sum(log(p.cl)*sumw.tmp))
   }


## ----------------------------------------------------------------------------------------------------------------
## Optimum calculation
## ----------------------------------------------------------------------------------------------------------------

## Function for computing the optimum of the multinomial part of the Q function.
## Argument:
##   w: the full weight array; dimensions: (nb of groups) x (nb of potential edges)  
Qmultinom.optimum.fun <- function(w)
   {
   	J <- ncol(w)
    apply(w, MAR = 1, sum) / J
   }


## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Q-alarm
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## The function can be optimized independently for each group.
## Arguments:
##   x.mat: observations of the state of edges over time; dimensions: (nb of potential edges) x (nb of times) 
##   w: weights for a single group; length: (nb of potential edges)
##   pars: vector of autoregressive coefficients and constant for the same group as w, in this order (to optimize
##       function for); length: (autoregressive order + 1)
Qbinom.singlecluster.fun <- function(pars, x.mat, w)
   {
   	# x.mat <- xx.tmp
	# w <- table.w.fun( pars.mat = pars.mat, p.cl = p.cl, x.mat = xx.tmp)[1,]
	# pars <- pars.mat[1,]
	## Compute the log-likelihood of each time series:
	ll.tmp <- t(apply(x.mat, MAR = 1, singlets.balarmll.fun, pars = pars) )
    sum(ll.tmp*w)
   }


## ----------------------------------------------------------------------------------------------------------------
## Single-data loglikelihood for an ALARM
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##   pars: time series parameters  for a single group (order: c(AR pars, constant))
##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
singledata.balarmll.fun <- function(pars, x, time)
   {
   	K <- length(pars) - 1
   	if(time < K+1) stop("Missing data in autoregressive sum for this time value (time < K+1)")
	b.tmp <- pars[-length(pars)]
	c.tmp <- pars[length(pars)]
	eta <- sum(x[(time-1):(time-K)]*b.tmp) + c.tmp
	return(x[time]*eta - log(1+exp(eta)))
   }


## ----------------------------------------------------------------------------------------------------------------
## Single-ts loglikelihood for an ALARM
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##   pars: time series parameters  for a single group (order: c(AR pars, constant))
##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
singlets.balarmll.fun <- function(pars, x)
   {
   	K <- length(pars) - 1
    sum(sapply((K+1):length(x), singledata.balarmll.fun, pars = pars, x = x))
   }


## ----------------------------------------------------------------------------------------------------------------
## Singledata probability density for an ALARM
## ----------------------------------------------------------------------------------------------------------------
  
singledata.balarmpdf.fun <- function(pars, x, time)
   {
   	K <- length(pars) - 1
   	if(time < K+1) stop("Missing data in autoregressive sum for this time value (time < K+1)")
	b.tmp <- pars[-length(pars)]
	c.tmp <- pars[length(pars)]
	eta <- sum(x[(time-1):(time-K)]*b.tmp) + c.tmp
	return(exp(x[time]*eta - log(1+exp(eta))))
   }
 
singledata.balarmpdf.fun2 <- function(pars, x, time)
   {
   	K <- length(pars) - 1
   	if(time < K+1) stop("Missing data in autoregressive sum for this time value (time < K+1)")
	b.tmp <- pars[-length(pars)]
	c.tmp <- pars[length(pars)]
	eta <- sum(x[(time-1):(time-K)]*b.tmp) + c.tmp
	fac1 <- (exp(eta) / (1+exp(eta)))^x[time]
	fac2 <- (1 / (1+exp(eta)))^(1-x[time])
	return(fac1*fac2)
   }


## ----------------------------------------------------------------------------------------------------------------
## Singlets probability density for an ALARM
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##   pars: time series parameters  for a single group (order: c(AR pars, constant))
##   x: the time series of a single edge 
singlets.balarmpdf.fun <- function(pars, x)
   {
   	K <- length(pars) - 1
    prod(sapply((K+1):length(x), singledata.balarmpdf.fun, pars = pars, x = x))
   }


## ----------------------------------------------------------------------------------------------------------------
## Weight calculation
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## For not too big data sets, not too small probabilities:

## Internal function to compute the weight for the observed time series of an edge
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x: data vector; length: (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components); length: (nb of clusters)
singlets.w.fun <- function(pars.mat, x, p.cl)
   {
   	v.tmp <- apply(pars.mat, 1, singlets.balarmpdf.fun, x = x) * p.cl
   	return(v.tmp / sum(v.tmp))
   }

## Function to compute the weight table for all clusters, for a whole observed data set
## (every edge).
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x.mat: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components)
## Value: 
##    a matrix of dimensions (nb of clusters) x (nb of potential edges)
table.w.fun <- function(pars.mat, x.mat, p.cl)
   {
	m.tmp <- apply(x.mat, MAR = 1, FUN = singlets.w.fun, pars.mat = pars.mat, p.cl = p.cl)
	return(m.tmp)
   }

## ----------------------------------------------------------------------------------------------------------------
## For big data sets or small probabilities:

## Alternative internal function to compute the weight for the observed time series of an edge, using the 
## log-likelihood (hopefully works when there are too many prob values to multiply together in 
## singlets.cbalarmpdf.fun)
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x: data vector; length: (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components); length: (nb of clusters)
singlets.w.ll.fun <- function(pars.mat, x, p.cl)
   {
   	## Compute the value of  log( [likelihoods of the ts in the cluster] * [prob of the cluster])
   	## for each cluster:
   	v.tmp <- apply(pars.mat, 1, singlets.balarmll.fun, x = x) + log(p.cl)
   	u.tmp <- sapply(v.tmp, function(ll1, all.ll) 
   	   {
   	   	sum(exp(all.ll - ll1))^(-1)
   	   }, all.ll = v.tmp)
   	return(u.tmp)
   }

## Alternative function to compute the weight table for all clusters, for a whole observed data set
## (every edge), using the log-likelihood (hopefully works when there are too many
## prob values to multiply together in singlets.cbalarmpdf.fun)
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x.mat: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components)
##    model.matrix: the matrix of the covariate values for the regression part
## Value: 
##    a matrix of dimensions (nb of clusters) x (nb of potential edges)
table.w.ll.fun <- function(pars.mat, x.mat, p.cl)
   {
	m.tmp <- apply(x.mat, MAR = 1, FUN = singlets.w.ll.fun, pars.mat = pars.mat, p.cl = p.cl)
	return(m.tmp)
   }

## table.w.ll.fun(pars.mat, x.mat[1:10,], p.cl = c(0.5,0.5))
## table.w.fun(pars.mat, x.mat[1:10,], p.cl = c(0.5,0.5))
## Looks alright now.


## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## EM algorithm
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

##     verbose: if 0: no output; 
##              if 1: only writes out the CPU spent at the end of each iteration;
##              if 2: outputs intermediate resulst for debugging
EM_BALARM.fun <- function(initprob, initpars, x, verbose = 1, tol.iter = 1e-6, max.iter = 250)	
   {
	# initprob <- rep(1/3, 3)
	# initpars <- matrix((1:6)/6, nrow = 3, ncol = 2)
	# x defined above, under Data
	G <- length(initprob)
	K <- ncol(initpars) - 1
	Nedges <- nrow(x)
	## Starting values stored:
	p.cl <- initprob
	if(verbose == 2)	cat("p.cl =", p.cl, "\n")
	pars.mat <- initpars  
	if(verbose == 2)	cat("pars.mat =", pars.mat, "\n")

	## Create storage for results during iteration:
	mll.tmp <- -5e8
	all.tmp <- -5e8
	ll.tmp <- -1e9
	clprob.tmp <- matrix(p.cl, ncol = 1, nrow = G)
	membprob.tmp <- array(NA, dim = c(1, G, Nedges))
	pars.tmp <- array(pars.mat, dim = c(1, nrow(pars.mat), ncol(pars.mat)))
	t.tmp <- NA
	
	## Start iteration
	p.old <- p.cl
	if(verbose == 2)	cat("p.old =", p.old, "\n")
	pars.old <- pars.mat
	if(verbose == 2)	 cat("pars.old =", pars.old, "\n")
	w.old <- table.w.ll.fun(pars.mat = pars.old, x.mat = x, p.cl = p.old)
	membprob.tmp[1,,] <- w.old
	if(verbose == 2)	cat("w.old[1,1:10] =", w.old[1,1:10], "\n")
	if(verbose == 2)	cat("w.old[2,1:10] =", w.old[2,1:10], "\n")
	ii <- 2
	relimp.tmp <- 1e10
	arelimp.tmp <- 1e10
	
	## Iteration
	## while(abs(relimp.tmp) > 1e-6)
	while(relimp.tmp > tol.iter & ii <= max.iter)
	## while(ii <= max.iter)
	   {
		s.tmp <- system.time({
			p.new <- Qmultinom.optimum.fun(w = w.old)
			if(verbose == 2)	 cat("p.new =", p.new, "\n")
			o.tmp <- list()
			for(gg in 1:G)
			    o.tmp[[gg]] <- optim(par = pars.old[gg,], fn = Qbinom.singlecluster.fun, hessian = F,
			        control = list(fnscale = -1, trace = 0), x.mat = x, w = w.old[gg,])
			pars.new <- t(sapply(o.tmp, function(lst) lst$par))
			w.new <- table.w.ll.fun(pars.mat = pars.new, x.mat = x, p.cl = p.new)
			if(verbose == 2)	 cat("pars.new =", pars.new, "\n")
			apply(w.old, 1, sum)
			apply(w.new, 1, sum)
			mll.tmp[ii] <- multinomll.fun(p.cl = p.old, w = w.old)
			all.tmp[ii] <- sum(sapply(o.tmp, function(lst) lst$value))
			if(verbose == 2)	 cat("all.tmp =", all.tmp, "\n")
			ll.tmp[ii] <- all.tmp[ii] + mll.tmp[ii]
			if(verbose == 2)	 cat("ll.tmp =", ll.tmp, "\n")
			relimp.tmp <- (ll.tmp[ii] - ll.tmp[ii-1]) / abs(ll.tmp[ii])
			arelimp.tmp <- (all.tmp[ii] - all.tmp[ii-1]) / abs(all.tmp[ii])
			if(verbose == 2)	 cat("Rel. improvement on total likelihood:", relimp.tmp, "\n")
			# if(verbose == 2)	 cat("Rel. improvement on ALARM likelihood:", arelimp.tmp, "\n")
			clprob.tmp <- cbind(clprob.tmp, as.numeric(p.new))
			p.tmp <- array(dim = c(ii, nrow(pars.mat), ncol(pars.mat)))
			p.tmp[-ii,,] <- pars.tmp
			p.tmp[ii,,] <- pars.new
			pars.tmp <- p.tmp
			rm(p.tmp)
			p.tmp <- array(dim = c(ii, G, Nedges))
			p.tmp[-ii,,] <- membprob.tmp
			p.tmp[ii,,] <- w.new
			membprob.tmp <- p.tmp
			rm(p.tmp)
		 	if(verbose == 2)	 cat("w.new[1,1:10] =", w.new[1,1:10], "\n")
			if(verbose == 2)	 cat("w.new[2,1:10] =", w.new[2,1:10], "\n")
			w.old <- w.new
			p.old <- p.new
			pars.old <- pars.new
		   })[3]
		if(verbose > 0) cat(paste0("Iteration ", ii, " is done, time: ", round(s.tmp, 3), "s \n"))
		t.tmp[ii] <- s.tmp
		ii <-  ii+1
		if(verbose == 2)	 cat("relimp.tmp = ", relimp.tmp, "\n")
		# if(verbose == 2)	 cat("arelimp.tmp = ", arelimp.tmp, "\n")
		if(verbose == 2)	 cat("tol.iter = ", tol.iter, "\n")
		if(verbose == 2)	 cat("New ii = ", ii, "\n")
		if(verbose == 2)	 cat("max.iter = ", max.iter, "\n")
		if(verbose == 2)	 cat("Condition (relimp.tmp > tol.iter & ii <= max.iter) =", relimp.tmp > tol.iter & ii <= max.iter, "\n")
		# if(verbose == 2)	 cat("Condition (arelimp.tmp > tol.iter & ii <= max.iter) =", arelimp.tmp > tol.iter & ii <= max.iter, "\n")
		if(verbose == 2)	 cat("\n")
	   }
	list(est = pars.tmp, cluster.prob = clprob.tmp, membership.prob = membprob.tmp, ll = ll.tmp,
	   multinom.ll = mll.tmp, alarm.ll = all.tmp)
   }  


## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Fitting functions for BALARM with covariates and AR terms
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## The multinomial part is the same as for BALARM.


## ----------------------------------------------------------------------------------------------------------------
## Singledata probability density for an ALARM with covariates
## ----------------------------------------------------------------------------------------------------------------
 
## First version, the likelihood expanded and simplified to two terms in the exponent.
## Probably faster than the direct implementation of the likelihood.
## Arguments:
##   pars: time series parameters  for a single group (necessary order: c(regr. pars, AR pars, constant))
##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
##   model.matrix: the matrix of the covariate values for the regression part. Dimensions: 
##      (number of time points) x (number of covariates)
##   time: the index in the time series for which the probability density is required
singledata.cbalarmpdf.fun <- function(pars, x, model.matrix, time)
   {
   	D <- ncol(model.matrix)
   	K <- length(pars) - D - 1   	
   	if(time < K+1) stop("Missing data in autoregressive sum for this time value (time < K+1)")
   	## Regression part:
   	a.tmp <- pars[1:D]
   	## Autoregressive part:
	b.tmp <- pars[-c(1:D, length(pars))]
	## Constant:
	c.tmp <- pars[length(pars)]
	eta <- sum(model.matrix[time, ]*a.tmp) + sum(x[(time-1):(time-K)]*b.tmp) + c.tmp
	return(exp(x[time]*eta - log(1+exp(eta))))
   }

## Second version, just using the logit expression in its unsimplified form. Easier to check.
## Arguments:
##   pars: time series parameters  for a single group (necessary order: c(regr. pars, AR pars, constant))
##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
##   model.matrix: the matrix of the covariate values for the regression part
##   time: the index in the time series for which the probability density is required
singledata.cbalarmpdf.fun2 <- function(pars, x, model.matrix, time)
   {
   	D <- ncol(model.matrix)
   	K <- length(pars) - D - 1   	
   	if(time < K+1) stop("Missing data in autoregressive sum for this time value (time < K+1)")
   	## Regression part:
   	a.tmp <- pars[1:D]
   	## Autoregressive part:
	b.tmp <- pars[-c(1:D, length(pars))]
	## Constant:
	c.tmp <- pars[length(pars)]
	eta <- sum(model.matrix[time, ]*a.tmp) + sum(x[(time-1):(time-K)]*b.tmp) + c.tmp
	fac1 <- (exp(eta) / (1+exp(eta)))^x[time]
	fac2 <- (1 / (1+exp(eta)))^(1-x[time])
	return(fac1*fac2)
   }

## Logarithm of this: 
## First version, the likelihood expanded and simplified to two terms in the exponent.
## Probably faster than the direct implementation of the likelihood.
## Arguments:
##   pars: time series parameters  for a single group (necessary order: c(regr. pars, AR pars, constant))
##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
##   model.matrix: the matrix of the covariate values for the regression part. Dimensions: 
##      (number of time points) x (number of covariates)
##   time: the index in the time series for which the probability density is required
singledata.cbalarmll.fun <- function(pars, x, model.matrix, time)
   {
   	D <- ncol(model.matrix)
   	K <- length(pars) - D - 1   	
   	if(time < K+1) stop("Missing data in autoregressive sum for this time value (time < K+1)")
   	## Regression part:
   	a.tmp <- pars[1:D]
   	## Autoregressive part:
	b.tmp <- pars[-c(1:D, length(pars))]
	## Constant:
	c.tmp <- pars[length(pars)]
	eta <- sum(model.matrix[time, ]*a.tmp) + sum(x[(time-1):(time-K)]*b.tmp) + c.tmp
	return(x[time]*eta - log(1+exp(eta)))
   }


## ----------------------------------------------------------------------------------------------------------------
## Singlets probability density and log-likelihood for an ALARM with covariates
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##   pars: time series parameters  for a single group (necessary order: c(regr. pars, AR pars, constant))
##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
##   model.matrix: the matrix of the covariate values for the regression part
singlets.cbalarmpdf.fun <- function(pars, x, model.matrix)
   {
   	K <- length(pars) - 1
    prod(sapply((K+1):length(x), singledata.cbalarmpdf.fun, pars = pars, x = x, model.matrix = model.matrix))
   }

## Arguments:
##   pars: time series parameters  for a single group (necessary order: c(regr. pars, AR pars, constant))
##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
##   model.matrix: the matrix of the covariate values for the regression part
singlets.cbalarmll.fun <- function(pars, x, model.matrix)
   {
   	K <- length(pars) - 1
   	l.tmp <- sapply((K+1):length(x), singledata.cbalarmll.fun, pars = pars, x = x, model.matrix = model.matrix)
    sum(l.tmp)
   }


## ----------------------------------------------------------------------------------------------------------------
## Weight calculation with covariates
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## For small data sets, not very small prob values:

## Internal function to compute the weight for the observed time series of an edge
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x: data vector; length: (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components); length: (nb of clusters)
##   model.matrix: the matrix of the covariate values for the regression part
singlets.cw.fun <- function(pars.mat, x, p.cl, model.matrix)
   {
   	v.tmp <- apply(pars.mat, 1, singlets.cbalarmpdf.fun, x = x, model.matrix = model.matrix) * p.cl
   	return(v.tmp / sum(v.tmp))
   }

## Function to compute the weight table for all clusters, for a whole observed data set
## (every edge).
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x.mat: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components)
##    model.matrix: the matrix of the covariate values for the regression part
## Value: 
##    a matrix of dimensions (nb of clusters) x (nb of potential edges)
table.cw.fun <- function(pars.mat, x.mat, p.cl, model.matrix)
   {
	m.tmp <- apply(x.mat, MAR = 1, FUN = singlets.cw.fun, pars.mat = pars.mat, p.cl = p.cl,
	   model.matrix = model.matrix)
	return(m.tmp)
   }

## ----------------------------------------------------------------------------------------------------------------
## For large data sets, or very small prob values:

## Internal function to compute the weight for the observed time series of an edge, using the 
## log-likelihood (hopefully works when there are too many prob values to multiply together in 
## singlets.cbalarmpdf.fun)
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x: data vector; length: (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components); length: (nb of clusters)
##   model.matrix: the matrix of the covariate values for the regression part
singlets.cw.ll.fun <- function(pars.mat, x, p.cl, model.matrix)
   {
   	## Compute the value of  log( [likelihoods of the ts in the cluster] * [prob of the cluster] )
   	## for each cluster:
   	v.tmp <- apply(pars.mat, 1, singlets.cbalarmll.fun, x = x, model.matrix = model.matrix) + log(p.cl)
   	u.tmp <- sapply(v.tmp, function(ll1, all.ll) 
   	   {
   	   	sum(exp(all.ll - ll1))^(-1)
   	   }, all.ll = v.tmp)
   	return(u.tmp)
   }
# singlets.cw.ll.fun(pars.mat = init.pars, x = x1[7,], p.cl = c(0.5, 0.5), model.matrix = model.matrix1)
# singlets.cw.ll.fun(pars.mat = init.pars, x = x1[8,], p.cl = c(0.5, 0.5), model.matrix = model.matrix1)
# singlets.cw.fun(pars.mat = init.pars, x = x1[8,], p.cl = c(0.5, 0.5), model.matrix = model.matrix1)
## Looks okay: the latter two gives the same value, while the first is now not NaN, but
## sg very close to (1,0).

## For large data sets, or very small prob values:
## Function to compute the weight table for all clusters, for a whole observed data set
## (every edge), using the log-likelihood (hopefully works when there are too many
## prob values to multiply together in singlets.cbalarmpdf.fun)
## Arguments:
##    pars.mat: matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x.mat: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    p.cl: cluster probabilities (weight of mixture components)
##    model.matrix: the matrix of the covariate values for the regression part
## Value: 
##    a matrix of dimensions (nb of clusters) x (nb of potential edges)
table.cw.ll.fun <- function(pars.mat, x.mat, p.cl, model.matrix)
   {
	m.tmp <- apply(x.mat, MAR = 1, FUN = singlets.cw.ll.fun, pars.mat = pars.mat, p.cl = p.cl,
	   model.matrix = model.matrix)
	return(m.tmp)
   }


## ----------------------------------------------------------------------------------------------------------------
## Q-cbalarm with use of sparsity
## ----------------------------------------------------------------------------------------------------------------

## The function can be optimized independently for each group.
## Arguments:
##   pars: vector of time covariate coefficients, autoregressive coefficients and constant for the 
##       same group as w, in this order (to optimize function for); 
##       length: (nb. of time covariates + autoregressive order + 1)
##   x.mat1: matrix of non-empty time series; dimensions: (nb of edges with at least one switched-on link) x 
##      (nb of times) 
##   model.matrix0: the matrix of the covariate values for the regression part. Dimensions: 
##      (nb of times) x (number of covariates)
##   model.array1: the matrices of the autoregressive covariates (the past values of the ts)
##      for all non-empty time series, arranged in an array. To the dimension of autoregressive covariates,
##      one extra column is added, an indicator variable which indicates whether the row contains 
##      at least one nonzero value. Also, this array here must not contain the regressors corresponding
##      to the first K elements of the time series, where the regressors in the AR term would be NA.
##      Dimensions therefore: 
##      (number of non-empty ts) x (nb of times - autoregressive order) x (autoregressive order + 1).
##   w0, w1: weights for a single cluster, separated for empty and non-empty ; 
##     length: (nb of potential edges)
singlecluster.cbalarmQ.fun2 <- function(pars, x.mat1, model.matrix0, model.array1, w0, w1)
   {
	D <- ncol(model.matrix0)
   	K <- dim(model.array1)[3] - 1
   	## Compute the sum_d (a_d * f_d(t_i) + c) that is common to everything, even for the empty time series:
    d.tmp <- apply(model.matrix0[-(1:K), ], 1, function(vec) sum(vec*pars[1:D]) + pars[length(pars)])
    ## Compute the likelihood contribution of the empty time series, which is equal to a fn of only d.tmp,
    ## and is weighted by the weight vector.
	ll.tmp0 <- sum(w0) * sum(log(1 + exp(d.tmp)))
	## Compute the single-ts log-likelihood for the non-empty time series.
	## First, compute the linear predictor eta for the non-empty time series:
	eta.tmp <- apply(model.array1, 1, function(mat)
	   {
	   	# mat <- model.array1[1,,]
	   	## Check the last column of the matrix for the nonzero AR terms and restrict
	   	## the matrix to only these:
	   	mat2 <- as.matrix(mat[mat[, ncol(mat)] == 1, -ncol(mat)])
	   	## Compute the contribution of the AR terms to the linear predictor at each time point
	   	## where it was nonzero:
	   	ar.tmp <- apply(mat2, 1, function(vec) sum(vec*pars[(D+1):(D+K)]))
	   	v.tmp <- d.tmp
	   	v.tmp[mat[, ncol(mat)] == 1] <- d.tmp[mat[, ncol(mat)] == 1] + ar.tmp
	   	v.tmp
	   })
    ## Then, first term of the log-likelihood contribution of the non-empty sequences:
	ll.tmp11 <- sum(w1 * sapply(1:nrow(x.mat1), function(ii) sum(eta.tmp[,ii]*x.mat1[ii,-(1:K)])))
    ## Second term of the log-likelihood contribution of the non-empty sequence:
	ll.tmp12 <- sum(w1 * apply(eta.tmp, 2, function(vec) sum(log(1 + exp(vec)))))
	## Finally the total likelihood:
    -  ll.tmp0 + ll.tmp11 - ll.tmp12
   }


## ----------------------------------------------------------------------------------------------------------------
## Single-cluster all-data loglikelihood gradient for an ALARM with covariates
## ----------------------------------------------------------------------------------------------------------------

## The function can be optimized independently for each group.
## Arguments:
##   pars: vector of time covariate coefficients, autoregressive coefficients and constant for the 
##       same group as w, in this order (to optimize function for); 
##       length: (nb. of time covariates + autoregressive order + 1)
##   x.mat1: matrix of non-empty time series; dimensions: (nb of edges with at least one switched-on link) x 
##      (nb of times) 
##   model.matrix0: the matrix of the covariate values for the regression part. Dimensions: 
##      (nb of times) x (number of covariates)
##   model.array1: the matrices of the autoregressive covariates (the past values of the ts)
##      for all non-empty time series, arranged in an array. To the dimension of autoregressive covariates,
##      one extra column is added, an indicator variable which indicates whether the row contains 
##      at least one nonzero value. Also, this array here must not contain the regressors corresponding
##      to the first K elements of the time series, where the regressors in the AR term would be NA.
##      Dimensions therefore: 
##      (number of non-empty ts) x (nb of times - autoregressive order) x (autoregressive order + 1).
##   w0, w1: weights for a single cluster, separated for empty and non-empty ; 
##     length: (nb of potential edges)
singlecluster.cbalarmQder.fun2 <- function(pars, x.mat1, model.matrix0, model.array1, w0, w1)
   {
	D <- ncol(model.matrix0)
   	K <- dim(model.array1)[3] - 1
   	## Compute the sum_d (a_d * f_d(t_i) + c) that is common to everything, even for the empty time series:
    d.tmp <- apply(model.matrix0[-(1:K), ], 1, function(vec) sum(vec*pars[1:D]) + pars[length(pars)])
    ## The common factor for all empty series, the derivative of the loglikelihood with respect 
    ## to eta occurring in the derivative:
	commonfac.tmp0 <- - 1/(1 + exp(-d.tmp))
	dercov.tmp0 <- apply(model.matrix0[-(1:K),], 2, function(vec) sum(w0) * sum(vec*commonfac.tmp0))
	derar.tmp0 <- rep(0, K)
	derc.tmp0 <- sum(w0) * sum(commonfac.tmp0)
	grad.tmp0 <- c(dercov.tmp0, derar.tmp0, derc.tmp0)
	## Compute the gradient single-ts log-likelihood for the non-empty time series
	## (this is the same as in the function for the loglikelihood).
	## First, compute the linear predictor eta for the non-empty time series:
	## (output: a matrix with nb of rows = nb of time points, nb of cols = nb of time series)
	eta.tmp <- apply(model.array1, 1, function(mat)
	   {
	   	# mat <- model.array1[1,,]
	   	## Check the last column of the matrix for the nonzero AR terms and restrict
	   	## the matrix to only these:
	   	mat2 <- as.matrix(mat[mat[, ncol(mat)] == 1, -ncol(mat)])
	   	## Compute the contribution of the AR terms to the linear predictor at each time point
	   	## where it was nonzero:
	   	ar.tmp <- apply(mat2, 1, function(vec) sum(vec*pars[(D+1):(D+K)]))
	   	v.tmp <- d.tmp
	   	v.tmp[mat[, ncol(mat)] == 1] <- d.tmp[mat[, ncol(mat)] == 1] + ar.tmp
	   	v.tmp
	   })
	## Second, compute the derivatives of the log-likelihood by eta now for non-empty time series
	## (now they will not be common for all non-empty ts). The first term of the product is the sum 
	## of two terms, one is the time series values (NOT the autoregressive terms), second is the 
	## - (1 + exp(-eta))^(-1) thing.
	term.tmp12 <- - 1/(1 + exp(-eta.tmp))
	## Find the derivative of the likelihood with respect to eta by summing the input matrix and
	## the above:
	commonfac.tmp1 <- t(x.mat1[,-(1:K)]) + term.tmp12
    ## Then, contribution of the non-empty sequences to the gradient of the log-likelihood:
    ## The gradient with respect to the time-covariate coefficients:
    dercov.tmp1 <- t(as.matrix(w1)) %*% t(commonfac.tmp1) %*% model.matrix0[-(1:K), ]
    ## The gradient with respect to the autoregressive coefficients:
    ## First the contributions to this of each time series:
	r.tmp <- apply(array(model.array1[,,-(K+1)], dim = c(dim(model.array1)[1],   
	   dim(model.array1)[2], dim(model.array1)[3]-1)), 3, function(mat)
		   {
		   	## dim of mat: (nb of time series) x (nb of time points)
		   	m.tmp <- mat * t(commonfac.tmp1)
		   	apply(m.tmp, 1, sum)
		   })
	## Then, the weighted sum of the contributions from all the time series:
	derar.tmp1 <- apply(r.tmp, 2, function(vec) sum(w1 * vec))
	## The derivative wrt the constant:
	derc.tmp1 <- sum(w1 * apply(commonfac.tmp1, 2, sum))
	grad.tmp1 <- c(dercov.tmp1, derar.tmp1, derc.tmp1)
    ## The total gradient is the sum of the two terms, the contribution from the empty series
    ## and the contribution from the non-empty ones:
    grad.tmp0 + grad.tmp1
    # list(grad.empty = grad.tmp0, grad.nonempty = grad.tmp1)
   }


## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## EM algorithm with covariates, with BFGS and with gradient
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##    initprob: initial guess for cluster probabilities (weight of mixture components)
##    initpars: initial guess for the matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    model.matrix: the matrix of the covariate values for the regression part
##    verbose: if 0: no output; 
##             if 1: only writes out the CPU spent at the end of each iteration;
##             if 2: outputs intermediate results for debugging
##    tol.iter: minimum relative improvement on the log-likelihood required for continuing the iteration
##    max.iter: the maximum number of iterations allowed
EM_cBALARM.fun <- function(initprob, initpars, model.matrix, x, verbose = 1, tol.iter = 1e-6, max.iter = 250)	
   {
	G <- length(initprob)
	D <- ncol(model.matrix)
	K <- ncol(initpars) - D - 1
	Ntimes <- ncol(x)
	Nedges <- nrow(x)
	## Starting values stored:
	p.cl <- initprob
	if(verbose == 2)	cat("p.cl =", p.cl, "\n")
	if(verbose == 2)	cat("Initial pars =", initpars, "\n")

	## Create storage for results during iteration:
	mll.tmp <- -5e8
	all.tmp <- -5e8
	ll.tmp <- -1e9
	clprob.tmp <- matrix(p.cl, ncol = 1, nrow = G)
	membprob.tmp <- array(NA, dim = c(1, G, Nedges))
	pars.tmp <- array(initpars, dim = c(1, nrow(initpars), ncol(initpars)))
	t.tmp <- NA
	
	## Partition the data set into sets of empty and non-empty time series
	i0 <- apply(x,1,sum) == 0
	i1 <- apply(x,1,sum) > 0
	x.mat1 <- x[i1, ]
	model.matrix0 <- model.matrix
	model.array1 <- array(dim = c(nrow(x.mat1), ncol(x.mat1)-K, K+1))
	for(ii in 1:nrow(x.mat1))
	   {
	   	model.array1[ii,,1:K] <- sapply(1:K, function(k)  x.mat1[ii, (K+1-k):(Ntimes-k)])
	   	model.array1[ii,,K+1] <- as.numeric(apply(as.matrix(model.array1[ii,,1:K]), 1, function(vec) any(vec == 1)))
	   }

	## Start iteration (computation of weights still with the old, slow procedure, but at least,
	## it's only once per iteration)
	p.old <- initprob
	if(verbose == 2)	cat("p.old =", p.old, "\n")
	pars.old <- initpars
	if(verbose == 2)	 cat("pars.old =", pars.old, "\n")
	w.old <- table.cw.ll.fun(pars.mat = pars.old, x.mat = x, p.cl = p.old, model.matrix = model.matrix)
	membprob.tmp[1,,] <- w.old
	w.old0 <- w.old[, i0]
	w.old1 <- w.old[, i1]
	if(verbose == 2)	cat("w.old0[1,1:10] =", w.old[1,1:10], "\n")
	if(verbose == 2)	cat("w.old1[1,1:10] =", w.old[2,1:10], "\n")
	ii <- 2
	relimp.tmp <- 1e10
	arelimp.tmp <- 1e10
	
	## Iteration
	while(relimp.tmp > tol.iter & ii <= max.iter)
	   {
		s.tmp <- system.time({
			p.new <- Qmultinom.optimum.fun(w = w.old)
			if(verbose == 2)	 cat("p.new =", p.new, "\n")
			o.tmp <- list()
			for(gg in 1:G)
			    o.tmp[[gg]] <- optim(par = pars.old[gg,], fn = singlecluster.cbalarmQ.fun2, hessian = F,
			        gr = singlecluster.cbalarmQder.fun2, control = list(fnscale = -1, trace = 0), 
			        x.mat1 = x.mat1, w0 = w.old0[gg,], w1 = w.old1[gg,], model.matrix0 = model.matrix, 
			        model.array1 = model.array1, method = "BFGS")
			pars.new <- t(sapply(o.tmp, function(lst) lst$par))
			w.new <- table.cw.ll.fun(pars.mat = pars.new, x.mat = x, p.cl = p.new, model.matrix = model.matrix)
			w.new0 <- w.new[, i0]
			w.new1 <- w.new[, i1]
			if(verbose == 2)	 cat("pars.new[1, ] =", pars.new[1,], "\n")
			mll.tmp[ii] <- multinomll.fun(p.cl = p.old, w = w.old)
			all.tmp[ii] <- sum(sapply(o.tmp, function(lst) lst$value))
			ll.tmp[ii] <- all.tmp[ii] + mll.tmp[ii]
			if(verbose == 2)	 cat("ll.tmp =", ll.tmp, "\n")
			relimp.tmp <- (ll.tmp[ii] - ll.tmp[ii-1]) / abs(ll.tmp[ii])
			# arelimp.tmp <- (all.tmp[ii] - all.tmp[ii-1]) / abs(all.tmp[ii])
			if(verbose == 2)	 cat("Rel. improvement on total likelihood:", relimp.tmp, "\n")
			# if(verbose == 2)	 cat("Rel. improvement on cALARM likelihood:", arelimp.tmp, "\n")
			clprob.tmp <- cbind(clprob.tmp, as.numeric(p.new))
			p.tmp <- array(dim = c(ii, nrow(initpars), ncol(initpars)))
			p.tmp[-ii,,] <- pars.tmp
			p.tmp[ii,,] <- pars.new
			pars.tmp <- p.tmp
			rm(p.tmp)
			p.tmp <- array(dim = c(ii, G, Nedges))
			p.tmp[-ii,,] <- membprob.tmp
			p.tmp[ii,,] <- w.new
			membprob.tmp <- p.tmp
			rm(p.tmp)
		 	if(verbose == 2)	 cat("w.new[1,1:10] =", w.new[1,1:10], "\n")
			if(verbose == 2)	 cat("w.new[2,1:10] =", w.new[2,1:10], "\n")
			w.old <- w.new
			w.old0 <- w.new0
			w.old1 <- w.new1
			p.old <- p.new
			pars.old <- pars.new
		   })[3]
		if(verbose > 0) cat(paste0("Iteration ", ii, " is done, time: ", round(s.tmp, 3), "s \n"))
		t.tmp[ii] <- s.tmp
		ii <-  ii+1
		if(verbose == 2)	 cat("relimp.tmp = ", relimp.tmp, "\n")
		# if(verbose == 2)	 cat("arelimp.tmp = ", arelimp.tmp, "\n")
		# if(verbose == 2)	 cat("tol.iter = ", tol.iter, "\n")
		if(verbose == 2)	 cat("New ii = ", ii, "\n")
		# if(verbose == 2)	 cat("max.iter = ", max.iter, "\n")
		if(verbose == 2)
			 cat("Condition (relimp.tmp > tol.iter & ii <= max.iter) =", relimp.tmp > tol.iter & ii <= max.iter, "\n")
		# if(verbose == 2)
		#	 cat("Condition (arelimp.tmp > tol.iter & ii <= max.iter) =", arelimp.tmp > tol.iter & ii <= max.iter, "\n")
		if(verbose == 2)	 cat("\n")
	   }
	list(est = pars.tmp[-(ii-1),,], cluster.prob = clprob.tmp[,-(ii-1)], membership.prob = membprob.tmp[-(ii-1),,],
	   ll = ll.tmp[-(ii-1)], multinom.ll = mll.tmp[-(ii-1)], alarm.ll = all.tmp[-(ii-1)])
   }  


## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Auxiliary functions
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## Function to get a reasonable initial value
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##    x.mat: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    model.matrix: the matrix of the covariate values 
##    K: autoregressive order of the fitted model
##    G: number of clusters for the k-means fit
## Value: a list with components...
##    prob: cluster probabilities (for initprob in EM_cBALARM.fun)
##    pars: matrix of parameters, cluster centers from the k-means algorithm;
##             dimensions: (nb of clusters) x (nb of AR parameters + 1)
##             (for initpars in EM_cBALARM.fun)

initval_arc.fun <- function(x.mat, model.matrix, K, G)
   {
   	glm.tmp <- t(apply(x.mat, 1, myglm.fun, model.matrix = model.matrix, K = K))
   	k.tmp <- kmeans(glm.tmp, centers = G)
	return(list(prob = k.tmp$size / sum(k.tmp$size), pars = k.tmp$centers))
   }


## ----------------------------------------------------------------------------------------------------------------
## Special GLM-fitting function
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##    x.mat: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    model.matrix: the matrix of the covariate values 
##    K: autoregressive order of the fitted model
## Value: a vector containing the fitted parameters, in the order of 
##    (covariate coefficients, ar1,...arK, const)

myglm.fun <- function(x.vec, model.matrix, K)
   {
   	if(sum(x.vec) == length(x.vec))
   	   {
   	   	c.tmp <- c(runif(1, min = 18, max = 20), rep(0, K + ncol(model.matrix)))
   	   }
   	if(sum(x.vec) == 0) 
   	   {
   	   	c.tmp <- c(runif(1, min = -20, max = -18), rep(0, K + ncol(model.matrix)))
   	   }
    if(sum(x.vec) < length(x.vec) & sum(x.vec) > 0)
       {
	   	N <- length(x.vec)
   	   	xmat <- sapply(0:K, function(k)  x.vec[(K+1-k):(N-k)])
		m.tmp <- as.data.frame(model.matrix)
		if(K > 0) m.tmp <- m.tmp[-(1:K),]
		m.tmp[,(ncol(m.tmp)+1):(ncol(m.tmp)+K+1)] <- xmat
		colnames(m.tmp)[(ncol(model.matrix)+1):(ncol(model.matrix)+K+1)] <- c("y", paste0("ar", 1:K))
	    g.tmp <- glm(y ~ ., data = m.tmp, family = "binomial")
	    c.tmp <- coef(g.tmp)
       }
    return(c.tmp[c(2:length(c.tmp), 1)])
   }
































## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Old (too slow but correct)
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------


## # ----------------------------------------------------------------------------------------------------------------
# ## ----------------------------------------------------------------------------------------------------------------
# ## Q-cbalarm
# ## ----------------------------------------------------------------------------------------------------------------
# ## ----------------------------------------------------------------------------------------------------------------

# ## The function can be optimized independently for each group.
# ## Arguments:
# ##   pars: vector of autoregressive coefficients and constant for the same group as w, in this order (to optimize
# ##       function for); length: (autoregressive order + 1)
# ##   x.mat: observations of the state of edges over time; dimensions: (nb of potential edges) x (nb of times) 
# ##   model.matrix: the matrix of the covariate values for the regression part. Dimensions: 
# ##      (number of time points) x (number of covariates)
# ##   w: weights for a single group; length: (nb of potential edges)
# singlecluster.cbalarmQ.fun <- function(pars, x.mat, model.matrix, w)
   # {
   	# # x.mat <- xx.tmp
	# # w <- table.w.fun( pars.mat = pars.mat, p.cl = p.cl, x.mat = xx.tmp)[1,]
	# # pars <- pars.mat[1,]
	# ## Compute the log-likelihood of each time series:
	# ll.tmp <- t(apply(x.mat, MAR = 1, singlets.cbalarmll.fun, pars = pars, model.matrix = model.matrix) )
    # sum(ll.tmp*w)
   # }


# ## ----------------------------------------------------------------------------------------------------------------
# ## Single-ts loglikelihood for an ALARM with covariates
# ## ----------------------------------------------------------------------------------------------------------------

# ## Arguments:
# ##   pars: time series parameters  for a single group (order: c(time covariate pars, AR pars, constant))
# ##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
# ##   model.matrix: the matrix of the covariate values for the regression part. Dimensions: 
# ##      (number of time points) x (number of covariates)
# singlets.cbalarmll.fun <- function(pars, x, model.matrix)
   # {
   	# D <- ncol(model.matrix)
   	# K <- length(pars) - D - 1   	
    # sum(sapply((K+1):length(x), singledata.cbalarmll.fun, pars = pars, model.matrix = model.matrix, x = x))
   # }


# ## ----------------------------------------------------------------------------------------------------------------
# ## Single-data loglikelihood for an ALARM with covariates
# ## ----------------------------------------------------------------------------------------------------------------

# ## Arguments:
# ##   pars: time series parameters  for a single group (necessary order: c(regr. pars, AR pars, constant))
# ##   x: the time series of a single edge (necessary for even a single time, because of autoregressivity) 
# ##   model.matrix: the matrix of the covariate values for the regression part. Dimensions: 
# ##      (number of time points) x (number of covariates)
# singledata.cbalarmll.fun <- function(pars, x, model.matrix, time)
   # {
   	# D <- ncol(model.matrix)
   	# K <- length(pars) - D - 1   	
   	# if(time < K+1) stop("Missing data in autoregressive sum for this time value (time < K+1)")
   	# ## Regression part:
   	# a.tmp <- pars[1:D]
   	# ## Autoregressive part:
	# b.tmp <- pars[-c(1:D, length(pars))]
	# ## Constant:
	# c.tmp <- pars[length(pars)]
	# eta <- sum(model.matrix[time, ]*a.tmp) + sum(x[(time-1):(time-K)]*b.tmp) + c.tmp
	# return(x[time]*eta - log(1+exp(eta)))
   # }


## ----------------------------------------------------------------------------------------------------------------
## EM algorithm for an ALARM with covariates
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##    initprob: initial guess for cluster probabilities (weight of mixture components)
##    initpars: initial guess for the matrix of parameters; dimensions: (nb of clusters) x (nb of AR parameters + 1)
##    x: data matrix, dimensions: (nb of possible edges) x (nb of observation times)
##    model.matrix: the matrix of the covariate values for the regression part
##    verbose: if 0: no output; 
##             if 1: only writes out the CPU spent at the end of each iteration;
##             if 2: outputs intermediate results for debugging
##    tol.iter: minimum relative improvement on the log-likelihood required for continuing the iteration
##    max.iter: the maximum number of iterations allowed
# EM_cBALARM.fun <- function(initprob, initpars, model.matrix, x, verbose = 1, tol.iter = 1e-6, max.iter = 250)	
   # {
	# G <- length(initprob)
	# D <- ncol(model.matrix)
	# K <- ncol(initpars) - D - 1
	# Nedges <- nrow(x)
	# ## Starting values stored:
	# p.cl <- initprob
	# if(verbose == 2)	cat("p.cl =", p.cl, "\n")
	# if(verbose == 2)	cat("Initial pars =", initpars, "\n")

	# ## Create storage for results during iteration:
	# mll.tmp <- -5e8
	# all.tmp <- -5e8
	# ll.tmp <- -1e9
	# clprob.tmp <- matrix(p.cl, ncol = 1, nrow = G)
	# membprob.tmp <- array(NA, dim = c(1, G, Nedges))
	# pars.tmp <- array(initpars, dim = c(1, nrow(initpars), ncol(initpars)))
	# t.tmp <- NA
	
	# ## Start iteration
	# p.old <- p.cl
	# if(verbose == 2)	cat("p.old =", p.old, "\n")
	# pars.old <- initpars
	# if(verbose == 2)	 cat("pars.old =", pars.old, "\n")
	# w.old <- table.cw.ll.fun(pars.mat = pars.old, x.mat = x, p.cl = p.old, model.matrix = model.matrix)
	# membprob.tmp[1,,] <- w.old
	# if(verbose == 2)	cat("w.old[1,1:10] =", w.old[1,1:10], "\n")
	# if(verbose == 2)	cat("w.old[2,1:10] =", w.old[2,1:10], "\n")
	# ii <- 2
	# relimp.tmp <- 1e10
	# arelimp.tmp <- 1e10
	
	# ## Iteration
	# while(relimp.tmp > tol.iter & ii <= max.iter)
	   # {
		# s.tmp <- system.time({
			# p.new <- Qmultinom.optimum.fun(w = w.old)
			# if(verbose == 2)	 cat("p.new =", p.new, "\n")
			# o.tmp <- list()
			# for(gg in 1:G)
			    # o.tmp[[gg]] <- optim(par = pars.old[gg,], fn = singlecluster.cbalarmQ.fun, hessian = F,
			        # control = list(fnscale = -1, trace = 0), x.mat = x, w = w.old[gg,], model.matrix = model.matrix)
			# pars.new <- t(sapply(o.tmp, function(lst) lst$par))
			# w.new <- table.cw.ll.fun(pars.mat = pars.new, x.mat = x, p.cl = p.new, model.matrix = model.matrix)
			# if(verbose == 2)	 cat("pars.new =", pars.new, "\n")
			# mll.tmp[ii] <- multinomll.fun(p.cl = p.old, w = w.old)
			# all.tmp[ii] <- sum(sapply(o.tmp, function(lst) lst$value))
			# if(verbose == 2)	 cat("all.tmp =", all.tmp, "\n")
			# ll.tmp[ii] <- all.tmp[ii] + mll.tmp[ii]
			# if(verbose == 2)	 cat("ll.tmp =", ll.tmp, "\n")
			# relimp.tmp <- (ll.tmp[ii] - ll.tmp[ii-1]) / abs(ll.tmp[ii])
			# arelimp.tmp <- (all.tmp[ii] - all.tmp[ii-1]) / abs(all.tmp[ii])
			# if(verbose == 2)	 cat("Rel. improvement on total likelihood:", relimp.tmp, "\n")
			# # if(verbose == 2)	 cat("Rel. improvement on cALARM likelihood:", arelimp.tmp, "\n")
			# clprob.tmp <- cbind(clprob.tmp, as.numeric(p.new))
			# p.tmp <- array(dim = c(ii, nrow(initpars), ncol(initpars)))
			# p.tmp[-ii,,] <- pars.tmp
			# p.tmp[ii,,] <- pars.new
			# pars.tmp <- p.tmp
			# rm(p.tmp)
			# p.tmp <- array(dim = c(ii, G, Nedges))
			# p.tmp[-ii,,] <- membprob.tmp
			# p.tmp[ii,,] <- w.new
			# membprob.tmp <- p.tmp
			# rm(p.tmp)
		 	# if(verbose == 2)	 cat("w.new[1,1:10] =", w.new[1,1:10], "\n")
			# if(verbose == 2)	 cat("w.new[2,1:10] =", w.new[2,1:10], "\n")
			# w.old <- w.new
			# p.old <- p.new
			# pars.old <- pars.new
		   # })[3]
		# if(verbose > 0) cat(paste0("Iteration ", ii, " is done, time: ", round(s.tmp, 3), "s \n"))
		# t.tmp[ii] <- s.tmp
		# ii <-  ii+1
		# if(verbose == 2)	 cat("relimp.tmp = ", relimp.tmp, "\n")
		# # if(verbose == 2)	 cat("arelimp.tmp = ", arelimp.tmp, "\n")
		# if(verbose == 2)	 cat("tol.iter = ", tol.iter, "\n")
		# if(verbose == 2)	 cat("New ii = ", ii, "\n")
		# if(verbose == 2)	 cat("max.iter = ", max.iter, "\n")
		# if(verbose == 2)
			 # cat("Condition (relimp.tmp > tol.iter & ii <= max.iter) =", relimp.tmp > tol.iter & ii <= max.iter, "\n")
		# # if(verbose == 2)
		# #	 cat("Condition (arelimp.tmp > tol.iter & ii <= max.iter) =", arelimp.tmp > tol.iter & ii <= max.iter, "\n")
		# if(verbose == 2)	 cat("\n")
	   # }
	# list(est = pars.tmp[-(ii-1),,], cluster.prob = clprob.tmp[,-(ii-1)], membership.prob = membprob.tmp[-(ii-1),,],
	   # ll = ll.tmp[-(ii-1)], multinom.ll = mll.tmp[-(ii-1)], alarm.ll = all.tmp[-(ii-1)])
   # }  













