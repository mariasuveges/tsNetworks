

## A couple more nice functions TO DO:
##  - covariance/correlation plot similar to the marginal probability plot



## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Simulations
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## ALARM simulation (Agaskar & Lu 2013)
## ----------------------------------------------------------------------------------------------------------------

## function to simulate from the ALARM model:
## Arguments: 
##   n: the required time series length.
##   model: the terms of the linear model producing the Bernoulli parameter, in list form; 
##       - component named ar contains the coeffs of the linear predictor (the coefficients
##       of the past Bernoulli values of the time series);
##       - component named const contains a single value, the constant of the linear predictor; 
##       - component named prob.start contains the Bernoulli parameter for the generation of 
##         starting values;
##   burn.in: the length of the burn-in period.
   
ALARM.fun <- function(n, model, burn.in = 5*length(model$ar), seed = NULL)
   {
   	if(!is.null(seed)) set.seed(seed)
   	p <- length(model$ar)
   	N <- n + burn.in
   	x <- numeric(N)
   	x[1:p] <- rbinom(n = p, size = 1, prob = model$prob.start)
   	for(ii in (p+1):N)
   	   {
   	   	eta <- model$ar*x[(ii-1):(ii-p)] + model$const
   	   	pr <- exp(eta) / (1 + exp(eta))
   	   	x[ii] <- rbinom(1, 1, prob = pr)
   	   }
    return(x[-(1:burn.in)])
   }


## ----------------------------------------------------------------------------------------------------------------
## Thresholding ARMA processes
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##    n: length of the generated time series
##    model: model parameters for the generating ARMA process, as arima.sim requires it.
##    rand.gen: the distribution of the noise in the ARMA
##    thr: threshold 
## Value: the time series x_t = I(y_t > thr), where y_t is the simulated ARMA model

thrARMA.fun <- function(n, model, rand.gen = rnorm, thr = 0, seed = NULL,  ...)
   {
   	if(!is.null(seed)) set.seed(seed)
   	pmax(0, sign(arima.sim(n = n, model = model, rand.gen = rand.gen, ...) - thr))
   }



## ----------------------------------------------------------------------------------------------------------------
## BAR (Katselis et al. 2016)
## ----------------------------------------------------------------------------------------------------------------

## function to simulate from the BAR model:
## Arguments: 
##   n: the required time series length.
##   model: the terms of the linear model producing the Bernoulli parameter, in list form; 
##       - component named ar contains the autoregressive parameters (the coefficients of the
##       past Bernoulli values of the time series);
##       - component named prob contains the Bernoulli parameter of the innovation process;
##       - component named ass contains the type of association: "pos" for positive, "neg" for
##       negative. 
##   burn.in: the length of the burn-in period.
##   
BAR.fun <- function(n, model, burn.in = 5*length(model$ar), seed = NULL)
   {
   	if(!is.null(seed)) set.seed(seed)
   	p <- length(model$ar)
   	ma <- 1 - sum(model$ar)
   	N <- n + burn.in
   	w <- rbinom(n = N, size = 1, prob = model$prob)
   	x <- numeric(N)
   	x[1:p] <- rbinom(n = p, size = 1, prob = model$prob)
   	if(model$ass == "pos")
   	   {
   	    for(ii in (p+1):N)
   	   	   {
   	   	   	x[ii] <- rbinom(1, 1, prob = model$ar*x[(ii-1):(ii-p)] + ma*w[ii])
   	   	   }
   	   }
   	if(model$ass == "neg")
   	   {
   	    for(ii in (p+1):N)
   	   	   {
   	   	   	x[ii] <- rbinom(1, 1, prob = model$ar*(1 - x[(ii-1):(ii-p)]) + ma*w[ii])
   	   	   }
   	   }
    return(x[-(1:burn.in)])
   }


## ----------------------------------------------------------------------------------------------------------------
## cALARM (covariate-ALARM) processes
## ----------------------------------------------------------------------------------------------------------------

## function to simulate from an ALARM model with functions of time as covariates:
## Arguments: 
##   n: the required time series length.
##   model: the terms of the linear model producing the Bernoulli parameter, in list form; 
##       - component named ar contains the autoregressive parameters (the coefficients of the
##       past Bernoulli values of the time series);
##       - component named prob.start contains the Bernoulli parameter of the innovation process;
##       - component named const contains the value of the constant in the linear predictor;
##       - component named regrcoefs contains the values of the regression parameters, in the order
##         in which the column of model.matrix contains the covariate values for each observation.
##   model.matrix: the matrix of the covariate values for the regression part. Dimensions: 
##       (number of time points) x (number of covariates). The order of the covariates must be the
##       same as in the parameter vector model$regrcoefs.
##   burn.in: the length of the burn-in period.
##   seed: an optional random seed value.
##   x.only: only the simulated time series is needed, or the variation of the edge probs over time
##       too?
cALARM.fun <- function(n, model, model.matrix, burn.in = 5*length(model$ar), seed = NULL,
    x.only = TRUE)
   {
   	if(!is.null(seed)) set.seed(seed)
   	p <- length(model$ar)
   	N <- n + burn.in
   	x <- numeric(N)
   	## For the first p values in the extended ts, we don't even have all the former values
   	## for the autoregr. contribution, and don't have the covariates:
   	x[1:p] <- rbinom(n = p, size = 1, prob = model$prob.start)
   	## For the burn-in period, we don't necessarily have the covariate values:
   	for(ii in (p+1):burn.in)
   	   {
   	   	eta <- model$ar*x[(ii-1):(ii-p)] + model$const
   	   	pr <- exp(eta) / (1 + exp(eta))
   	   	x[ii] <- rbinom(1, 1, prob = pr)
   	   }
   	## Then for the rest of time, we have everything we need:
   	covar.term <- model.matrix %*% t(t(model$regrcoefs))
   	probs <- numeric(n)
   	for(ii in (burn.in+1):N)
   	   {
   	   	eta <- model$ar*x[(ii-1):(ii-p)] + covar.term[ii-burn.in,] + model$const
   	   	pr <- exp(eta) / (1 + exp(eta))
   	   	probs[ii-burn.in] <- pr
   	   	x[ii] <- rbinom(1, 1, prob = pr)
   	   }
    if(x.only) {return(x[-(1:burn.in)])} else {return(list(x = x[-(1:burn.in)], probs = probs))}
   }



## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Diagnostic functions
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## Marginal probabilities in adjacency matrix form
## ----------------------------------------------------------------------------------------------------------------

## To plot the marginal probabilities of each time series for dynamic network data.
## Arguments:
##    mp.mat: the matrix of the marginal probabilities.
##    filename: character string, the filename (.png extension is necessary) if a file should be created.
##         Must be NULL (the default) if a plotting device should be opened in R.
##    breaks: the breakpoints for the grayscale.
## Output: either a figure in R or a .png stored, a greyscale image of the marginal probabilities. A
##    colour key is added to the right border.

mp.plot.fun <- function(mp.mat, filename = NULL, breaks = seq(min(mp.mat), max(mp.mat), length.out = 51))
   {
   	nnodes <- nrow(mp.mat)
   	if(is.null(filename)) 
   	   {
		quartz(height = 7.5, width = 11/10*7.5)
		par(oma = c(0,0,2,1.5))
		layout(matrix(c(1,2), ncol = 2), widths = c(10, 1))
		par(mar = c(3.5,3.5,0.1,0.1), mgp = c(1.6,0.6,0))
		image(mp.mat, axes = F, xlab = "", ylab = "", breaks = breaks,  
		   col = colorRampPalette(c("grey20","white"), space = "rgb")(length(breaks)-1))
		axis(1, labels = rownames(mp.mat), at = seq(0,1, length.out = nnodes), las = 2, tcl = 0, cex.axis = 0.4)
		axis(2, labels = colnames(mp.mat), at = seq(0,1, length.out = nnodes), las = 2, tcl = 0, cex.axis = 0.4)
		box()
		par(mar = c(3.5,0.5,0.1,2), mgp = c(1.6,0.6,0))
		image(1, breaks, t(matrix(breaks, ncol = 1)), breaks = breaks, axes = F, xlab = "", ylab = "",
		   col = colorRampPalette(c("grey20","white"), space = "rgb")(length(breaks)-1))
		box()
		axis(4, las = 1, at = pretty(breaks))
		mtext("Marginal probability of edge", side = 4, outer = T, line = -0.5)
       } else {
        png(filename = filename, width = 11/10 * 1200, height = 1200) 
		par(oma = c(0,0,2,1.5))
		layout(matrix(c(1,2), ncol = 2), widths = c(10, 1))
		par(mar = c(3.5,3.5,0.1,0.1), mgp = c(1.6,0.6,0))
		image(mp.mat, axes = F, xlab = "", ylab = "", breaks = breaks,  
		   col = colorRampPalette(c("grey20","white"), space = "rgb")(length(breaks)-1))
		axis(1, labels = rownames(mp.mat), at = seq(0,1, length.out = nnodes), las = 2, tcl = 0, cex.axis = 0.9)
		axis(2, labels = colnames(mp.mat), at = seq(0,1, length.out = nnodes), las = 2, tcl = 0, cex.axis = 0.9)
		box()
		par(mar = c(3.5,1,0.1,4), mgp = c(1.6,0.6,0))
		image(1, breaks, t(matrix(breaks, ncol = 1)), breaks = breaks, axes = F, xlab = "", ylab = "",
		   col = colorRampPalette(c("grey20","white"), space = "rgb")(length(breaks)-1))
		box()
		axis(4, las = 1, at = pretty(breaks))
		mtext("Marginal probability of edge", side = 4, outer = T, line = -0.5)
		dev.off()
       }
   }


## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------
## Geometric qq plots
## ----------------------------------------------------------------------------------------------------------------
## ----------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------
## Computation of run lengths
## ----------------------------------------------------------------------------------------------------------------

## Function to compute the run lengths of both 0s and 1s in a binary time series.
## Arguments:
##   ts: the time series, composed of ones and zeros
## Value: a list with two components named 'ones' and 'zeros', containing the run lengths of the ones and
##        the zeroes, resp.

runlengths.fun <- function(ts)
   {
   	#mpr.tmp <- sum(ts) / length(ts)
   	if(sum(ts) == length(ts) | sum(1-ts) == length(ts))
   	   {
   	   	if(sum(ts) == length(ts)) {rl.tmp <- list(ones = length(ts), zeros = NA)}
   	   	if(sum(1-ts) == length(ts)) {rl.tmp <- list(ones = NA, zeros = length(ts))}
   	   } else {
		v.tmp <- strsplit(paste0(ts, collapse = ""), split = "0")[[1]]
		ones.tmp <- as.numeric(sapply(v.tmp, nchar))
		v.tmp <- strsplit(paste0(ts, collapse = ""), split = "1")[[1]]
		zeros.tmp <- as.numeric(sapply(v.tmp, nchar))
		rl.tmp <- list(ones = ones.tmp, zeros = zeros.tmp)	
	   }	
	return(rl.tmp)    
   }

## ----------------------------------------------------------------------------------------------------------------
## Computation of a pointwise CI band for a geometric qq plot
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##   ndata: number of data points, which should be checked against the geometric distribution
##   prob: the probability of the geometric distribution
##   N: number of repetitions 
##   alpha: confidence level

cigeom.fun <- function(ndata, prob, N = 10000, alpha = 0.05)
   {
   	d.tmp <- replicate(N, sort(rgeom(n = ndata, prob = prob)))
   	p.tmp <- qgeom(ppoints(ndata), prob = prob)
   	d.tmp1 <- split(d.tmp, p.tmp)
   	sapply(d.tmp1, quantile, probs = c(alpha/2, 1-alpha/2))
   }


## ----------------------------------------------------------------------------------------------------------------
## Geometric qq plot
## ----------------------------------------------------------------------------------------------------------------

## Arguments:
##   rl: list with two components, named 'ones' and 'zeros' (output of runlengths.fun)  
##   mp: the probability of the geometric distribution
##   nname1, nname2: names of the nodes to add to the plot (to specify the edge figuring on the plot)
##   which.rl: which component of the list to use, 'ones' or 'zeros'? Possible values: 'ones', 'zeros' or 'both'.
##      Defaults to 'both'

geomQQplot.fun <- function(rl, mp, nname1, nname2, which.rl = "both")
   {
   	n.tmp0 <- length(rl$ones)
	n.tmp1 <- length(rl$zeros)
	q.tmp0 <- qgeom(ppoints(n.tmp0), prob = 1-mp)
	q.tmp1 <- qgeom(ppoints(n.tmp1), prob = mp)
	if(mp == 0 | mp == 1)
	   {
	   	plot(1,1, type = "n", xlab = "", ylab = "", axes = F)
	   } else {
		if(which.rl == "both") 
		   {
			ci.tmp0 <- cigeom.fun(ndata = n.tmp0, prob = mean(1-mp, na.rm = T))
			ci.tmp1 <- cigeom.fun(ndata = n.tmp1, prob = mean(mp, na.rm = T))
		   	r.tmp <- range(unlist(rl))
		   } else {
		   	if(which.rl == "ones") 
		   	   {
				ci.tmp0 <- cigeom.fun(ndata = n.tmp0, prob = mean(1-mp, na.rm = T))
		   	   	r.tmp <- range(rl$ones)
		   	   }
		   	if(which.rl == "zeros") 
		   	   {
				ci.tmp1 <- cigeom.fun(ndata = n.tmp1, prob = mean(mp, na.rm = T))
		   	   	r.tmp <- range(rl$zeros)
		   	   }
		   }
		plot(q.tmp0, sort(rl$ones),
		   type = "n", xlab = "", ylab = "", axes = F, xlim = range(q.tmp0, q.tmp1),
		   ylim = r.tmp)
		if(which.rl == "both" | which.rl == "ones")
		   {
			points(q.tmp0, sort(rl$ones), pch = 15, cex = 1, col = "violetred")
			lines(as.numeric(colnames(ci.tmp0)), ci.tmp0[1,], col = "violetred")
			lines(as.numeric(colnames(ci.tmp0)), ci.tmp0[2,], col = "violetred")
		   }
		if(which.rl == "both" | which.rl == "ones")
		   {
			points(q.tmp1, sort(rl$zeros), pch = 17, cex = 0.9, col = "deepskyblue")
			lines(as.numeric(colnames(ci.tmp1)), ci.tmp1[1,], col = "deepskyblue")
			lines(as.numeric(colnames(ci.tmp1)), ci.tmp1[2,], col = "deepskyblue")
		   }
		abline(c(0,1))
		box()
		text(1.1*min(c(q.tmp0, q.tmp1)), 0.95*max(unlist(rl)), 
		   labels = paste0(nname1, "--", nname2), pos = 4)
	   }
   }




















