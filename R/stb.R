#' Simultaneous Tolerance Bands (STB).
#' 
#' Compute And/Or Plot Simultaneous Tolerance Bands for numeric vectors.
#' 
#' Function takes a numeric vector 'vec' and computes the 100(1-alpha)\%-simultaneous tolerance band (STB) for
#' the (DEFAULT )Null-hypothesis H0: vec~N(my, sigma^2) distributed, which is equal to checking whether the residuals of 
#' the simplest linear model y = mu + e (y~1) are normally distributed, i.e. 'e ~ N(0, sigma^2)'. By specification of argument \code{rand.func} 
#' other null-distributions can be specified. One has to specify a function with a single argument 'n',  which returns a random sample with 'n' 
#' observations, randomly drawn from the desired null-distribution (see description argument \code{rand.func} below). 
#' Note that all random samples as well as vector \code{vec} will be centered to mean=0 and scaled to sd=1.
#' 
#' One can choose between three methods for constructing the 100(1-alpha)\% STB. There are two implementations of the quantile-based algorithm
#' ("C", "R" see 1st reference) and one of the rank-based algorithm (see 2nd reference). Methods "C" and "R" can be run in parallel. The rank-based
#' algorithm does not benefit form parallel processing, at least not in the current implementation. It is still the default and recommended for small to medium 
#' sized vectors and 10000 <= N <= 20000 simulations, which should be sufficiently accurate reflect the null-distribution. The "C" and "R" options refer
#' to implementations of the quantile-based algorithm. The "C" implementation benefits most from using multiple cores, i.e. the larger 'Ncpu' the better,
#' and should be used for large problems, i.e. rather large number of elements and large number of simulations.
#' 
#' The table below gives an impression how these algorithms perform. Runtimes were measured under Windows 7 on a Intel Xeon E5-2687W 3.1 GHz workstation with 16
#' logical cores and 16 GB RAM. The runtime of the C-implementation of the quantile-based algorithm is denoted as "t_qC12" applied parallely with 12 cores.
#' Each setting was repeated 10 times and the overall run time was then divided by 10 providing sufficiently robust simulation results.
#' Note, that for smaller problem sizes a large proportion of the overall runtime is due to simulating, i.e. drawing from the null-distribution.
#' 
#' \tabular{rrrr}{
#' 	\strong{_____N_obs} \tab \strong{_____N_sim} \tab \strong{____t_rank} \tab \strong{____t_qC12} \cr
#'  25    \tab 5000  \tab  0.4s \tab   0.5s \cr
#'  25	  \tab 10000 \tab  0.8s \tab   1.3s \cr
#'  50	  \tab 10000 \tab  1.0s \tab   3.2s \cr 
#'  100	  \tab 10000 \tab  1.7s \tab   2.9s \cr
#'  100	  \tab 20000 \tab  3.0s \tab   4.8s \cr
#'  225	  \tab 20000 \tab  5.1s \tab   8.3s \cr
#'  300	  \tab 30000 \tab  9.6s \tab  17.2s \cr
#'  300	  \tab 50000 \tab 16.1s \tab  24.9s \cr
#'  1000  \tab 50000 \tab 47.8s \tab 123.5s \cr
#' }                                      
#'
#' @param obj 			(numeric) vector, which is supposed to be N(my, sigma^2)-distributed
#' @param N 			(integer) value specifying the number of random samples to be used for constructing the STB
#' @param alpha 		(numeric) value specifying the simultaneous tolerance level, i.e. 100(1-alpha)\% of all 'N'
#'              		random samples have to be completely enclosed by the bounds of the STB
#' @param rand.func 	(function) a function which generates random samples, e.g. \code{rand.func=rnorm} which corresponds to random
#'                  	sampling from the standard normal distribution. Another example is defining func=function(n){rchisq(n=n, df=3, ncp=2)}
#'                  	and using \code{rand.func=func}. See examples for further examples.
#' @param tol 			(numeric) value specifying the max. acceptable deviation from 'alpha' used in the bisection algorithm
#' @param max.iter 		(integer) value specifying the max. number of iteration for finding the bounds of the STB
#' @param algo			(character) (string) specifying the method to be used for constructing a 100(1-alpha)\% STB,
#' 						choose "rank" for the rank-based, "C" for a C-implementation of the quantile-based, and
#' 						"R" for an R-implentation of the quantile-based algorithm (see details). "C" uses SAS PCTLDEF5 definition of quantiles,
#' 						whereas "R" can use any of the built-in R types of quantiles (see \code{\link{quantile}}.
#' @param Ncpu			(integer) specifying the number cores/CPUs to be used, for N>1 multi-processing is applied 
#' @param q.type 		(integer) the quantile-type used if \code{algo="R"}, see ? quantile for details.
#' @param stb.col 		(character) string, a valid specification of a color to be used for the STB
#' @param main 			(characer) string for a main title appearing over the plot
#' @param add.pwb 		(logical) should the point-wise tolerance band be plotted for comparison?
#' @param quiet 		(logical) TRUE = no additional output ist printed (progress bar etc.)
#' @param add 			(logical) TRUE = the 100(1-alpha)\% STB is added to an existing plot
#' @param col.points 	(character) color for the points in the QQ-plot
#' @param col.out 		(character) color for points outsied of the 100(1-alpha)\% STB
#' @param col.pwb 		(character) color for the point-wise STB (not adjusted for multiplicity), defaults to "#0000FF40" which is "blue"
#'                		with 80\% transparency
#' @param plot 			(logical) TRUE = either a QQ-plot with STB (add=FALSE) or a STB added to an existing plot (add=TRUE) is plotted.
#'                      FALSE = only computations are carried out without plotting, an object of class 'STB' is returned which
#'                              can be stored an plotted later on, e.g. to avoid computing an STB every time a Sweave/mWeave report is updated
#' @param legend 		(logical) TRUE a legend is plotted "topleft"
#' @param timer 		(logical) TRUE = the time spent for computing the STB will be printed
#' @param pch 			(integer) plotting symbols for the QQ-plot
#' @param pch.out 		(integer) plotting symbols for outlying points
#' @param seed 			(numeric) value interpreted as integer, setting the random number generator (RNG) to a defined state
#' @param ... further graphical parameters passed on
#' 
#' @return invisibly returns a list-type object of class \code{STB}, which comprises all arguments accepted by this function.
#'  
#' @author  Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \link{plot.STB}
#' 
#' @S3method stb default
#'
#' @references 
#' 
#' Schuetzenmeister, A., Jensen, U., Piepho, H.P. (2011), Checking assumptions of normality and homoscedasticity in the general linear model. 
#'             Communications in Statistics - Simulation and Computation; S. 141-154
#' 
#' Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
#' Computational Statistics and Data Analysis, 56, 1405-1416
#' 
#' @examples 
#' ### log-normal vector to be checked for normality
#' \dontrun{
#' set.seed(111)
#' stb(exp(rnorm(30)), col.out="red", legend=TRUE)
#' 
#' ### uniformly distributed sample checked for Chi-Squared Distribution with DF=1, degrees of freedom
#' set.seed(707)
#' stb(runif(25, -5, 5), rand.func=function(n){rchisq(n=n, df=1)}, 
#'     col.out="red", legend=TRUE, main="Chi-Squared with DF=1")
#' 
#' ### check whether an Chi-Squared (DF=1) random sample better fits 
#' stb(rchisq(25, df=1), rand.func=function(n){rchisq(n=n, df=1)}, 
#'     col.out="red", legend=TRUE, main="Chi-Squared with DF=1")
#' 
#' ### add STB to an existing plot
#' plot(sort(rnorm(30)), sort(rnorm(30)))
#' stb(rnorm(30), add=TRUE)
#' 
#' ### compute STB for later use and prevent plotting
#' STB <- stb(rnorm(30), plot=FALSE)
#' }

stb.default <- function(obj, N=10000L, alpha=.05, rand.func=rnorm, tol=.0001, max.iter=100L, algo=c("rank", "C", "R"), 
				Ncpu=1, q.type=2L, stb.col="#0000FF40", col.points="black", col.out="red", col.pwb="#0000FF40",
                main=NULL, add.pwb=FALSE, quiet=FALSE, add=FALSE, plot=TRUE, legend=FALSE, timer=FALSE, pch=16, pch.out=16, seed=NULL, ... )
{
    ARGS <- list(...)
    stopifnot(alpha>0 && alpha<1)
    algo <- toupper(algo)[1]
    algo <- match.arg(algo, c("RANK", "C", "R"))	
	stopifnot(q.type %in% 1:9)
	vec <- obj
    if(any(is.na(vec)))
        vec <- na.omit(vec)
    vec <- c(scale(vec, center=TRUE, scale=TRUE))                                       # center and scale 'vec' to have mean=0, and sd=1
    vec <- sort(vec)
    Nobs <- length(vec)   
	
	if(is.null(seed))
		seed <- runif(1, -2^30, 2^30)
	
	set.seed(seed)
    
	Sample <- eval(call(name="rand.func", n=N*Nobs))                                    # random sampling according to 'rand.func'
		
	if(!isTRUE(all.equal(rand.func, rnorm)))
		fun <- function(x){sort(scale(x, scale=TRUE, center=TRUE))}                     # mean-centering, scaling to sd=1 and sorting if 'rand.func' is not N(0,1)
	else
		fun <- sort                                                                     # only sorting if N(0,1)    
	
	mc.mat <- t(apply(matrix(Sample, ncol=Nobs), 1, fun))                               # 'N' random samples from N(0,1)
    
    if(algo=="C")
	{
		stb <- fastSTB(mc.mat, alpha=alpha, tol=tol, max.iter=max.iter, timer=timer, Ncpu=Ncpu )
	}
    else if(algo=="R")
	{
		stb <- getSTB(	mc.mat, alpha=alpha, output=quiet, tol=tol, max.iter=max.iter,   # compute bounds of the STB
						q.type=q.type, timer=timer, Ncpu=Ncpu)#, q.Rtype=q.Rtype)
	}
	else
	{
		stb <- rankSTB(mc.mat, alpha=alpha)
	}
	
    means <- apply(mc.mat, 2, mean)                                                     # compute means for each order-statistic, i.e. X-values for the QQ-plot
    Q <- stb$Q                                                                          # lower and upper bound of the simultaneous tolerance band (STB)
    
    # extend STB beyond smallest and largest observation for a nicer look
    
    ll <- -10 * ( (Q[1,2]-Q[1,1])/ (means[2]-means[1]) )                                # Y-value for x=-10 corresponding to the left-lower bound
    lu <- -10 * ( (Q[2,2]-Q[2,1])/ (means[2]-means[1]) )                                # left-upper
    rl <-  10 * ( (Q[1,Nobs]-Q[1,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )             # right-lower
    ru <-  10 * ( (Q[2,Nobs]-Q[2,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )             # right-upper
    
    if(add.pwb)
    {
        pw.Q <- apply(mc.mat, 2, quantile, probs=c(alpha/2, 1-alpha/2))
        ll2 <- -10 * ( (pw.Q[1,2]-pw.Q[1,1])/ (means[2]-means[1]) )                     # Y-value for x=-10 corresponding to the left-lower bound
        lu2 <- -10 * ( (pw.Q[2,2]-pw.Q[2,1])/ (means[2]-means[1]) )                     # left-upper
        rl2 <- 10 * ( (pw.Q[1,Nobs]-pw.Q[1,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )  # right-lower
        ru2 <- 10 * ( (pw.Q[2,Nobs]-pw.Q[2,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )  # right-upper
        pw.Q <- cbind( c(ll2, lu2), pw.Q, c(rl2, ru2) )
    }
    if(!add && plot)
        plot(means, means, ylim=c(min(Q[1,]), max(Q[2,])), type="n", main=main, xlab="Theoretical Quantiles", ylab="Sample Quantiles" )
    
    Q <- cbind(c(ll,lu) ,Q)
    Q <- cbind(Q, c(rl,ru))
    
    means <- c(-10, means, 10)
        
    if(plot)
    {       
        polygon(x=c(means, rev(means)), y=c(Q[1,],rev(Q[2,])), border=NA, col=stb.col )
        if(!add)
        {
            abline(0,1)
            box()
            points( means[-c(1,length(means))], vec, pch=pch, col=col.points, ... ) 
            out.high <- which(vec > Q[2,-c(1,ncol(Q))])
            out.low <- which(vec < Q[1, -c(1,ncol(Q))])
            out <- unique(c(out.high, out.low))
            if(length(out > 0))
                points(means[out+1], vec[out], pch=pch.out, col=col.out)     # means[1] and means[N] have no counterpart in vec
            if(legend)
                legend("topleft", fill=stb.col, legend=paste("STB (", 100*round(stb$cov,5),"% coverage)",sep=""),
                        paste(N," Simulations", sep=""), box.lty=0, border=NA)
        }
    }
    
    if(add.pwb && plot)
    {
        polygon(x=c(means, rev(means)), y=c(pw.Q[1,],rev(pw.Q[2,])), border=NA, col=col.pwb )
        if(legend)
            legend("topleft", fill=c(stb.col, rgb(t(col2rgb("blue")), maxColorValue=255, alpha=120), NA), legend=c(paste("STB (", 100*round(stb$cov,5),"% coverage)",sep=""),
                            paste("PWB (", 100*(1-alpha),"% coverage)",sep=""), paste(N," Simulations", sep="")), box.lty=0, border=c("black","black", NA))
    }
    res <- list(vec=vec, means=means, Q=Q, N=N, alpha=alpha, stb.col=stb.col, col.pwb=col.pwb, main=main, add.pwb=add.pwb, 
				quiet=quiet, add=add, col=col.points, col.out=col.out, coverage=stb$cov, pch=pch, pch.out=pch.out, 
                rand.func=rand.func, quantile.type=q.type, stb.algo=algo, legend=legend, ARGS=ARGS)
    class(res) <- "STB"
    invisible(res)
}


#' Plot Objects of Class 'STB'.
#' 
#' Standard plotting method for objects of class 'STB'.
#' 
#' This function plots objects of class 'STB' as generated by function \code{\link{stb}}.
#' Objects of S3-class 'STB' are list-type objects storing all the information
#' needed to plot QQ-plots with simultaneous tolerance bounds.
#' 
#' @param x     (object) of class 'STB' as generated by function \code{getSTB}
#' @param ...   arguments passed to other methods
#' 
#' @author  Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @method plot STB
#' @S3method plot STB
#' 
#' @seealso \link{stb}
#' 
#' @examples 
#' \dontrun{
#' ### generate an 'STB' object without plotting
#' obj <- stb(rnorm(30), plot=FALSE)
#' plot(obj)
#' 
#' ### manipulate the 'STB' object for plotting
#' obj$legend=TRUE
#' plot(obj)
#' 
#' ### add a previously generated STB-ocject to an existing plot
#' plot(sort(rnorm(30)), sort(rnorm(30)))
#' obj$add <- TRUE
#' plot(obj) 
#' }

plot.STB <- function(x, ...)
{
    obj <- x
    stopifnot(class(obj) == "STB")
    N <- length(obj$means)
    if(!obj$add)
        plot(obj$means[-c(1,N)], obj$means[-c(1,N)], ylim=c(min(obj$Q[1,-1]), max(obj$Q[2,-N])), 
                type="n", main=obj$main, xlab="Theoretical Quantiles", ylab="Sample Quantiles")
    
    polygon(x=c(obj$means, rev(obj$means)), y=c(obj$Q[1,], rev(obj$Q[2,])), border=NA, col=obj$stb.col )
    
    if(!obj$add)
    {
        abline(0,1)
        box()
        points( obj$means[-c(1,length(obj$means))], obj$vec, pch=obj$pch, col=obj$col ) 
        out.high <- which(obj$vec > obj$Q[2,-c(1,ncol(obj$Q))])
        out.low <- which(obj$vec < obj$Q[1, -c(1,ncol(obj$Q))])
        out <- unique(c(out.high, out.low))
        if(length(out > 0))
            points(obj$means[out+1], obj$vec[out], pch=obj$pch.out, col=obj$col.out)     # means[1] and means[N] have no counterpart in vec
        if(obj$legend)
            legend("topleft", fill=obj$stb.col, legend=paste("STB (", 100*round(obj$coverage,5),"% coverage)",sep=""),
                    paste(obj$N," Simulations", sep=""), box.lty=0, border=NA)
    }
}




#' Simultaneous Tolerance Bands Using a Fast C-Implementation.
#' 
#' This function represents the interface to a C-implementation of the bisection algorithm for computing
#' simultaneous tolerance bounds as described in Schuetzenmeister et al. 2012.
#' 
#' Quantiles will be computed according to SAS PCTLDEF5 definition of quantiles, which is identical to 'type=2' 
#' in function \code{\link{quantile}}. This is also the default-option throughout this package. Function \code{\link{SASquantile}}
#' is a R-implementation identical to the C-implementation used here.
#'
#' @param mat 			(numeric) matrix with rows representing simulations of 'ncol' values, rows are supposed to be sorted
#' @param alpha 		(numeric) 100(1-alpha)\% simultaneous tolerance level (coverage)
#' @param tol 			(numeric) convergence tolerance level for the bisection algorithm
#' @param max.iter 		(integer) number of bisection-steps for the algorithm
#' @param Ncpu			(integer) specifying the number of processors (CPUs, cores, ...) to be used
#' @param timer 		(logical) TRUE = the time spent for computing the STB will be printed
#' 
#' @return (list) with elements:
#' \item{mat}{(numeric) matrix with sorted random vector stored in rows}
#' \item{nCol}{(integer) number of columns of 'mat'}
#' \item{nRow}{(integer) number of rows of 'mat'}
#' \item{alpha}{(numeric) values used for the 100(1-alpha)\% STB}
#' \item{tol}{(numeric) the tolerance value used as stopping criterion for the bisection algorithm}
#' \item{max.iter}{(integer) the max. number of iteration of the bisection algorithm}
#' \item{Q}{(numeric) matrix with 2 rows, corresponding to the bounds of the STB (1st row = lower bounds, 2nd row = upper bounds)}
#' \item{coverage}{(numeric) value corresponding to the actual coverage of the STB}
#'
#' @author  Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \link{getSTB} \link{stb}
#'
#' @references Schuetzenmeister, A., Jensen, U., Piepho, H.P. (2011), Checking assumptions of normality and homoscedasticity in the general linear model. 
#'             Communications in Statistics - Simulation and Computation; S. 141-154
#'
#' @examples
#' \dontrun{
#' ### example for 10000 x 30 matrix
#'  set.seed(333)
#'  mat <- t(apply(matrix(rnorm(10000*30), ncol=30), 1, sort))
#'  stb.obj1 <- fastSTB(mat, timer=TRUE)
#'  stb.obj1$coverage
#' 	stb.obj2 <- fastSTB(mat, timer=TRUE, Ncpu=4)
#'  stb.obj3 <- fastSTB(mat, timer=TRUE, Ncpu=6)
#'  stb.obj4 <- fastSTB(mat, timer=TRUE, Ncpu=8)
#' }


fastSTB <- function(mat, alpha=.05, tol=0.0001, max.iter=100L, Ncpu=2, timer=FALSE)
{
  t1 <- Sys.time()
  stopifnot(alpha>0 && alpha<1)
  stopifnot(is.matrix(mat) && is.numeric(mat))
  stopifnot(is.numeric(max.iter) && max.iter > 1)
  stopifnot(is.numeric(tol) && tol > 0)
  stopifnot(is.numeric(Ncpu) && Ncpu > 0)
  stopifnot(is.logical(timer))
  
  NumCores <- detectCores()
  
  if(Ncpu <= 0)
	  Ncpu <- 1									# sequential processing
  if(Ncpu > NumCores)
  {
	  warning("There are only ", NumCores," cores available but Ncpu set to ",Ncpu,"! Use ", NumCores," cores instead!")
	  Ncpu <- NumCores				# rule of thumb for parallel-processing 
  }

  res <- .C("getSTB", mat=as.double(mat), nCol=as.integer(ncol(mat)), nRow=as.integer(nrow(mat)), alpha=as.double(alpha),
                      tol=as.double(tol), max.iter=as.integer(max.iter), Ncpu=as.integer(Ncpu),
					  Q=double(2*ncol(mat)), coverage=double(1), PACKAGE="STB")

  res$Q <- matrix(res$Q, nrow=2)
  rownames(res$Q) <- c("lower", "upper")
  res$mat <- matrix(res$mat, ncol=ncol(mat))
  if(timer)
    print(Sys.time()-t1)
  return(res)
}


#' Load/unload C-lib.

.onLoad <- function(lib, pkg)
{
	library.dynam(chname="STB", package=pkg, lib.loc=lib)
}

.onUnload <- function(lib)
{
	library.dynam.unload(chname="STB", libpath=lib)
}



#' Simultaneous Tolerance Bands Using R-Implementation.
#' 
#' Compute simultaneous tolerance bounds (STB) from a matrix.
#' 
#' Function computes 100(1-alpha)\% simultaneous tolerance bounds (intervals, bands) from a matrix, where rows correspond
#' to order statistics (sorted) of random samples of a, e.g. normal distribution. It starts by computing joint coverage 
#' of the naive unadjusted (point-wise) alpha-level  intervals, computed as percentiles across each order statistic (columns).
#' Alpha is decreased using a bisection search until the joint coverage becomes at least 100(1-alpha)\% depending on arguments
#' \code{max.iter} and \code{tol}. If the number of simulated random samples goes to infinity, the STB becomes exact.
#' Note that checking whether a numeric vector is normally distributed is equal to checking whether the residuals of the simplest
#' linear model y = mu + e (y~1) are normally distributed [e ~ N(0, sigma^2)].
#'
#' @param mat 		(numeric) matrix with nrow=N_simulation, ncol=N_points, where each row is sorted in increasing order
#' @param alpha 	(numeric) 100(1-alpha)\% simultaneous tolerance-level will be computed
#' @param tol 		(numeric) value controlling the termination of the bisection algorithm, i.e. the condition \code{coverage-(1-alpha)<=tol}
#'            		has to be TRUE.
#' @param max.iter 	(integer) maximum number of iterations of the bisection algorithm. If this number is reached the algorithm terminates
#'                 	and returns the best approximation of lower and upper bounds reached up to this iteration.
#' @param q.type 	(integer) the quantile-type used if \code{algo="R"}, see ? quantile for details.
#' @param logfile 	(character) string specifying the name of the (optional) log-file, storing the iteration history 
#' @param output 	(logical) TRUE a txtProgressBar will be printed as well as additional status information
#' @param timer 	(logical) TRUE = the time spent for computing the STB will be printed
#' @param Ncpu		(integer) specifying the number cores/CPUs to be used in computing the coverage on C-level,
#' 							  for N>1 multi-processing is applied 
#' 
#' @return (list) with elements:
#' \item{Q}{(numeric) matrix with two rows, where the 1st row = lower bounds of the STB, 2nd row = upper bounds of the STB}
#' \item{cov}{(numeric) value indicating the actual coverage of the STB}
#'
#' @author  Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \link{fastSTB} \link{stb}
#'
#' @references Schuetzenmeister, A., Jensen, U., Piepho, H.P. (2011), Checking assumptions of normality and homoscedasticity in the general linear model. 
#'             Communications in Statistics - Simulation and Computation; S. 141-154
#' 
#' @examples 
#' \dontrun{
#' set.seed(333)
#' mat <- t(apply(matrix(rnorm(10000*30), ncol=30), 1, sort))
#' stb.obj <- getSTB(mat, timer=TRUE)
#' stb.obj$cov
#' }

getSTB <- function( mat, alpha=.05, tol=.0001, max.iter=100L, q.type=2L, logfile=NULL, output=TRUE, timer=FALSE, Ncpu=2)
{
    t1 <- Sys.time()
		
	stopifnot(alpha>0 && alpha<1)
	stopifnot(is.matrix(mat) && is.numeric(mat))
	stopifnot(is.numeric(max.iter) && max.iter > 1)
	stopifnot(is.numeric(tol) && tol > 0)
	stopifnot(is.numeric(q.type) && q.type %in% 1:9)
	stopifnot(is.numeric(Ncpu) && Ncpu > 0)
	stopifnot(is.logical(timer))
	stopifnot(is.logical(output))
	
    if( !is.null(logfile) && class(logfile) == "character" && length(logfile) > 0)
        file.remove( logfile )
    check <- TRUE
    include <- 1-alpha                                                                          # requested confidence limit
    alpha.old <- alpha                                                                          # for the 1st iteration
    alpha <- alpha/2                                                                            # actually 1st alpha value to be used
    iter <- 0
    if( output ){
        cat("\n\nConstruct simultaneous tolerance bounds - Progress:\n\n")                      # progress bar
        pb <- txtProgressBar(min=0, max=max.iter, width=getOption("width")*.25, char="*", style=3)
    }

	nc <- ncol(mat)
	nr <- nrow(mat)	
    best.cov <- 1
	
    while( check )
	{
		iter <- iter + 1			
		Qs <- apply( mat, 2, quantile, probs=c(alpha,1-alpha), type=q.type )  				# point-wise confidence-interval bounds

		coverage <- .C( "getCoverage", mat=as.double(mat), lower=as.double(Qs[1,]),				# coverage computed as C-function
						upper=as.double(Qs[2,]), nCol=as.integer(ncol(mat)), 
						nRow=as.integer(nrow(mat)), nCPU=as.integer(Ncpu), 
						cov=double(1), PACKAGE="STB")$cov

		if(iter == 1)
			best.Qs <- Qs
		
		if( coverage >= include ){                                                                # at least (1-alpha)*100% coverage
            if( (coverage < best.cov) ){                                                          # new best coverage
                best.alpha <- alpha
                best.Qs <- Qs
                best.cov <- coverage
            }
        }
        if( !is.null(logfile) && class(logfile) == "character" && length(logfile) > 0){           # create or update logfile
            if( file.exists(logfile) ){                                                           # update
                line <- paste(iter, alpha, alpha.old, abs(alpha-alpha.old)/2, (coverage * 100), coverage - include, sep="\t" )
                cat( line, file=logfile, append=TRUE, sep="\n" )
            }
            else{                                                                                   # create log-file
                cat( paste("Iter", "alpha_i", "alpha_i-1", "stride", "%coverage", "Diff", sep="\t"), file=logfile, append=FALSE, sep="\n" )
                line <- paste(iter, alpha, alpha.old, abs(alpha-alpha.old)/2, (coverage * 100), coverage - include, sep="\t" )
                cat( line, file=logfile, append=TRUE, sep="\n" )
            }
        }
        if( abs(alpha-alpha.old)/2 == 0 ){                                            ########## difference in local alpha of current iteration and previous iteration is 0
            if(output){
                setTxtProgressBar(pb,max.iter)
                close(pb)
                cat(paste(best.cov*100,"%", sep=""),"of",nrow(mat),"simulations covered simultaneously.\n\n")
            }
            if(!is.null(logfile) && class(logfile) == "character")
                cat("\nAlgorithm stopped because the value abs(alpha-alpha.old)/2 == 0!\n", file=logfile, append=TRUE, sep="\n" )
            if(timer)
                print(Sys.time()-t1)
            return( list( mat=mat, Q=best.Qs, coverage=best.cov ) )
        }
        if( iter == max.iter ){                                                       ######### max. number of iterations reached --> use best coverage so far obtained
            if(output){
                setTxtProgressBar(pb,max.iter)
                close(pb)
                cat("\n\nMaximum number of iterations reached (max.iter =",max.iter,").\n\n")
                cat(paste(best.cov*100,"%", sep=""),"of",nrow(mat),"simulated samples simultaneously covered.\n\n")
            }
            if(timer)
                print(Sys.time()-t1)
            return( list( mat=mat, Q=best.Qs, coverage=best.cov ) )
        }
        if( (abs(coverage - include) <= tol) && (coverage - include >= 0) ){          ######### tolerance level reached AND difference is positive
            if(output){
                setTxtProgressBar(pb,max.iter)
                close(pb)
                cat("\n\n",(coverage*100),"% of",nrow(mat),"simulated samples simultaneously covered.\n\n")
            }
            if(!is.null(logfile) && class(logfile) == "character")
                cat(paste("\n\nAlgorithm terminated because the tolerance ",tol," was reached.", sep=""), file=logfile, append=TRUE, sep="\n" )
            if(timer)
                print(Sys.time()-t1)
            return( list(mat=mat, Q=best.Qs, coverage=best.cov ) )
        }
        else{                                             ######### bisection step ##########
            if( coverage - include < 0 ){                   # alpha too large
                tmp <- alpha
                alpha <- alpha - abs(alpha.old-alpha)/2
                alpha.old <- tmp
            }
            else{                                           # alpha too small
                tmp <- alpha
                alpha <- alpha + abs(alpha.old-alpha)/2
                alpha.old <- tmp
            }
        }
        if(output)
            setTxtProgressBar(pb,iter)
    }
}



#' Implements SAS-quantile Calculation (PCTLDEF=5) as Described by SAS-Help.
#' 
#' This implementation seems to be identical with 'type=2' in function \code{\link{quantile}}
#' but less efficiently implemented, i.e. for large vectors x it is much slower than the built-in
#' quantile-function.
#' 
#' @param x     (numeric) vector
#' @param prob  (numeric) value 0 < prob < 1
#' @param tol	(numeric) value used for checking numerical equivalence
#' @param type  (character) "R" = uses the R-implementation, "C" = calls the C-implementation
#' n
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' @examples
#' \dontrun{
#' SASquantile(1:100, .75)
#' 
#' ### compare to R-default
#' quantile(1:100, .75) 
#' 
#' ### or to what R calls SAS-definition
#' quantile(1:100, .75, type=3)
#' 
#' # should work for any vector (no seed)
#' v <- rnorm(50000,20,7)
#' Q.R2    <- quantile(v, probs=c(.01, .025, .05, .1, .25, .5, .75, .9, .95, .975, .99), type=2)
#' Q.SAS.R <- SASquantile(v, prob=c(.01, .025, .05, .1, .25, .5, .75, .9, .95, .975, .99), type="R")
#' Q.SAS.C <- SASquantile(v, prob=c(.01, .025, .05, .1, .25, .5, .75, .9, .95, .975, .99), type="C")
#' 
#' Q.R2
#' Q.SAS.R
#' Q.SAS.C
#' }

SASquantile <- function(x, prob, tol=1e-12, type=c("R", "C"))
{
	v <- x
	
	type <- match.arg(toupper(type[1]), c("R", "C"))
	
	if(type == "C")
	{ 
		tmp.res <- .C(	"getSASquantile", v=as.double(v), p=as.double(prob), len=as.integer(length(v)),
				  		numQ=as.integer(length(prob)), tol=as.double(tol), Q=double(length(prob)), PACKAGE="STB")
				
		res <- tmp.res$Q
		names(res) <- paste(100*prob, "%", sep="")
	}
	else
	{	
		
		if( length(prob) > 1)
			return(sapply(prob, function(x) return(SASquantile(v, x, type=type))))
		
		v <- na.omit(v)
		n <- length(v)
		vo <- sort(v)                             			# sort increasing order
		np <- n * prob
		j <- floor(np)                            			# integer part
		g <- np - j                                			# fractional part
		res <- ifelse(g < tol, mean(vo[c(j,j+1)]), vo[j+1])	# can g be treated as zero?
		names(res) <- paste(100*prob, "%", sep="")
	}
	return(res)
}



#' Simultaneous Tolerance Bounds on Residuals and Random Effects for 'VCA' Objects.
#' 
#' Simulate \eqn{N}-times data incorporating the estimated variance-covariance
#' matrix of observations \eqn{y} and construct a 100(1-alpha)\% simultaneous tolerance band. 
#' 
#' A Linear Mixed Models, noted in standard matrix notation, can be written as \eqn{y = Xb + Zg + e}, where
#' \eqn{y} is the column vector of observations, \eqn{X} and \eqn{Z} are design matrices assigning fixed (\eqn{b}),
#' respectively, random (\eqn{g}) effects to observations, and \eqn{e} is the column vector of residual errors.
#' 
#' Here, simulation is performed incorporating the variance-covariance matrix \eqn{V = ZGZ^{T}+R}{V = ZGZ'+R} of observations \eqn{y}. 
#' There are two types of random variates in a mixed model, random effects \eqn{g} and residual errors \eqn{e}. These
#' follow a multivariate normal distribution with expectation zero and covariance matrices \eqn{G} and \eqn{R}. See the 1st 
#' reference for a detailed description of the properties. 
#' Following Schuetzenmeister and Piepho (2012), a Cholesky decomposition \eqn{V = CC'} is applied to \eqn{V}, 
#' yielding the upper triangular matrix \eqn{C}, which can be used to simulate a new set of observations 
#' \eqn{y_{sim}=Cz}{y_sim = Cz}, where \eqn{z} is a vector of independent standard normal deviates of the same size as \eqn{y}.
#' Actually, \eqn{y_sim = C'z} is used here, because the Cholesky decomposition in \code{R} is defined as \eqn{V=C'C}. 
#' For each simulated dataset, random variates of interest ('term') are extracted, possibly transformed ('mode') and
#' stored in ordered form (order statistics) in a \eqn{N x n} matrix, \eqn{N} being the number of simulated datasets and
#' \eqn{n} being the number of random variates of interest. For each column of this matrix tolerance intervals 
#' are computed iteratively untill the joint coverage is as close to but >= 100(1-alpha)/% or the max. number of 
#' iterations is reached. This quantile-based algorithm is exact for \eqn{N --> Infinity}.
#' 
#' SAS-quantile definition PCTLDEF=5 is used in the fast C-implementation of the STB-algorithm (\code{\link{SASquantile}}),
#' i.e. in case \code{algo="C"}. One can compute and plot two types of plots (see argument 'type'). Simultaneous tolerance
#' bands (STB) are helpful in assessing the general distribution of a random variate, i.e. checking for departure from
#' the normality assumption. Outlying observations may also be detected using STB. Simultaneous tolerance intervals (STI)
#' are taylored for identification of extreme values (possible outliers). STI are a simplification of STB, where simultaneous
#' coverage is only required for extreme values of each simulation, i.e. an STB is constructed from the min and max values 
#' from all N simulations. This results in lower and upper bounds, which can be used in residuals plots for assessing outliers.
#' 
#' One can choose between 3 methods for constructing the 100(1-alpha)\% STB. The fastest one is the rank-based algorithm ("rank"), which
#' should only be applied for reasonably large number of simulations (rule of thumb: N>5000). For fewer simulations,
#' the quantile-based algorithm is recommended. It exists in two flavours, a native R-implementation ("R") and a pure C-implementation ("C").
#' Both can be applied using parallel-processing (see arguments 'parallel' and 'Ncpu'). Only the R-implementation allows to specify
#' a specific quantile-definition other than \code{type=2} of function \code{\link{quantile}}.
#' 
#' @param obj			(VCA) object
#' @param term			(character, integer) specifying a type of residuals if one of c("conditional",
#' 						"marginal"), or, the name of a random term (one of obj$re.assign$terms). If 'term'
#' 						is a integer, it is interpreted as the i-th random term in 'obj$re.assign$terms'. 
#' @param mode			(character) string specifying a possible transformation of random effects or 
#'                      residuals (see \code{\link{residuals.VCA}} and \code{\link{ranef.VCA}}for details)
#' @param N				(integer) specifying the number of simulated datasets \eqn{y_sim}
#' @param alpha			(numeric) value 0 < alpha < 1 specifying the min. 100(1-alpha)\% coverage of the
#'                      simultaneous tolerance band (STB)
#' @param algo			(character) (string) specifying the method to be used for constructing a 100(1-alpha)\% STB,
#' 						choose "rank" for the rank-based, "C" for a C-implementation of the quantile-based, and
#' 						"R" for an R-implentation of the quantile-based algorithm (see details). 
#' @param q.type		(integer) value specifying the quantile type to be used as offered in \code{\link{quantile}}
#' 						in case 'algo="R"'.	Whenever 'algo="C"', quantiles are computed according to SAS PCTLDEF5, 
#' 						which is identical to type 2. The rank-based algorithm does not employ quantiles.						
#' @param plot			(logical) TRUE = create 'stbVCA' object and plot it, FALSE = only create the 'stbVCA' object
#' @param legend		(logical) TRUE = add legend to the plot(s)
#' @param main1			(character) string specifying an user-defined main-title of the 'type=1' plot (STB)
#' @param main2			(character) string specifying an user-defined main-title of the 'type=2' plot (STI)
#' @param seed			(integer) value used as seed for the RNG 
#' @param orient		(integer) in QQ-plot, specifying whether to plot expected values vs. observed values (1)
#'                      or observed vs. expected (2)
#' @param type			(integer) 1 = QQ-plot with simultaneous tolerance band (STB), 2 = residual plot with simultaneous 
#' 						tolerance interval (STI), 3 = both plot at once
#' @param pb			(logical) TRUE = a text-based progress bar will display the simulation progress
#' @param parallel		(logical) TRUE = parallel processing will be attempted on 'Ncpu' cores of the local machine.
#' 						FALSE = no parallel processing applied, only in this case, 'pb' will have an effect, since this 
#'                      is only available for non-parallel processing.
#' @param Ncpu			(integer) specifying the number of CPUs on which the parallelization will be carried out.
#' 						In case that 'Ncup' is larger than the number of existing CPUs, the max. number of CPUs will be used
#' 						instead. Note, that setting 'Ncpu' to the max. number available may not result in the min. time 
#'                      spent on computing.
#' @param ...			additional arguments passed to other methods
#' 
#' @return (stbVCA) object, a list with all information needed to create the QQ-plot with ~100(1-alpha)\% STB.
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \code{\link{getSTB}}, \code{\link{fastSTB}}, \code{\link{rankSTB}}
#' 
#' @examples 
#' \dontrun{
#' library(VCA)
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1)
#' fit
#' 
#' # use studentized conditional residuals
#' stb.obj1 <- stb(fit, term="cond", mode="student", N=1000)
#' 
#' # plot it again
#' plot(stb.obj1)
#' 
#' # now request also plotting the corresponding residual plot
#' # capture additional computation results which are invisibly 
#' # returned
#' stb.obj1 <- plot(stb.obj1, type=3)
#' 
#' # use other type of legend in QQ-plot
#' plot(stb.obj1, stb.lpos="topleft")
#' 
#' # use random effects "day" and apply standardization
#' stb.obj2 <- stb(fit, term="day", mode="stand", N=1000)
#' # plot it again
#' plot(stb.obj2)
#' 
#' # more complex example
#' data(Orthodont)
#' Ortho <- Orthodont
#' Ortho$age2 <- Ortho$age - 11
#' Ortho$Subject <- factor(as.character(Ortho$Subject))
#' fit.Ortho <- anovaMM(distance~Sex+Sex:age2+(Subject)+(Subject):age2-1, Ortho)
#' 
#' # studentized conditional residuals
#' stb.cr.stud <- stb(fit.Ortho, term="cond", mode="stud", N=1000)
#' 
#' # same model fitted via REML (same covariance structure of random effects by
#' # constraining it to be diagonal)
#' fit.Ortho.reml1 <- remlMM(distance~Sex*age2+(Subject)*age2, Ortho, cov=FALSE)
#' 
#' # allow block-diagonal covariance structure of random effects due to non-zero
#' # correlation between intercept and slope of random regression part,
#' # not 'cov=TRUE' is the default
#' fit.Ortho.reml2 <- remlMM(distance~Sex*age2+(Subject)*age2, Ortho)
#' fit.Ortho.reml1
#' fit.Ortho.reml2
#' 
#' # covariance matrices of random effects 'G' differ
#' getMat(fit.Ortho.reml1, "G")[1:10, 1:10]
#' getMat(fit.Ortho.reml2, "G")[1:10, 1:10]
#' 
#' # therefore, (conditional) residuals differ
#' resid(fit.Ortho.reml1)
#' resid(fit.Ortho.reml2)
#' 
#' # therefore, STB differ
#' 
#' # studentized conditional residuals
#' system.time({
#' stb.cr.stud.reml1 <- stb(fit.Ortho.reml1, term="cond", mode="stud", 
#'                          N=5000, Ncpu=2, seed=11) })
#' system.time({
#' stb.cr.stud.reml2 <- stb(fit.Ortho.reml2, term="cond", mode="stud", 
#'                          N=5000, Ncpu=4, seed=11) })
#' 
#' # same seed-value should yield identical results
#' system.time({
#' stb.cr.stud.reml3 <- stb(fit.Ortho.reml2, term="cond", mode="stud", 
#'                          N=5000, Ncpu=4, seed=11) })
#' 
#' par(mfrow=c(1,2))
#' plot(stb.cr.stud.reml2)
#' plot(stb.cr.stud.reml3)
#' 
#' # both type of plots side-by-side
#' plot(stb.cr.stud.reml2, type=3)
#' 
#' # and enabling identification of points
#' # identified elements in the 1st plot will
#' # be automatically added to the 2nd one
#' plot(stb.cr.stud.reml2, type=3, pick=TRUE)
#' 
#' # raw "day" random effects
#' stb.re.subj <- stb(fit.Ortho, term="Subject", N=1000)
#' 
#' # identify points using the mouse
#' stb.re.subj <- plot(stb.re.subj, pick=TRUE, type=3)
#' 
#' # now click on points
#' }
#' 
#' @S3method stb VCA
#' 
#' @seealso \code{\link{fastSTB}}
#' 
#' @references 
#' 
#' Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
#' Computational Statistics and Data Analysis, 56, 1405-1416
#' 
#' Schuetzenmeister, A., Jensen, U., Piepho, H.P., 2012. Checking the assumptions of normality and homoscedasticity in 
#' the general linear model using diagnostic plots. Communications in Statistics-Simulation and Computation 41, 141-154.

stb.VCA <- function(obj, term=NULL, mode=c("raw", "student", "standard", "pearson"), 
					N=5000, alpha=.05, algo=c("rank", "R", "C"), q.type=2L, plot=TRUE, legend=TRUE, orient=1, main1=NULL, 
					main2=NULL, seed=NULL, type=1, pb=TRUE, parallel=TRUE, Ncpu=2, ...)
{
	stopifnot(class(obj) == "VCA")
	stopifnot(orient %in% 1:2)
	stopifnot(!is.null(term))
	stopifnot(type %in% 1:3)
	stopifnot(q.type %in% 1:9)
	stopifnot(alpha>0 && alpha<1)
	stopifnot(is.logical(plot))
	stopifnot(is.logical(legend))
	stopifnot(is.logical(parallel))
	stopifnot(is.numeric(Ncpu) && Ncpu > 0)
	stopifnot(is.logical(pb))
	stopifnot(is.null(main1) || is.character(main1))
	stopifnot(is.null(main2) || is.character(main2))
	stopifnot(is.numeric(N) && N > 0)
	
	algo <- toupper(algo)[1]							# use 1st option if multiple specified
	algo <- match.arg(algo, c("RANK", "R", "C"))
	
	if(algo == "RANK" && N < 5000)
		warning("The rank-based algorithm is most robustly used with larger number of simulations!")
	
	if(is.null(obj$re.assign))							# solve mixed model equations
		obj <- solveMME(obj)
	
	if(length(mode) > 1)
		mode <- mode[1]
	
	if(is.null(seed))
		seed <- runif(1, -2^30, 2^30)

	if(!is.character(term))
	{
		term <- obj$re.assign$terms[as.integer(term)]
		RE <- TRUE
		mode.opts <- c("raw", "student", "standard")		
	}
	else
	{
		term <- match.arg(term, choices=c("conditional", "marginal", obj$re.assign$terms))
		
		if(term %in% c("conditional", "marginal"))
		{
			RE <- FALSE
			mode.opts <- c("raw", "student", "standard", "pearson")
		}
		else
		{
			RE <- TRUE
			mode.opts <- c("raw", "student", "standard")			
		}
	}

	mode <- match.arg(mode, choices=mode.opts)
	m.names <- c(raw="raw", student="studentized", standard="standardized", pearson="Pearson-type")
	
	if(is.null(main1))
		main1 <- paste("QQ-Plot of ", m.names[mode], " ", ifelse(RE, "", term), ifelse(RE, " Random Effects for", " Residuals"), ifelse(RE, paste(" '", term, "'", sep=""), ""), sep="")
	if(is.null(main2))
		main2 <- "Residual Plot"
	
	sti.ylab <- paste(term, " - ", m.names[mode], ifelse(RE, "' Random Effects", "' Residuals"), sep="")

	V <- getMat(obj, "V")
	C <- Matrix(t(chol(V)))  					# Cholesky decomposition of V = C'C, transpose C since in the reference it is V = CC'
	
	if(RE)
	{
		vec  <- ranef(obj, term=term, mode=mode)
		nam <- rownames(vec)
		vec  <- vec[,1]
		names(vec) <- nam
	}
	else
	{
		vec  <- resid(obj, type=term, mode=mode)
	}
	Nvec <- length(vec)
	
	G  <- getMat(obj, "G")
	Z  <- getMat(obj, "Z")
	Vi <- getMat(obj, "Vi")
	
	REasgn <- obj$re.assign
	
	set.seed(seed)
	seeds <- runif(N, -2^30, 2^30)
	
	if(parallel)
	{
		# function to be called when parallel processing is activated
	
		Ncores <- detectCores()						# number of available cores
		if(is.null(Ncpu) || Ncpu > Ncores)
			Ncpu <- Ncores

		loopFunc <- function(iter, obj, Z, C, Vi, G, RE, term, mode, seeds)					# applied for each simulation cycle
		{
			set.seed(seeds[iter])
			obj <- obj
			#y <- C %*% Matrix:::Matrix(rnorm(obj$Nobs), ncol=1)
			y <- C %*% Matrix(rnorm(obj$Nobs), ncol=1)
			
			#simRE <- G %*% Matrix:::t(Z) %*% Vi %*% y
			simRE <- G %*% t(Z) %*% Vi %*% y
			
			obj$RandomEffects <- simRE
			obj$FixedEffects[,"Estimate"] <- 0
			obj$Matrices$y <- y
			
#			if(RE)
#				return(sort(as.numeric(VCA:::ranef.VCA(obj, term=term, mode=mode)[,1])))
#			else
#				return(sort(as.numeric(VCA:::residuals.VCA(obj, type=term, mode=mode))))
			if(RE)
				return(sort(as.numeric(ranef.VCA(obj, term=term, mode=mode)[,1])))
			else
				return(sort(as.numeric(residuals.VCA(obj, type=term, mode=mode))))
		}
		
		cluster <- makePSOCKcluster(Ncpu)
		
		vecs <- t(parSapply(cluster, 1:N, loopFunc, obj, Z, C, Vi, G, RE, term, mode, seeds))		# actually do parallel-processing of N-loops

		stopCluster(cluster)
	}	
	else
	{			
		if(pb)
		{
			cat("\n   Simulation Progress\n")	
			PB <- txtProgressBar(max=N, style=3, width=.5*unlist(options("width")), char=">")
		}
		
		for(i in 1:N)
		{
			y <- C %*% Matrix(rnorm(obj$Nobs), ncol=1)
	
			#simRE <- G %*% Matrix:::t(Z) %*% Vi %*% y
			simRE <- G %*% t(Z) %*% Vi %*% y
			
			obj$RandomEffects <- simRE
			obj$FixedEffects[,"Estimate"] <- 0
			obj$Matrices$y <- y
			
			if(RE)
				sim <- sort(ranef(obj, term=term, mode=mode))
			else
				sim <- sort(resid(obj, type=term, mode=mode))
			
			if(i == 1)
			{
				vecs <- matrix(sim, nrow=1)
			}
			else
			{
				vecs <- rbind(vecs, sim)
			}
			
			if(pb)
				setTxtProgressBar(PB, i)
		}	
		
		if(pb)
			close(PB)
	}
		
	res <- list()
	res$VCAobj <- obj
	res$mat <- vecs
	res$alpha <- alpha
	res$vec <- vec 
	res$term <- term
	res$mode <- mode
	res$plot$pch <- 16
	res$plot$col <- "black"
	res$plot$col.out <- "red"
	res$plot$pch.out <- 16
	res$plot$sti.lty=2
	res$plot$sti.lwd=1
	res$plot$sti.col="red"
	res$plot$sti.ylab <- sti.ylab
	res$plot$sti.lpos <- "title"
	res$plot$sti.main <- main2
	res$plot$stb.col <- "#0000FF40"
	res$plot$stb.lpos <- "title"
	res$plot$stb.main <- main1
	res$legend <- TRUE
	res$RE <- RE
	res$seed <- seed
	
	if(type %in% c(1,3))
	{
		if(algo == "C")
			stbObj <- fastSTB(vecs, alpha=alpha, Ncpu=Ncpu)
		else if(algo == "R")
			stbObj <- getSTB(vecs, alpha=alpha, Ncpu=Ncpu, q.type=q.type)
		else
		{
			stbObj <- rankSTB(vecs, alpha=alpha)
			stbObj$mat <- vecs
		}
		STB <- unclass(stbObj)

		means <- apply(STB$mat, 2, mean)
		res$N <- N
		STB$add <- FALSE
		STB$legend <- legend
		
		# extend STB beyond smallest and largest observation for a nicer look
		
		Q    <- STB$Q
		Nobs <- length(vec)
		
		ll <- -10 * ( (Q[1,2]-Q[1,1])/ (means[2]-means[1]) )                                # Y-value for x=-10 corvecsponding to the left-lower bound
		lu <- -10 * ( (Q[2,2]-Q[2,1])/ (means[2]-means[1]) )                                # left-upper
		rl <-  10 * ( (Q[1,Nobs]-Q[1,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )             # right-lower
		ru <-  10 * ( (Q[2,Nobs]-Q[2,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )             # right-upper
		
		Q <- cbind(c(ll,lu) ,Q)
		Q <- cbind(Q, c(rl,ru))
		
		means <- c(-10, means, 10)
		
		STB$means <- means
		STB$Q     <- Q
		STB$mat   <- NULL
		res$STB   <- STB
	}	
	if(type %in% c(2,3))
	{
		if(algo == "C")
			stbObj <- fastSTB(vecs[,c(1, ncol(vecs))], alpha=alpha, Ncpu=Ncpu)
		else if(algo == "R")
			stbObj <- getSTB(vecs[,c(1, ncol(vecs))], alpha=alpha, Ncpu=Ncpu)
		else
			stbObj <- rankSTB(vecs[,c(1, ncol(vecs))], alpha=alpha)
		STI <- unclass(stbObj)
		STI$mat <- NULL
		res$STI <- STI 
	}	
	
	if(type == 3)
	{
		old.par <- par(mfrow=c(1,2))
	}
	else
		old.par <- par

	res$type <- type
	class(res) <- "stbVCA"
	
	if(plot)
		plot(res, orient=orient)

	par(old.par)
	invisible(res)
}


#' Plot Objects of Class 'stbVCA'.
#' 
#' Standard plotting method for objects of Class 'stbVCA'.
#' 
#' This function plots objects of class 'stbVCA' as generated by function \code{\link{stb.VCA}}.
#' Objects of S3-class 'stbVCA' are list-type objects storing all the information
#' needed to plot QQ-plots with simultaneous tolerance bounds. Additionally to the information
#' contained in ordinary 'STB' objects, a copy of the 'VCA' object is stored as well as the
#' type of random variate and the mode, i.e. the type of transformation applied.
#'
#' One can specify additional parameters for changing the appearance of the plot(s). Any of the following
#' parameters can be set to a value different from the default: \cr
#' \tabular{lcl}{
#' \code{legend}   \tab ... \tab (logical) TRUE = will add legend to plot(s) (is default) \cr 
#' \code{pch}      \tab ... \tab plotting symbol for non-outlying points (default=16) \cr 
#' \code{col}      \tab ... \tab point color for non-outlying points (default="black") \cr 
#' \code{col.out}  \tab ... \tab point color for outlying points (default="red") \cr 
#' \code{pch.out}  \tab ... \tab plotting symbold for outlying points (default=16) \cr 
#' \code{stb.col}  \tab ... \tab color of the STB in the QQ-plot (default="#0000FF40") \cr 
#' \code{stb.lpos} \tab ... \tab position placement as done in function 'legend', or, character string "title" indicating\cr
#'                 \tab     \tab that legend information should be displayed as (sub-main) title (default="title") \cr 
#' \code{stb.main} \tab ... \tab character string used as main title in the QQ-plot with STB \cr 
#' \code{sti.lty}  \tab ... \tab line type for STI-bounds in the residual plot (default=2) \cr 
#' \code{sti.lwd}  \tab ... \tab line width for STI-bounds (default=1) \cr 
#' \code{sti.col}  \tab ... \tab line color for STI bounds (default="red") \cr 
#' \code{sti.ylab} \tab ... \tab character string specifying the Y-axis label in the resiudal plot with STI \cr 
#' \code{sti.lpos} \tab ... \tab position placement as done in function 'legend' (default="topright") \cr 
#' \code{sti.main} \tab ... \tab character string used as main title in the residual plot with STI \cr
#' }
#' 
#' @param x     		(object) of class 'STB' as generated by function \code{getSTB}
#' @param orient		(integer) in 'type=1' plots, 1 = expected vs. observed values, 2 = observed vs. expected values
#' @param pick			(logical, integer) TRUE = triggers function \code{\link{identify}} for identification
#' 						of points using the mouse. Stop interactive identification by pressing 'Esc' or
#' 						"Right-Mousclick --> Stop". If 'TRUE' and 'type=3', identified points will be added to the residual
#'                      plot with STI automatically. If 'pick=1', only in the QQ-plot identification will be toggled, and if
#' 						'pick=2' only in the residual plot identification will be possible. 
#' @param type			(integer) 1 = plot simultaneous tolerance band (STB), 2 = plot simultaneous tolerance interval (STI),
#' 						3 = plot STB and STI
#' @param ...   		additional arguments changing the visual appearance of the plot, e.g.
#' 						'pch', 'pch.out', 'col', 'col.out', 'stb.col', 'stb.main', "'sti.main', 'legend', ...
#' 						in fact all elements of 'x'
#' 
#' @return (stbVCA) object is invisibly returned, any additional compution results are added 
#' 
#' @author  Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @seealso \link{stb.VCA}
#' 
#' @method plot stbVCA
#' @S3method plot stbVCA
#' 
#' @examples 
#' \dontrun{
#' library(VCA)
#' data(dataEP05A2_1)
#' fit <- anovaVCA(y~day/run, dataEP05A2_1)
#' fit
#' 
#' # use studentized conditional residuals
#' stb.obj1 <- stb.VCA(fit, term="cond", mode="student", N=1000)
#' 
#' # plot it again
#' plot(stb.obj1)
#' 
#' # use random effects "day" and apply standardization
#' stb.obj2 <- stb.VCA(fit, term="day", mode="stand", N=1000)
#' 
#' # plot it again
#' plot(stb.obj2)
#' 
#' # initially, request QQ-plot with STB
#' stb.obj3 <- stb.VCA(fit, term="day", mode="stand", N=1000, type=1)
#'
#' # now request plotting of the residual plot as well
#' # catch computation result which are invisibly returned
#' stb.obj4 <- plot(stb.obj3, type=3)
#' 
#' # individualize the appearance of the plot
#' plot(stb.obj4, sti.lpos="top", col="darkblue", out.pch=17, out.col="green") 
#' }

plot.stbVCA <- function(x, orient=1, pick=FALSE, type=NULL, ...)
{
	if(!is.null(type))
		stopifnot(type %in% 1:3)
	else
		type <- x$type
	
	args <- list(...)
	x$plot[names(args)] <- args
	
	XLIM <- YLIM <- NULL
	
	if("ylim" %in% names(args))
		YLIM <- args[["ylim"]]
	if("xlim" %in% names(args))
		XLIM <- args[["xlim"]]
	
	obj <- x
	stopifnot(class(obj) == "stbVCA")
	
	if(type == 3)
		old.par <- par(mfrow=c(1,2))
	else
		old.par <- par("mfrow")
	
	m.names <- c(raw="raw", student="studentized", standard="standardized", pearson="Pearson-type")
	
	stbYlim <- Ident <- NULL
	
	if(type %in% c(1,3))
	{
		STB <- x$STB
		
		if(is.null(STB))					# STB results no yet computed
		{
			STB <- unclass(fastSTB(x$mat, alpha=x$alpha))
			means <- apply(STB$mat, 2, mean)
			STB$N <- N
			STB$add <- FALSE
			STB$stb.col <- "#0000FF40"
			STB$legend <- legend
			
			# extend STB beyond smallest and largest observation for a nicer look
			
			Q    <- STB$Q
			Nobs <- length(vec)
			
			ll <- -10 * ( (Q[1,2]-Q[1,1])/ (means[2]-means[1]) )                                # Y-value for x=-10 corvecsponding to the left-lower bound
			lu <- -10 * ( (Q[2,2]-Q[2,1])/ (means[2]-means[1]) )                                # left-upper
			rl <- 10 * ( (Q[1,Nobs]-Q[1,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )             # right-lower
			ru <- 10 * ( (Q[2,Nobs]-Q[2,(Nobs-1)])/ (means[Nobs]-means[(Nobs-1)]) )             # right-upper
			
			Q <- cbind(c(ll,lu) ,Q)
			Q <- cbind(Q, c(rl,ru))
			
			means <- c(-10, means, 10)
			
			STB$means <- means
			STB$Q     <- Q
			STB$mat <- NULL						# simulated data stored in upper most level
			x$STB <- STB
		}
		
		N <- length(STB$means)
		
		nam <- names(x$vec)
		ord <- order(x$vec)
		vec <- x$vec[ord]
		names(vec) <- nam[ord]
		
		out.high <- which(vec > STB$Q[2,-c(1,ncol(STB$Q))])
		out.low <- which(vec < STB$Q[1, -c(1,ncol(STB$Q))])
		out  <- unique(c(out.high, out.low))
		
		if(orient[1] == 1)
		{
			Xval <- STB$means[-c(1,length(STB$means))]
			Yval <- vec
			Xpol <- c(STB$means, rev(STB$means))
			Ypol <- c(STB$Q[1,], rev(STB$Q[2,]))
			xlim <- range(Xval)
			ylim <- c(min(c(STB$Q[1,-1], Yval)), max(c(STB$Q[2,-N], Yval)))
			xlab <- "Expected Values"
			ylab <- "Observed Values"
		}
		else
		{
			Yval <- STB$means[-c(1,length(STB$means))]
			Xval <- vec
			Ypol <- c(STB$means, rev(STB$means))
			Xpol <- c(STB$Q[1,], rev(STB$Q[2,]))
			ylim <- range(Yval)
			xlim <- c(min(c(STB$Q[1,-1], Xval)), max(c(STB$Q[2,-N], Xval)))
			ylab <- "Expected Values"
			xlab <- "Observed Values"
		}	
		if(is.null(x$stb.main))
			x$plot$stb.main <- paste("QQ-Plot of ", m.names[x$mode], " ", ifelse(x$RE, "", x$term), ifelse(x$RE, " Random Effects for", " Residuals"), ifelse(x$RE, paste(" '", x$term, "'", sep=""), ""), sep="")
		
		if(!is.null(XLIM))
			xlim <- XLIM
		if(!is.null(YLIM))
			ylim <- YLIM
		
		plot(Xval, Yval, ylim=ylim, xlim=xlim,
			 type="n", main=x$plot$stb.main, xlab=xlab, ylab=ylab)
		
		polygon(x=Xpol, y=Ypol, border=NA, col=x$plot$stb.col )
		
		if(!STB$add)
		{
			abline(0,1)
			box()
						
			if(length(out > 0))
			{
				points( Xval[-out], Yval[-out], pch=x$plot$pch, col=x$plot$col ) 			
				points(Xval[out], Yval[out], pch=x$plot$pch.out, col=x$plot$col.out)     # means[1] and means[N] have no counterpart in vec
			}
			else
				points( Xval, Yval, pch=x$plot$pch, col=x$plot$col ) 
			
			if(x$legend)
			{
				if(x$plot$stb.lpos == "title")
				{
					mtext(text=bquote(paste("Simultaneous Tolerance Band (coverage = ", .(100*round(STB$coverage,5)),"% ;", N[sim]==.(x$N), ")", sep="")),
						  side=3, line=.25, cex=.9)
				}
				else
				{
					legend(x$plot$stb.lpos, fill=c(x$plot$stb.col, par("bg")), 
						   legend=c(paste("STB (", 100*round(STB$coverage,5),"% coverage)",sep=""),
						   			paste("N=", x$N," Simulations", sep="")), 
						   box.lty=0, border=NA)
			    }
		   }
		}
		
		if(orient == 1 && type == 3)
			stbYlim <- ylim
		
		if(pick && pick != 2)
		{
			Ident <- identify(x=Xval, y=Yval, labels=names(vec), pos=TRUE)
			nam <- names(vec[Ident$ind])
			Ident$ind <- which(names(x$vec) %in% nam)
		}
	}
	
	if(type %in% c(2,3))
	{
		STI <- x$STI
		
		if(is.null(STI))
		{
			STI <- unclass(fastSTB(x$mat[,c(1, ncol(x$mat))], alpha=x$alpha))
			STI$mat <- NULL
			x$STI <- STI 
			if(is.null(x$plot$sti.main))
				x$plot$sti.main <- "Residual Plot" 
		}
		else
		{
			STI <- x$STI
		}
		
		Xval <- 1:length(x$vec)
		Yval <- x$vec

		if(!is.null(stbYlim))
			ylim <- stbYlim
		else
			ylim <- range(c(x$vec, c(STI$Q)))

		plot(Xval, Yval, type="n", main=x$plot$sti.main, ylim=ylim, 
			 ylab=x$plot$sti.ylab, xlab="i-th random variate")
	 
	 	abline(h=STI$Q[1,1], col=x$plot$sti.col, lty=x$plot$sti.lty, lwd=x$plot$sti.lwd)
		abline(h=STI$Q[2,2], col=x$plot$sti.col, lty=x$plot$sti.lty, lwd=x$plot$sti.lwd)
		
		out.high <- which(x$vec > STI$Q[2,2])
		out.low  <- which(x$vec < STI$Q[1,1])
		out  <- unique(c(out.high, out.low))
		
		if(length(out > 0))
		{
			points( Xval[-out], Yval[-out], pch=x$plot$pch, col=x$plot$col ) 			
			points(Xval[out], Yval[out], pch=x$plot$pch.out, col=x$plot$col.out)     # means[1] and means[N] have no counterpart in vec
		}
		else
			points( Xval, Yval, pch=x$plot$pch, col=x$plot$col ) 
		
		if(x$legend)
		{
			if(x$plot$sti.lpos=="title")
			{
				mtext(text=bquote(paste("Simultaneous Tolerance Interval (coverage = ", .(100*round(STI$coverage,5)),"% ;", N[sim]==.(x$N), ")", sep="")),
					  side=3, line=.25, cex=.9)
			}
			else
			{
				legend(x$plot$sti.lpos, lty=c(x$plot$sti.lty, 0), col=c(x$plot$sti.col, par("bg")), 
						legend=c(paste("STB (", 100*round(STI$coverage,5),"% coverage)",sep=""),
								paste("N=",x$N," Simulations", sep="")), 
						box.lty=1, bg="white")
			}
		}
		
		if(pick)
		{			
			if(is.logical(pick))
			{
				if(is.null(Ident))
					identify(x=Xval, y=Yval, labels=names(x$vec))
				else
				{
					if(length(Ident$ind) > 0)
						text(Xval[Ident$ind], Yval[Ident$ind], labels=names(x$vec[Ident$ind]), pos=Ident$pos)
				}
			}
			else
			{
				if(pick == 2)
					identify(x=Xval, y=Yval, labels=names(x$vec))
			}
		}
	}
	x$type <- type
	par(old.par)
	invisible(x)
}



#' Rank-Based Algorithm for Computing 100(1-alpha)\% Simulataneous Tolerance Bounds.
#' 
#' Implementation of a rank-based algorithm for constructing 100(1-alpha)\% STBs as outlined in the reference.
#' 
#' This function is a performance optimized version of the original rank-based algorithm avoiding the time-consuming
#' iteration. In principle it sorts out simulation results which have at least one extreme order statistic untill
#' exactly 100(1-alpha)\% of all simulation results remain. From these, bounds of the STB are constructed determining
#' extreme-values per order-statistic (column).
#' 
#' This implementation also corrects step 4) of the published algorithm, which has to find those indices of elements being
#' equal to "min(c)" OR being equal to "N-min(c)+1". This reflects the construction of vector "c", where max. rank-values
#' are transformed to min. rank-values. In step 6) the "N_{k-1}-(1-alpha)*N largest" elements of "d_{l}^{theta}" have to
#' be selected, which needs also correction.
#' 
#' Parallel processing did not minimize the computation time in contrast to the algorithms for computing the quantile-based algorithm.
#' Thus, parallel processing is not supported for this algorithm.
#' 
#' @param mat			(numeric) matrix of dimension (N, n), where i-th row corrsponds to ordered values 
#'                      of the i-th simulation
#' @param alpha			(numeric) value defining the desired coverage as 100(1-alpha)\% 
#' 
#' @return (list) with two elements:\cr
#' \item{Q}{(matrix) 1st row stores lower bounds, 2nd row upper bounds}
#' \item{cov}{(numeric) value corresponding the coverage}
#' 
#' @author Andre Schuetzenmeister \email{andre.schuetzenmeister@@roche.com}
#' 
#' @references
#' 
#' Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
#' Computational Statistics and Data Analysis, 56, 1405-1416
#' 
#' @examples 
#' \dontrun{
#' # for following problem size the rank-based algo
#' # outperforms the quantile based one, although,
#' # ran serially
#' mat <- matrix(rnorm(10000*100), ncol=100)
#' mat <- t(apply(mat, 1, sort))
#' system.time(stb.rank <- rankSTB(mat))
#' system.time(stb.q.R  <- getSTB(mat))
#' system.time(stb.q.C  <- fastSTB(mat))
#' x <- apply(mat, 2, mean)
#' plot(x,x, ylim=c(-5,5))
#' lines(x, stb.q.R$Q[1,], col="blue", lwd=2)
#' lines(x, stb.q.R$Q[2,], col="blue", lwd=2)
#' lines(x, stb.q.C$Q[1,], col="red",  lwd=2)
#' lines(x, stb.q.C$Q[2,], col="red",  lwd=2)
#' lines(x, stb.rank$Q[1,],  col="cyan", lwd=2)
#' lines(x, stb.rank$Q[2,],  col="cyan", lwd=2)
#' legend("top", legend=c("R-quantile", "C-quantile", "rank-based"), 
#'        fill=c("blue", "red", "cyan"))
#' 
#' # varying Ncpu for the C-implementation of the quantile-based algo
#' system.time(stb.q.C  <- fastSTB(mat, Ncpu=4))
#' system.time(stb.q.C  <- fastSTB(mat, Ncpu=6))
#' system.time(stb.q.C  <- fastSTB(mat, Ncpu=8))
#' system.time(stb.q.C  <- fastSTB(mat, Ncpu=10))
#' system.time(stb.q.C  <- fastSTB(mat, Ncpu=12))
#' }

rankSTB <- function(mat, alpha=0.05)
{
	stopifnot(is.matrix(mat) && is.numeric(mat))
	stopifnot(alpha>0 && alpha<1)
	
	ind <- min_v <- NULL
	N   <- nrow(mat)
	N0  <- N1 <- N													# is N_{k-1} in the paper 
	n <- ncol(mat)

	D <- apply(mat, 2, function(x) abs(scale(x)))	# scaled absolute values per column
	C <- apply(mat, 2, rank)						# per-column rank values of mat
	v <- apply(C, 1, function(x){
				return(min(min(x), N-max(x)+1))})

	v2 <- sort(v)												# this code replaces the loop in the paper

	thresh <- N*alpha											
	cs <- 0
	for(i in 1:length(v2))										# determines the rank-value needed for exact coverage
	{
		cs <- cs + 1
		if(cs == thresh)
		{
			min_v <- v2[i]-1
			break
		}
	}

	if(min_v > 0)												# case: number of rows with rank = 1 > N*(1-alpha)
	{
		ind <- which(v <= min_v)		
		mat <- mat[-ind,]
		C   <- C[-ind,]
		D   <- D[-ind,]
		v   <- v[-ind]
		N0  <- nrow(mat)
	}
	
	ind <- which(v == (min_v+1))								# if reached exact coverage N0-N1 == N*(1-alpha)
	min_v <- min_v + 1
	N1  <- length(ind)											# loop-replacing code ends here

	if(N0-N1 < N*(1-alpha))										# apply post-processing
	{
		D0 <- D[ind,,drop=F]									# is D^{theta} in the paper
		C0 <- C[ind,,drop=F]									# is C^{theta} in the paper
		ind0 <- apply(C0, 1, function(x) which(x %in% c(min_v, N-min_v+1)))		# column-indices

		if(is.list(ind0))										# multiple elements in row corresponding to min_v
		{
			d0 <- NULL
			for(i in 1:length(ind0))
			{
				d0 <- c(d0, max(D0[i,ind0[[i]]]))
			}
		}
		else
			d0   <- D0[,ind0]

		num  <- N0 - N*(1-alpha)								# that many items need to be removed
		ind1 <- order(d0)[1:num]
		mat <- mat[-ind1,]										# remove as many rows as needed to have exact coverage
	}
	else
		mat <- mat[-ind,]										# exact coverage without post-processing required

	lower <- apply(mat, 2, min)
	upper <- apply( mat, 2, max)
	res <- list(Q=rbind(lower, upper), coverage=nrow(mat)/N)			
	return(res)									
}

#' Generic Method for computing Simultaneous Tolerance Bounds (bands and intervals) 
#' for fitted models or numeric vectors.
#' 
#' If the generic method is applied to an object of class \code{VCA} function \code{\link{stb.VCA}}
#' is called and all its arguments may be specified. Otherwise, function \code{\link{stb.default}}
#' is called and all its arguments may be specified. The latter was developed for numeric vectors
#' (see ?stb.default). Therefore this method will most likely produce strange results on other
#' types of objects than \code{VCA} and numeric vectors.
#' 
#' @param obj		(object) passed on
#' @param ...		additional parameters

stb <- function(obj, ...)
	UseMethod("stb")
