% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stb.R
\name{stb.default}
\alias{stb.default}
\title{Simultaneous Tolerance Bands (STB).}
\usage{
\method{stb}{default}(
  obj,
  N = 10000L,
  alpha = 0.05,
  rand.func = rnorm,
  tol = 1e-04,
  max.iter = 100L,
  algo = c("rank", "C", "R"),
  Ncpu = 1,
  q.type = 2L,
  stb.col = "#0000FF40",
  col.points = "black",
  col.out = "red",
  col.pwb = "#0000FF40",
  main = NULL,
  add.pwb = FALSE,
  quiet = FALSE,
  add = FALSE,
  plot = TRUE,
  legend = FALSE,
  timer = FALSE,
  pch = 16,
  pch.out = 16,
  seed = NULL,
  ...
)
}
\arguments{
\item{obj}{(numeric) vector, which is supposed to be N(my, sigma^2)-distributed}

\item{N}{(integer) value specifying the number of random samples to be used for constructing the STB}

\item{alpha}{(numeric) value specifying the simultaneous tolerance level, i.e. 100(1-alpha)\% of all 'N'
random samples have to be completely enclosed by the bounds of the STB}

\item{rand.func}{(function) a function which generates random samples, e.g. \code{rand.func=rnorm} which corresponds to random
sampling from the standard normal distribution. Another example is defining 'sfunc=function(n) rchisq(n=n, df=3, ncp=2)'
and using \code{rand.func=func}. See examples for further examples.}

\item{tol}{(numeric) value specifying the max. acceptable deviation from 'alpha' used in the bisection algorithm}

\item{max.iter}{(integer) value specifying the max. number of iteration for finding the bounds of the STB}

\item{algo}{(character) (string) specifying the method to be used for constructing a 100(1-alpha)\% STB,
choose "rank" for the rank-based, "C" for a C-implementation of the quantile-based, and
"R" for an R-implentation of the quantile-based algorithm (see details). "C" uses SAS PCTLDEF5 definition of quantiles,
whereas "R" can use any of the built-in R types of quantiles (see \code{\link{quantile}}.}

\item{Ncpu}{(integer) specifying the number cores/CPUs to be used, for N>1 multi-processing is applied}

\item{q.type}{(integer) the quantile-type used if \code{algo="R"}, see ? quantile for details.}

\item{stb.col}{(character) string, a valid specification of a color to be used for the STB}

\item{col.points}{(character) color for the points in the QQ-plot}

\item{col.out}{(character) color for points outsied of the 100(1-alpha)\% STB}

\item{col.pwb}{(character) color for the point-wise STB (not adjusted for multiplicity), defaults to "#0000FF40" which is "blue"
with 80\% transparency}

\item{main}{(characer) string for a main title appearing over the plot}

\item{add.pwb}{(logical) should the point-wise tolerance band be plotted for comparison?}

\item{quiet}{(logical) TRUE = no additional output ist printed (progress bar etc.)}

\item{add}{(logical) TRUE = the 100(1-alpha)\% STB is added to an existing plot}

\item{plot}{(logical) TRUE = either a QQ-plot with STB (add=FALSE) or a STB added to an existing plot (add=TRUE) is plotted.
FALSE = only computations are carried out without plotting, an object of class 'STB' is returned which
        can be stored an plotted later on, e.g. to avoid computing an STB every time a Sweave/mWeave report is updated}

\item{legend}{(logical) TRUE a legend is plotted "topleft"}

\item{timer}{(logical) TRUE = the time spent for computing the STB will be printed}

\item{pch}{(integer) plotting symbols for the QQ-plot}

\item{pch.out}{(integer) plotting symbols for outlying points}

\item{seed}{(numeric) value interpreted as integer, setting the random number generator (RNG) to a defined state}

\item{...}{further graphical parameters passed on}
}
\value{
invisibly returns a list-type object of class \code{STB}, which comprises all arguments accepted by this function.
}
\description{
Compute And/Or Plot Simultaneous Tolerance Bands for numeric vectors.
}
\details{
Function takes a numeric vector 'vec' and computes the 100(1-alpha)\%-simultaneous tolerance band (STB) for
the (DEFAULT )Null-hypothesis H0: vec~N(my, sigma^2) distributed, which is equal to checking whether the residuals of 
the simplest linear model y = mu + e (y~1) are normally distributed, i.e. 'e ~ N(0, sigma^2)'. By specification of argument \code{rand.func} 
other null-distributions can be specified. One has to specify a function with a single argument 'n',  which returns a random sample with 'n' 
observations, randomly drawn from the desired null-distribution (see description argument \code{rand.func} below). 
Note that all random samples as well as vector \code{vec} will be centered to mean=0 and scaled to sd=1.

One can choose between three methods for constructing the 100(1-alpha)\% STB. There are two implementations of the quantile-based algorithm
("C", "R" see 1st reference) and one of the rank-based algorithm (see 2nd reference). Methods "C" and "R" can be run in parallel. The rank-based
algorithm does not benefit form parallel processing, at least not in the current implementation. It is still the default and recommended for small to medium 
sized vectors and 10000 <= N <= 20000 simulations, which should be sufficiently accurate reflect the null-distribution. The "C" and "R" options refer
to implementations of the quantile-based algorithm. The "C" implementation benefits most from using multiple cores, i.e. the larger 'Ncpu' the better,
and should be used for large problems, i.e. rather large number of elements and large number of simulations.

The table below gives an impression how these algorithms perform. Runtimes were measured under Windows 7 on a Intel Xeon E5-2687W 3.1 GHz workstation with 16
logical cores and 16 GB RAM. The runtime of the C-implementation of the quantile-based algorithm is denoted as "t_qC12" applied parallely with 12 cores.
Each setting was repeated 10 times and the overall run time was then divided by 10 providing sufficiently robust simulation results.
Note, that for smaller problem sizes a large proportion of the overall runtime is due to simulating, i.e. drawing from the null-distribution.

\tabular{rrrr}{
	\strong{_____N_obs} \tab \strong{_____N_sim} \tab \strong{____t_rank} \tab \strong{____t_qC12} \cr
 25    \tab 5000  \tab  0.4s \tab   0.5s \cr
 25	  \tab 10000 \tab  0.8s \tab   1.3s \cr
 50	  \tab 10000 \tab  1.0s \tab   3.2s \cr 
 100	  \tab 10000 \tab  1.7s \tab   2.9s \cr
 100	  \tab 20000 \tab  3.0s \tab   4.8s \cr
 225	  \tab 20000 \tab  5.1s \tab   8.3s \cr
 300	  \tab 30000 \tab  9.6s \tab  17.2s \cr
 300	  \tab 50000 \tab 16.1s \tab  24.9s \cr
 1000  \tab 50000 \tab 47.8s \tab 123.5s \cr
}
}
\examples{
### log-normal vector to be checked for normality
\dontrun{
set.seed(111)
stb(exp(rnorm(30)), col.out="red", legend=TRUE)

### uniformly distributed sample checked for Chi-Squared Distribution with DF=1, degrees of freedom
set.seed(707)
stb(runif(25, -5, 5), rand.func=function(n){rchisq(n=n, df=1)}, 
    col.out="red", legend=TRUE, main="Chi-Squared with DF=1")

### check whether an Chi-Squared (DF=1) random sample better fits 
stb(rchisq(25, df=1), rand.func=function(n){rchisq(n=n, df=1)}, 
    col.out="red", legend=TRUE, main="Chi-Squared with DF=1")

### add STB to an existing plot
plot(sort(rnorm(30)), sort(rnorm(30)))
stb(rnorm(30), add=TRUE)

### compute STB for later use and prevent plotting
STB <- stb(rnorm(30), plot=FALSE)
}
}
\references{
Schuetzenmeister, A., Jensen, U., Piepho, H.P. (2011), Checking assumptions of normality and homoscedasticity in the general linear model. 
            Communications in Statistics - Simulation and Computation; S. 141-154

Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach.
Computational Statistics and Data Analysis, 56, 1405-1416
}
\seealso{
\link{plot.STB}
}
\author{
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}
}
