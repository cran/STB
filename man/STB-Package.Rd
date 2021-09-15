\docType{package}
\name{STB-package}
\alias{STB}
\alias{STB-package}
\title{(S)imultaneous (T)olerance (B)ounds}
\description{
	Compute simultaneous tolerance bounds for arbitrary null-distributions or random variates of linear mixed models (LMM).
	A common problem is to check whether a numeric values follow a distinct distribution. Formal significance tests lack power
	when its most needed, when sample size is small and have sort of too much power for large sample sizes, i.e. rejecting
	e.g. normality even if departure from normality is negligible for drawing inference.\cr
	
	Graphical methods are usually preferred over formal testing but lack some objectivity, i.e.
	the same plot is differently interpreted by two different persons. Simultaneous tolerance bounds (intervals or bands)
	add some objectivity to plots like the QQ-plot, a plot frequently used to check for specific distributions. See the first reference
	for details about simultaneous tolerance bounds (there, used for checking normality in the general linear model). 
	Various other null-distributions can be checked, since resampling can be applied to all sorts of distributions (see  argument 'rand.func' of \code{\link{stb}}).\cr
	
	Besides checking distributional assumptions of numeric vectors, graphical methods are available for linear mixed models fitted using R-package \code{\link{VCA}}.
	Various types of random variates exist in this framework, i.e. random effects and at least two types of residuals. All types of random variates of 
	LMM need to be checked for their particular distributional assumptions. Departure from these may indicate that a transformation of the response variable 
	is required or that extreme values (possible outliers) have a negative effect on the distribution of random variates. See the 2nd reference for
	a discussion of this topic and for practical examples.\cr
	
	Generic function \code{\link{stb}} applies to either numeric vectors or to objects of class \code{VCA}. The former generates a quantile-quantile (QQ) plot 
	following a specified null-distribution with simultaneous tolerance band (see \code{\link{stb.default}}). The latter checks residuals or random effects,
	possibly transformed (e.g. studentized) to be normally distributed, incorporating the covariance-structure of the fitted model (see \code{\link{stb.VCA}}). 
	Plotting methods exists for the resulting objects (see \code{\link{plot.STB}} and \code{\link{plot.stbVCA}} for details). 
	There are two graphical methods available for random variates of LMM, the usual QQ-plot with simultaneous tolerance band (argument \code{type} set to 1 or 3)
	and a plot of the random variates with a simultaneous tolerance interval (STI) (argument \code{type} set to 2 or 3), where the STI is derived from the extreme
	values (first and last order statistics).
}
\details{\tabular{ll}{
Package: \tab STB\cr
Type: \tab Package\cr
Version: \tab 0.6.5\cr
Date: \tab 2021-09-14\cr
License: \tab GPL (>= 3)\cr
LazyLoad: \tab yes\cr
}}
\author{Andre Schuetzenmeister <andre.schuetzenmeister@roche.com>}
\references{
	Schuetzenmeister, A., Jensen, U., Piepho, H.P. (2011), Checking assumptions of normality and homoscedasticity in the general linear model. Communications in Statistics - Simulation and Computation; S. 141-154
	
	Schuetzenmeister, A. and Piepho, H.P. (2012). Residual analysis of linear mixed models using a simulation approach. Computational Statistics and Data Analysis, 56, 1405-1416
}
\keyword{package}

