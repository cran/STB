% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stb.R
\name{plot.STB}
\alias{plot.STB}
\title{Plot Objects of Class 'STB'.}
\usage{
\method{plot}{STB}(x, ...)
}
\arguments{
\item{x}{(object) of class 'STB' as generated by function \code{getSTB}}

\item{...}{arguments passed to other methods}
}
\description{
Standard plotting method for objects of class 'STB'.
}
\details{
This function plots objects of class 'STB' as generated by function \code{\link{stb}}.
Objects of S3-class 'STB' are list-type objects storing all the information
needed to plot QQ-plots with simultaneous tolerance bounds.
}
\examples{
\dontrun{
### generate an 'STB' object without plotting
obj <- stb(rnorm(30), plot=FALSE)
plot(obj)

### manipulate the 'STB' object for plotting
obj$legend=TRUE
plot(obj)

### add a previously generated STB-ocject to an existing plot
plot(sort(rnorm(30)), sort(rnorm(30)))
obj$add <- TRUE
plot(obj) 
}
}
\seealso{
\link{stb}
}
\author{
Andre Schuetzenmeister \email{andre.schuetzenmeister@roche.com}
}
