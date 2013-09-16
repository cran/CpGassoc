\name{cpg.everything}
\alias{cpg.everything}
\alias{cpg.everything.character}
\alias{cpg.everything.complex}
\alias{cpg.everything.logical}
\alias{cpg.everything.matrix}
\alias{cpg.everything.numeric}


\title{
Multi-Task function   
}
\description{
A function designed to do a group of smaller functions required for the \code{cpg.assoc}
}
\usage{
cpg.everything(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Just a generic object
}
  \item{\dots}{
Arguments to be passed to methods.
}

}
\details{
A function created to do a bunch of much smaller tasks within \pkg{CpGassoc} based on the class of x.
}


\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
}
\note{
\code{cpg.everything} is designed to perform a multitude of smaller tasks for \code{cpg.assoc} that do not warrant a full
help page. These include:warnings, getting the name of the independent variable, designing the random function
to be used, and getting the names for the values returned.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cpg.assoc}}
\code{\link{cpg.perm}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Has four methods:character,complex, numeric/matrix, and logical
#They correspon to getting the indep variable name, warnings, getting the random function,
#and getting the names for the values returned. For the design of these functions see
#the R code. 
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
