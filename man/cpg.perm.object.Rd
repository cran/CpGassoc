\name{Class cpg.perm}
\alias{plot.cpg.perm}
\alias{summary.cpg.perm}
\alias{print.cpg.perm}
\alias{sort.cpg.perm}

\title{
Methods for object of class \code{"cpg.perm"}.
}
\usage{
  \method{plot}{cpg.perm}(x, save.plot = NULL, file.type = "pdf", popup.pdf = FALSE,
  main.title = NULL, eps.size = c(5, 5), tplot = FALSE, perm.ci = TRUE, classic = TRUE,
  gc.p.val = FALSE, gcdisplay = FALSE, ...)

  \method{summary}{cpg.perm}(object,\dots)

  \method{print}{cpg.perm}(x,\dots)

  \method{sort}{cpg.perm}(x,decreasing,\dots)
}

\arguments{
  \item{x}{
Output from \code{cpg.perm}. Of class \code{"cpg.perm"}.
}
  \item{save.plot}{
Name of the file for the plot to be saved to. If not specified, plot will not be saved.
}
  \item{file.type}{
Type of file to be saved. Can either be \code{"pdf"} or \code{"eps"}. Selecting \code{file.type="eps"} will
result in publication quality editable postscript files that can be opened by Adobe Illustrator or Photoshop.
  }
  \item{popup.pdf}{
\code{TRUE} or \code{FALSE}. If creating a pdf file, this indicates if the plot should appear in a popup window as well. If running in a 
cluster-like environment, best to leave \code{FALSE}.
}
  \item{main.title}{
Main title to be put on the graph. If \code{NULL} one based on the analysis will be used.
}
  \item{eps.size}{
Vector indicating the size of .eps file (if creating one). Corresponds to the options horizontal and height in the
\code{postscript} function.
}
  \item{tplot}{
Logical. If \code{TRUE}, ordered t-statistics will be plotted against their expected quantities. If \code{FALSE} (default), -log(p) will be plotted. If indep is a class variable this option will be
ignored.
}
  \item{perm.ci}{
Logical. If \code{TRUE}, the confidence intervals computed will be from the permutated values, otherwise will be based on the theoretical values.
}
  \item{classic}{
Logical. If \code{TRUE}, a classic qq-plot will be generated, with all p-values plotted against predicted values (including significant).
If \code{FALSE} Holm-significant CpG sites will not be used to compute expected quantiles and will be plotted separately.
}
  \item{gc.p.val}{
Logical. If true, plot will use the genomic control adjusted p-values.
}
  \item{gcdisplay}{
Logical.If true, plot will display the genomic control value in the legend.
}
  \item{object}{
Output of class \code{"cpg.perm"} from \code{"cpg.perm"}.
  }
  \item{decreasing}{
logical. Should the sort be increasing or decreasing? Not available for partial sorting.
}
  \item{\dots}{
Arguments to be passed to methods, such as graphical parameters.
}
}

\description{
Methods and extra functions for class \code{"cpg.perm"}.
\code{plot.cpg.perm} creates a QQ plot based on the association p-values or t-statistics from the function \code{cpg.perm}.
}

\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <barfieldrichard8@gmail.com>
}
\note{
Empirical confidence intervals will be computed only if there are a hundred or more permutations. Otherwise the theoretical confidence intervals will be plotted.
}


\seealso{
\code{\link{cpg.perm}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{manhattan}}
\code{\link{cpg.assoc}}
}
\examples{
data(samplecpg,samplepheno,package="CpGassoc")

#The qq plot:
Testperm<-cpg.perm(samplecpg[1:300,],samplepheno$weight,seed=2314,nperm=10,large.data=FALSE)
plot(Testperm)
#The t-statistic plot from cpg.perm has confidence intervals since we were allowed 
#to perform permutations on the T-values.
plot(Testperm,tplot=TRUE)
#If there was 100 or more permutations, there would be emperical confidence intervals.

#Getting an example of the non classic QQ plot
plot(Testperm,classic=FALSE)


###Now for Sort
head(sort(Testperm)$results)
head(Testperm$results)



}
