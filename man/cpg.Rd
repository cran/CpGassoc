\name{cpg}
\alias{plot.cpg}

\alias{summary.cpg}
\alias{print.cpg}
\alias{sort.cpg}


\title{Methods for object of class \code{"cpg"}}

\usage{

  \method{plot}{cpg}(x, save.plot = NULL, file.type = "pdf", popup.pdf = FALSE, 
                     tplot = FALSE, classic = TRUE,main.title = NULL, eps.size = c(5, 5),
                     gc.p.val = FALSE, gcdisplay = FALSE, \dots)

  \method{summary}{cpg}(object,\dots)

  \method{print}{cpg}(x,\dots)

  \method{sort}{cpg}(x,decreasing,\dots)
}

\arguments{
   \item{x}{
Output of class \code{"cpg"} from cpg.assoc or cpg.work.
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
  \item{tplot}{
Logical. If \code{TRUE}, t-statistics will be plotted vs. their expected quantiles. If \code{FALSE} (default), -log(p) will be
plotted. (Note: if \code{class(x$indep)=='factor'} this option will be ignored.)
}
  \item{classic}{
Logical. If \code{TRUE}, a classic qq-plot will be generated, with all p-values plotted against predicted values (including significant).
If \code{FALSE} Holm-significant CpG sites will not be used to compute expected quantiles and will be plotted separately.
}
  \item{main.title}{
Main title to be put on the graph. If \code{NULL} one based on the analysis will be used.
}
  \item{eps.size}{
Vector indicating the size of .eps file (if creating one). Corresponds to the options horizontal and height in the
\code{postscript} function.
}
  \item{gc.p.val}{
Logical. If true, plot will use the genomic control adjusted p-values.
}
  \item{gcdisplay}{
Logical.If true, plot will display the genomic control value in the legend.
}
  \item{object}{
Output of class \code{"cpg"} from \code{cpg.assoc} or \code{cpg.work}.
  }
  \item{decreasing}{
logical. Should the sort be increasing or decreasing? Not available for partial sorting.
}
  \item{\dots}{
Arguments to be passed to methods, such as graphical parameters.
}
}
\description{
Methods and extra functions for class \code{"cpg"}.
\code{plot.cpg} creates a QQ plot based on the association p-values or t-statistics from the function \code{cpg.assoc}.
}

\value{
\code{sort.cpg} returns an item of class \code{"cpg"} that is sorted by p-value.
\code{summary.cpg} creates a qq-plot based on the data, and scatterplots or boxplots for the top sites.
}

\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <barfieldrichard8@gmail.com>
}
\note{
Plots with empirical confidence intervals based on permutation tests can be obtained from \code{cpg.perm}.
See \code{\link{plot.cpg.perm}} for more info.
}



\seealso{
\code{\link{cpg.perm}}
\code{\link{cpg.assoc}}
\code{\link{scatterplot}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}

}
\examples{
##QQ Plot:
data(samplecpg,samplepheno,package="CpGassoc")
test<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
plot(test)
##t-statistic plot:
plot(test,tplot=TRUE)


#Getting our plot:
plot(test,classic=FALSE)


##Now an example of sort
head(sort(test)$results)

##Summary
summary(test)
}
