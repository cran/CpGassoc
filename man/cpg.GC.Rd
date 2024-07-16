\name{cpg.GC}
\alias{cpg.GC}

\alias{print.cpg.gc}
\alias{print.cpg.perm.gc}
\title{
cpg.GC and methods for output of function
}
\description{
cpg.GC accepts an object of class \code{"cpg.perm"} or \code{"cpg"} and returns information
regarding Holm and FDR-significance of the GC (genomic control) adjusted test statistics. For \code{"cpg.perm"}
will return permutation p-values based on the GC-adjusted values from each permutation.
}
\usage{
   cpg.GC(x)

  \method{print}{cpg.gc}(x,\dots)
  
  \method{print}{cpg.perm.gc}(x,\dots)
 
}

\arguments{
  \item{x}{
Object of class \code{"cpg.perm"} or \code{"cpg"}. 
}
  \item{\dots}{
Arguments to be passed to methods, such as graphical parameters.
}
}
\details{
\code{cpg.GC} will display the number of Holm and FDR-significant sites using the genomic control adjusted p-values
test statistics. It will also display the estimated genomic control inflation factor.
}
\value{
\code{cpg.GC} returns an object of class \code{"cpg.gc"} or \code{"cpg.perm.gc"}
  \item{gc.results}{
Matrix consisting of GC-adjusted test statistics for each CpG site. Similar to the results output of \code{cpg.assoc}. 
      }
  \item{gc.info}{
Data frame with information on the number of Holm and FDR significant sites. Will also have the genomic control
inflation estimate. Objects from \code{"cpg.perm"} will also have information concerning the permutation p-values.
 }
}
\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <barfieldrichard8@gmail.com>
}


\seealso{
\code{\link{cpg.work}}
\code{\link{cpg.perm}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{cpg.combine}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}
\code{\link{sort.cpg.perm}}
\code{\link{sort.cpg}}
\code{\link{cpg.qc}}
\code{\link{cpg.assoc}}
}
\examples{
data(samplecpg,samplepheno,package="CpGassoc")
results<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)

cpg.GC(results)
##If the genomic inflation factor is less than one there is no need for adjustment

}
