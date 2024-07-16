\name{cpg.combine}
\alias{cpg.combine}

\title{
Combine various objects of class \code{"cpg"}
}
\description{
Takes a list containing objects of class \code{"cpg"} and combines them into one cpg item. Assumes that there are no repeated CpG sites between the various objects (i.e. analysis wasn't performed on the same sites twice).

}
\usage{
cpg.combine(allvalues, fdr.method="BH",fdr.cutoff=.05,return.data=FALSE)
}

\arguments{
  \item{allvalues}{
A list containing the \code{"cpg"} objects that are desired to be consolidated.
}
  \item{fdr.method}{
  FDR method that user wants to use. For options see the \code{cpg.assoc} help page.
  }
  \item{fdr.cutoff}{
  The desired FDR threshold. The default setting is .05. The set of CpG sites with FDR < fdr.cutoff will be labeled as significant.
    }
  \item{return.data}{
Logical. cpg.assoc can return dataframes containing the the variable of interest, covariates, and the chip id (if present). Defaults to FALSE. Set to TRUE if plan on using the downstream scaterrplot functions).
}
}

\value{
\item{info.data}{An object of class \code{"cpg"} that is the consolidated version of the objects of class cpg that were passed in.
  }
}

\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <barfieldrichard8@gmail.com>
}
\note{
This is designed to be used by \code{cpg.assoc} when it does analysis on large data sets or by the user if they split up the analysis by chromosome or some other such partition.
}



\seealso{
\code{\link{cpg.assoc}}
\code{\link{cpg.perm}}
\code{\link{cpg.work}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}
}
\examples{
data(samplecpg,samplepheno,package="CpGassoc")
test1<-cpg.assoc(samplecpg[1:100,],samplepheno$weight,large.data=FALSE)
test2<-cpg.assoc(samplecpg[101:200,],samplepheno$weight,large.data=FALSE)
bigtest<-list(test1,test2)
overall<-cpg.combine(bigtest)
overall
      
}
