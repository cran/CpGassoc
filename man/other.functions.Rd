\name{Other CpGassoc Functions}
\alias{cpg.length}
\alias{pointsizefunction}
\alias{cpgassocsummary}


\title{
Information on miscellaneous other functions
}
\description{
\code{cpg.length} compares the dimensions of the covariates, the independent phenotype, chip, and the matrix of beta values from
\code{cpg.assoc} or \code{cpg.perm}. If the number of individuals does not match up. Stops the code.

\code{pointsizefunction} simply scales the size of the points for the qq-plot for \code{plot.cpg} and \code{plot.cpg.perm}

\code{cpgassocsummary} is used by \code{cpg.work} when the phenotypes contains such small amount of variance such that matrix methods
can not be used. Takes an object of class \code{"aov"} or \code{"mlm"} and gets the test statistics from these objects while avoiding doing a loop
or the list structures which is typically done if one does summary on one of these objects. 
}
\usage{
cpg.length(indep, numpatients, covariates, chip.id)

pointsizefunction(x)

cpgassocsummary(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indep}{
object \code{indep} from \code{cpg.assoc}, \code{cpg.work} or \code{cpg.perm}.
}
  \item{numpatients}{
The number of column of the object beta.values in \code{cpg.work} or \code{cpg.perm} or number of columns of beta.val in \code{cpg.assoc}

}
  \item{covariates}{
object covariate from \code{cpg.assoc}, \code{cpg.work} or \code{cpg.perm}.
}
  \item{chip.id}{
object \code{chip.id} from \code{cpg.assoc}, \code{cpg.work} or \code{cpg.perm}.
}

  \item{x}{
the p-values from the object \code{"cpg"} or \code{"cpg.perm"} 
  }
  \item{object}{
An object of \code{"aov"} or \code{"mlm"}
  }
  
}
\details{
\code{cpg.length} stops the functions if the dimensions do not add up. Used Internally by \pkg{CpGassoc}

\code{pointsizefunction} simply returns a vector of equal length of x to be used in plot for the cex option.

\code{cpgassocsummary} returns a matrix of columns two and three of the \code{results} value of an \code{"cpg"} or \code{"cpg.perm"} object. This is for the cpg sites with non-missing data
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cpg.perm}}
\code{\link{cpg.assoc}}
\code{\link{scatterplot}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}

}