\name{cpg.work}
\alias{cpg.work}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Does the analysis between the CpG sites and phenotype of interest
}
\description{
Association Analysis Between Methylation Beta Values and Phenotype of Interest. This function contains the code that does the brunt of the work for \code{cpg.assoc} and \code{cpg.perm}.
}
\usage{
cpg.work(beta.values, indep, covariates = NULL, data = NULL, logit.transform = FALSE,
chip.id = NULL, subset = NULL, random = FALSE, fdr.cutoff = 0.05, callarge = FALSE,
fdr.method = "BH", logitperm = FALSE,big.split=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta.values}{
A vector, matrix, or data frame containing the beta values of interest (1 row per CpG site, 1 column per individual).
}
  \item{indep}{
A vector containing the main variable of interest.  \code{cpg.work} will evaluate the association between indep and the beta values. 
}
  \item{covariates}{
A data frame consisting of the covariates of interest. covariates can also be a matrix if it is a model matrix minus the intercept column. 
It can also be a vector if there is only one covariate of interest. Can also be a formula (e.g. \code{~cov1+cov2}).
}
  \item{data}{
an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from the environment from which \code{cpg.work} is called.
}
  \item{logit.transform}{
logical. If \code{TRUE}, the logit transform of the beta values log(beta.val/(1-beta.val)) will be used. Any values equal to zero or one will be set to the 
next smallest or  next largest value respectively; values <0 or >1 will be set to NA.
}
  \item{chip.id}{
An optional vector containing chip or batch identities. If specified, chip id will be included as a factor in the model.
}
  \item{subset}{
an optional logical vector specifying a subset of observations to be used in the fitting process.
}
  \item{random}{
logical. If \code{TRUE}, the \option{chip.id} will be included in the model as a random effect, and a random intercept model will be fitted. If \code{FALSE}, \option{chip.id} will be included
in the model as an ordinary categorical covariate, for a much faster analysis.
}
  \item{fdr.cutoff}{
The threshold at which to compare the FDR values. The default setting is .05. Any FDR values less than .05 will be considered significant.
}
  \item{callarge}{
logical. Used by \code{cpg.assoc} when it calls \code{cpg.work}. If \code{TRUE} it means that beta.values is actually split up from a larger data set and that \code{memory.limit} may be a problem.
This tells \code{cpg.work} to perform more \code{rm()} and \code{gc()} to clear up space.
}
  \item{fdr.method}{
Character.Method used to calculate False Discovery Rate. Can be any of the methods listed in \code{\link{p.adjust}} or "qvalue" 
for John Storey's qvalue method (required to have \pkg{qvalue} package installed). The default method is "BH" for the Benjamini & Hochberg method.
}
  \item{logitperm}{
Passes from \code{cpg.perm} when permutation test is performed. Stops from future checks involving the logistic transformation.
}
  \item{big.split}{
Passes from \code{cpg.assoc}. Internal flag to inform \code{cpg.work} that the large data did not need to be split up.
}
}
\details{
\code{cpg.work} does the analysis between the methylation and the phenotype of interest. It is called by \code{cpg.assoc} to do the brunt of the work. It can be called itself with the same input as \code{cpg.assoc}, it just cannot handle large data sets.
}
\value{
\code{cpg.work} will return an object of class \code{"cpg"}.
The functions summary and plot can be called to get a summary of results and to create QQ plots. The output is in the same order as the original input. To sort it by p-value, use the \code{sort} function.

  \item{results}{
A data frame consisting of the statistics and P-values for each CpG site. Also has the adjusted p-value based on the fdr.method and whether the site was Holm significant.
          }
  \item{Holm.sig}{
A list of sites that met criteria for Holm significance.
            }
  \item{FDR.sig}{
A data.frame of the sites that were FDR significant by the fdr method. 
  }
  \item{info}{
A data frame consisting of the minimum P-value observed, the fdr method used, what the phenotype of interest was, and the number of covariates in the model.
  }
  \item{indep}{
The main phenotype of interest.
  }
  \item{covariates}{
If covariates was non \code{NULL}, the covariates will be included.Otherwise will be \code{NULL}.
  }
  \item{chip}{
If chip.id was non \code{NULL}, the chip will be included. Otherwise will be \code{NULL}.
  }

  \item{coefficients}{
A data frame consisting of the degrees of freedom, and if object is continous the intercept effect adjusted for possible covariates in the model, the estimated effect size, and the standard error. 
The degrees of freedom is used in \code{\link{plot.cpg}} to compute the genomic inflation factors.\cr
 }
}

\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
}



\seealso{
\code{\link{cpg.perm}}
\code{\link{cpg.assoc}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{cpg.combine}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}
\code{\link{sort.cpg.perm}}
\code{\link{sort.cpg}}
\code{\link{cpg.qc}}
\code{\link{cpg.GC}}
}
\examples{
##See the examples listed in cpg.assoc for ways in which to use cpg.work.
##Just change the cpg.assoc to cpg.work.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
