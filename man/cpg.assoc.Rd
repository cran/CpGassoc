\name{cpg.assoc}
\alias{cpg.assoc}

\title{
Association Analysis Between Methylation Beta Values and Phenotype of Interest
}
\description{
Association Analysis Between Methylation Beta Values and Phenotype of Interest. For more detail see Tutorial.
}
\usage{
cpg.assoc(beta.val, indep, covariates = NULL, data = NULL, logit.transform = FALSE, 
chip.id = NULL, subset = NULL, random = FALSE, fdr.cutoff = 0.05, large.data = TRUE,
fdr.method = "BH", logitperm = FALSE)
}

\arguments{
  \item{beta.val}{
A vector, matrix, or data frame containing the beta values of interest (1 row per CpG site, 1 column per individual).
}
  \item{indep}{
A vector containing the variable to be tested for association. \code{cpg.assoc} will evaluate the association between the beta values (dependent
variable) and indep (independent variable).
 }
  \item{covariates}{
A data frame consisting of additional covariates to be included in the model. covariates can also be specified as a matrix 
if it takes the form of a model matrix with no intercept column, or can be specified as a vector if there is
only one covariate of interest. Can also be a formula(e.g. \code{~cov1+cov2}).
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the
variables in the model. If not found in data, the variables are taken from the environment from which cpg.assoc is
called.
}
\item{logit.transform}{
logical. If \code{TRUE}, the logit transform of the beta values log(beta.val/(1-beta.val)) will be used. Any values equal to zero or one will be set to the 
next smallest or  next largest value respectively; values <0 or >1 will be set to NA.
}
  \item{chip.id}{
An optional vector containing chip, batch identities, or other categorical factor of interest to the researcher. If specified, chip id will be included as a factor in the model.
}
  \item{subset}{
An optional logical vector specifying a subset of observations to be used in the fitting process.
}
  \item{random}{
Logical. If \code{TRUE}, \option{chip.id} will be included in the model as a random effect, and a random intercept model will be
fitted. If \code{FALSE}, \option{chip.id} will be included in the model as an ordinary categorical covariate, for a much faster analysis.
}
  \item{fdr.cutoff}{
The desired FDR threshold. The default setting is .05. The set of CpG sites with FDR < \option{fdr.cutoff} will be labeled as significant.
}
  \item{large.data}{
Logical.  Enables analyses of large datasets.  When \code{large.data=TRUE}, \code{cpg.assoc} avoids memory problems by performing the analysis in chunks. 
}
  \item{fdr.method}{
Character. Method used to calculate False Discovery Rate. Choices include any of the methods available in \code{p.adjust()} or
"qvalue" for John Storey's qvalue method (requires that \pkg{qvalue} package is installed). The default method is "BH" for the Benjamini & Hochberg method.
}
  \item{logitperm}{
Logical. For internal use only.
}
}
\details{
\code{cpg.assoc} is designed to test for association between an independent variable and methylation at a number of CpG sites, with the option to include additional covariates and factors. \cr
\code{cpg.assoc} assesses significance with the Holm (step-down Bonferroni) and FDR methods. \cr


If \code{class(indep)='factor'}, \code{cpg.assoc} will perform an ANOVA test of the variable conditional on the covariates specified. Covariates, if entered, should be in the form of a data 
frame, matrix, or vector. For example, \code{covariates=data.frame(weight,age,factor(city))}. The data frame can also be specified prior to calling \code{cpg.assoc}.  The covariates 
should either be vectors or columns of a matrix or data.frame.\cr
\code{cpg.assoc} is also designed to deal with large data sets. Setting \code{large.data=TRUE} will make \code{cpg.assoc} split up the data to enable efficient analysis of large datasets.

}
\value{
\code{cpg.assoc} will return an object of class \code{"cpg"}. The functions \code{summary} and \code{plot} can be called to get a summary of results and to create QQ plots.
  \item{results}{
A data frame consisting of the t or F statistics and P-values for each CpG site, as well as indicators of Holm and FDR significance. CpG sites will be in the same order as the original
input, but the \code{sort()} function can be used directly on the \code{cpg.assoc} object to sort CpG sites by p-value.
          }
  \item{Holm.sig}{
A list of sites that met criteria for Holm significance.
            }
  \item{FDR.sig}{
A data.frame of the CpG sites that were significant by the FDR method specified.
  }
  \item{info}{
A data frame consisting of the minimum P-value observed, the FDR method that was used, the phenotype of interest, the number of covariates in the model, the name of the matrix or data frame the methylation beta values were
taken from, the FDR cutoff value and whether a mixed effects analysis was performed.
  }
  \item{indep}{
The independent variable that was tested for association.
  }
  \item{covariates}{
Data.frame or matrix of covariates, if specified (otherwise \code{NULL}).
  }
  \item{chip}{
chip.id vector, if specified (otherwise \code{NULL}).
  }
  \item{coefficients}{
A data frame consisting of the degrees of freedom, and if object is continous the intercept effect adjusted for possible covariates in the model, the estimated effect size, and the standard error. 
The degrees of freedom is used in \code{\link{plot.cpg}} to compute the genomic inflation factors.\cr

  }
}

\author{
Barfield, R.; Conneely, K.; Kilaru,V.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
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
\code{\link{cpg.GC}}
}
\examples{
 # Sample output from CpGassoc
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. 
###The default option is true.
##This will involve partitioning up the data and performing more gc() to clear up space
data(samplecpg,samplepheno,package="CpGassoc")
results<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
results
#Analysis with covariates. There are multiple ways to do this. One can define the
#dataframe prior or do it in the function call.
test<-cpg.assoc(samplecpg,samplepheno$weight,data.frame(samplepheno$Distance,
samplepheno$Dose),large.data=FALSE)
# or
covar<-data.frame(samplepheno$Distance,samplepheno$Dose)
test2<-cpg.assoc(samplecpg,samplepheno$weight,covar,large.data=FALSE)


#Doing a mixed effects model. This does take more time, so we will do a subset of
#the samplecpg
randtest<-cpg.assoc(samplecpg[1:10,],samplepheno$weight,chip.id=samplepheno$chip,
random=TRUE,large.data=FALSE)
          }

