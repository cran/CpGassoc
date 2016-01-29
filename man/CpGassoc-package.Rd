\name{CpGassoc-package}
\alias{CpGassoc-package}
\alias{CpGassoc}
\docType{package}
\title{
Association Between Methylation and a Phenotype of Interest 
}
\description{
Is designed to test for association between methylation at CpG sites across the genome and a phenotype of interest, 
adjusting for any relevant covariates. The package can perform standard analyses of large datasets very quickly with no 
need to impute the data. It can also handle mixed effects models with chip or batch entering the model as a random intercept. 
Also includes tools to apply quality control filters, perform permutation tests, and create QQ plots, manhattan plots, 
and scatterplots for individual CpG sites.

There is a tutorial for CpGassoc available at the website <genetics.emory.edu/conneely/>
}
\details{
\tabular{ll}{
Package: CpGassoc
Type: Package
Title: Association between Methylation and a phenotype of interest
Version: 2.55
Date: 2016-1-28
Author: Barfield, R., Conneely, K., Kilaru,V
Maintainer: R Barfield <rbarfield01@fas.harvard.edu>
Description: CpGassoc is designed to test for association between methylation at CpG sites across the genome and a phenotype of interest, adjusting for any relevant covariates. The package can perform standard analyses of large datasets very quickly with no need to impute the data. It can also handle mixed effects models with chip or batch entering the model as a random intercept. CpGassoc also includes tools to apply quality control filters, perform permutation tests, and create QQ plots, manhattan plots, and scatterplots for individual CpG sites.
Depends:nlme,methods
Suggests:qvalue
License: GPL (>= 2)
}
CpGassoc is a suite of R functions designed to perform flexible analyses of methylation array data. 
The two main functions are  \code{cpg.assoc} and  \code{cpg.perm}.  \code{cpg.assoc} will perform an association test 
between the CpG sites and the phenotype of interest. Covariates can be added to the model, and can be 
continuous or categorical in nature.  \code{cpg.assoc} allows users to set their own false discovery rate threshold,
to transform the beta values to log(beta/(1-beta)), and to subset if required.  cpg.assoc can also fit a 
linear mixed effects model with a single random effect to control for possible technical difference due to 
batch or chip.  \code{cpg.assoc} uses the Holm method to determine significance. The user can also specify 
an FDR method to determine significance based on the function  \code{p.adjust} or  \code{qvalue}.  \code{cpg.perm} 
performs the same tasks as \code{cpg.assoc} followed by a permutation test on the data, repeating the analysis 
multiple times after randomly permuting the main phenotype of interest. The user can 
specify the seed and the number of permutations. If over one hundred permutations are performed 
QQ plots can be created with empirical confidence intervals based on the permuted t-statistics. 
For more information see \code{\link{plot.cpg.perm}}. For more information on how to perform  \code{cpg.assoc} or 
 \code{cpg.perm} see their corresponding help pages. CpGassoc can also perform quality control (see \code{\link{cpg.qc}}).
}
\author{
Barfield, R.; Kilaru,V.; Conneely, K.\cr
Maintainer: R. Barfield: <rbarfield01@fas.harvard.edu>
}


\keyword{ package }
\seealso{
\code{\link{cpg.assoc}}
\code{\link{cpg.combine}}
\code{\link{cpg.perm}}
\code{\link{cpg.work}}
\code{\link{plot.cpg}}
\code{\link{scatterplot}}
\code{\link{manhattan}}
\code{\link{plot.cpg.perm}}
\code{\link{cpg.qc}}
}
\examples{
###NOTE: If you are dealing with large data, do not specify large.data=FALSE. 
###The default option is true
##This will involve partitioning up the data and performing more gc() to clear up space
##Using cpg.assoc:
data(samplecpg,samplepheno,package="CpGassoc")
results<-cpg.assoc(samplecpg,samplepheno$weight,large.data=FALSE)
results

##Using cpg.perm:
Testperm<-cpg.perm(samplecpg[1:200,],samplepheno$weight,data.frame(samplepheno$Dose),
                seed=2314,nperm=10,large.data=FALSE)
Testperm


#For more examples go to those two pages main help pages.
}
