cpg.qc <-
function(beta.orig,siga,sigb,pval,p.cutoff=.001,cpg.miss=NULL,sample.miss=NULL,constant100=FALSE) {

#Get measures of mean intensity, use this to flag low-signal individuals for removal
    avg_sig=colMeans(siga+sigb,na.rm=TRUE)
    experimentwide_median_sig=median(avg_sig)
    flag_to_remove=(avg_sig < 0.5*experimentwide_median_sig | avg_sig < 2000)
    gc()

#Remove low-signal individuals
    print(paste("Removed",sum(flag_to_remove),"samples with low signal"))
    siga=siga[,!flag_to_remove]
    sigb=sigb[,!flag_to_remove]

#Create new beta value with no 100 in denominator (unless the 100 is requested by the user)
    const<-ifelse(constant100,100,0)
    beta.new = sigb/(siga+sigb+const)

#Clean up large files to make way for other files
    rm(siga,sigb)
    gc()

#Remove low-signal individuals from original beta value file
    beta.orig=beta.orig[,!flag_to_remove]

    beta.new<-as.matrix(beta.new)
#Use beta.orig values to set bad data points to missing, then discard
#(Bad data points are those set to missing by GenomeStudio, usually due to 0 signal or Nbead < 3)
    beta.new[is.na(beta.orig)] = NA
    rm(beta.orig)
    gc()

#Remove low-signal individuals
    pval=pval[,!flag_to_remove]

#Set to missing any values with detection p-values > cutoff
    beta.new[pval>p.cutoff] = NA
    rm(pval)
    gc()

#Optional checks removing any sample or sites with missing rate above a user-specified cutoff
    if(!is.null(sample.miss) | !is.null(cpg.miss)) {
      missing.location<-which(is.na(beta.new))
      sample.missing<- (missing.location %/% nrow(beta.new)) +1
      cpg.missing<-missing.location %% nrow(beta.new)
      cpg.missing[which(cpg.missing==0)]<-nrow(beta.new)
      cpg.missing.table<-table(cpg.missing)/ncol(beta.new)
      sample.missing.table<-table(sample.missing)/nrow(beta.new)

      if(!is.null(cpg.miss)){
      	remove.cpg<-which(cpg.missing.table>cpg.miss)
      	remove.cpg<-as.numeric(names(remove.cpg))
   	    if(length(remove.cpg)>0) {
		      beta.new<-beta.new[-remove.cpg,]
      	}
        print(paste("Removed",length(remove.cpg),"CpG sites with missing data for >",cpg.miss,"of samples"))
      	}
      if(!is.null(sample.miss)) {
        remove.samp<-which(sample.missing.table>sample.miss)
      	remove.samp<-as.numeric(names(remove.samp))
      	if(length(remove.samp)>0) {
      		beta.new<-beta.new[,-remove.samp]
      	}
        print(paste("Removed",length(remove.samp),"samples with missing data for >",sample.miss,"of CpG sites"))
        }
    }
    beta.new
    }
