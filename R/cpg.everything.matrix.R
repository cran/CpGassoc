cpg.everything.matrix <-
function(x,i.p,levin,chip.id,...) {
      p<-x
      if(levin){
        if(is.null(i.p)) {
          ranmixed<-function(z) {
              zz<-z
              temp<-try(anova(lme(zz~p,random=~1|factor(chip.id),na.action=na.omit)),silent=TRUE)
              temp<-as.matrix(temp)
              if(length(temp)<2) {
                info<-matrix(NA,3)}
              else {
                if(temp[nrow(temp),4]==0) {
                  info<-c(temp[nrow(temp),3],df(temp[nrow(temp),3],temp[nrow(temp),1],temp[nrow(temp),2]))
                      }
                else { 
                   info<-temp[nrow(temp),3:4]
                      }
                   info<-c(info,temp[nrow(temp),2])
                   }
              info
              }}
        else {
          ranmixed<-function(z) {
              zz<-z
              temp<-try(anova(lme(zz~p+i.p,random=~1|factor(chip.id),na.action=na.omit)),silent=TRUE)
              temp<-as.matrix(temp)
              if(length(temp)<2) {
                info<-matrix(NA,3)}
              else {
                if(temp[nrow(temp),4]==0) {
                  info<-c(temp[nrow(temp),3],df(temp[nrow(temp),3],temp[nrow(temp),1],temp[nrow(temp),2]))
                }
                else { 
                   info<-temp[nrow(temp),3:4]
                      }
                info<-c(info,temp[nrow(temp),2])
          }
          info
                } }}
    else {
       ranmixed<-function(z) {
       zz<-z
       temp<-try(summary(lme(zz~p,random=~1|factor(chip.id),na.action=na.omit))[[19]],silent=TRUE)
      if(length(temp)>1) {
        info <- temp[2,c(4:5,3)]
     if(is.null(ncol(p))) {
       info<-c(info,temp[1:2,1])
        }
      else{ 
        betas<-temp[-2,1]
        altdesign<-as.matrix(p[!is.na(zz),2:ncol(p)])
        thecolmeans<-colMeans(altdesign)
        thecolmeans<-thecolmeans[colSums(altdesign)!=0]
        newint<-betas[1]+sum(thecolmeans*betas[2:length(betas)])
        info<-c(info,newint,temp[2,1])
                    } }
      else {
           info<-matrix(NA,5)}
      info
        }}
      return(ranmixed)}
