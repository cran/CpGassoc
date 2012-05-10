plot.cpg <-
function(x,save.plot=NULL,file.type="pdf",popup.pdf=FALSE,tplot=FALSE,classic=TRUE,main.title=NULL,eps.size=c(5,5),...) {
 nam.ind<-as.character(x$info$Phenotype)
  ob=x$results[which(!is.na(x$results$P.value)),]
  u=(1:nrow(ob)-1/2)/nrow(ob)


  gcvalue<-format(median(-log(ob$P.value,base=10))/.30103,digits=3)
  sig<-nrow(x$FDR.sig)      

    if(sig>0) {
      nonsig<-ob$P.value[-(which(ob$P.value %in% x$FDR.sig$P.value))]
      }
    else {
      nonsig<-ob$P.value
      }

   u2=(1:length(nonsig)-1/2)/length(nonsig)
  if(tplot & is.factor(x$indep)) {
    warning("Can not do t-statistic plot with a factor variable\n")
    tplot=FALSE }
  if(!is.null(save.plot)){
    if(!(file.type %in% c("eps","pdf"))) {
       stop("Incorrect file type. Must be pdf or eps\n")
              }
     if(file.type=="eps"){
       postscript(paste(save.plot[1],".eps",sep=""), horizontal = FALSE, 
                 onefile = FALSE, paper = "special",width=eps.size[1],height=eps.size[2])
             }
     if(file.type=="pdf" & !popup.pdf) {
       pdf(file=paste(save.plot[1],".pdf",sep=""))
           }
       }
  if(!classic) {
    u=(1:length(nonsig)-1/2)/length(nonsig)
    pvalues<-nonsig 
              }
  else {
    pvalues<-ob$P.value
        }
  adjusp<-p.adjust(ob$P.value,"holm")
  lu=length(u)
  uu=1:lu
  upper=qbeta(.95,uu,lu-uu+1)
  lower=qbeta(.05,uu,lu-uu+1)
  if(!tplot) {
    ob<-ob$P.value
    if(is.null(main.title)) {
        main.title<-paste("QQ plot for association\nbetween methylation and",nam.ind)
          }
   
    plot(-log(u,base=10),-log(sort(pvalues),base=10),xlab=expression(paste("Expected -log", scriptstyle(10), "(P-values)",sep="")),
            ylab=expression(paste("Observed -log ", scriptstyle(10), "(P-values)",sep="")),main=main.title,ylim=c(0,max(-log(ob,base=10))),cex=pointsizefunction(sort(pvalues)),...)
    legend(0,-log(min(ob),base=10),c("Holm-significant",paste("FDR-significant (",as.character(x$info$FDR.method),")"),"95% confidence interval",paste("Genomic control factor = ", gcvalue,sep="")),lty=c(-1,-1,2,NA),pch=c(19,1,-1,NA),col=c("red","red","black","black"))
    
    abline(a=0,b=1)

    if(sig>0) {
      holm<-sort(ob[which(adjusp<.05)])
      if(nrow(x$Holm.sig) >0){
        if(!classic) {
          points(rep(-log(min(u),base=10),length(holm)),-log(holm,base=10),col="red",pch=19)
            }
        if(classic) {
          points(-log(u[1:length(holm)],base=10),-log(holm,base=10),col="red",pch=19)
            }
        }
      if(nrow(x$Holm.sig)==0 | nrow(x$Holm.sig) < nrow(x$FDR.sig)) {
        fdrsig<-sort(ob[which(adjusp > .05 & ob %in% x$FDR.sig$P.value)])
        if(!classic) {
          points(rep(-log(min(u),base=10),length(fdrsig)),-log(fdrsig,base=10),col="red")
            }
        else {
           points(-log(u[(length(holm)+1):(length(holm)+length(fdrsig))],base=10),-log(fdrsig,base=10),col="red")
              }
      }}
    
    points(-log(u,base=10),-sort(log(upper,base=10)),type='l',lty=2)
    points(-log(u,base=10),-sort(log(lower,base=10),na.last=FALSE),type='l',lty=2)

     }

  if(tplot) {
  k=order(ob$T.statistic)
  df_use<-x$coefficients[k,1]
  t.val<-qt(u2,df_use)
  index3<-ob$T.statistic %in% x$FDR.sig$T.statistic
  adjusp<-adjusp[k]
  tstatistic<-ob$T.statistic[k]
  if(sig>0) { nonsig<-ob$T.statistic[-(which(ob$T.statistic %in% x$FDR.sig$T.statistic))]}
    else {nonsig<-ob$T.statistic}
  index3<-index3[k]
   if(is.null(main.title)) {
        main.title<-paste("Expected vs. Observed T-statistics",
            "for association\nbetween methylation and",x$info$Phenotype)
          }
  k2<-order(nonsig)
  remsig<-which(ob$T.statistic %in% x$FDR.sig$T.statistic)
  if(!classic) {
        tstatistic <-nonsig[k2]
        if(length(remsig)>0) {
          df_use<-x$coefficients[-(remsig),1][k2]        
                 }}
   
  y.val<-ob$T.statistic[k]
  t.val2<-qt(u,df_use)
  upper=qt(upper,df_use)
  lower=qt(lower,df_use)
  
  qqplot(t.val2,tstatistic,xlab="Expected",ylab="Observed",main=main.title,xlim=c(min(t.val),max(t.val)),ylim=c(min(y.val),max(y.val)),...)
  
  legend(min(t.val),max(y.val),c("Holm-significant","Holm-significant",paste("FDR-significant (",as.character(x$info$FDR.method),")"),paste("FDR-significant (",as.character(x$info$FDR.method),")"),
                                    "95% confidence interval",paste("Genomic control factor = ", gcvalue,sep="")),lty=c(-1,-1,-1,-1,2,NA),pch=c(19,19,1,1,-1,NA),col=c("red","green","red","green","black","black"))
  
  if(sig>0) {
      holmreds<-which(adjusp<.05 & (y.val>0) )
      holmgreen<-which(adjusp<.05 & (y.val<0))
      if(nrow(x$Holm.sig) >0){
        if(!classic) {
          points(rep(max(t.val2),length(holmreds)),y.val[holmreds], col="red",pch=19)
          points(rep(min(t.val2),length(holmgreen)),y.val[holmgreen], col="green",pch=19)
           }
        else {
          if(length(holmreds)>0) {points(t.val2[(length(t.val2)-length(holmreds)+1):length(t.val2)],y.val[holmreds], col="red",pch=19)}
          if(length(holmgreen)>0) {points(t.val2[1:length(holmgreen)],y.val[holmgreen], col="green",pch=19)}
           }
        }
      if(nrow(x$Holm.sig)==0 | nrow(x$Holm.sig) < nrow(x$FDR.sig)) {
        sorted<-ob$P.value[k]
        holmlered<-which(adjusp > .05 & (sorted %in% x$FDR.sig$P.value) & y.val >0)
        holmlegreen<-which(adjusp > .05 & sorted %in% x$FDR.sig$P.value & y.val<0)
        if(!classic){
          points(rep(max(t.val2),length(holmlered)),y.val[holmlered], col="red")
          points(rep(min(t.val2),length(holmlegreen)),y.val[holmlegreen], col="green")
          }
        else{
         if(length(holmlered)>0) { points(t.val2[(length(t.val2)-length(holmreds)-length(holmlered)+1):(length(t.val2)-length(holmreds))],y.val[holmlered],col="red") }
         if(length(holmlegreen)>0) { points(t.val2[(1+length(holmgreen)):(length(holmgreen)+length(holmlegreen))],y.val[holmlegreen],col="green")}
            }
      }
       }
  reds<-which((index3==TRUE)&(y.val>0))
  greens<-which((index3==TRUE)&(y.val<0))
  points(t.val2,upper,type='l',lty=2)
  points(t.val2,lower,type='l',lty=2)
  
        }
      if(!is.null(save.plot))  {
      if(file.type=="eps") {
            dev.off()
            }
      else{
        if(popup.pdf) {
          dev.copy2pdf(file=paste(save.plot[1],".pdf",sep=""))
            }
        else{
            dev.off()
            }
      
      }}
      }
