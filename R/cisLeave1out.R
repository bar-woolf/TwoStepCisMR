
#' @title cisLeave1out
#' @description This function implements a leave one out analysis using IVW with correlated variants. Requires installation of the meta R package.
#' @param betaXG SNP-exposure associations.
#' @param betaYG SNP-outcome associations.
#' @param sebetaXG Standard error of the SNP-exposure associations.
#' @param sebetaYG Standard error of the SNP-outcome associations
#' @param matrix SNP correlation matrix
#' @param sm Estimator type (see meta R package documentation)
#' @param method fixed (FE) or random effects (RE, default) IVW
#' @param gene colname specifying which gene SNPs come from if stratifying by locus will be useful.
#' @param FirstOrderSE Use first or second order approximations of the standard error. Default is to use first order (i.e. SE = sebetaYG/betaXG)
#' @export





cisLeave1out<-function(SNP,beta.exposure,se.exposure, beta.outcome,se.outcome,  method="RE",sm,gene=NA,matrix){

  dat2<-as.data.frame(cbind(SNP,gene,beta.exposure,se.exposure,beta.outcome,se.outcome))
  colnames(dat2)<-c("SNP","gene","beta.exposure","se.exposure","beta.outcome","se.outcome")
  dat2$beta.exposure<-as.numeric(dat2$beta.exposure)
  dat2$se.exposure<-as.numeric(dat2$se.exposure)
  dat2$beta.outcome<-as.numeric(dat2$beta.outcome)
  dat2$se.outcome<-as.numeric(dat2$se.outcome)
  # matrix<-as.matrix(matrix)

  if (is.na(dat2$gene[1])){
    for (j in 1:nrow(dat2)){
      dat4<-dat2[-j,]
      IVW<-IVWcorrel( method=method, betaYG=dat4$beta.outcome,      sebetaYG=dat4$se.outcome,      betaXG=dat4$beta.exposure,      sebetaXG=dat4$se.exposure,      rho=matrix[dat4$SNP,dat4$SNP]    )
      dat2$LOO_beta[dat2$SNP==dat2$SNP[j]]<-IVW[1]
      dat2$LOO_se[dat2$SNP==dat2$SNP[j]]<-IVW[2]

    }

    dat2$LOO_beta<-as.numeric(dat2$LOO_beta)
    dat2$LOO_se<-as.numeric(dat2$LOO_se)
    a<-meta::forest(meta::metagen(data=dat2,TE=LOO_beta,overall=F, seTE =LOO_se ,  random = F,fixed=F ,studlab = SNP,  sm=sm, tau.common = FALSE),    leftcols = c("studlab"), leftlabs = c(""),smlab="",lwd=1.5,squaresize = .2,  weight.study = "same",col.square="black")#,plotwidth="12cm", rightlabs = c("Odds ratio", "[ 95% CI ]")

  }
  if (!is.na(dat2$gene[1])){

    for (i in unique(dat2$gene)){
      dat3<-dat2[dat2$gene==i,]
      for (j in 1:nrow(dat3)){
        dat4<-dat3[-j,]
        IVW<-IVWcorrel(method=method, betaYG=dat4$beta.outcome,      sebetaYG=dat4$se.outcome,      betaXG=dat4$beta.exposure,      sebetaXG=dat4$se.exposure,      rho=matrix[dat4$SNP,dat4$SNP]    )
        dat2$LOO_beta[dat2$SNP==dat3$SNP[j]]<-IVW[2,1]
        dat2$LOO_se[dat2$SNP==dat3$SNP[j]]<-IVW[2,2]

      }}
    dat2$LOO_beta<-as.numeric(dat2$LOO_beta)
    dat2$LOO_se<-as.numeric(dat2$LOO_se)
    a<-meta::forest(meta::metagen(data=dat2,TE=LOO_beta,overall=F, seTE =LOO_se ,  test.subgroup=F, random = F,fixed=F ,studlab = SNP,  sm=sm, subgroup = gene, tau.common = FALSE),    leftcols = c("studlab"), leftlabs = c(""),smlab="",sep.subgroup = "",subgroup.name ="",lwd=1.5,squaresize = .2,  weight.study = "same",col.square="black")#,plotwidth="12cm", rightlabs = c("Odds ratio", "[ 95% CI ]")
  }
  #return(a)
}

