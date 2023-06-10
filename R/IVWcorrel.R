
# this function does regular IVW accounting for the correlation matrix
# code mostly taken from DOI: 10.1002/gepi.22077. Eleanor Sanderson  checked F stat


#' @title IVWcorrel
#' @description This function implements IVW with correlated variants. Code taken from DOI: 10.1002/gepi.22077
#' @param betaXG SNP-exposure associations.
#' @param betaYG SNP-outcome associations.
#' @param sebetaXG Standard error of the SNP-exposure associations.
#' @param sebetaYG Standard error of the SNP-outcome associations
#' @param rho SNP correlation matrix
#' @param method fixed (FE) or random effects (RE, default)
#' @param FirstOrderSE Use first or second order approximations of the standard error. Default is to use first order (i.e. SE = sebetaYG/betaXG)
#' @export

IVWcorrel<-function(betaYG,sebetaYG,betaXG,sebetaXG,rho, method="RE", FirstOrderSE=T){

  if(FirstOrderSE ==F ){ # if you want to use 2nd order standard error
    sebetaYG<-(sebetaXG^2 + sebetaYG^2)^0.5
  }

  if (length(betaYG)==1) {
    out<-matrix( c(betaYG/betaXG, abs(sebetaYG/betaXG),  2*pnorm(abs((betaYG/betaXG)/(sebetaYG/betaXG)),low=F), 1,   (betaXG^2/sebetaXG^2)      ),nrow=1,ncol=5)
    colnames(out)<-c("beta_IVWcorrel","se_IVWcorrel.random","p_IVWcorrel.random","n_snp","F_stat")
    return(out)

  }else if (length(betaYG)>1) {
    if (method == "FE"|method == "Fixed"|method == "fixed") {
      #IVW accouning for the correlation in a GLM in a generalised weigthed lm
      Omega = sebetaYG%o%sebetaYG*rho
      beta_IVWcorrel = solve(t(betaXG)%*%solve(Omega)%*%betaXG)*t(betaXG)%*%solve(Omega)%*%betaYG
      se_IVWcorrel.fixed = sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))
      out<-matrix(c(beta_IVWcorrel,  se_IVWcorrel.fixed,  2*pnorm(abs(beta_IVWcorrel/se_IVWcorrel.fixed),low=F), length(betaYG), (solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2)    ),nrow=1,ncol=5)
      colnames(out)<-c("beta_IVWcorrel","se_IVWcorrel.fixed","p_IVWcorrel.fixed","n_snp","F_stat")
      return(out)
    }

    if (method == "RE"|method =="random"|method =="Random") {
      Omega = sebetaYG%o%sebetaYG*rho
      beta_IVWcorrel = solve(t(betaXG)%*%solve(Omega)%*%betaXG)*t(betaXG)%*%solve(Omega)%*%betaYG
      resid = betaYG-beta_IVWcorrel[1,1]*betaXG
      se_IVWcorrel.random = sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))*  max(sqrt(t(resid)%*%solve(Omega)%*%resid/(length(betaXG)-1)),1)
      out<- matrix(c(beta_IVWcorrel,se_IVWcorrel.random,  2*pnorm(abs(beta_IVWcorrel/se_IVWcorrel.random),low=F), length(betaYG), (solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2)   ),nrow=1,ncol=5)
      colnames(out)<-c("beta_IVWcorrel","se_IVWcorrel.random","p_IVWcorrel.random", "n_snp","F_stat")
      return(out)
    }
    if (method == "PCA") {
      # #IVW estimate (accounting for correlation) using principal components:
      Phi = (betaXG/sebetaYG)%o%(betaXG/sebetaYG)*rho
      summary(prcomp(Phi, scale=FALSE))
      K = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2/sum((prcomp(Phi, scale=FALSE)$sdev^2)))>0.99)[1]
      # K is number of principal components to include in analysis
      # this code includes principal components to explain 99% of variance in the risk factor
      betaXG0 = as.numeric(betaXG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
      betaYG0 = as.numeric(betaYG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
      Omega = sebetaYG%o%sebetaYG*rho
      pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]
      beta_IVWcorrel.pc = solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0
      se_IVWcorrel.fixed.pc = sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0))

      out<-matrix(c(beta_IVWcorrel.pc,se_IVWcorrel.fixed.pc, 2*pnorm(abs(beta_IVWcorrel.pc/se_IVWcorrel.fixed.pc),low=F),length(betaYG),(solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2) ),nrow=1,ncol=5)
      colnames(out)<-c("beta_IVWcorrel","se_IVWcorrel.fixed", "p_IVWcorrel.fixed", "n_snp",   "F_stat")
      return(out)
    } } }

