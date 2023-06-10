
#this runs MR Egger while accounting for SNP correlations
#method taken from supplement of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5506233/

#' @title EGGERcorrel
#' @description This function implements MR Egger with correlated variants. Code taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5506233/
#' @param betaXG SNP-exposure associations.
#' @param betaYG SNP-outcome associations.
#' @param sebetaXG Standard error of the SNP-exposure associations.
#' @param sebetaYG Standard error of the SNP-outcome associations
#' @param rho SNP correlation matrix
#' @param method method for accounting for SNP correlations. Can be set to "GWLR" (generalized weighted linear regression, the default) or "LRwCd"  (linear regression command after weighting the data by a Cholesky decomposition)
#' @param FirstOrderSE Use first or second order approximations of the standard error. Default is to use first order (i.e. SE = sebetaYG/betaXG)
#' @export


EGGERcorrel<-function(betaYG,sebetaYG,betaXG,sebetaXG,rho, method="GWLR", FirstOrderSE=T){

  # if(FirstOrderSE ==F ){ # if you want to use 2nd order standard error
  #   sebetaYG<-(sebetaXG^2 + sebetaYG^2)^0.5
  # }

  if (length(betaYG)<4) {
    return(warning("Error: MR Egger needs more than 3 SNPs"))

  }else if (length(betaYG)>3) {
    if (method == "GWLR") {# generalized weighted linear regression

      Omega = sebetaYG%o%sebetaYG*rho
      betaIVW.correl = solve(t(betaXG)%*%solve(Omega)%*%betaXG)*t(betaXG)%*%solve(Omega)%*%betaYG
      seIVW.correl.fixed = sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))
      residIVW = betaYG-as.vector(betaIVW.correl)*betaXG
      sebetaIVW.correl.random = sqrt(solve(t(betaXG)%*%solve(Omega)%*%betaXG))*
        max(sqrt(t(residIVW)%*%solve(Omega)%*%residIVW/(length(betaXG)-1)),1)
      betaEGGER.correl = solve(t(cbind(rep(1, length(betaXG)), betaXG))%*%solve(Omega)%*%
                                 cbind(rep(1, length(betaXG)), betaXG))%*%
        t(cbind(rep(1, length(betaXG)), betaXG))%*%solve(Omega)%*%betaYG
      # first component is intercept term, second component is slope term (causal estimate)
      residEGGER = betaYG-betaEGGER.correl[1]-betaEGGER.correl[2]*betaXG
      varEGGER.correl.random = solve(t(cbind(rep(1, length(betaXG)), betaXG))%*%solve(Omega)%*%
                                       cbind(rep(1, length(betaXG)), betaXG))*
        max(sqrt(t(residEGGER)%*%solve(Omega)%*%residEGGER/(length(betaXG)-2)),1)
      seinterEGGER.correl.random = sqrt(varEGGER.correl.random[1,1])
      sebetaEGGER.correl.random = sqrt(varEGGER.correl.random[2,2])

      out<-matrix(c(betaEGGER.correl[2],betaEGGER.correl[1],
                    sebetaEGGER.correl.random, seinterEGGER.correl.random,
                    2*pnorm(abs(betaEGGER.correl[2]/sebetaEGGER.correl.random),low=F), 2*pnorm(abs(betaEGGER.correl[1]/seinterEGGER.correl.random),low=F),
                    length(betaYG), NA,
                    (solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2),NA
      ),nrow=2,ncol=5)
      colnames(out)<- c("beta.random", "se.random",  "p.random",  "n_snp",   "F_stat")
      rownames(out)<-c("effect estimate","intercept")
      return(out)

    }

    if (method == "LRwCd") {#  linear regression command after weighting the data by a Cholesky decomposition:

      Omega = sebetaYG%o%sebetaYG*rho
      c_betaXG = solve(t(chol(Omega)))%*%betaXG
      c_betaYG = solve(t(chol(Omega)))%*%betaYG
      c_inter = solve(t(chol(Omega)))%*%rep(1, length(betaXG))
      betaIVW.correl = lm(c_betaYG~c_betaXG-1)$coef[1]
      sebetaIVW.correl.fixed = sqrt(1/(t(betaXG)%*%solve(Omega)%*%betaXG))
      sebetaIVW.correl.random = sqrt(1/(t(betaXG)%*%solve(Omega)%*%betaXG))*
        max(summary(lm(c_betaYG~c_betaXG-1))$sigma,1)
      interEGGER.correl = lm(c_betaYG~c_inter+c_betaXG-1)$coef[1]
      betaEGGER.correl = lm(c_betaYG~c_inter+c_betaXG-1)$coef[2]
      seinterEGGER.correl.random = sqrt(solve(t(cbind(rep(1, length(betaXG)), betaXG))%*%solve(Omega)%*%
                                                cbind(rep(1, length(betaXG)), betaXG))[1,1])*
        max(summary(lm(c_betaYG~c_inter+c_betaXG-1))$sigma,1)
      sebetaEGGER.correl.random = sqrt(solve(t(cbind(rep(1, length(betaXG)), betaXG))%*%solve(Omega)%*%
                                               cbind(rep(1, length(betaXG)), betaXG))[2,2])*
        max(summary(lm(c_betaYG~c_inter+c_betaXG-1))$sigma,1)



      out<-matrix(c(betaEGGER.correl,interEGGER.correl,
                    sebetaEGGER.correl.random, seinterEGGER.correl.random,
                    2*pnorm(abs(betaEGGER.correl/sebetaEGGER.correl.random),low=F), 2*pnorm(abs(interEGGER.correl/seinterEGGER.correl.random),low=F),
                    length(betaYG), NA,
                    (solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2),NA
      ),nrow=2,ncol=5)
      colnames(out)<- c("beta.random", "se.random",  "p.random",  "n_snp",   "F_stat")
      rownames(out)<-c("effect estimate","intercept")
      return(out)

    }


  } }

