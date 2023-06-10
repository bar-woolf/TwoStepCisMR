
#' @title F_cor
#' @description This function estimates the F stat with correlated variants. If N is specified then it will also return the r^2
#' @param betaXG SNP-exposure associations.
#' @param sebetaXG Standard error of the SNP-exposure associations.
#' @param rho SNP correlation matrix
#' @param N  The sample size of the GWAS
#' @export

#code for F stat (not R2) checked by E Salndlerson
F_cor<-function(betaXG,sebetaXG,rho,N=NA){ #N is sample size of GWAS. If this is provided then the function returns the r^2 of the instrument too
  if (is.na(N)){
    return(solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2)}
  if (!is.na(NA)){
    f<-solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2
    r2<-f/(N-2+f)
    matrix<-matrix(c(f,r2),1,2)
    colnames(matrix)<-c("F", "r2" )
    return(matrix)
  }
}

