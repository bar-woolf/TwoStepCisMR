

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


# #The IVW method accounting for correlation can also be performed using the standard linear regression command after weighting the data by the Cholesky decomposition:
# #n.b. I've not made this function properlly 
# #n.b. work out how to convert it to MR-Egger
# IVWcorrel<-function(betaYG,sebetaYG,betaXG,sebetaXG,rho){
#   
# 
# Omega = sebetaYG%o%sebetaYG*rho
# c_betaXG = solve(t(chol(Omega)))%*%betaXG
# c_betaYG = solve(t(chol(Omega)))%*%betaYG
# beta_IVWcorrel = lm(c_betaYG~c_betaXG-1)$coef[1]
# se_IVWcorrel.fixed = sqrt(1/(t(betaXG)%*%solve(Omega)%*%betaXG))
# se_IVWcorrel.random = sqrt(1/(t(betaXG)%*%solve(Omega)%*%betaXG))*max(summary(lm(c_betaYG~c_betaXG-1))$sigma,1)
# 
# return(matrix(
#   c("beta_IVWcorrel",beta_IVWcorrel,
#     "se_IVWcorrel.fixed",se_IVWcorrel.fixed,
#     "p_IVWcorrel.fixed",  2*pnorm(abs(beta_IVWcorrel/se_IVWcorrel.fixed),low=F),
#     "n_snp", length(betaYG),
#     "F_stat", (solve(t(sebetaXG^2)%*%solve(rho)%*%sebetaXG^2)*t(sebetaXG^2)%*%solve(rho)%*%betaXG^2)
#   ),nrow=2,ncol=5))
# }



# #IVW estimate (accounting for correlation) using principal components:
# Phi = (betaXG/sebetaYG)%o%(betaXG/sebetaYG)*rho
# summary(prcomp(Phi, scale=FALSE))
# K = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2/sum((prcomp(Phi, scale=FALSE)$sdev^2)))>0.99)[1]
# # K is number of principal components to include in analysis
# # this code includes principal components to explain 99% of variance in the risk factor
# betaXG0 = as.numeric(betaXG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
# betaYG0 = as.numeric(betaYG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K])
# Omega = sebetaYG%o%sebetaYG*rho
# pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]
# beta_IVWcorrel.pc = solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0
# se_IVWcorrel.fixed.pc = sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0))




# theis function performs ld check, see https://www.nature.com/articles/s41588-020-0682-6#Sec30 to support robustness of coloc when there are fewer than 50 snps
# n.b. this function assumes you have the data formated using the df function 
#follwoing Chris' advice, either implemetning using top 30 snps or snps with p less than 1e-5

# ldcheck<-function(dat, snp_EA_OA,SNPrank=NA,ldr2=0.8){ 
#   
#     if (is.na(SNPrank)==T){
#     eqtl2<-dat[dat$pvalues<1e-5,c("snp","pvalues")] #ordering eqtl
#     names(eqtl2)<-c("SNP","pval.exposure")
#     } else if (is.na(SNPrank)==F){
#       eqtl2<-dat[order(dat$pvalues),c("snp","pvalues")] #ordering eqtl
#       names(eqtl2)<-c("SNP","pval.exposure")
#       eqtl2$rank<-1:nrow(eqtl2)
#       eqtl2<-eqtl2[eqtl2$rank<(SNPrank+1),]
#     }
#   
# 
# # nsps<-nrow(eqtl2)
#   
#   if (hyprcoloc_results$results$candidate_snp %in% eqtl2$snp){
#     eqtl2[eqtl2$snp==hyprcoloc_results$results$candidate_snp]$pval.exposure<-0
#   } else if (!(hyprcoloc_results$results$candidate_snp %in% eqtl2$snp)) {
#     eqtl2<-rbind(eqtl2,c(hyprcoloc_results$results$candidate_snp,0,0) )
# }
#     ld<-ld_matrix(eqtl2$SNP)
#     ld<-ld[,snp_EA_OA]
#     nsps2<-length(ld[abs(ld)<ldr2])
#     return( matrix(c("n_ld_SNPs",length(ld)-nsps2,"n_inc_SNPs",length(ld), "prop_SNPs",((length(ld)-nsps2)/length(ld))),2,3) )
# } #returns the number of snps in LD with teh cuasal SNP, the number of SNPs included in the anlaysis, and the proportion of these that are in LD

#' @title ldcheck_2SMR
#' @description this function performs ld check, see https://www.nature.com/articles/s41588-020-0682-6#Sec30 
#' @param dat a two-sample MR data frame
#' @param candidate_snp candidate lead SNP rsid
#' @param snp_EA_OA candidate SNP formatted as rsid_effectallele_otherallele
#' @param SNPrank  The number of lead SNPs to include in the analysis. Put as NA to instead check for snpss p < 1e-5
#' @param ldr2 minimum LD r2 needed to pass LD check. 
#' @export

ldcheck_2SMR<-function(dat,candidate_snp,snp_EA_OA,SNPrank=30,ldr2=0.8){ 
  eqtl2<-dat[order(dat[,"pval.exposure"]),c("SNP","pval.exposure")] #ordering eqtl and subsettign to only the cols used
  if (SNPrank==0){
    eqtl2<-eqtl2[eqtl2[,"pval.exposure"]<1e-5,] #ordering eqtl
  } 
  if (SNPrank>0){
    eqtl2<-eqtl2[1:SNPrank+1,] #removing rows bellow rank of interest
  }  
  
 # nsps<-nrow(eqtl2)
  
  eqtl2<-eqtl2$SNP
  inc<-"F"
if (!(candidate_snp %in% eqtl2)) {
    eqtl2<-c(eqtl2,candidate_snp) #c(candidate_snp,0,0) ) )
    inc<-"T"
  }
  ld<-ld_matrix(eqtl2)#[,"SNP"])
  ld<-ld[,snp_EA_OA]
  nsps2<-length(ld[abs(ld)<ldr2])
  return( matrix(c("n_ld_SNPs",length(ld)-nsps2,"n_inc_SNPs",length(ld), "prop_SNPs",((length(ld)-nsps2)/length(ld)), "inc",inc),2,4) )
}


# #this function creates plots for gwasglue
# 
# plots_gwasglue<-function(hyper){
#   dat2<-df(hyper[2])
#  dat1<-df(hyper[1])
#     p1 <- ggplot(dat1, aes(x = pos, y = -log10(pvalues))) +
#       geom_point() +
#       labs(x = "", y = bquote(-log[10](italic(p))),
#            title = id[1]) +
#       theme_bw() +
#       theme(plot.title = element_text(hjust = 0.5)) +
#       scale_x_continuous(                       breaks = c( gene_start - window,gene_start - window/2,gene_start, (gene_start+gene_end)/2, gene_end, gene_end + window/2, gene_end, gene_end + window),
#                                                 labels = c( gene_start - window,gene_start - window/2,gene_start, (gene_start+gene_end)/2, gene_end, gene_end + window/2, gene_end, gene_end + window))
# 
# 
#   p4 <- ggplot(dat2, aes(x =  pos, y = -log10(pvalues))) +
#     geom_point() +
#     labs(x = "", y = bquote(-log[10](italic(p))),
#          title = lab2) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     scale_x_continuous(name = paste("Chromosome ",chr[1]," Position",sep=""),
#                        breaks = c( gene_start - window,gene_start - window/2,gene_start, (gene_start+gene_end)/2, gene_end, gene_end + window/2, gene_end, gene_end + window),
#                        labels = c( gene_start - window,gene_start - window/2,gene_start, (gene_start+gene_end)/2, gene_end, gene_end + window/2, gene_end, gene_end + window))
# 
#   fig<-ggpubr::ggarrange(p1,  p4,
#                          heights = c(1, 1), nrow = 2,
#                          ncol = 1, align = "hv")
# 
#   return(fig)
# }
















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


















#this rusn MR Egger while accounting for SNP correlations 
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