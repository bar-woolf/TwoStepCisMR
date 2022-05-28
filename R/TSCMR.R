#' @title Two-Step Cis MR: A novel method and R package for attenuating bias in cis-MR in the absence of uncorrelated variants
#' @description This function implements a simple adjsutment to the snp-outcome association to remove bias from cis-MR estimates in the absence of uncorrelated variants.
#' @param Bgc SNP-confounder associations.
#' @param Bgo SNP-outcome associations.
#' @param SEgc Standard error of the SNP-confounder associations.
#' @param SEgo Standard error of the SNP-outcome associations.
#' @param Bco confounder-outcome association.
#' @param SEco Standard error of the confounder-outcome association.
#' @param straps Number of iterations for the bootstrap standard error
#' @keywords 2SMR
#' @export



TSCMR<-function(Bgo, SEgo,Bgc, SEgc,Bco, SEco,straps=1e6){

  #adjusted gene-outcome association
  Bgo2<-Bgo-(Bgc*Bco)

  #boostrap se
  set.seed(496)
  Bgo_boot = rnorm(straps, Bgo, SEgo)
  Bgc_boot = rnorm(straps, Bgc, SEgc)
  Bco_boot = rnorm(straps, Bco, SEco)
  indirect = Bgo_boot- (Bgc_boot*Bco_boot)
  #quantile(indirect, prob=c(0.025, 0.975)) #this gives you the 95% CI
  BSSE<-sd(indirect) # this gives you the boot strapped SE

  #error propgagated SE
  EPSE<- (SEgo^2 + ((Bgc*Bco)^2)*(  (SEgc/Bgc)^2 + (SEco/Bco)^2    )   )^0.5

  out<-as.data.frame(1)
  out$Bgo<-Bgo2
  out$BSSE<-BSSE
  out$EPSE<-EPSE
  out<-out[,-1]
  return(out)
}
