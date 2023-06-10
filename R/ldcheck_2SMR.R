
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


