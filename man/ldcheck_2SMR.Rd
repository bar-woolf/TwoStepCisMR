% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ldcheck_2SMR.R
\name{ldcheck_2SMR}
\alias{ldcheck_2SMR}
\title{ldcheck_2SMR}
\usage{
ldcheck_2SMR(dat, candidate_snp, snp_EA_OA, SNPrank = 30, ldr2 = 0.8)
}
\arguments{
\item{dat}{a two-sample MR data frame}

\item{candidate_snp}{candidate lead SNP rsid}

\item{snp_EA_OA}{candidate SNP formatted as rsid_effectallele_otherallele}

\item{SNPrank}{The number of lead SNPs to include in the analysis. Put as NA to instead check for snpss p < 1e-5}

\item{ldr2}{minimum LD r2 needed to pass LD check.}
}
\description{
this function performs ld check, see https://www.nature.com/articles/s41588-020-0682-6#Sec30
}
