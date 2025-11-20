#' Gene Expression Data for Lymphoblastoid Cell Lines of 210 Unrelated HapMap
#' individuals from four populations
#'
#' A 210 times 3221 matrix with indviduals along the rows and expression data along the columns
#'
#' The first column corresponds to column number 24266 (with name GI_6005726-S) in the original data.
#' Column names give the name of the variables, row names the "name" of the individuals.
#' This is a subset of SangerData where the 3220 last rows are select among all original rows following the same
#' pre-processing procedure as in (section 1.6.1). See also the file  Read_sanger_data.R
#'
#'
#' @docType data
#' @keywords datasets
#' @name SangerData2
#' @usage data(SangerData2)
#' @format A data frame  with 210 rows and 3221 variables
#' @source Dataset downloaded from 
#' \url{https://ftp.sanger.ac.uk/pub/genevar/}
#'
#' References:
#'
#' Stranger, BE et al (2007): Relative impact of nucleotide and copy number variation on gene expression phenotypes
#' Science, 2007•science.org
#' 
#' Bogdan et al (2020): Handbook of Multiple Comparisons, \url{https://arxiv.org/pdf/2011.12154}
#'
NULL
