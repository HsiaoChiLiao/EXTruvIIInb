#' The Triana Dataset
#'
#' A dataset containing a subset of QC-ed matched single-cell transcriptomics and proteomics data from Triana et al. (2021, Nature Immunology).
#'
#' @format
#' An R list object consisting of an ordered collection of six components. 
#' \describe{
#'   \item{rna_raw}{The count matrix of 9663 cells on the rows and 422 targeted mRNA features on the columns.}
#'   \item{adt_raw}{The count matrix of 9663 cells on the rows and 97 protein (ADT) features on the columns.}
#'   \item{metadata}{The metadata of each cell with sample information (sample_class and sample_batch columns), and cell type annotation (cellType_broad and cellType_refined).}
#'   \item{nc1}{The negative control set for inferring W1.}
#'   \item{nc2}{The negative control set for inferring W2.}
#'   \item{nc3}{The negative control set for inferring W3.}
#' }
#' @source <https://ega-archive.org/datasets/EGAD00001008188>
"triana_data"

