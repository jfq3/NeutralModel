#' rplB gene count data from a pouchitis study
#'
#' An experiment level phyloseq object including rplB gene count data from a pouchitis study
#'
#' @format A phyloseq object containing OTU, sample data, and taxonomy tables.
#'
"rplb"

#' rslt
#'
#' Neutral Model result calculated from the Schmidt data
#'
#' @format A list of 8 objects returned by the function \code{neutral_model}
#' \describe{
#'   \item{stats}{A list (Ntm, mean.Ntm, m, and Spearman's rank correlation with p-value)}
#'   \item{obs}{A list (target.data.freq.occur, target.data.mean.abundance, source data freq.ocurr, source.data.mean.abund, target.data.mean.abund.se, source.data.mean.abund.se)}
#'   \item{bestneutralmatrix}{A numerical matrix}
#'   \item{neutral_matrix}{a data frame giving target.data.freq.occur, target.data.mean.abund,
#'     source.data.freq.occur, source.data.mean.abund, target.data.mean.abund.se, and
#'     source.data.mean.abund.se for the neutral OTUs}
#'   \item{for_matrix}{a data frame giving target.data.freq.occur, target.data.mean.abund,
#'     source.data.freq.occur, source.data.mean.abund, target.data.mean.abund.se, and
#'     source.data.mean.abund.se for the OTUs selected for}
#'   \item{against_matrix}{a data frame giving target.data.freq.occur, target.data.mean.abund,
#'     source.data.freq.occur, source.data.mean.abund, target.data.mean.abund.se, and
#'     source.data.mean.abund.se for OTUs selected against}
#'   \item{sum1}{a matrix giving counts of OTUs in 22 categories (SpearmanRank, No.OTUsInTarget, No.OTUsInSource, No.OTUsShared, No.Grey, No.Green, No.Red, SharedAbundance, UniqueAbundanceToTarget, UniqueAbundanceToSource, GreyAbundance, GreenAbundance, RedAbundance, TargetSeqs, SharedSeqs, NeutralSeqs, GreenSeqs, RedSeqs, SourceAnalyzed, SourceNeutral, SourceSelectedFor, SourceSelectedAgainst)}
#'   \item{sum2}{a matrix giving counts of OTUs in 10 categories (TotalSeqsInSource, TotalSeqsInTarget, SharedSeqsInSource, SharedSeqsInTarget, NeutralSeqsInSource, NeutralSeqsInTarget, SelectedForSeqsInSource, SelectedForSeqsInTarget, SelectedAgainstSeqsInSource, SelectedAgainstSeqsInTarget)}
#' }
"rslt"

#' schmidt data
#'
#' Data from Schmidt paper as a phyloseq object containing OTU, sample data and taxonomy tables
#'
#' @format An experiment level \code{phyloseq} object.
#'
"schmidt"
