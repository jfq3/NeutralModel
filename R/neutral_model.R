#' Neutral Model
#'
#' Fits OTU frequencies in the target community to OTU relative abundances in the source
#'   community. Returns an S3 object of class NCM with print, summary, and plot functions.
#'
#' @aliases neutral_mode
#' @usage neutral_model(physeq, source.otu, target.otu, add_taxa = FALSE)
#' @param physeq A phyloseq object containing samples for the source and target communities.
#' @param source.otu A matrix of counts with samples in rows and OTUs in columns.
#' @param target.otu A matrix of counts with samples in rows and OTUs in columns.
#' @param add_taxa A logical. If TRUE, taxa from the phyloseq object are included in the output.
#'
#' @return An S3 object of class NCM.
#' @export
#' @details Returns an S3 object of class NCM with print, summary, and plot functions.
#'   The object consists of a list of 8 items:
#'   \enumerate{
#'   \item stats: Ntm, mean.Nt, amd m.
#'   \item obs: Statistics used by the plot function.
#'   \item bestneutralmatrix: tatistics used by the plot function.
#'   \item neutral_matrix: For each OTU, statistics on frequencies of occurrence and relative abundances
#'      for OTUs occurring in the target community as expected according to the neutral model. May also
#'      include taxonomic assignment of each OTU.
#'   \item against_matrix: Same data for OTUs occurring in the target community less frequently than
#'      expected.
#'   \item for_matrix: Same data for OTUs occurring in the target community more frequently than expected.
#'   \item sum1: Summary statistics included by the summary function.
#'   \item sum2: Summary statistics included by the summary function.
#'   }
#'
#' @importFrom Hmisc binconf
#' @author John Quensen, Tom Schmidt & Arvind Venkataraman
#' @examples
#' data(schmidt)
#' source.com <- veganotu(subset_samples(schmidt, Site=="Mouth"))
#' target.com <- veganotu(subset_samples(schmidt, Site=="Lung"))
#' test <- neutral_model(schmidt, source.com, target.com)

neutral_model <-
function(physeq, source.otu, target.otu, add_taxa=FALSE) {

  # For source and target, calculate frequency of occurrence of each OTU,
  # mean relative abundance, and se of relative abundance.
  source.data <- freq_occur_abund(source.otu)
  target.data <- freq_occur_abund(target.otu)

  # Assemble all neutral data table.
  neutral.data.table <- data.frame(target.data$freq.occur,
                                   target.data$mean.abund,
                                   source.data$freq.occur,
                                   source.data$mean.abund,
                                   target.data$mean.abund.se,
                                   source.data$mean.abund.se)

  # Sort on decreasing relative abundance in source group.
  neutral.data.table <- neutral.data.table[order(neutral.data.table[,4],
                                                 decreasing=TRUE), ]

  # Get lists of OTUs in each category.
  "%w/o%" <- function(x, y) x[!x %in% y] #--  x without y
  source.otus <- colnames(source.otu[ , colSums(source.otu)>0])
  target.otus <- colnames(target.otu[ , colSums(target.otu)>0])
  shared.otus <- intersect(source.otus, target.otus)
  unique.source.otus <- source.otus %w/o% target.otus
  unique.target.otus <- target.otus %w/o% source.otus

  # Determine the detection limit - basically the lowest abundance in the source community.
  obs <- neutral.data.table  #[shared.otu.names, ]
  obs <- obs[shared.otus, ]
  obs <- obs[order(obs[ , 4], decreasing=TRUE), ]
  detlim <- obs[nrow(obs), 4]

  # Calculate bestNtm
  bestNtm <- calc_Ntm(obs)

  # Calculate all the parameters for the best neutral model
  bestneutralmatrix <- matrix(nrow=nrow(obs),ncol=7)
  conf <- matrix(nrow=nrow(obs), ncol=3)
  diff <- matrix(nrow(obs), ncol=1)
  bestneutralmatrix[ , 1] <- bestNtm$bestNtm*obs[ , 4]
  bestneutralmatrix[ , 2] <- bestNtm$bestNtm*(1-obs[ , 4])
  bestneutralmatrix[ , 3] <- pbeta(detlim, bestneutralmatrix[ , 1], bestneutralmatrix[ , 2])
  bestneutralmatrix[ , 4] <- 1-bestneutralmatrix[ , 3]
  bestneutralmatrix[ , 5] <- (obs[ , 1]-bestneutralmatrix[ , 4])^2
  # Determine the 95% confidence intervals using a binomial distribution
  conf[ , ] <- unname(binconf(x=nrow(target.otu)*bestneutralmatrix[ , 4],
                    n=nrow(target.otu),alpha=0.05,method="wilson"))
  # Lower limit
  bestneutralmatrix[ , 6] <- conf[ , 2]
  # Upper limit
  bestneutralmatrix[ , 7] <- conf[ , 3]
  rownames(bestneutralmatrix) <- rownames(obs)

# Get OTUs that are over represented in the target community.
  OTUs.for=c()
  j <- 1
  for(i in 1:nrow(bestneutralmatrix)){
    j <- j
    diff <- (obs[i,1]/bestneutralmatrix[i, 4]) # Ratio between observed frequency
    # and probability it is not (sic) observed in target community.
    # If the OTU falls above the upper confidence interval, deem it as over-represented.
    # if (diff >=1.5 & obs[i,1]>= bestneutralmatrix[i,7]){ # Upper limit
    if (obs[i,1] >= bestneutralmatrix[i,7]){ # Upper limit
      OTUs.for[j] <- rownames(bestneutralmatrix)[i]
      j <- j+1
    }
  }

  # If they exist, subset obs matrix to for_matrix.
  n.for <- length(OTUs.for)
  if(n.for>0) {
    for_matrix <- obs[OTUs.for, ]
    if(add_taxa==TRUE) {
      for_taxa <- phyloseq::tax_table(physeq)[OTUs.for, ]
      for_matrix <- cbind(for_matrix, for_taxa)
    }
  } else {
    for_matrix <- NULL
  }

  # Get OTUs that are under represented in target community.
  OTUs.against=c()
  j <- 1
  for(i in 1:nrow(bestneutralmatrix)){
    j <- j
    diff <- (obs[i,1]/bestneutralmatrix[i, 4]) # Ratio between observed frequency
    # and probability it is not (sic) observed in target community.
    # If the OTU falls below the lower confidence interval, deem it as under-represented.
    # if (diff <=0.75 & obs[i,1]<=bestneutralmatrix[i,6]){
    if (obs[i,1]<=bestneutralmatrix[i,6]){
      OTUs.against[j]=rownames(bestneutralmatrix)[i]
      j=j+1
    }
  }
  # If such exist, subset obs to against_matrix.
  n.against <- length(OTUs.against)
  if(n.against>0) {
    against_matrix <- obs[OTUs.against, ]
    if(add_taxa==TRUE) {
      against_taxa <- phyloseq::tax_table(physeq)[OTUs.against, ]
      against_matrix <- cbind(against_matrix, against_taxa)
    }
  } else {
    against_matrix <- NULL
  }

  # Get OTUs that are neutral.
  OTUs.sel <- c(OTUs.for, OTUs.against)
  OTUs.neutral <- rownames(obs) %w/o% OTUs.sel
  n.neutral <- length(OTUs.neutral)

  # If such exist, subset obs matrix to neutral_matrix.
  if(n.neutral>0) {
    if(n.neutral>0) {
      neutral_matrix <- obs[OTUs.neutral, ]
      if(add_taxa==TRUE) {
        neutral_taxa <- phyloseq::tax_table(physeq)[OTUs.neutral, ]
        neutral_matrix <- cbind(neutral_matrix, neutral_taxa)
      }
    } else {
      neutral_matrix <- NULL
    }
  }

  # Determine the Spearman correlation coefficient
  # correlcoeff <- cor.test(obs[ , 1],bestneutralmatrix[ , 4], method="spearman")

  # # Determine the Pearson correlation coefficient
  # correlcoeff <- cor.test(obs[ , 1],bestneutralmatrix[ , 4], method="pearson")

  # Parameter summary
  mean.Ntm <- mean(rowSums(source.otu))
  mean.m <- bestNtm$bestNtm/mean.Ntm
  stats <- list(Ntm=bestNtm$bestNtm, mean.Nt=mean.Ntm, m=mean.m)
  # stats <- list(Ntm=bestNtm, mean.Ntm=mean.Ntm, m=mean.m, correlcoef=correlcoeff)

  # OTU names
  otu.names.all <- colnames(source.otu)
  otu.names.source <- colnames(source.otu[,colSums(source.otu)>0])
  otu.names.target <- colnames(target.otu[,colSums(target.otu)>0])
  otu.names.unique.to.source <- otu.names.source %w/o% otu.names.target
  otu.names.unique.to.target <- otu.names.target %w/o% otu.names.source
  otu.names.shared <- intersect(otu.names.source, otu.names.target)
  otu.names.neutral <- rownames(neutral_matrix)
  otu.names.for <- rownames(for_matrix)
  otu.names.against <- rownames(against_matrix)

  otu.names <- list(all=otu.names.all,
                    source=otu.names.source,
                    target=otu.names.target,
                    shared=otu.names.shared,
                    unique.to.source=otu.names.unique.to.source,
                    unique.to.target=otu.names.unique.to.target,
                    neutral=otu.names.neutral,
                    selected.for=otu.names.for,
                    selected.against=otu.names.against)

  # Integer output
  OTUsInTarget <- length(otu.names.target)
  OTUsInSource <- length(otu.names.source)
  SharedOTUs <- length(otu.names.shared)
  NeutralOTUs <- length(otu.names.neutral)
  OTUsSelctedAgainst <- length(otu.names.against)
  OTUsSelctedFor <- length(otu.names.for)
  SeqsInTarget <- sum(target.otu)
  SeqsInSource <- sum(source.otu)
  SharedSeqsInSource <- sum(source.otu[, otu.names.shared])
  SharedSeqsInTarget <- sum(target.otu[, otu.names.shared])
  NeutralSeqsInSource <- sum(source.otu[ , otu.names.neutral])
  NeutralSeqsInTarget <- sum(target.otu[, otu.names.neutral])
  SelectedForSeqsInSource <- sum(source.otu[ , otu.names.for])
  SelectedForSeqsInTarget <- sum(target.otu[, otu.names.for])
  SelectedAgainstSeqsInSource  <- sum(source.otu[ , otu.names.against])
  SelectedAgainstSeqsInTarget <- sum(target.otu[, otu.names.against])
  counts.output <- as.matrix(t(data.frame(
    OTUsInTarget,
    OTUsInSource,
    SharedOTUs,
    NeutralOTUs,
    OTUsSelctedAgainst,
    OTUsSelctedFor,
    SeqsInTarget,
    SeqsInSource,
    SharedSeqsInSource,
    SharedSeqsInTarget,
    NeutralSeqsInSource,
    NeutralSeqsInTarget,
    SelectedForSeqsInSource,
    SelectedForSeqsInTarget,
    SelectedAgainstSeqsInSource,
    SelectedAgainstSeqsInTarget
  )))
  colnames(counts.output) <- "Counts"

  # Real number output
  # SpearmanRank <- unname(correlcoeff$estimate) # SpearmanRank
  # Pearson <- unname(correlcoeff$estimate) # Pearson correlation coefficeint
  SharedAbundance <- sum(neutral_matrix[otu.names.shared, 2], na.rm=TRUE) #8
  UniqueAbundanceToTarget <- sum(neutral.data.table[otu.names.unique.to.target, 2], na.rm=TRUE) #9
  UniqueAbundanceToSource <- sum(neutral.data.table[otu.names.shared, 4], na.rm=TRUE) #10
  NeutralAbundance <- sum(neutral_matrix[ , 2], na.rm=TRUE) #11
  ForAbundance <- sum(for_matrix[ , 2], na.rm=TRUE) #12
  AgainstAbundance <- sum(against_matrix[ , 2], na.rm=TRUE) #13
  SourceAnalyzed <- sum(source.otu[ , otu.names.shared])/sum(source.otu) # SourceAnalyzed
  SourceNeutral <- sum(source.otu[ , otu.names.neutral])/sum(source.otu) #SourceNeutral
  SourceSelectedFor <- sum(source.otu[ , otu.names.for])/sum(source.otu) #SourceEnvSelFor
  SourceSelectedAgainst <- sum(source.otu[ , otu.names.against])/sum(source.otu) #SourceEnvSelAgainst
  prop.output <- data.frame(SharedAbundance, UniqueAbundanceToTarget,
                            UniqueAbundanceToSource, NeutralAbundance, ForAbundance,
                            AgainstAbundance, SourceAnalyzed, SourceNeutral,
                            SourceSelectedFor, SourceSelectedAgainst)
  # prop.output <- data.frame(SpearmanRank, SharedAbundance, UniqueAbundanceToTarget,
  #                           UniqueAbundanceToSource, NeutralAbundance, ForAbundance,
  #                           AgainstAbundance, SourceAnalyzed, SourceNeutral,
  #                           SourceSelectedFor, SourceSelectedAgainst)
  # prop.output <- data.frame(Pearson, SharedAbundance, UniqueAbundanceToTarget,
  #                           UniqueAbundanceToSource, NeutralAbundance, ForAbundance,
  #                           AgainstAbundance, SourceAnalyzed, SourceNeutral,
  #                           SourceSelectedFor, SourceSelectedAgainst)
  prop.output <- t(prop.output)
  colnames(prop.output) <- "Proportions"

  rslt <- list(stats=stats, otu.names=otu.names, obs=obs,
               bestneutralmatrix=bestneutralmatrix, neutral_matrix=neutral_matrix,
               against_matrix=against_matrix, for_matrix=for_matrix, sum1=counts.output,
               sum2=prop.output, bestNtm = bestNtm)
  class(rslt) <- "NCM"
  return(rslt)
}
