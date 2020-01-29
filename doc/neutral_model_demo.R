## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(phyloseq)
library(vegan)
library(Hmisc)
library(NeutralModel)
data(schmidt)

## -----------------------------------------------------------------------------
schmidt
range(sample_sums(schmidt))
table(sample_data(schmidt)[,1])
rank_names(schmidt)

## -----------------------------------------------------------------------------
sum(taxa_sums(schmidt)==0)

## -----------------------------------------------------------------------------
source.com <- veganotu(subset_samples(schmidt, Site=="Mouth"))
target.com <- veganotu(subset_samples(schmidt, Site=="Lung"))

## -----------------------------------------------------------------------------
test <- neutral_model(physeq=schmidt, source.otu=source.com,
                      target.otu=target.com, add_taxa=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  test <- neutral_model(schmidt, source.com, target.com, add_taxa=TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  test <- neutral_model(NULL, source.com, target.com)

## -----------------------------------------------------------------------------
plot(test)

## -----------------------------------------------------------------------------
plot(test, xlab="Relative Abundance in Mouth", ylab="Detected Frequency in Lung",
     main="My NCM Plot")

## -----------------------------------------------------------------------------
test

## -----------------------------------------------------------------------------
p.val <- pval_NeutralModel(test)
p.val

## -----------------------------------------------------------------------------
summary(test)

## -----------------------------------------------------------------------------
test$for_matrix[,8:12]

## -----------------------------------------------------------------------------
test$for_matrix[ , c("Phylum", "Genus")]

## -----------------------------------------------------------------------------
data(rplb)
sample_variables(rplb)
unique(sample_data(rplb)[,"Diagnosed"])
unique(sample_data(rplb)[,"Group"])

## -----------------------------------------------------------------------------
temp <- subset_samples(rplb, Diagnosed=="No")

## -----------------------------------------------------------------------------
temp <- subset_samples(temp, Group!="IPAA.Healthy")
table(sample_data(temp)[,"Site"])

## -----------------------------------------------------------------------------
temp <- prune_taxa(taxa_sums(temp)>9, temp)
range(sample_sums(temp))

## -----------------------------------------------------------------------------
set.seed(123) # To make reproducible.
temp <- rarefy_even_depth(temp, rngseed=TRUE)
range(sample_sums(temp))

## -----------------------------------------------------------------------------
source.com <- veganotu(subset_samples(temp, Site=="Colon"))
target.com <- veganotu(subset_samples(temp, Site=="IPAA"))

## -----------------------------------------------------------------------------
rslt <- neutral_model(temp, source.com, target.com)
rslt
summary(rslt)

## -----------------------------------------------------------------------------
plot(rslt, main="rplB Counts", xlab="Relative Abundance in Colon",
     ylab="Detection Frequency in Pouch")

## -----------------------------------------------------------------------------
p.val <- pval_NeutralModel(rslt, n.iter=999, plot=FALSE)
p.val

