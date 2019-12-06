#' veganotu
#'
#' Extracts a community matrix from a phyloseq object.
#'
#' @param physeq A phyloseq object containing an OTU table.
#'
#' @return A community matrix with samples as rows and taxa as columns.
#' @export
#' @importFrom phyloseq taxa_are_rows
#' @author John Quensen
#' @examples
#' data(rplb)
#' veganotu(rplb)
#'
veganotu <-
function(physeq) {
    OTU <- phyloseq::otu_table(physeq)
    if (phyloseq::taxa_are_rows(OTU)) {
      OTU <- t(OTU)
    }
    OTU <- as(OTU, "matrix")
    return(OTU)
  }
