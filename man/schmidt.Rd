\name{schmidt}
\alias{schmidt}
\docType{data}
\title{
schmidt
}
\description{
A phyloseq object containing 16S rRNA DNA sequence data from a human microbiome study.
}
\usage{data("schmidt")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : int [1:85, 1:475] 127 84 76 77 60 80 50 95 86 76 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:85] "Sam_1" "Sam_2" "Sam_3" "Sam_4" ...
  .. .. .. .. ..$ : chr [1:475] "OTU0001" "OTU0002" "OTU0003" "OTU0004" ...
  .. .. ..@ taxa_are_rows: logi FALSE
  ..@ tax_table:Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
  .. .. ..@ .Data: chr [1:475, 1:6] "Bacteria" "Bacteria" "Bacteria" "Bacteria" ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:475] "OTU0001" "OTU0002" "OTU0003" "OTU0004" ...
  .. .. .. .. ..$ : chr [1:6] "Domain" "Phylum" "Class" "Order" ...
  ..@ sam_data :'data.frame':	85 obs. of  1 variable:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 1
  .. .. .. ..$ : Factor w/ 2 levels "Lung","Mouth": 2 2 2 2 2 2 2 2 2 2 ...
  .. .. ..@ names    : chr "Site"
  .. .. ..@ row.names: chr [1:85] "Sam_1" "Sam_2" "Sam_3" "Sam_4" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree : NULL
  ..@ refseq   : NULL
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Venkataraman, A., C. M. Bassis, J. M. Beck, V. B. Young, J. L. Curtis, G. B. Huffnagle, and T. M. Schmidt. 2015. Application of a neutral community model to assess structuring of the human lung microbiome. MBio 6:e02284-14.
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(schmidt)
## maybe str(schmidt) ; plot(schmidt) ...
}
\keyword{datasets}
