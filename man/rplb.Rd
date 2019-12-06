\name{rplb}
\alias{rplb}
\docType{data}
\title{
rplB count data
}
\description{
rplB gene count data from a pouchitis study.
}
\usage{data("rplb")}
\format{
  The format is:
Formal class 'phyloseq' [package "phyloseq"] with 5 slots
  ..@ otu_table:Formal class 'otu_table' [package "phyloseq"] with 2 slots
  .. .. ..@ .Data        : int [1:377, 1:96] 250 59 45 49 131 86 25 119 0 1 ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:377] "OTU_1" "OTU_2" "OTU_3" "OTU_4" ...
  .. .. .. .. ..$ : chr [1:96] "A1" "A2" "B1" "B2" ...
  .. .. ..@ taxa_are_rows: logi TRUE
  ..@ tax_table:Formal class 'taxonomyTable' [package "phyloseq"] with 1 slot
  .. .. ..@ .Data: chr [1:377, 1] "Bacteroides " "unclassified " "Ruminococcus " "Eubacterium " ...
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:377] "OTU_1" "OTU_2" "OTU_3" "OTU_4" ...
  .. .. .. .. ..$ : chr "Genus"
  ..@ sam_data :'data.frame':	96 obs. of  6 variables:
Formal class 'sample_data' [package "phyloseq"] with 4 slots
  .. .. ..@ .Data    :List of 6
  .. .. .. ..$ : Factor w/ 34 levels "A","B","C","D",..: 1 1 2 2 3 3 4 4 5 5 ...
  .. .. .. ..$ : Factor w/ 2 levels "Healthy","Pouchitis": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ : Factor w/ 2 levels "Colon","IPAA": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ : Factor w/ 9 levels "M1","M10","M2",..: 1 3 1 3 1 3 1 3 1 3 ...
  .. .. .. ..$ : Factor w/ 3 levels "Colon.Healthy",..: 1 1 1 1 1 1 1 1 1 1 ...
  .. .. .. ..$ : Factor w/ 2 levels "No","Yes": 1 1 1 1 1 1 1 1 1 1 ...
  .. .. ..@ names    : chr [1:6] "Patient" "Condition" "Site" "Time" ...
  .. .. ..@ row.names: chr [1:96] "A1" "A2" "B1" "B2" ...
  .. .. ..@ .S3Class : chr "data.frame"
  ..@ phy_tree : NULL
  ..@ refseq   : NULL
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(rplb)
## maybe str(rplb) ; plot(rplb) ...
}
\keyword{datasets}
