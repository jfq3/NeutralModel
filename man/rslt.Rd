\name{rslt}
\alias{rslt}
\docType{data}
\title{
rslt
}
\description{
Neutral Model result of class NCM calculated from the Schmidt data.
}
\usage{data("rslt")}
\format{
  The format is:
List of 8
 $ stats            :List of 4
  ..$ Ntm       : num 3.15
  ..$ mean.Ntm  : num 1415
  ..$ m         : num 0.00223
  ..$ correlcoef:List of 8
  .. ..$ statistic  : Named num 24922
  .. .. ..- attr(*, "names")= chr "S"
  .. ..$ parameter  : NULL
  .. ..$ p.value    : num 0.0779
  .. ..$ estimate   : Named num 0.233
  .. .. ..- attr(*, "names")= chr "rho"
  .. ..$ null.value : Named num 0
  .. .. ..- attr(*, "names")= chr "rho"
  .. ..$ alternative: chr "two.sided"
  .. ..$ method     : chr "Spearman's rank correlation rho"
  .. ..$ data.name  : chr "obs[, 1] and bestneutralmatrix[, 4]"
  .. ..- attr(*, "class")= chr "htest"
 $ obs              :'data.frame':	58 obs. of  6 variables:
  ..$ target.data.freq.occur   : num [1:58] 0.938 0.938 0.812 0.312 0.188 ...
  ..$ target.data.mean.abund   : num [1:58] 0.1686 0.0695 0.0642 0.0127 0.0011 ...
  ..$ source.data.freq.occur   : num [1:58] 1 1 0.956 0.956 0.689 ...
  ..$ source.data.mean.abund   : num [1:58] 0.3797 0.0609 0.0608 0.0533 0.0469 ...
  ..$ target.data.mean.abund.se: num [1:58] 0.040117 0.018449 0.022953 0.00845 0.000968 ...
  ..$ source.data.mean.abund.se: num [1:58] 0.0327 0.01339 0.00848 0.0086 0.01739 ...
 $ bestneutralmatrix: num [1:58, 1:7] 1.196 0.192 0.191 0.168 0.148 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:58] "OTU_1" "OTU_3" "OTU_4" "OTU_8" ...
  .. ..$ : NULL
 $ neutral_matrix   :'data.frame':	33 obs. of  6 variables:
  ..$ target.data.freq.occur   : num [1:33] 0.938 0.938 0.812 1 0.375 ...
  ..$ target.data.mean.abund   : num [1:33] 0.16864 0.06952 0.06422 0.12244 0.00389 ...
  ..$ source.data.freq.occur   : num [1:33] 1 1 0.956 1 0.844 ...
  ..$ source.data.mean.abund   : num [1:33] 0.3797 0.0609 0.0608 0.0369 0.0143 ...
  ..$ target.data.mean.abund.se: num [1:33] 0.04012 0.01845 0.02295 0.03274 0.00288 ...
  ..$ source.data.mean.abund.se: num [1:33] 0.0327 0.01339 0.00848 0.00295 0.00598 ...
 $ against_matrix   :'data.frame':	4 obs. of  6 variables:
  ..$ target.data.freq.occur   : num [1:4] 0.312 0.188 0.125 0.188
  ..$ target.data.mean.abund   : num [1:4] 0.01268 0.0011 0.00115 0.00358
  ..$ source.data.freq.occur   : num [1:4] 0.956 0.689 1 0.911
  ..$ source.data.mean.abund   : num [1:4] 0.0533 0.0469 0.0393 0.0184
  ..$ target.data.mean.abund.se: num [1:4] 0.00845 0.000968 0.001102 0.002765
  ..$ source.data.mean.abund.se: num [1:4] 0.0086 0.01739 0.00525 0.00344
 $ for_matrix       :'data.frame':	21 obs. of  6 variables:
  ..$ target.data.freq.occur   : num [1:21] 1 0.688 0.938 0.438 0.312 ...
  ..$ target.data.mean.abund   : num [1:21] 0.046908 0.056272 0.081095 0.000795 0.020583 ...
  ..$ source.data.freq.occur   : num [1:21] 1 0.933 0.756 0.622 0.556 ...
  ..$ source.data.mean.abund   : num [1:21] 0.03625 0.01046 0.00319 0.00217 0.00203 ...
  ..$ target.data.mean.abund.se: num [1:21] 0.015179 0.052055 0.029879 0.000309 0.010729 ...
  ..$ source.data.mean.abund.se: num [1:21] 0.004722 0.001191 0.000619 0.000408 0.000578 ...
 $ sum1             : num [1:22, 1] 0.233 63 142 58 33 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:22] "SpearmanRank" "No.OTUsInTarget" "No.OTUsInSource" "No.OTUsShared" ...
  .. ..$ : NULL
 $ sum2             : num [1:10, 1] 63675 22640 52350 22536 38315 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:10] "TotalSeqsInSource" "TotalSeqsInTarget" "SharedSeqsInSource" "SharedSeqsInTarget" ...
  .. ..$ : NULL
 - attr(*, "class")= chr "NCM"
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
data(rslt)
## maybe str(rslt) ; plot(rslt) ...
}
\keyword{datasets}
