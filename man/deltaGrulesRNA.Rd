\name{deltaGrulesRNA}
\alias{deltaGrulesRNA}
\docType{data}
\title{
Pseudoenergy Parameters for RNA Quadruplets
}
\description{
An 8D array with four adjacent base pairs of the RNA duplex.  Each dimension has five elements defining the nucleotide at that position ("A", "C", "G", "U", or "-").  The array contains the pseudoenergy of duplex formation for every quadruple base pairing.
}
\usage{data("deltaGrulesRNA")}
\format{
  The format is:
 num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.776 -0.197 -0.197 -0.291 0 ...
 - attr(*, "dimnames")=List of 8
  ..$ : chr [1:5] "A" "C" "G" "U" ...
  ..$ : chr [1:5] "A" "C" "G" "U" ...
  ..$ : chr [1:5] "A" "C" "G" "U" ...
  ..$ : chr [1:5] "A" "C" "G" "U" ...
  ..$ : chr [1:5] "A" "C" "G" "U" ...
  ..$ : chr [1:5] "A" "C" "G" "U" ...
  ..$ : chr [1:5] "A" "C" "G" "U" ...
  ..$ : chr [1:5] "A" "C" "G" "U" ...
}
\details{
The first four dimensions correspond to the four top strand positions from 5' to 3'.  The fifth to eighth dimensions correspond to the four bottom strand positions from 5' to 3'.
}
\source{
Psuedoenergy values of ungapped quadruplets are inferred from their log-odds of being found in palindromes within hairpin regions across thousands of non-coding RNA alignments.  Each value represents the log-odds of \emph{in vivo} folding relative to chance.
}
\examples{
data(deltaGrulesRNA)
# dG of the duplex AGCU / ACCU pairing (1 mismatch)
deltaGrulesRNA["A", "G", "C", "U", "A", "C", "C", "U"]
}
\keyword{datasets}
