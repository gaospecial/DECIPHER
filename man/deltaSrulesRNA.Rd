\name{deltaSrulesRNA}
\alias{deltaSrulesRNA}
\docType{data}
\title{
Change in Entropy of Hybridization of RNA/RNA Quadruplets in Solution
}
\description{
An 8D array with four adjacent base pairs of the RNA duplex.  Each dimension has five elements defining the nucleotide at that position ("A", "C", "G", "T", or "-").  The array contains the standard entropy change of probe binding (dS, [kcal/mol]) for every quadruple base pairing.
}
\usage{data(deltaSrulesRNA)}
\format{
  The format is:
 num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.0182 0 0 0 0 ...
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
Data from a variety of publications by SantaLucia \emph{et al}.
}
\references{
SantaLucia, J., Jr., & Hicks, D. (2004) The Thermodynamics of DNA Structural Motifs. Annual Review of Biophysics and Biomolecular Structure, 33(1), 415-440. doi:10.1146/annurev.biophys.32.110601.141800.
}
\examples{
data(deltaSrulesRNA)
# dS of the duplex AGCU / A-CU pairing
deltaSrulesRNA["A", "G", "C", "U", "A", "-", "C", "U"]
}
\keyword{datasets}
