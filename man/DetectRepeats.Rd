\name{DetectRepeats}
\alias{DetectRepeats}
\title{
Detect Repeats in a Sequence
}
\description{
Detects approximate copies of sequence patterns that likely arose from duplication events and therefore share a common ancestor.
}
\usage{
DetectRepeats(myXStringSet,
              type = "tandem",
              minScore = 10,
              allScores = FALSE,
              maxPeriod = 10000,
              maxFailures = 2,
              maxShifts = 5,
              alphabet = AA_REDUCED[[125]],
              processors = 1,
              verbose = TRUE,
              \dots)
}
\arguments{
  \item{myXStringSet}{
An \code{AAStringSet}, \code{DNAStringSet}, or \code{RNAStringSet} object of unaligned sequences.
}
  \item{type}{
Character string indicating the type of repeats to detect.  This should be (an abbreviation of) one of \code{"tandem"}, \code{"interspersed"}, or \code{"both"}.  (See details section below.)
}
  \item{minScore}{
Numeric giving the minimum log-odds score of repeats in \code{myXStringSet} to report.
}
  \item{allScores}{
Logical specifying whether all repeats should be returned (\code{TRUE}) or only the top scoring repeat when there are multiple overlapping matches in the same region.
}
  \item{maxPeriod}{
Numeric indicating the maximum periodicity of tandem repeats to consider.  Interspersed repeats will only be detected that are at least \code{maxPeriod} nucleotides apart.
}
  \item{maxFailures}{
Numeric determining the maximum number of failing attempts to extend a repeat that are permitted.  Numbers greater than zero may increase accuracy at the expense of speed, with decreasing marginal returns as \code{maxFailures} gets higher and higher.
}
  \item{maxShifts}{
Numeric determining the maximum number of failing attempts to shift a repeat left or right that are permitted.  Numbers greater than zero may increase accuracy at the expense of speed, with decreasing marginal returns as \code{maxShifts} gets higher and higher.
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between sequences.  A non-reduced amino acid alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
  \item{\dots}{
Further arguments to be passed directly to \code{\link{FindSynteny}} if \code{type} is \code{"interspersed"} or \code{"both"}.
}
}
\details{
Many sequences are composed of a substantial fraction of repetitive sequence.  Two main forms of repetition are tandem repeats and interspersed repeats, which can be caused by duplication events followed by divergence.  Detecting duplications is challenging because of variability in repeat length and composition due to evolution.  The significance of a repeat can be quantified by its time since divergence from a common ancestor.  \code{DetectRepeats} uses a seed-and-extend approach to identify candidate repeats, and tests whether a set of repeats is statistically distinct from a background distribution using a goodness of fit multinomial test.  The background distribution is derived from input sequences (\code{myXStringSet}) and represents the distribution of characters that would be expected if repeats had diverged infinitely long ago (see Schaper et al. (2012) for examples).  The reported score is the product of statistical significance (i.e., log-odds) and average similarity among the repeats.  Therefore, a high-score implies the repeats diverged within a finite amount of time from a common ancestor.

Two possible types of repeats are detectable:  (1)  Tandem repeats are contiguous approximate copies of a nucleotide or amino acid sequence.  Once a k-mer seed is identified, repeated attempts are made to optimize the beginning and ending positions, as well as attempting to extend the repeat to the left and right.  (2)  Interspersed repeats are dispersed approximate copies of a nucleotide sequence on the same strand or opposite strands.  These are identified with \code{\link{FindSynteny}}, aligned with \code{\link{AlignSynteny}}, and then scored using the same statistical framework as tandem repeats.  In both cases, the highest scoring repeat in each region is returned, unless \code{allScores} is \code{TRUE}, in which case overlapping repeats are permitted in the result.
}
\value{
If \code{type} is \code{"tandem"}, a \code{data.frame} giving the \code{"Index"} of the sequence in \code{myXStringSet}, \code{"Begin"} and \code{"End"} positions of tandem repeats, \code{"Left"} and \code{"Right"} positions of each repeat, and its \code{"Score"}.

If \code{type} is \code{"interspersed"}, a \code{data.frame} similar to the matrix in the lower diagonal of Synteny objects (see \code{\link{Synteny-class}}).

If \code{type} is \code{"both"}, a \code{list} with the above two elements.
}
\references{
Schaper, E., et al. (2012). Repeat or not repeat?-Statistical validation of tandem repeat prediction in genomic sequences. Nucleic Acids Research, \bold{40(20)}, 10005-17.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\examples{
data(yeastSEQCHR1)
dna <- DNAStringSet(yeastSEQCHR1)

x <- DetectRepeats(dna)
x

# number of tandem repeats
lengths(x[, "Left"])

# average periodicity of tandem repeats
per <- mapply(function(a, b) b - a + 1,
	x[, "Left"],
	x[, "Right"])
sapply(per, mean)

# extract a tandem repeat
i <- 1
reps <- extractAt(dna[[x[i, "Index"]]],
	IRanges(x[[i, "Left"]], x[[i, "Right"]]))
reps
reps <- AlignSeqs(reps) # align the repeats
reps
BrowseSeqs(reps)

y <- DetectRepeats(dna, type="interspersed")
y
}
