\name{IdTaxa}
\alias{IdTaxa}
\title{
Assign Sequences a Taxonomic Classification
}
\description{
Classifies sequences according to a training set by assigning a confidence to taxonomic labels for each taxonomic rank.
}
\usage{
IdTaxa(test,
       trainingSet,
       type = "extended",
       strand = "both",
       threshold = 60,
       bootstraps = 100,
       samples = L^0.47,
       minDescend = 0.98,
       processors = 1,
       verbose = TRUE)
}
\arguments{
  \item{test}{
A \code{DNAStringSet} or \code{RNAStringSet} of unaligned sequences.
}
  \item{trainingSet}{
An object of class \code{Taxa} and subclass {Train}.
}
  \item{type}{
Character string indicating the type of output desired.  This should be (an abbreviation of) one of \code{"extended"} or \code{"collapsed"}.  (See value section below.)
}
  \item{strand}{
Character string indicating the orientation of the \code{test} sequences relative to the \code{trainingSet}.  This should be (an abbreviation of) one of \code{"both"}, \code{"top"}, or \code{"bottom"}.  The top strand is defined as the input \code{test} sequences being in the same orientation as the \code{trainingSet}, and the bottom strand is its reverse complement orientation.  The default of \code{"both"} will classify using both orientations and choose the result with highest confidence.
}
  \item{threshold}{
Numeric specifying the confidence at which to truncate the output taxonomic classifications.  Lower values of \code{threshold} will classify deeper into the taxonomic tree at the expense of accuracy, and vise-versa for higher values of \code{threshold}.
}
  \item{bootstraps}{
Integer giving the number of bootstrap replicates to perform for each sequence.
}
  \item{samples}{
A function or call written as a function of `L', which will evaluate to a numeric vector the same length as `L'.  Typically of the form ``\code{A + B*L^C}'', where `A', `B', and `C' are constants.
}
  \item{minDescend}{
Numeric giving the minimum fraction of \code{bootstraps} required to descend the tree during the initial tree descend phase of the algorithm.  Higher values are less likely to descend the tree, causing direct comparison against more sequences in the \code{trainingSet}.  Lower values may increase classification speed at the expense of accuracy.  Suggested values are between \code{1.0} and \code{0.9}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Sequences in \code{test} are each assigned a taxonomic classification based on the \code{trainingSet} created with \code{\link{LearnTaxa}}.  Each taxonomic level is given a confidence between 0\% and 100\%, and the taxonomy is truncated where confidence drops below \code{threshold}.  If the taxonomic classification was truncated, the last group is labeled with ``unclassified_'' followed by the final taxon's name.  Note that the reported confidence is not a p-value but does directly relate to a given classification's probability of being wrong.  The default \code{threshold} of \code{60\%} is intended to minimize the rate of incorrect classifications.  Lower values of \code{threshold} (e.g., \code{50\%}) may be preferred to increase the taxonomic depth of classifications.
}
\value{
If \code{type} is \code{"extended"} (the default) then an object of class \code{Taxa} and subclass {Train} is returned.  This is stored as a list with elements corresponding to their respective sequence in \code{test}.  Each list element contains components:
  \item{taxon}{
A character vector containing the taxa to which the sequence was assigned.
}
  \item{confidence}{
A numeric vector giving the corresponding percent confidence for each taxon.
}
  \item{rank}{
If the classifier was trained with a set of \code{rank}s, a character vector containing the rank name of each taxon.
}

If \code{type} is \code{"collapsed"} then a character vector is returned with the taxonomic assignment for each sequence.  This takes the repeating form ``Taxon name [rank, confidence\%]; ...'' if \code{rank}s were supplied during training, or ``Taxon name [confidence\%]; ...'' otherwise.
}
\references{
Murali, A., et al. (2018). IDTAXA: a novel approach for accurate taxonomic classification of microbiome sequences. Microbiome, 6, 140. https://doi.org/10.1186/s40168-018-0521-5
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{LearnTaxa}}, \code{\link{Taxa-class}}
}
\examples{
data("TrainingSet_16S")

# import test sequences
fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)

# remove any gaps in the sequences
dna <- RemoveGaps(dna)

# classify the test sequences
ids <- IdTaxa(dna, TrainingSet_16S, strand="top")
ids

# view the results
plot(ids, TrainingSet_16S)
}