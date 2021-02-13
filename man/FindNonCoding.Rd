\name{FindNonCoding}
\alias{FindNonCoding}
\title{
Find Non-Coding RNAs in a Genome
}
\description{
Searches for conserved patterns representing a family of non-coding RNAs.  Returns the start and end boundaries of potential matches along with their log-odds score.
}
\usage{
FindNonCoding(x,
              myXStringSet,
              minScore = 16,
              allScores = FALSE,
              verbose = TRUE)
}
\arguments{
  \item{x}{
A \code{NonCoding} object or a list of \code{NonCoding} objects for searching.
}
  \item{myXStringSet}{
A \code{DNAStringSet} or \code{RNAStringSet} object of unaligned sequences, typically representing a genome.
}
  \item{minScore}{
Numeric giving the minimum log-odds score of matches to \code{x} in \code{myXStringSet} to report.  The maximum false discovery rate is approximately \code{exp(-minScore)} per nucleotide per object in \code{x}.
}
  \item{allScores}{
Logical specifying whether all matches should be returned (\code{TRUE}) or only the top matches when there are multiple matches in the same region.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
Non-coding RNAs are identified by the location of representative sequence patterns relative to the beginning and end of the non-coding RNA.  Potential matches to each \code{NonCoding} object in \code{x} are scored based on their log-odds relative to a background that is derived from the input sequence (\code{myXStringSet}).  Matches of at least \code{minScore} are returned as a \code{Genes} object with the \code{"Gene"} column set to the negative index of the list element of \code{x} that was identified.
}
\value{
An object of class \code{Genes}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{LearnNonCoding}}, \code{\link{NonCoding-class}}, \code{\link{ExtractGenes}}, \code{\link{Genes-class}}, \code{\link{WriteGenes}}
}
\examples{
data(NonCodingRNA_Bacteria)
x <- NonCodingRNA_Bacteria
names(x)

# import a test genome
fas <- system.file("extdata",
	"Chlamydia_trachomatis_NC_000117.fas.gz",
	package="DECIPHER")
genome <- readDNAStringSet(fas)

z <- FindNonCoding(x, genome)
z
sort(table(names(x)[-z[, "Gene"]]))

genes <- ExtractGenes(z, genome, type="RNAStringSet")
genes
}
