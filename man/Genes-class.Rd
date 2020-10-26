\name{Genes}
\alias{Genes-class}% the class
\alias{[.Genes}
\alias{print.Genes}
\alias{plot.Genes}
\title{
Genes objects and accessors
}
\description{
Gene prediction consist of delimiting the boundaries of regions that function as genes within a genome.  Class \code{Genes} provides objects and functions for storing the boundaries of genes and associated information resulting from gene prediction.
}
\usage{
\method{plot}{Genes}(x,
    xlim = c(1, 1e4),
    ylim = c(-100, 500),
    interact = FALSE,
    \dots)

\method{print}{Genes}(x, \dots)

\method{[}{Genes}(x, i, j, \dots)
}
\arguments{
  \item{x}{
An object of class \code{Genes}.
}
  \item{xlim}{
Numeric vector of length 2 specifying the x-axis limits for plotting.
}
  \item{ylim}{
Numeric vector of length 2 specifying the y-axis limits for plotting.
}
  \item{interact}{
Logical determining whether the plot is interactive.  If \code{TRUE}, clicking the plot on the right or left side will scroll one frame in that direction.  To end interaction, either right-click, press the escape key, or press the stop button depending on the graphics device in use.
}
  \item{i}{
Numeric or character vector of row indices to extract from \code{x}.
}
  \item{j}{
Numeric or character vector of column indices to extract from \code{x}.  If \code{j} is missing, all columns are included and the returned object will also belong to class \code{Genes}.
}
  \item{\dots}{
Other optional parameters.
}
}
\details{
Objects of class \code{Genes} are stored as numeric matrices containing information pertaining to gene predictions.  The matrix columns include the index (\code{"Index"}) of the corresponding sequence in the original genome, the strand (\code{"Strand"}) where the gene is located (either \code{"+"} (\code{0}) or \code{"-"} (\code{1}), the beginning (\code{"Begin"}) and ending (\code{"End"}) positions of the gene, several scores used for prediction, the fraction of replicates (\code{"FractionReps"}) where the gene was included in the prediction, and whether (\code{> 0}) or not (\code{0}) the region was predicted to be a gene.  Note that the start codon is at the beginning position when the strand is \code{"+"} and end when the strand is \code{"-"}.

The \code{print} method displays summary information about the genes, including their length range, coding density, number of different initiation codons, and the range of correlations among the different scores used in prediction.  Note that subsetting a \code{Genes} object will result in a lower coding density being displayed because the density is calculated relative to the entire genome originally used for prediction.

The \code{plot} method will show the total score of each prediction along the genome.  This is most useful when displaying the result of setting \code{allScores} to \code{TRUE} in \code{FindGenes}.  Here, possible genes on either strand will be shown, with the predicted genes highlighted.  The start (solid) and stop (dashed) positions are denoted by vertical lines.  Note that the x-axis is cumulative genome position and changes between genome sequences indices are demarcated by dashed vertical lines.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{ExtractGenes}}, \code{\link{FindGenes}}, \code{\link{WriteGenes}}
}
\examples{
# import a test genome
fas <- system.file("extdata",
	"Chlamydia_trachomatis_NC_000117.fas.gz",
	package="DECIPHER")
genome <- readDNAStringSet(fas)

x <- FindGenes(genome, allScores=TRUE)
x
head(unclass(x)) # the underlying structure

plot(x)
}
