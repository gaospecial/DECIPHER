\name{Clusterize}
\alias{Clusterize}
\title{
Cluster Sequences By Similarity
}
\description{
Groups the sequences into approximate clusters of similarity.
}
\usage{
Clusterize(myXStringSet,
           cutoff = 0,
           method = "overlap",
           includeTerminalGaps = FALSE,
           penalizeGapLetterMatches = NA,
           minCoverage = 0.5,
           maxReps = 1000,
           avgComparisons = 5000,
           maxAlignments = 100,
           invertCenters = FALSE,
           processors = 1,
           verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
The \code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet} to cluster.
}
  \item{cutoff}{
A vector of maximum distances (approximately) separating sequences in the same cluster (i.e., 1 - similarities).  Multiple cutoffs may be provided in ascending or descending order.  (See details section below.)
}
  \item{method}{
Character string determining the region in which distance is calculated.  This should be (an unambiguous abbreviation of) one of \code{"overlap"}, \code{"shortest"}, or \code{"longest"}.  The default \code{method} (\code{"overlap"}) calculates distance from the overlapping region between terminal gaps when \code{includeTerminalGaps} is \code{FALSE} and the entire alignment otherwise.  Setting \code{method} to \code{"shortest"} or \code{"longest"} will use the region between the start and end of the shortest or longest sequence, respectively, for each pairwise distance.  The \code{method} is only applicable when \code{includeTerminalGaps} is \code{TRUE}.
}
  \item{includeTerminalGaps}{
Logical specifying whether or not to include terminal gaps ("-" characters) in the pairwise alignments into the calculation of distance.
}
  \item{penalizeGapLetterMatches}{
Logical specifying whether or not to consider gap-to-letter matches as mismatches.  If \code{FALSE}, then gap-to-letter matches are not included in the total length used to calculate distance, and if \code{TRUE} all gaps-to-letter pairs are considered mismatches.  The default (\code{NA}) is to penalize gap-to-letter mismatches once per insertion or deletion, which treats runs of gaps (i.e., indels) as equivalent to a single mismatch.
}
  \item{minCoverage}{
Numeric between zero and one giving the minimum fraction of sequence positions (not gap or mask) in the shortest sequence that must be overlapping with the longer sequence for the sequences to be clustered.
}
  \item{maxReps}{
An integer specifying the maximum number of replicates to perform when ordering sequences based on their k-mer similarity.
}
  \item{avgComparisons}{
An integer determining the average number of comparisons per sequence to perform when attempting to find cluster centers.
}
  \item{maxAlignments}{
An integer giving the maximum number of alignments to perform when attempting to assign a sequence to an existing cluster.
}
  \item{invertCenters}{
Logical controlling whether the cluster center is inverted (i.e., multiplied by \code{-1}), which allows the centers to be determined from the results.  The default (\code{FALSE}) only returns positive cluster numbers.  If \code{TRUE}, the center sequence(s) of each cluster are negative.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{Clusterize} groups the input sequences into approximate clusters using a heuristic algorithm.  First, the sequences are grouped (if needed) using spherical k-means of their 3-mer frequencies.  Second, the sequences are ordered by k-mer similarity with a stochastic partitioning approach.  Third, the sequences are iteratively clustered in this order by their similarity in a pairwise alignment.  That is, the first sequence becomes the representative of cluster #1.  If the second sequence is within \code{cutoff} distance then it is added to the cluster, otherwise it becomes a new cluster representative.  The remaining sequences are matched to cluster representatives in a similar fashion until all sequences belong to a cluster.  In the majority of cases, this process results in clusters with members separated by less than \code{cutoff} distance, and all cluster members must be within \code{cutoff} distance of their cluster representative.

The calculation of distance can be controlled in multiple ways, with each parameterization of distance having advantages and disadvantages.  By default, distance is the fraction of positions that are different, including gaps, within the overlapping region in a pairwise alignment.  The defaults will handle partial-length sequences well, but also cluster sequences with high similarity between their opposite ends.  For this reason, it is important to set \code{minimumCoverage} (by default \code{0.5}) such that distances are based off of considerable overlap between sequences in the pairwise alignment.  This distance parameterization works well, but there are reasonable alternatives.

If \code{penalizeGapLetterMatches} is \code{FALSE}, the distance will exclude gap regions.  If \code{includeTerminalGaps} is \code{TRUE}, the calculation of distance will use the entire (global) alignment.  If \code{method} is \code{"shortest"} and \code{includeTerminalGaps} is set to \code{TRUE}, then the distance is calculated for the region encompassed by the shorter sequence in each pair, which is the common definition of distance used by other clustering programs.  This common definition of distance will more sometimes separate partly overlapping sequences, which is why it is not the default.

The algorithm requires time proportional to the number of input sequences in \code{myXStringSet}.  The first step, spherical k-means requires linear time.   Next, the sequences are compared with up to \code{maxReps} passes that each take linear time.  Ordering of the sequences is performed in linear time using radix sorting.  Finally, each sequence is compared with up to \code{maxComparisons} previous cluster representatives.  This is possible because the sequences are in an order that is related to their similarity, such that nearby sequences are more similar.  Hence, the complete algorithm scales in linear time asymptotically and returns clusters of sequences within \code{cutoff} distance of their center sequence.

Multiple cutoffs can be provided in sorted order.  If the \code{cutoff}s are provided in \emph{descending} order then clustering at each new value of \code{cutoff} is continued within the prior \code{cutoff}'s clusters.  In this way clusters at lower values of \code{cutoff} are completely contained within their ``umbrella'' clusters at higher values of \code{cutoff}.  This slightly accelerates the clustering process, because each subsequent group is only clustered within the previous group.  If multiple cutoffs are provided in \emph{ascending} order then clustering at each level of \code{cutoff} is independent of the prior level.

Note, the algorithm stochastically selects partitions for ordering the sequences.  Hence, clusters can vary from run-to-run unless the random number seed is set for repeatability (i.e., with \code{set.seed}).  Also, \code{invertCenters} can be used to determine the center sequence of each cluster from the output.  Since identical sequences will always be assigned the same cluster numbers, it is possible for more than one input sequence in \code{myXStringSet} to be assigned as the center of a cluster.
}
\value{
A data.frame is returned with dimensions \eqn{N*M}, where each one of \eqn{N} sequences is assigned to a cluster at the \eqn{M}-level of cutoff.  The row.names of the data.frame correspond to the \emph{names} of \code{myXStingSet}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\seealso{
\code{\link{AA_REDUCED}}, \code{\link{DistanceMatrix}}, \code{\link{TreeLine}}
}
\examples{
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
aa <- translate(dna)

# cluster each cutoff within the previous cluster (faster)
clusters <- Clusterize(aa, cutoff=seq(0.7, 0, -0.1))
head(clusters)
apply(clusters, 2, max) # number of clusters per cutoff

# cluster each cutoff independently (fewer clusters per cutoff)
clusters <- Clusterize(aa, cutoff=seq(0, 0.7, 0.1))
head(clusters)
apply(clusters, 2, max) # number of clusters per cutoff

# make cluster center(s) negative for tracking
clusters <- Clusterize(aa, cutoff=0.5, invertCenters=TRUE)
head(clusters)
clusters[clusters$cluster < 0,, drop=FALSE]
unique(aa[clusters$cluster < 0]) # unique cluster centers
apply(clusters, 2, function(x) max(abs(x))) # number of clusters

# cluster nucleotide sequences
clusters <- Clusterize(dna, cutoff=0.5, invertCenters=TRUE)
head(clusters)
apply(clusters, 2, function(x) max(abs(x))) # number of clusters
}
