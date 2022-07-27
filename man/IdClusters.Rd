\name{IdClusters}
\alias{IdClusters}
\title{
Cluster Sequences By Similarity
}
\description{
Groups the sequences into approximate clusters of similarity.
}
\usage{
IdClusters(myXStringSet,
           cutoff = 0,
           maxReps = 1000,
           maxAttempts = 100,
           alphabet = AA_REDUCED[[30]],
           processors = 1,
           verbose = TRUE)
}
\arguments{
  \item{myXStringSet}{
The \code{DNAStringSet}, \code{RNAStringSet}, or \code{AAStringSet} to cluster.
}
  \item{cutoff}{
A vector of maximum distances (approximately) separating sequences in the same cluster.  Multiple cutoffs may be provided in ascending or descending order.  (See details section below.)
}
  \item{maxReps}{
An integer specifying the maximum number of replicates to perform in the first phase of the algorithm when ordering sequences based on their k-mer similarity.
}
  \item{maxAttempts}{
An integer giving the maximum number of alignments to perform in the second phase of the algorithm when attempting to assign a sequence to an existing cluster.
}
  \item{alphabet}{
Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups.  Alphabet reduction helps to find more distant homologies between sequences.  A non-reduced amino acid alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.  Only applicable if \code{myXStringSet} is an \code{AAStringSet}.
}
  \item{processors}{
The number of processors to use, or \code{NULL} to automatically detect and use all available processors.
}
  \item{verbose}{
Logical indicating whether to display progress.
}
}
\details{
\code{IdClusters} groups the input sequences into approximate clusters using a heuristic algorithm.  First the sequences are ordered by k-mer similarity using a stochastic partitioning approach.  Next, sequences are iteratively clustered in this order by the distance between their overlapping region in a pairwise alignment.  That is, the first sequence becomes the representative of cluster #1.  If the second sequence is within \code{cutoff} distance then it is added to the cluster, otherwise it becomes a new cluster representative.  The remaining sequences are matched to cluster representatives in a similar fashion until all sequences belong to a cluster.  In the majority of cases, this process results in clusters with members separated by less than \code{cutoff} distance.

The algorithm requires time proportional to the number of unique input sequences in \code{myXStringSet}.  The first phase (i.e., ordering) performs up to \code{maxReps} passes that each take linear time.  Ordering of the sequences is performed in linear time using radix sorting.  The second phase (i.e., clustering) requires linear time because comparisons are only made for up to \code{maxAttempts} previous cluster representatives.  This is possible because the sequences are in an order that is related to their similarity, such that nearby sequences are more similar.  Hence, the complete algorithm scales in linear time and returns clusters of sequences within \code{cutoff} distance of their representative.

Multiple cutoffs may be provided in sorted order.  If the \code{cutoff}s are provided in \emph{descending} order then clustering at each new value of \code{cutoff} is continued within the prior \code{cutoff}'s clusters.  In this way clusters at lower values of \code{cutoff} are completely contained within their umbrella clusters at higher values of \code{cutoff}.  This is useful for accelerating the clustering process, because each subsequent group is only clustered within the previous group.  If multiple cutoffs are provided in \emph{ascending} order then clustering at each level of \code{cutoff} is independent of the prior level.

Note that first phase of the algorithm stochastically selects partitions for ordering the sequences.  Hence, clusters can vary from run-to-run unless the random number seed is set for repeatability.
}
\value{
A data.frame is returned with dimensions \eqn{N*M}, where each one of \eqn{N} sequences is assigned to a cluster at the \eqn{M}-level of cutoff.  The row.names of the data.frame correspond to the \emph{names} of \code{myXStingSet}.
}
\author{
Erik Wright \email{eswright@pitt.edu}
}
\examples{
fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
dna <- readDNAStringSet(fas)
aa <- translate(dna)
clusters <- IdClusters(aa, cutoff=seq(0.7, 0.1, -0.1))
head(clusters)
apply(clusters, 2, max) # number of clusters per cutoff
}