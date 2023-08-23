#' Predict Protein Secondary Structure as Helix, Beta-Sheet, or Coil
#' 
#' Predicts 3-state protein secondary structure based on the primary (amino
#' acid) sequence using the GOR IV method (Garnier et al., 1996).
#' 
#' The GOR (Garnier-Osguthorpe-Robson) method is an information-theory method
#' for prediction of secondary structure based on the primary sequence of a
#' protein.  Version IV of the method makes 3-state predictions based on the
#' mutual information contained in single residues and pairs of residues within
#' \code{windowSize} residues of the position being assigned.  This approach is
#' about 65\% accurate, and is one of the most accurate methods for assigning
#' secondary structure that only use a single sequence.  This implementation of
#' GOR IV does not use decision constants or the number of contiguous states
#' when assigning the final state.  Note that characters other than the
#' standard 20 amino acids are not assigned a state.
#' 
#' @name PredictHEC
#' @param myAAStringSet An \code{AAStringSet} object of sequences.
#' @param type Character string indicating the type of results desired.  This
#' should be (an unambiguous abbreviation of) one of \code{"states"},
#' \code{"scores"}, or \code{"probabilities"}.
#' @param windowSize Numeric specifying the number of residues to the left or
#' right of the center position to use in the prediction.
#' @param background Numeric vector with the background ``scores'' for each of
#' the three states (H, E, and C).
#' @param HEC_MI1 An array of dimensions 20 x 21 x 3 giving the mutual
#' information for single residues.
#' @param HEC_MI2 An array of dimensions 20 x 20 x 21 x 21 x 3 giving the
#' mutual information for pairs of residues.
#' @return If \code{type} is \code{"states"} (the default), then the output is
#' a character vector with the secondary structure assignment ("H", "E", or
#' "C") for each residue in \code{myAAStringSet}.
#' 
#' Otherwise, the output is a list with one element for each sequence in
#' \code{myAAStringSet}.  Each list element contains a matrix of dimension 3
#' (H, E, or C) by the number of residues in the sequence.  If \code{type} is
#' \code{"scores"}, then values in the matrix represent log-odds ``scores''.
#' If \code{type} is \code{"probabilities"} then the values represent the
#' normalized probabilities of the three states at a position.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{HEC_MI1}}, \code{\link{HEC_MI2}},
#' \code{\link{PredictDBN}}
#' @references Garnier, J., Gibrat, J. F., & Robson, B. (1996). GOR method for
#' predicting protein secondary structure from amino acid sequence.
#' \emph{Methods in Enzymology}, \bold{266}, 540-553.
#' @examples
#' 
#' fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' aa <- translate(dna)
#' hec <- PredictHEC(aa)
#' head(hec)
#' 
#' @export PredictHEC
PredictHEC <- function(myAAStringSet,
	type="states",
	windowSize=7,
	background=c(H=-0.12, E=-0.25, C=0.23),
	HEC_MI1=NULL,
	HEC_MI2=NULL) {
	
	# error checking
	if (!is(myAAStringSet, "AAStringSet"))
		stop("myAAStringSet must be an AAStringSet.")
	TYPES <- c("states", "scores", "probabilities")
	if (length(type)==0)
		stop("No type specified.")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (!is.numeric(windowSize))
		stop("windowSize must be a numeric.")
	if (floor(windowSize) != windowSize)
		stop("windowSize must be an integer number.")
	if (windowSize < 1)
		stop("windowSize must be at least 1.")
	if (!is.double(background))
		stop("background must be a numeric.")
	if (length(background) != 3)
		stop("background must be length 3.")
	if (is.null(HEC_MI1)) {
		data("HEC_MI1", envir=environment(), package="DECIPHER")
	} else {
		if (!is.double(HEC_MI1))
			stop("HEC_MI1 must be an array of numerics.")
		if (length(dim(HEC_MI1)) != 3)
			stop("HEC_MI1 must be a three dimensional array.")
		if (dim(HEC_MI1)[1] != 20 ||
			windowSize > ((dim(HEC_MI1)[2] - 1)/2) ||
			(dim(HEC_MI1)[2] %% 2) != 1 ||
			dim(HEC_MI1)[3] != 3)
			stop("HEC_MI1 must have dimensions 20 x (2*windowSize + 1) x 3.")
	}
	if (is.null(HEC_MI2)) {
		data("HEC_MI2", envir=environment(), package="DECIPHER")
	} else {
		if (!is.double(HEC_MI2))
			stop("HEC_MI2 must be an array of numerics.")
		if (length(dim(HEC_MI2)) != 5)
			stop("HEC_MI2 must be a three dimensional array.")
		if (dim(HEC_MI2)[1] != 20 ||
			dim(HEC_MI2)[2] != 20 ||
			windowSize > ((dim(HEC_MI2)[3] - 1)/2) ||
			dim(HEC_MI2)[3] != dim(HEC_MI2)[4] ||
			(dim(HEC_MI2)[3] %% 2) != 1 ||
			dim(HEC_MI2)[5] != 3)
			stop("HEC_MI2 must have dimensions 20 x 20 x (2*windowSize + 1) x(2*windowSize + 1) x 3.")
	}
	
	states <- .Call("predictHEC",
		myAAStringSet,
		windowSize,
		background,
		HEC_MI1,
		HEC_MI2,
		type)
	
	if (type > 1) {
		states <- lapply(states, function(x) {
			rownames(x) <- c("H", "E", "C")
			return(x)
		})
	}
	
	return(states)
}
