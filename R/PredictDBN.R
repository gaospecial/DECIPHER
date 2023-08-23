.getWeights <- function(weight,
	myXStringSet,
	verbose,
	processors) {
	
	l <- length(myXStringSet)
	
	if (length(weight)==1) {
		if (is.numeric(weight)) {
			weight <- rep(1, l)
		} else if (is.na(weight)) {
			if (verbose) {
				cat("Determining distance matrix based on alignment:\n")
				flush.console()
			}
			d <- DistanceMatrix(myXStringSet,
				type="dist",
				penalizeGapLetterMatches=FALSE,
				verbose=verbose,
				processors=processors)
			
			if (verbose) {
				cat("Determining sequence weights:\n")
				flush.console()
			}
			suppressWarnings(guideTree <- TreeLine(myDistMatrix=d,
				method="UPGMA",
				verbose=verbose,
				processors=processors))
			
			.weights <- function(guideTree) {
				# initialize a stack of maximum length (l)
				stack <- vector("list", l)
				visit <- logical(l) # node already visited
				parent <- integer(l) # index of parent node
				index <- integer(l) # index in parent node
				pos <- 1L # current position in the stack
				stack[[pos]] <- guideTree
				while (pos > 0L) { # more nodes to visit
					if (visit[pos]) { # ascending tree
						visit[pos] <- FALSE # reset visit
						dend <- stack[[pos]]
						
						# align subtrees
						treeLengths <- numeric(length(dend))
						inherit <- attr(dend, "inherit")
						if (length(dend) > 1) {
							for (i in seq_len(length(dend))) {
								h <- attr(dend[[i]], "treeLength")
								if (!is.null(h))
									treeLengths[i] <- h
							}
							
							h <- attr(dend, "height")
							members <- vector("list", length(dend))
							for (i in seq_len(length(dend))) {
								m <- unlist(dend[i])
								treeLengths[i] <- treeLengths[i] + h - height[m[1]]
								weight[m] <<- weight[m] + (h - height[m])/length(m)
								height[m] <<- h
								members[[i]] <- m
							}
						} else if (is.leaf(dend)) {
							height[unlist(dend)] <- attr(dend, "height")
						} else { # inherit from subtree
							treeLengths[1] <- attr(dend[[1]], "treeLength")
						}
						
						attr(stack[[pos]], "treeLength") <- sum(treeLengths)
						# replace self in parent
						if (parent[pos] > 0)
							stack[[parent[pos]]][[index[pos]]] <- stack[[pos]]
						pos <- pos - 1L # pop off of stack
					} else { # descending tree
						visit[pos] <- TRUE
						p <- pos
						for (i in seq_along(stack[[p]])) {
							if (!is.leaf(stack[[p]][[i]])) {
								# push subtree onto stack
								pos <- pos + 1L
								stack[[pos]] <- stack[[p]][[i]]
								parent[[pos]] <- p
								index[[pos]] <- i
							}
						}
					}
				}
				
				return(attr(stack[[1]], "treeLength"))
			}
			
			weight <- height <- numeric(l)
			.weights(guideTree)
			weight <- weight/mean(weight)
		} else {
			stop("weight must be a numeric or NA.")
		}
	} else {
		if (length(weight) != l)
			stop("Length of weight must equal one or the length of the myXStringSet.")
		if (any(is.na(weight)))
			stop("weight must be a numeric.")
		if (!isTRUE(all.equal(1, mean(weight))))
			stop("The mean of weight must be 1.")
	}
	
	return(weight)
}

.parseDBN <- function(dbn) {
	s <- strsplit(dbn, "")[[1]]
	n <- sum(s != "." & s != "-")
	c1 <- c2 <- c3 <- integer(n/2)
	i <- 0L
	parens <- square <- curly <- straight <- integer()
	for (j in seq_along(s)) {
		if (s[j] != "." && s[j] != "-") {
			if (s[j]=="(") {
				parens <- c(parens, j)
			} else if (s[j]==")") {
				i <- i + 1L
				c1[i] <- parens[length(parens)]
				c2[i] <- j
				length(parens) <- length(parens) - 1L
			} else if (s[j]=="[") {
				square <- c(square, j)
			} else if (s[j]=="]") {
				i <- i + 1L
				c1[i] <- square[length(square)]
				c2[i] <- j
				c3[i] <- 1L
				length(square) <- length(square) - 1L
			} else if (s[j]=="{") {
				curly <- c(curly, j)
			} else if (s[j]=="}") {
				i <- i + 1L
				c1[i] <- curly[length(curly)]
				c2[i] <- j
				c3[i] <- 2L
				length(curly) <- length(curly) - 1L
			} else if (s[j]=="<") {
				straight <- c(straight, j)
			} else if (s[j]==">") {
				i <- i + 1L
				c1[i] <- straight[length(straight)]
				c2[i] <- j
				c3[i] <- 3L
				length(straight) <- length(straight) - 1L
			} else {
				stop("Unrecognized character in structure: ", s[j])
			}
		}
	}
	
	o <- order(c3, c1, c2)
	matrix(c(c1[o], c2[o], c3[o]),
		ncol=3,
		dimnames=list(NULL, c("(", ")", "order")))
}



















#' Predict RNA Secondary Structure in Dot-Bracket Notation
#' 
#' Predicts a consensus RNA secondary structure from a multiple sequence
#' alignment using mutual information.
#' 
#' \code{PredictDBN} employs an extension of the method described by Freyhult
#' et al. (2005) for determining a consensus RNA secondary structure.  It uses
#' the mutual information (\eqn{H}) measure to find covarying positions in a
#' multiple sequence alignment.  The original method is modified by the
#' addition of different weights for each type of base pairing and each input
#' sequence.  The formula for mutual information between positions \eqn{i} and
#' \eqn{j} then becomes:
#' 
#' \deqn{H(i,j) = \sum_{XY \in bp}^{} \left( impact(XY) \cdot f_{i,j}(XY) \cdot
#' \log_2 \left( \frac{f_{i,j}(XY)}{f_{i}(X) \cdot f_{j}(Y)} \right) \right)}
#' 
#' where, \eqn{bp} denotes the base pairings A/U, C/G, and G/U; \code{impact}
#' is their weight; \eqn{f} is the frequency of single bases or pairs weighted
#' by the corresponding \code{weight} of each sequence.
#' 
#' A penalty is then added for bases that are inconsistent with pairing:
#' 
#' \deqn{H_{mod}(i,j) = H(i,j) + \sum_{XY \notin bp}^{} \Big( impact(XY) \cdot
#' f_{i,j}(XY) \Big)}
#' 
#' Next an average product correction (Buslje et al., 2009) is applied to the
#' matrix \eqn{H}:
#' 
#' \deqn{H_{APC}(i,j) = H_{mod}(i,j) - avgProdCorr \cdot
#' \frac{\overline{H_{mod}(i,.)} \cdot
#' \overline{H_{mod}(.,j)}}{\overline{H_{mod}(.,.)}}}
#' 
#' The mutual information values are then rescaled between \code{0} and
#' \code{1} by applying a sigmoidal transformation, which is controlled by
#' \code{shift} and \code{slope}:
#' 
#' \deqn{H_{final}(i,j) = \left( 1 + \exp \left( slope \cdot log_e \left(
#' \frac{H_{APC}(i,j)}{shift \cdot H_{APC}[n]} \right) \right) \right)^{-1}}
#' 
#' where, \eqn{n} is the number of positions having \code{minOccupancy} divided
#' by two (i.e., the maximum possible number of paired positions) and
#' \eqn{H_{APC}[n]} denotes the \eqn{n^{th}} highest value in the matrix
#' \eqn{H_{APC}}.
#' 
#' If \code{useFreeEnergies} is \code{TRUE}, mutual information is supplemented
#' with a probabalistic model of folding based on \code{deltaGrules}.  That is,
#' palindromes in each sequence are ranked by their free energy, and converted
#' to probabilities of base pairing by assuming an exponential distribution of
#' free energies.  This tends to improve predictive accuracy when the aligned
#' sequences are insufficiently diverse for considerable evidence of
#' compensatory mutations.
#' 
#' If \code{type} is \code{"states"} or \code{"pairs"}, the secondary structure
#' is determined using a variant of the Nussinov algorithm similar to that
#' described by Venkatachalam et al. (2014).  Pairings with a score below
#' \code{threshold} are not considered during the traceback.  If
#' \code{psuedoknots} is greater than \code{0}, paired positions are removed
#' from consideration and the method is applied again to find pseudoknots.
#' 
#' In practice the secondary structure prediction is most accurate when the
#' input alignment is of high quality, contains a wide diversity of sequences,
#' the number of sequences is large, no regions are completely conserved across
#' all sequences, and most of the sequences span the entire alignment (i.e.,
#' there are few partial/incomplete sequences).
#' 
#' @name PredictDBN
#' @param myXStringSet A \code{DNAStringSet} or \code{RNAStringSet} object
#' containing aligned sequences.
#' @param type Character string indicating the type of results desired.  This
#' should be (an unambiguous abbreviation of) one of \code{"states"},
#' \code{"pairs"}, \code{"evidence"}, \code{"scores"}, \code{"structures"}, or
#' \code{"search"}.  (See value section below.)
#' @param minOccupancy Numeric specifying the minimum occupancy (1 - fraction
#' of gaps) required to include a column of the alignment in the prediction.
#' @param impact A vector with four elements giving the weights of A/U, G/C,
#' G/U, and other pairings, respectively.  The last element of \code{impact} is
#' the penalty for pairings that are inconsistent with two positions being
#' paired (e.g., A/- or A/C).
#' @param avgProdCorr Numeric specifying the weight of the average product
#' correction (APC) term, as described in Buslje et al. (2009).
#' @param slope Numeric giving the slope of the sigmoid used to convert mutual
#' information values to scores ranging from zero to one.
#' @param shift Numeric giving the relative shift of the sigmoid used to
#' convert mutual information values to scores ranging from zero to one.
#' @param threshold Numeric specifying the score threshold at which to consider
#' positions for pairing.  Only applicable if \code{type} is \code{"states"} or
#' \code{"pairs"}.
#' @param pseudoknots Integer indicating the maximum order of pseudoknots that
#' are acceptable.  A value of \code{0} will prevent pseudoknots in the
#' structure, whereas \code{1} (the default) will search for first-order
#' psuedoknots.  Only used if \code{type} is \code{"states"} or \code{"pairs"}.
#' @param weight Either a numeric vector of weights for each sequence, a single
#' number implying equal weights, or \code{NA} (the default) to automatically
#' calculate sequence weights based on \code{myXStringSet}.
#' @param useFreeEnergy Logical determining whether RNA free energy predictions
#' should be incorporated along with mutual information into the secondary
#' structure prediction.
#' @param deltaGrules Free energy rules for all possible base pairings in
#' quadruplets.  If NULL, defaults to pseudoenergies
#' (\code{\link{deltaGrulesRNA}}).  Only applicable if \code{useFreeEnergies}
#' is \code{TRUE}.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @return If \code{type} is \code{"states"} (the default), then the output is
#' a character vector with the predicted secondary structure assignment for
#' each position in \code{myXStringSet}.  Standard dot-bracket notation (DBN)
#' is used, where ``.'' signifies an unpaired position, ``('' and ``)'' a
#' paired position, and successive ``[]'', ``{}'', and ``<>'' indicate
#' increasing order pseudoknots.  Columns below \code{minOccupancy} are denoted
#' by the ``-'' character to indicate that they contained too many gaps to be
#' included in the consensus structure.
#' 
#' If \code{type} is \code{"pairs"}, then a matrix is returned with one row for
#' each base pairing and three columns giving the positions of the paired bases
#' and their \code{pseudoknot} order.
#' 
#' If \code{type} is \code{"evidence"}, then a matrix is returned with one row
#' for each base pairing and three columns giving the positions of the paired
#' bases and their respective scores (greater than or equal to
#' \code{threshold}).  This differs from \code{type} \code{"pairs"} in that
#' \code{"evidence"} does not perform a traceback.  Therefore, it is possible
#' to have conflicting evidence where a single base has evidence for pairing
#' with multiple other bases.
#' 
#' If \code{type} is \code{"scores"}, then a matrix of three rows is returned,
#' where the values in a column represent the maximum score for a state in each
#' position.  Columns sum to \code{1} if the position was above
#' \code{minOccupancy} and \code{0} otherwise.
#' 
#' If \code{type} is \code{"structures"}, then the output is a list with one
#' element for each sequence in \code{myXStringSet}.  Each list element
#' contains a matrix of dimension 3 (each state) by the number of nucleotides
#' in the sequence.  Columns of the matrix sum to zero where the nucleotide was
#' located in a position that was below \code{minOccupancy}.  Otherwise,
#' positions are considered paired if they are consistent with pairing (i.e.,
#' A/U, C/G, or G/U) in the consensus secondary structure.
#' 
#' If \code{type} is \code{"search"} then an attempt is made to find additional
#' secondary structure beyond positions exhibiting covariation.  First, anchors
#' are identified as pairs of covarying positions with their score above
#' \code{threshold}.  Next, the regions between anchors are searched for
#' previously unidentified stem loops.  Finally, any helices are assigned a
#' score according to their length, i.e. one minus the probability of finding
#' that many consecutive pairs within the anchor boundaries by chance.  Hence,
#' output \code{type} \code{"search"} will find secondary structure outside of
#' the consensus structure shared by most sequences, and can identify secondary
#' structure in conserved alignment regions.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{PredictHEC}}
#' @references Buslje, C., et al. (2009). Correction for phylogeny, small
#' number of observations and data redundancy improves the identification of
#' coevolving amino acid pairs using mutual information. Bioinformatics,
#' \bold{25(9)}, 1125-1131.
#' 
#' Freyhult, E., et al. (2005). Predicting RNA Structure Using Mutual
#' Information. Applied Bioinformatics, \bold{4(1)}, 53-59.
#' 
#' Venkatachalam, B., et al. (2014). Faster algorithms for RNA-folding using
#' the Four-Russians method. Algorithms for Molecular Biology : AMB,
#' \bold{9(1)}, 1-12.
#' 
#' Wright, E. S. (2020). RNAconTest: comparing tools for noncoding RNA multiple
#' sequence alignment based on structural consistency. RNA 2020, 26, 531-540.
#' @examples
#' 
#' # load the example non-coding RNA sequences
#' db <- system.file("extdata", "Bacteria_175seqs.sqlite", package="DECIPHER")
#' rna <- SearchDB(db, type="RNAStringSet")
#' 
#' \dontrun{
#' # predict the secondary structure in dot-bracket notation (dbn)
#' p <- PredictDBN(rna, "states") # predict the secondary structure in dbn
#' p # pairs are denoted by (), and (optionally) pseudoknots by [], {}, and <>
#' 
#' # convert the dot-bracket notation into pairs of positions within the alignment
#' p <- PredictDBN(rna, "pairs") # paired positions in the alignment
#' head(p) # matrix giving the pairs and their pseudoknot order (when > 0)
#' 
#' # plot an arc diagram with the base pairings
#' plot(NA, xlim=c(0, 1), ylim=c(0, 1),
#' 	xaxs="i", yaxs="i",
#' 	xlab="Alignment position", ylab="",
#' 	bty="n", xaxt="n", yaxt="n")
#' ticks <- pretty(seq_len(width(rna)[1]))
#' axis(1, ticks/width(rna)[1], ticks)
#' rs <- c(seq(0, pi, len=100), NA)
#' r <- (p[, 2] - p[, 1] + 1)/width(rna)[1]/2
#' r <- rep(r, each=101)
#' x <- (p[, 1] + p[, 2])/2/width(rna)[1]
#' x <- rep(x, each=101) + r*cos(rs)
#' y <- r*sin(rs)/max(r, na.rm=TRUE)
#' lines(x, y, xpd=TRUE)
#' 
#' # show all available evidence of base pairing
#' p <- PredictDBN(rna, "evidence") # all pairs with scores >= threshold
#' head(p) # matrix giving the pairs and their scores
#' 
#' # determine the score at every alignment position
#' p <- PredictDBN(rna, "scores") # score in the alignment
#' p["(", 122] # score for left-pairing at alignment position 122
#' p[")", 260] # score for right-pairing at alignment position 260
#' 
#' # find the scores individually for every sequence in the alignment
#' p <- PredictDBN(rna, "structures") # scores per sequence
#' p[[1]][, 1] # the scores for the first position in the first sequence
#' p[[2]][, 10] # the scores for the tenth position in the second sequence
#' # these positional scores can be used as shades of red, green, and blue:
#' BrowseSeqs(rna, patterns=p) # red = unpaired, green = left-pairing, blue = right
#' # positions in black are not part of the consensus secondary structure
#' 
#' # search for additional secondary structure between the consensus pairs
#' p <- PredictDBN(rna, "search") # scores per sequence after searching
#' BrowseSeqs(rna, patterns=p) # red = unpaired, green = left-pairing, blue = right
#' # note that "search" identified many more pairings than "structures"
#' }
#' 
#' @export PredictDBN
PredictDBN <- function(myXStringSet,
	type="states",
	minOccupancy=0.4,
	impact=c(1, 1.2, 0.4, -1),
	avgProdCorr=1,
	slope=2,
	shift=1.3,
	threshold=0.4,
	pseudoknots=1,
	weight=NA,
	useFreeEnergy=TRUE,
	deltaGrules=NULL,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!(is(myXStringSet, "DNAStringSet") ||
		is(myXStringSet, "RNAStringSet")))
		stop("myXStringSet must be a DNAStringSet or RNAStringSet.")
	l <- length(myXStringSet)
	if (l < 2)
		stop("myXStringSet must contain at least two sequences.")
	u <- unique(width(myXStringSet))
	if (length(u)!=1)
		stop("Sequences in myXStringSet must be the same width (aligned).")
	TYPES <- c("states", "pairs", "scores", "structures", "search", "evidence")
	if (length(type)==0)
		stop("No type specified.")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (u==0) { # zero width sequences
		if (type==1L) {
			return("")
		} else if (type==2L) {
			return(matrix(nrow=0, ncol=3, dimnames=list(NULL, c("(", ")", "order"))))
		} else if (type==3L) {
			return(matrix(nrow=3, ncol=0, dimnames=list(c(".", "(", ")"), NULL)))
		} else if (type==6L) {
			return(matrix(nrow=0, ncol=3, dimnames=list(NULL, c("(", ")", "score"))))
		} else {
			x <- matrix(0, nrow=3, ncol=0, dimnames=list(c(".", "(", ")"), NULL))
			x <- lapply(seq_len(l),
				function(...)
					return(x))
			return(x)
		}
	}
	if (!is.numeric(minOccupancy))
		stop("minOccupancy must be a numeric.")
	if (minOccupancy < 0 || minOccupancy > 1)
		stop("minOccupancy must be between zero and one.")
	if (!is.numeric(impact) || length(impact) != 4)
		stop("impact must be a numeric vector of length four.")
	if (impact[4] > 0)
		stop("The impact of inconsistent pairs (element four) can be at most zero.")
	if (any(impact[1:2] < 0))
		stop("The impact of A/T and G/C base pairs must be at least zero.")
	if (!is.numeric(avgProdCorr))
		stop("avgProdCorr must be a numeric.")
	if (avgProdCorr < 0)
		stop("avgProdCorr must be at least zero.")
	if (!is.numeric(slope))
		stop("slope must be a numeric.")
	if (slope <= 0)
		stop("slope must be greater than zero.")
	if (!is.numeric(shift))
		stop("shift must be a numeric.")
	if (!is.numeric(threshold))
		stop("threshold must be a numeric.")
	if (threshold < 0 || threshold > 1)
		stop("threshold must be between zero and one.")
	if (!is.numeric(pseudoknots))
		stop("pseudoknots must be a numeric.")
	if (pseudoknots < 0 || pseudoknots > 3)
		stop("pseudoknots must be between zero and three.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (!is.null(processors) && !is.numeric(processors))
		stop("processors must be a numeric.")
	if (!is.null(processors) && floor(processors)!=processors)
		stop("processors must be a whole number.")
	if (!is.null(processors) && processors < 1)
		stop("processors must be at least 1.")
	if (is.null(processors)) {
		processors <- detectCores()
	} else {
		processors <- as.integer(processors)
	}
	
	weight <- .getWeights(weight,
		myXStringSet,
		verbose,
		processors)
	
	if (useFreeEnergy) {
		# initialize a progress bar
		if (verbose) {
			time.1 <- Sys.time()
			cat("Computing Free Energies:\n")
			pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
		}
		
		if (is.null(deltaGrules)) {
			data("deltaGrulesRNA", envir=environment(), package="DECIPHER")
			deltaGrules <- deltaGrulesRNA
		} else {
			if (!is.numeric(deltaGrules))
				stop("deltaGrules must be numeric.")
			if (length(deltaGrules) != 390625L)
				stop("deltaGrules must be of dimensions 5 x 5 x 5 x 5 x 5 x 5 x 5 x 5.")
		}
		
		pos <- as.matrix(myXStringSet)
		pos <- lapply(seq_len(nrow(pos)),
			function(x)
				which(!(pos[x,] %in% c("-", "."))))
		
		noGaps <- RemoveGaps(myXStringSet)
		J <- W <- vector( "list", u)
		for (i in seq_along(myXStringSet)) {
			pals <- findPalindromes(noGaps[[i]],
				min.armlength=4,
				max.looplength=500,
				min.looplength=3,
				max.mismatch=1,
				allow.wobble=TRUE)
			
			dna <- DNAStringSet(pals)
			s1 <- as.character(subseq(dna, 1L, 1L))
			e1 <- as.character(subseq(dna, -1L))
			w <- which((s1 == "A" & e1 == "T") |
				(s1 == "T" & e1 == "A") |
				(s1 == "C" & e1 == "G") |
				(s1 == "G" & e1 == "C") |
				(s1 == "G" & e1 == "T") |
				(s1 == "T" & e1 == "G"))
			pals <- pals[w]
			dna <- dna[w]
			
			arms <- palindromeArmLength(dna,
				max.mismatch=0,
				allow.wobble=TRUE)
			arms1 <- palindromeArmLength(dna,
				max.mismatch=1,
				allow.wobble=TRUE)
			w <- which(arms + 1L != arms1)
			if (length(w) > 0)
				arms[w] <- arms1[w] # mismatch not at end
			max_arms <- as.integer((width(pals) - 3)/2)
			w <- which(arms > max_arms)
			if (length(w) > 0)
				arms[w] <- max_arms[w]
			
			s1 <- start(pals)
			s2 <- start(pals) + arms - 1L
			e1 <- end(pals) - arms + 1L
			e2 <- end(pals)
			p <- pos[[i]]
			
			w <- which(width(dna) - 3L > 2L*arms)
			arms[w] <- arms[w] + 1L # include end bases
			
			dG <- .Call("calculateHairpinDeltaG",
				dna,
				arms,
				deltaGrulesRNA,
				PACKAGE="DECIPHER")
			o <- order(dG)
			o <- o[dG[o] < 0]
			dG <- dG/mean(dG[o][seq_len(sum(cumsum(arms[o]) < length(p)))])
			dG <- 1 - exp(-dG) # convert to probability
			
			for (j in seq_along(o)) {
				p1 <- s1[o[j]]:s2[o[j]]
				p2 <- e2[o[j]]:e1[o[j]]
				p1 <- p[p1]
				p2 <- p[p2]
				for (k in seq_along(p1)) {
					m <- match(p2[k], J[[p1[k]]])
					if (is.na(m)) {
						J[[p1[k]]] <- c(J[[p1[k]]], p2[k])
						W[[p1[k]]] <- c(W[[p1[k]]], dG[o[j]]*weight[i])
					} else {
						W[[p1[k]]][m] <- W[[p1[k]]][m] + dG[o[j]]*weight[i]
					}
				}
			}
			
			if (verbose)
				setTxtProgressBar(pBar, i/l)
		}
		
		t <- TerminalChar(myXStringSet)
		w <- numeric(u)
		for (i in seq_along(myXStringSet))
			if (t[i, 3L] > 0)
				w[(t[i, 1L] + 1):(u - t[i, 2L])] <- w[(t[i, 1L] + 1):(u - t[i, 2L])] + weight[i]
		for (i in seq_along(J))
			W[[i]] <- W[[i]]/pmin(w[J[[i]]], w[i])
		
		for (i in seq_along(J)) {
			if (length(J[[i]]) > 0) {
				o <- order(J[[i]])
				J[[i]] <- J[[i]][o]
				W[[i]] <- W[[i]][o]
			}
		}
		patterns <- list(J, W)
		
		if (verbose) {
			close(pBar)
			time.2 <- Sys.time()
			cat("\n")
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\n")
		}
	} else {
		patterns <- NULL
	}
	
	# initialize a progress bar
	if (verbose) {
		time.1 <- Sys.time()
		cat("Predicting RNA Secondary Structures:\n")
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
	} else {
		pBar <- NULL
	}
	
	ans <- .Call("predictDBN",
		myXStringSet,
		type,
		minOccupancy,
		impact,
		avgProdCorr,
		slope,
		shift,
		weight,
		pseudoknots,
		threshold,
		patterns,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (type==2L) {
		ans <- .parseDBN(ans)
	} else if (type==3L) {
		rownames(ans) <- c(".", "(", ")")
	} else if (type==4L) {
		ans <- lapply(ans,
			function(x) {
				rownames(x) <- c(".", "(", ")")
				return(x)
			})
	} else if (type==6L) {
		colnames(ans) <- c("(", ")", "score")
	}
	
	if (verbose) {
		setTxtProgressBar(pBar, 100)
		close(pBar)
		time.2 <- Sys.time()
		cat("\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
	}
	
	return(ans)
}
