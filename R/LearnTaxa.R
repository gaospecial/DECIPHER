#' Train a Classifier for Assigning Taxonomy
#' 
#' Trains a classifier based on a reference taxonomy containing sequence
#' representatives assigned to taxonomic groups.
#' 
#' Learning about the training data is a two part process consisting of (i)
#' forming a taxonomic tree and then (ii) ensuring that the \code{train}ing
#' sequences can be correctly reclassified.  The latter step relies on
#' reclassifying the sequences in \code{train} by descending the taxonomic
#' tree, a process termed ``tree descent''.  Ultimately, the goal of tree
#' descent is to quickly and accurately narrow the selection of groups where a
#' sequence may belong.  During the learning process, tree descent is tuned so
#' that it performs well when classifying new sequences.
#' 
#' The process of training the classifier first involves learning the taxonomic
#' tree spanning all of the reference sequences in \code{train}.  Typically,
#' reference taxonomic classifications are provided by an authoritative source,
#' oftentimes along with a ``taxid'' file containing taxonomic \code{rank}
#' information.  The taxonomic tree may contain any number of levels (e.g.,
#' Root, Phylum, Class, Order, Family, Genus) as long as they are
#' hierarchically nested and always begin with ``Root''.
#' 
#' The second phase of training the classifier, tree descent, involves learning
#' the optimal set of k-mers for discerning between the different sub-groups
#' under each edge.  Here a fraction of the k-mers with the greatest discerning
#' power are matched to a training sequence, and this process is repeated with
#' 100 random subsamples to decide on the set of possible taxonomic groups to
#' which a training sequence may belong.
#' 
#' The learning process works by attempting to correctly re-classify each
#' training sequence in the taxonomy.  Initially, \code{maxFraction} of
#' informative k-mers are repeatedly sampled at each edge during tree descent.
#' Training sequences that are incorrectly classified at an edge will lower the
#' fraction of k-mers that are sampled by an amount that is proportional to
#' \code{multiplier}.  As the fraction of sampled k-mers decreases, the tree
#' descent process terminates at higher rank levels.
#' 
#' A major advantage of tree descent is that it both speeds up the
#' classification process and indicates where the training set likely contains
#' mislabeled sequences or incorrectly-placed taxonomic groups.  Training
#' sequences that are not correctly classified within \code{maxIterations} are
#' marked as ``problem sequences'', because it is likely that they are
#' mislabeled.  If enough sequences have difficulty being correctly classified
#' at an edge that the fraction drops below \code{minFraction}, then the edge
#' is recorded as a ``problem group''.
#' 
#' The final result is an object that can be used for classification with
#' \code{\link{IdTaxa}}, as well as information about \code{train} that could
#' be used to help correct any errors in the taxonomy.
#' 
#' @name LearnTaxa
#' @param train An \code{AAStringSet}, \code{DNAStringSet}, or
#' \code{RNAStringSet} of unaligned sequences.
#' @param taxonomy Character string providing the reference taxonomic
#' assignment for each sequence in \code{train}.  Taxonomic ranks are separated
#' by semicolons (``;'') beginning with ``Root''.
#' @param rank Optionally, a \code{data.frame} with 5 named columns giving the
#' ``Index'' (i.e., 0 to the number of unique taxa), ``Name'' (i.e., taxon
#' name), ``Parent'' (i.e., ``Index'' of the parent taxon), ``Level'' (i.e.,
#' integer rank level), and ``Rank'' (e.g., ``genus'') of each taxonomic rank.
#' This information is often provided in a separate ``taxid'' file along with
#' publicly available training sequence sets.
#' @param K Integer specifying the k-mer size or \code{NULL} (the default) to
#' calculate the k-mer size automatically.  The default value of \code{K} is
#' such that matches between sequences are found by chance every \code{N}
#' k-mers.
#' @param N Numeric indicating the approximate number of k-mers that can be
#' randomly selected before one is found by chance on average. For example, the
#' default value of \code{500} will set \code{K} (when \code{K} is unspecified)
#' so that every 500th k-mer is expected to match by chance.
#' @param minFraction Numeric giving the minimum fraction of k-mers to sample
#' during the initial tree descent phase of the classification algorithm.  (See
#' details section below.)
#' @param maxFraction Numeric giving the maximum fraction of k-mers to sample
#' during the initial tree descent phase of the classification algorithm.  (See
#' details section below.)
#' @param maxIterations Integer specifying the maximum number of iterations to
#' attempt re-classification of a training sequence before declaring it a
#' ``problem sequence''.  (See details section below.)
#' @param multiplier Numeric indicating the degree to which individual
#' sequences have control over the fraction of k-mers sampled at any edge
#' during the initial tree descent phase of the classification algorithm.  (See
#' details section below.)
#' @param maxChildren Integer giving the maximum number of child taxa of any
#' taxon at which to consider further descending the taxonomic tree.  A value
#' of \code{1} will prevent use of the tree descent algorithm altogether.
#' Lower values may decrease classification speed, but result in output objects
#' that require less memory.
#' @param alphabet Character vector of amino acid groupings used to reduce the
#' 20 standard amino acids into smaller groups.  Alphabet reduction helps to
#' find more distant homologies between sequences.  A non-reduced amino acid
#' alphabet can be used by setting \code{alphabet} equal to \code{AA_STANDARD}.
#' @param verbose Logical indicating whether to display progress.
#' @return An object of class \code{Taxa} and subclass Train, which is stored
#' as a list with components: \item{taxonomy}{ A character vector containing
#' all possible groups in the taxonomy. } \item{taxa}{ A character vector
#' containing the basal taxon in each taxonomy. } \item{ranks}{ A character
#' vector of rank names for each taxon, or \code{NULL} if \code{rank}
#' information was not supplied. } \item{levels}{ Integer giving the rank level
#' of each taxon. } \item{children}{ A list containing the index of all
#' children in the taxonomy for each taxon. } \item{parents}{ An integer
#' providing the index of the parent for each taxon. } \item{fraction}{ A
#' numeric between \code{minFraction} and \code{maxFraction} that represents
#' the learned fraction of informative k-mers to sample for each taxon during
#' the initial tree descent phase of the classification algorithm.  Problem
#' groups are marked by a fraction of \code{NA}. } \item{sequences}{ List
#' containing the integer indices of sequences in \code{train} belonging to
#' each taxon. } \item{kmers}{ List containing the unique sorted k-mers
#' (converted to integers) belonging to each sequence in \code{train}. }
#' \item{crossIndex}{ Integer indicating the index in taxonomy of each
#' sequence's taxonomic label. } \item{K}{ The value of \code{K} provided as
#' input. } \item{IDFweights}{ Numeric vector of length \code{4^K} providing
#' the inverse document frequency weight for each k-mer. }
#' \item{decisionKmers}{ List of informative k-mers and their associated
#' relative frequencies for each internal edge in the taxonomy. }
#' \item{problemSequences}{ A \code{data.frame} providing the ``Index'',
#' ``Expected'' label, and ``Predicted'' taxon for sequences that could not be
#' correctly classified during the initial tree descent phase of the algorithm.
#' } \item{problemGroups}{ Character vector containing any taxonomic groups
#' that repeatedly had problems with correctly re-classifying sequences in
#' \code{train} during the initial tree descent phase of the classification
#' algorithm.  Problem groups likely indicate that a number of the sequences
#' (or an entire group of sequences) assigned to the problem group are
#' incorrectly placed in the taxonomic tree. } \item{alphabet}{ The
#' \code{alphabet} when \code{train} is an ``AAStringSet'', otherwise
#' \code{NULL}. }
#' @note If \code{K} is \code{NULL}, the automatically determined value of
#' \code{K} might be too large for some computers, resulting in an error.  In
#' such cases it is recommended that \code{K} be manually set to a smaller
#' value.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{IdTaxa}}, \code{\link{Taxa-class}}
#' @references Murali, A., et al. (2018). IDTAXA: a novel approach for accurate
#' taxonomic classification of microbiome sequences. Microbiome, 6, 140.
#' https://doi.org/10.1186/s40168-018-0521-5
#' @examples
#' 
#' # import training sequences
#' fas <- system.file("extdata", "50S_ribosomal_protein_L2.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' 
#' # parse the headers to obtain a taxonomy
#' s <- strsplit(names(dna), " ")
#' genus <- sapply(s, `[`, 1)
#' species <- sapply(s, `[`, 2)
#' taxonomy <- paste("Root", genus, species, sep="; ")
#' head(taxonomy)
#' 
#' # train the classifier
#' trainingSet <- LearnTaxa(dna, taxonomy)
#' trainingSet
#' 
#' # view information about the classifier
#' plot(trainingSet)
#' 
#' \dontrun{
#' # train the classifier with amino acid sequences
#' aa <- translate(dna)
#' trainingSetAA <- LearnTaxa(aa, taxonomy)
#' trainingSetAA
#' }
#' 
#' @export LearnTaxa
LearnTaxa <- function(train,
	taxonomy,
	rank=NULL,
	K=NULL,
	N=500,
	minFraction=0.01,
	maxFraction=0.06,
	maxIterations=10,
	multiplier=100,
	maxChildren=200,
	alphabet=AA_REDUCED[[139]],
	verbose=TRUE) {
	
	# error checking
	if (!is(train, "DNAStringSet") && !is(train, "RNAStringSet") && !is(train, "AAStringSet"))
		stop("train must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	l <- length(train)
	if (l < 2)
		stop("At least two training sequences are required.")
	if (!is.character(taxonomy))
		stop("taxonomy must be a character vector.")
	if (length(taxonomy) != l)
		stop("taxonomy must be the same length as train.")
	if (!is.null(rank)) {
		if (!is.data.frame(rank))
			stop("rank must be a data.frame.")
		ns <- c('Index', 'Name', 'Parent', 'Level', 'Rank')
		w <- which(!(ns %in% names(rank)))
		if (length(w) > 0)
			stop("rank is missing columns named: ",
				paste(ns[w], collapse=", "))
		if (any(duplicated(rank$Index)))
			stop("All values of Index within rank must be unique.")
	}
	if (any(!grepl("Root;", taxonomy, fixed=TRUE)))
		stop("All elements of taxonomy must contain 'Root;'.")
	if (is(train, "AAStringSet")) {
		if (!is.character(alphabet))
			stop("alphabet must be a character vector.")
		if (any(alphabet==""))
			stop("No elements of alphabet can be empty.")
		r <- strsplit(alphabet, "", fixed=TRUE)
		alphabet <- setNames(rep(0L, 20),
			AA_STANDARD)
		for (i in seq_along(r)) {
			w <- which(!(r[[i]] %in% AA_STANDARD))
			if (length(w) > 0)
				stop("Unrecognized letter(s) found in alphabet:  ",
					paste(r[[i]][w], collapse=", "),
					".")
			w <- which(alphabet[r[[i]]] != 0L)
			if (length(w) > 0)
				stop("Repeated amino acids found in alphabet:  ",
					paste(r[[i]][w], collapse=", "),
					".")
			alphabet[r[[i]]] <- i
		}
		w <- which(alphabet==0L)
		if (length(w) > 0)
			stop("Standard amino acids missing from alphabet:  ",
				paste(names(w), collapse=", "),
				".")
		size <- max(alphabet)
		if (size==1)
			stop("More than one grouping of amino acids is required in the alphabet.")
		alphabet <- alphabet - 1L
	} else {
		alphabet <- NULL
		size <- 4
	}
	
	if (is.null(K)) {
		if (!is.numeric(N))
			stop("N must be a numeric.")
		if (N <= 1)
			stop("N must be greater than one.")
		quant <- quantile(width(train), 0.99)
		if (is(train, "AAStringSet")) {
			K <- floor(log(N*quant,
				.Call("alphabetSizeReducedAA",
					train,
					alphabet,
					PACKAGE="DECIPHER")))
		} else {
			K <- floor(log(N*quant,
				.Call("alphabetSize",
					train,
					PACKAGE="DECIPHER")))
		}
		if (is(train, "AAStringSet")) {
			maxK <- as.integer(floor(log(1e8, size))) # theoretically 4294967295, but may encounter memory errors
		} else {
			maxK <- 13 # theoretically 15, but may encounter memory errors
		}
		if (K < 1) {
			K <- 1
		} else if (K > maxK) {
			K <- maxK
			warning("K may be smaller than ideal.")
		}
	} else {
		if (!is.numeric(K))
			stop("K must be a numeric.")
		if (length(K) != 1)
			stop("K must be a single numeric.")
		if (K != floor(K))
			stop("K must be a whole number.")
		if (K < 1)
			stop("K must be at least one.")
		if (is(train, "AAStringSet")) {
			maxK <- as.integer(floor(log(4294967295, size)))
		} else {
			maxK <- 15
		}
		if (K > maxK)
			stop(paste("K can be at most ", maxK, ".", sep=""))
	}
	if (!is.numeric(minFraction))
		stop("minFraction must be a numeric.")
	if (minFraction <= 0)
		stop("minFraction must be greater than zero.")
	if (!is.numeric(maxFraction))
		stop("maxFraction must be a numeric.")
	if (maxFraction <= minFraction)
		stop("maxFraction must be greater than minFraction.")
	if (maxFraction >= 1)
		stop("maxFraction must be less than one.")
	if (!is.numeric(maxIterations))
		stop("maxIterations must be a numeric.")
	if (maxIterations < 1)
		stop("maxIterations must be at least one.")
	if (maxIterations != floor(maxIterations))
		stop("maxIterations must be a whole number.")
	if (!is.numeric(multiplier))
		stop("multiplier must be a numeric.")
	if (multiplier <= 0)
		stop("multiplier must be greater than zero.")
	if (!is.numeric(maxChildren))
		stop("maxChildren must be a numeric.")
	if (maxChildren < 1)
		stop("maxChildren must be at least one.")
	if (maxChildren != floor(maxChildren))
		stop("maxChildren must be a whole number.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- vcountPattern("-", train)
	if (any(a > 0))
		stop("Gap characters ('-') are not permitted in train.")
	a <- vcountPattern(".", train)
	if (any(a > 0))
		stop("Unknown characters ('.') are not permitted in train.")
	
	if (verbose) {
		num <- 0L
		time.1 <- Sys.time()
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	}
	
	# set default values
	nKmers <- size^K
	B <- 100L
	
	# get the full classification beginning from "Root"
	classes <- gsub(" *; *", ";", taxonomy)
	classes <- sapply(strsplit(classes, "Root;"),
		`[`,
		2L)
	classes <- gsub(" +$", "", classes)
	classes <- paste(classes,
		ifelse(endsWith(classes, ";"),
			"",
			";"),
		sep="")
	
	# index and sort all unique (non-ambiguous) k-mers
	if (is(train, "AAStringSet")) {
		kmers <- .Call("enumerateSequenceReducedAA",
			train,
			K,
			alphabet,
			FALSE, # mask repeats
			PACKAGE="DECIPHER")
	} else {
		kmers <- .Call("enumerateSequence",
			train,
			K,
			FALSE, # mask repeats
			PACKAGE="DECIPHER")
	}
	kmers <- lapply(kmers,
		function(x)
			sort(unique(x + 1L), na.last=NA))
	if (any(lengths(kmers)==0))
		stop("All training sequences must have at least one k-mer.")
	
	# create lists of:
	# taxonomy = all distinct taxa ranks
	# taxa = the taxon added in each taxonomy
	# crossIndex = index of taxonomy for each sequence
	# children = index (in taxonomy) of immediate children for each taxa
	# parents = index (in taxonomy) of the parent level
	# sequences = index (in classes/train/kmers) of sequences belonging to each taxa
	u_classes <- unique(classes)
	if (length(u_classes) <= 1)
		stop("More than one taxonomic group is required.")
	taxonomy <- strsplit(u_classes, ";", fixed=TRUE)
	taxonomy <- lapply(taxonomy,
		function(x) {
			sapply(seq_along(x),
				function(n)
					paste(x[seq_len(n)],
						";",
						sep="",
						collapse=""))
		})
	taxonomy <- unique(unlist(taxonomy))
	crossIndex <- match(classes, taxonomy) + 1L # offset by 1 to add "Root"
	taxonomy <- c("Root;", paste("Root;", taxonomy, sep=""))
	s <- strsplit(taxonomy, ";", fixed=TRUE)
	taxa <- unlist(lapply(s,
		function(x) x[length(x)]))
	levels <- lengths(s)
	
	maxL <- max(levels) - 1L # max level with children
	levs <- vector("list", maxL)
	for (i in seq_len(maxL))
		levs[[i]] <- which(levels==(i + 1L))
	starts <- rep(1L, maxL)
	children <- vector("list", length(taxonomy))
	for (i in seq_along(children)) {
		j <- levels[i]
		if (j <= maxL) {
			while (starts[j] <= length(levs[[j]]) &&
				levs[[j]][starts[j]] <= i)
				starts[j] <- starts[j] + 1L
			if (starts[j] <= length(levs[[j]])) {
				w <- levs[[j]][starts[j]:length(levs[[j]])]
				w <- w[startsWith(taxonomy[w], taxonomy[i])]
			} else {
				w <- integer()
			}
		} else {
			w <- integer()
		}
		children[[i]] <- w
	}
	
	parents <- numeric(length(taxonomy))
	for (i in seq_along(children))
		parents[children[[i]]] <- i
	# record the rank of each taxon
	if (is.null(rank)) {
		ranks <- NULL
	} else {
		ranks <- character(length(taxa))
		for (i in seq_along(ranks)) {
			w <- which(rank$Name==taxa[i])
			if (length(w)==0) {
				stop("rank is missing the name: ",
					taxa[i])
			} else if (length(w) > 1) {
				w <- w[which(rank$Level[w]==(levels[i] - 1L))]
				if (length(w)==0) {
					stop("rank is missing the name: ",
						taxa[i],
						" at level ",
						levels[i] - 1L)
				} else if (length(w) > 1) {
					stop("rank contains a duplicated name: ",
						taxa[i],
						" at level ",
						levels[i] - 1L)
				}
			} else {
				if (rank$Level[w] != (levels[i] - 1L))
					stop("rank matches a name at the wrong level: ",
						taxa[i])
			}
			ranks[i] <- rank$Rank[w]
			
			# confirm the correct parent
			w <- rank$Parent[w] # Index of Parent
			if (w > 0) { # not Root
				w <- which(rank$Index==w) # Index matching Parent
				if (length(w)==0)
					stop("All values of Parent in rank must have a corresponding Index.")
				if (rank$Name[w] != taxa[parents[i]])
					stop("rank contains an unexpected Parent for Name: ",
						taxa[i])
			}
		}
	}
	
	# find the sequences corresponding to each taxonomy
	end_taxonomy <- substring(taxonomy, 6) # without "Root;"
	sequences <- vector("list", length(end_taxonomy))
	nchars <- nchar(end_taxonomy)
	u_nchars <- unique(nchars)
	nchars <- match(nchars, u_nchars)
	nchars <- tapply(seq_along(nchars), nchars, c, simplify=FALSE)
	for (i in seq_along(u_nchars)) {
		m <- match(substring(classes, 0, u_nchars[i]),
			end_taxonomy[nchars[[i]]])
		m <- tapply(seq_along(m), m, c, simplify=FALSE)
		sequences[nchars[[i]]] <- m
	}
	nSeqs <- lengths(sequences)
	
	# record the most distinctive k-mers at each node in the taxonomic tree
	.createTree <- function(k) {
		l <- length(children[[k]])
		if (l > 0 && l <= maxChildren) {
			# obtain relative k-mer frequency by group
			profile <- vector("list", l)
			descendants <- integer(l)
			for (i in seq_len(l)) {
				results <- .createTree(children[[k]][i])
				profile[[i]] <- results[[1]]
				descendants[i] <- results[[2]]
			}
			
			profile <- matrix(unlist(profile),
				ncol=l)
			
			# calculate cross-entropy between parent and children
			q <- colSums(t(profile)*descendants) # weighted k-mer frequencies
			descendants <- sum(descendants) # total descendant groups
			q <- q/descendants # relative k-mer frequency overall
			H <- -profile*log(q) # entropy_i = -p_i*log(q_i)
			
			# select the most distinguishing k-mers for each group
			o <- apply(H, 2, order, decreasing=TRUE)
			keep <- logical(length(q))
			count <- 0L # number of k-mers
			i <- 1L # k-mer
			j <- 0L # group
			recordKmers <- ceiling(max(colSums(profile > 0)))*0.10
			while (count < recordKmers) {
				j <- j + 1L
				if (j > l) {
					j <- 1L
					i <- i + 1L
				}
				
				if (!keep[o[i, j]]) {
					count <- count + 1L
					keep[o[i, j]] <- TRUE
				}
			}
			keep <- which(keep)
			
			decision_kmers[[k]] <<- list(keep, t(profile[keep,]))
			invisible(list(q, descendants))
		} else {
			w <- sequences[[k]]
			
			profile <- tabulate(unlist(kmers[w]),
				nKmers)
			profile <- profile/sum(profile) # normalize
			invisible(list(profile, 1L))
		}
	}
	decision_kmers <- vector("list", length(taxonomy))
	.createTree(1) # traverse the taxonomic tree
	
	# learn the appropriate fraction for sampling
	fraction <- rep(maxFraction, length(taxonomy))
	incorrect <- rep(TRUE, l)
	predicted <- character(l)
	# allow each sequence maxIterations*delta control over fraction
	delta <- (maxFraction - minFraction)*multiplier
	for (it in seq_len(maxIterations)) {
		remainingSeqs <- which(incorrect)
		if (length(remainingSeqs)==0)
			break
		
		# test whether the sequences can be correctly classified
		for (i in remainingSeqs) {
			if (length(kmers[[i]])==0)
				next # no k-mers to test
			
			k <- 1L
			correct <- TRUE
			repeat {
				subtrees <- children[[k]]
				n <- length(decision_kmers[[k]][[1]])
				
				if (n==0) { # no decision k-mers
					break
				} else if (length(subtrees) > 1) {
					# set number of k-mers to choose each bootstrap replicate
					if (is.na(fraction[k])) {
						s <- ceiling(n*minFraction)
					} else {
						s <- ceiling(n*fraction[k])
					}
					
					# sample the decision k-mers
					sampling <- matrix(sample(n, s*B, replace=TRUE), B, s)
					
					hits <- matrix(0,
						nrow=length(subtrees),
						ncol=B)
					matches <- .Call("intMatch",
						decision_kmers[[k]][[1]],
						kmers[[i]],
						PACKAGE="DECIPHER")
					myweights <- decision_kmers[[k]][[2]]
					for (j in seq_along(subtrees)) {
						hits[j,] <- .Call("vectorSum",
							matches,
							myweights[j,],
							sampling,
							B,
							PACKAGE="DECIPHER")
					}
					
					maxes <- apply(hits, 2, max) # max hits per bootstrap replicate
					hits <- colSums(t(hits)==maxes & maxes > 0)
					w <- which.max(hits)
					if (hits[w] < B*0.8) # require 80% confidence to further descend the taxonomic tree
						break
				} else { # only one child
					w <- 1
				}
				
				if (!startsWith(classes[i], end_taxonomy[subtrees[w]])) {
					correct <- FALSE
					break
				}
				
				if (length(children[[subtrees[w]]])==0)
					break
				
				k <- subtrees[w]
			}
			
			if (correct) { # correct group
				incorrect[i] <- FALSE
				if (verbose) {
					num <- num + 1L
					setTxtProgressBar(pBar, num/l)
				}
			} else { # incorrect group
				predicted[i] <- taxonomy[subtrees[w]]
				if (is.na(fraction[k])) {
					incorrect[i] <- NA
					if (verbose) {
						num <- num + 1L
						setTxtProgressBar(pBar, num/l)
					}
				} else {
					fraction[k] <- fraction[k] - delta/nSeqs[k]
					if (fraction[k] <= minFraction) {
						incorrect[i] <- NA
						fraction[k] <- NA
						if (verbose) {
							num <- num + 1L
							setTxtProgressBar(pBar, num/l)
						}
					}
				}
			}
		}
	}
	
	# record problematic taxonomic groups and sequences
	w <- which(incorrect | is.na(incorrect))
	if (length(w) > 0) {
		problemSequences <- data.frame(Index=w,
			Expected=paste("Root;", classes[w], sep=""),
			Predicted=predicted[w],
			#row.names=names(train)[w],
			stringsAsFactors=FALSE)
	} else {
		problemSequences <- data.frame(Index=w,
			Expected=character(),
			Predicted=predicted[w],
			stringsAsFactors=FALSE)
	}
	
	w <- which(is.na(fraction))
	if (length(w) > 0) {
		problemGroups <- taxonomy[w]
	} else {
		problemGroups <- character()
	}
	
	# count the number of representatives per class
	classTable <- table(classes)
	weight <- 1/classTable[classes] # weigh each class equally
	# compute Inverse Document Frequency (IDF) weights
	counts <- numeric(nKmers)
	for (i in seq_along(classes)) {
		counts[kmers[[i]]] <- counts[kmers[[i]]] + weight[i]
	}
	counts <- length(classTable)/(1 + counts) # Max/(1 + F)
	counts <- log(counts) # IDF weight is log(Max/(1 + F))
	
	# prepare the result to return
	result <- list(taxonomy=taxonomy,
		taxa=taxa,
		ranks=ranks,
		levels=levels,
		children=children,
		parents=parents,
		fraction=fraction,
		sequences=sequences,
		kmers=kmers,
		crossIndex=crossIndex,
		K=K,
		IDFweights=counts,
		decisionKmers=decision_kmers,
		problemSequences=problemSequences,
		problemGroups=problemGroups,
		alphabet=alphabet)
	class(result) <- c("Taxa", "Train")
	
	if (verbose) {
		setTxtProgressBar(pBar, 1)
		close(pBar)
		cat("\n")
		time.2 <- Sys.time()
		print(round(difftime(time.2,
				time.1,
				units='secs'),
			digits=2))
		cat("\n")
	}
	
	invisible(result)
}
