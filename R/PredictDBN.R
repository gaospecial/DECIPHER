PredictDBN <- function(myXStringSet,
	type="states",
	minOccupancy=0.5,
	impact=c(1, 1.2, 0.4, -1),
	avgProdCorr=1,
	slope=2,
	shift=1.3,
	threshold=0.7,
	pseudoknots=1,
	weight=NA,
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
	if (length(weight)==1) {
		if (is.numeric(weight)) {
			weight <- rep(1, length(myXStringSet))
		} else if (is.na(weight)) {
			if (verbose) {
				cat("Determining distance matrix based on alignment:\n")
				flush.console()
			}
			d <- DistanceMatrix(myXStringSet,
				type="dist",
				verbose=verbose,
				processors=processors)
			
			if (verbose) {
				cat("Determining sequence weights:\n")
				flush.console()
			}
			suppressWarnings(guideTree <- IdClusters(d,
				method="UPGMA",
				type="dendrogram",
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
		if (length(weight)!=length(myXStringSet))
			stop("Length of weight must equal one or the length of the myXStringSet.")
		if (any(is.na(weight)))
			stop("weight must be a numeric.")
		if (!isTRUE(all.equal(1, mean(weight))))
			stop("The mean of weight must be 1.")
	}
	
	# initialize a progress bar
	if (verbose) {
		cat("Computing RNA Secondary Structures:\n")
		pBar <- txtProgressBar(min=0, max=100, initial=0, style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
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
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	if (type==2L) {
		s <- strsplit(ans, "")[[1]]
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
				}
			}
		}
		
		o <- order(c3, c1, c2)
		ans <- matrix(c(c1[o], c2[o], c3[o]),
			ncol=3,
			dimnames=list(NULL, c("(", ")", "order")))
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
