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

PredictDBN <- function(myXStringSet,
	type="states",
	minOccupancy=0.5,
	impact=c(1, 1.2, 0.4, -1),
	avgProdCorr=1,
	slope=2,
	shift=1.3,
	threshold=0.3,
	pseudoknots=1,
	weight=NA,
	useFreeEnergy=TRUE,
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
		
		ions <- 1
		temp <- 37
		data("deltaHrulesRNA", envir=environment(), package="DECIPHER")
		data("deltaSrulesRNA", envir=environment(), package="DECIPHER")
		deltaSrulesRNA <- deltaSrulesRNA + 0.368*log(ions)/1000
		deltaGrulesRNA <- deltaHrulesRNA - (273.15 + temp)*deltaSrulesRNA
		max_dG <- 4.065225 # 3.6 - (273.15 + temp)*(-1.5/1000 + 0.368*log(ions)/1000)
		
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
			arms <- palindromeArmLength(dna,
				max.mismatch=1,
				allow.wobble=TRUE)
			max_arms <- as.integer((width(pals) - 3)/2)
			arms <- ifelse(arms > max_arms,
				max_arms,
				arms)
			
			s1 <- start(pals)
			s2 <- start(pals) + arms - 1L
			e1 <- end(pals) - arms + 1L
			e2 <- end(pals)
			
			dG <- .Call("calculateHairpinDeltaG",
				dna,
				arms + 1L, # include end bases
				deltaGrulesRNA,
				PACKAGE="DECIPHER")
			o <- order(dG)
			o <- o[dG[o] <= max_dG]
			
			if (length(o) > 0) {
				p <- pos[[i]]
				row <- col <- integer(100)
				count <- 0L
				r <- c <- logical(length(p))
				for (j in o) {
					cont <- .Call("allZero",
						r,
						c,
						s1[j],
						s2[j],
						e2[j],
						e1[j],
						PACKAGE="DECIPHER")
					if (cont) {
						p1 <- s1[j]:s2[j]
						p2 <- e2[j]:e1[j]
						r[p1] <- TRUE
						c[p2] <- TRUE
						p1 <- p[p1]
						p2 <- p[p2]
						index <- (count + 1L):(count + length(p1))
						delta <- index[length(index)] - length(row)
						if (delta > 0) {
							row <- c(row, integer(max(100, delta)))
							col <- c(col, integer(max(100, delta)))
						}
						row[index] <- p1
						col[index] <- p2
						count <- count + length(p1)
					}
				}
				if (count > 0L) {
					for (j in seq_len(count)) {
						m <- match(col[j], J[[row[j]]])
						if (is.na(m)) {
							J[[row[j]]] <- c(J[[row[j]]],
								col[j])
							W[[row[j]]] <- c(W[[row[j]]],
								weight[i])
						} else {
							W[[row[j]]][m] <- W[[row[j]]][m] + weight[i]
						}
					}
				}
			}
			if (verbose)
				setTxtProgressBar(pBar, i/l)
		}
		W <- lapply(W, `/`, l)
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
