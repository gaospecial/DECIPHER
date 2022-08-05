IdClusters <- function(myXStringSet,
	cutoff=0,
	maxReps=10000,
	maxComparisons=1000,
	maxAlignments=100,
	invertCenters=FALSE,
	alphabet=AA_REDUCED[[30]],
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.numeric(cutoff))
		stop("cutoff must be a numeric.")
	if (is.integer(cutoff))
		cutoff <- as.numeric(cutoff)
	if (any(cutoff < 0))
		stop("cutoff must be at least zero.")
	if (any(cutoff >= 1))
		stop("cutoff must be less than one.")
	if (any(duplicated(cutoff)))
		stop("cutoff cannot contain duplicated values.")
	ASC <- TRUE
	if (is.unsorted(cutoff)) {
		if (is.unsorted(rev(cutoff))) {
			stop("cutoff must be sorted.")
		} else {
			ASC <- FALSE
		}
	}
	if (!is.numeric(maxReps))
		stop("maxReps must be a numeric.")
	if (length(maxReps) != 1)
		stop("maxReps must only be a single number.")
	if (maxReps < 1)
		stop("maxReps must be at least 1.")
	if (floor(maxReps) != maxReps)
		stop("maxReps must be a whole number.")
	if (!is.numeric(maxComparisons))
		stop("maxComparisons must be a numeric.")
	if (length(maxComparisons) != 1)
		stop("maxComparisons must only be a single number.")
	if (maxComparisons < 2)
		stop("maxComparisons must be at least 2.")
	if (floor(maxComparisons) != maxComparisons)
		stop("maxComparisons must be a whole number.")
	if (!is.numeric(maxAlignments))
		stop("maxAlignments must be a numeric.")
	if (length(maxAlignments) != 1)
		stop("maxAlignments must only be a single number.")
	if (maxAlignments > maxComparisons)
		stop("maxAlignments must be less than or equal to maxComparisons.")
	if (floor(maxAlignments) != maxAlignments)
		stop("maxAlignments must be a whole number.")
	if (!is.logical(invertCenters))
		stop("invertCenters must be a logical.")
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
	
	if (is(myXStringSet, "DNAStringSet")) {
		typeX <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		typeX <- 2L
	} else if (is(myXStringSet, "AAStringSet")) {
		typeX <- 3L
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
		sizeAA <- max(alphabet)
		if (sizeAA == 1L)
			stop("More than one grouping of amino acids is required in the alphabet.")
		maxK <- as.integer(floor(log(4294967295, sizeAA)))
		alphabet <- alphabet - 1L
	} else {
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') are not allowed in myXStringSet.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') are not allowed in myXStringSet.")
	widths <- width(myXStringSet)
	if (all(widths == 0L))
		stop("All sequences in myXStringSet are zero width.")
	l <- length(myXStringSet)
	if (l==0)
		stop("myXStringSet contains no sequences.")
	
	# initialize parameters
	alpha <- 0.1 # weight of exponential moving average
	attempts <- 10L # comparison attempts before possibly skipping
	minAttempts <- min(10L, maxComparisons) # minimum number of comparison attempts
	buffer <- 0.01 # minimum (approximate) distance buffer
	
	if (typeX==3L) { # AAStringSet
		wordSize <- ceiling(log(100*quantile(widths, 0.99),
			.Call("alphabetSizeReducedAA",
				myXStringSet,
				alphabet,
				PACKAGE="DECIPHER")))
		if (wordSize > maxK)
			wordSize <- maxK
		if (wordSize < 1)
			wordSize <- 1
		words <- 20^wordSize
	} else { # DNAStringSet or RNAStringSet
		wordSize <- ceiling(log(100*quantile(widths, 0.99),
			.Call("alphabetSize",
				myXStringSet,
				PACKAGE="DECIPHER")))
		if (wordSize > 15)
			wordSize <- 15
		if (wordSize < 1)
			wordSize <- 1
		words <- 4^wordSize
	}
	
	if (verbose) {
		time.1 <- Sys.time()
		cat("Ordering sequences by ",
			wordSize,
			"-mer similarity:\n",
			sep="")
		flush.console()
	}
	
	lc <- length(cutoff)
	if (is.null(names(myXStringSet))) {
		dNames <- 1:l
	} else {
		dNames <- names(myXStringSet)
		w <- which(duplicated(dNames))
		if (length(w) > 0) {
			warning("Duplicated names of myXStringSet appended with index.")
			dNames[w] <- paste(dNames[w],
				w,
				sep="_")
		}
	}
	if (lc > 1) {
		cNames <- paste("cluster",
			gsub(".", "_", prettyNum(cutoff), fixed=TRUE),
			sep="_")
	} else {
		cNames <- "cluster"
	}
	c <- matrix(0L,
		nrow=l,
		ncol=lc,
		dimnames=list(dNames,
			cNames))
	
	# identify duplicated sequences
	x <- selfmatch(myXStringSet)
	u <- unique(x)
	l <- length(u)
	t <- tabulate(x, length(x))[u]
	
	overlap <- function(x, y) {
		.Call("matchOverlap",
			x,
			y,
			v,
			z,
			wordSize,
			PACKAGE="DECIPHER")
	}
	
	similarity <- function(x, w1, w2) {
		n <- ncol(x)
		if (n > 0) {
			s <- sum(x[2L,] - x[1L,] + 1)
			p1 <- x[1L, 1L]
			p2 <- x[3L, 1L]
			if (p1 <= p2 && w1 <= w2) { # 1 within 2
				ov <- w1
			} else if (p2 <= p1 && w2 <= w1) { # 2 within 1
				ov <- w2
			} else if (p1 > p2) { # end of 1 overlaps start of 2
				ov <- w1 - p1 + p2
			} else { # start of 2 overlaps end of 1
				ov <- w2 - p2 + p1
			}
			s/ov
		} else {
			0
		}
	}
	
	if (typeX==3L) { # AAStringSet
		V <- v <- .Call("enumerateSequenceReducedAA",
			.subset(myXStringSet, u),
			wordSize,
			alphabet,
			TRUE, # mask repeats
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		V <- v <- .Call("enumerateSequence",
			.subset(myXStringSet, u),
			wordSize,
			TRUE, # mask repeats
			PACKAGE="DECIPHER")
	}
	Z <- z <- lapply(V, order)
	wu <- sapply(seq_along(V),
		function(x) {
			res <- overlap(x, x)[[1L]]
			sum(res[2L,] - res[1L,] + 1)
		})
	
	# order sequences by approximate similarity
	var <- rep(length(v), length(v))
	avg_cor <- 0
	for (i in seq_len(maxReps)) {
		r <- sample(length(v), 2L, prob=var)
		res1 <- overlap(r[1L], seq_along(v))
		res2 <- overlap(r[2L], seq_along(v))
		m1 <- mapply(similarity, res1, wu[r[1L]], wu)
		m2 <- mapply(similarity, res2, wu[r[2L]], wu)
		d <- m1 - m2
		
		if (i > 1) {
			if (cor(rS, d) < 0) # negatively correlated
				d <- -d
			rS <- rS + d
			
			o <- order(rS,
				wu,
				t, # frequency
				decreasing=TRUE,
				method="radix")
			r <- order(o,
				method="radix")
			
			var <- alpha*abs(R - r) + (1 - alpha)*var # exponential moving average of change in rank order
			avg_cor <- alpha*sum(var < maxComparisons/2)/length(var) + (1 - alpha)*avg_cor # exponential moving average of fraction stabilized in rank order
			
			O <- o
			R <- r
		} else {
			rS <- d
			
			O <- order(rS,
				wu,
				t, # frequency
				decreasing=TRUE,
				method="radix")
			R <- order(O,
				method="radix")
		}
		
		if (verbose && interactive()) {
			cat("\riteration = ",
				i,
				" of up to ",
				maxReps,
				" (stability = ",
				formatC(100*avg_cor, digits=1, format="f"),
				"%) ",
				sep="")
		}
		if (avg_cor >= 0.9995) # rounds to 100%
			break # reached rank order stability
	}
	rS <- rS/i # normalize per iteration
	
	if (verbose) {
		time.2 <- Sys.time()
		cat("\n\n")
		print(round(difftime(time.2,
			time.1,
			units='secs'),
			digits=2))
		cat("\n")
		
		time.1 <- time.2
		lastValue <- 0
		cat("Clustering sequences by similarity:\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	}
	
	# attempt to find cluster centers
	S <- smooth.spline(rep(rS[O], t[O]),
		w=rep(wu[O], t[O]),
		spar=0.2)
	S <- predict(S,
		x=cumsum(t[O]) - (t[O] - 1L)/2,
		deriv=1L)$y
	S <- abs(S)
	P <- order(S,
		method="radix")
	P <- O[P]
	Q <- R[P]
	
	buffer <- buffer + 2*cutoff
	cutoff <- 1 - cutoff
	for (i in seq_len(lc)) {
		if (!ASC && i > 1) {
			P <- order(c[u, i - 1][O],
				S,
				method="radix")
			P <- O[P]
			Q <- R[P]
		}
		
		bL <- bR <- integer(l)
		k <- 1L
		for (j in seq_len(l)) {
			while (j - k >= maxComparisons ||
				(j - k > minAttempts &&
				rS[O[k]] > rS[O[j]] + buffer[i]))
				k <- k + 1L
			bL[j] <- k
		}
		k <- l
		for (j in l:1) {
			while (k - j >= maxComparisons ||
				(k - j > minAttempts &&
				rS[O[j]] > rS[O[k]] + buffer[i]))
				k <- k - 1L
			bR[j] <- k
		}
		
		v <- V[P] # reorder k-mers
		z <- Z[P] # reorder k-mer ordering
		p <- u[P] # index of sequence
		
		cluster_num <- 1L
		offset <- 0L
		if (invertCenters) {
			c[p[1L], i] <- -1L
		} else {
			c[p[1L], i] <- 1L
		}
		seeds.index <- integer(l)
		seeds.index[Q[1L]] <- 1L
		
		for (j in seq_along(P)[-1L]) {
			if (verbose) {
				value <- round(((i - 1)*l + j)/lc/l, 2)
				if (value > lastValue) {
					lastValue <- value
					setTxtProgressBar(pBar, value)
				}
			}
			
			if (!ASC && i > 1L) {
				if (c[p[j], i - 1L] != c[p[j - 1L], i - 1L]) {
					# different clusters in last cutoff
					offset <- offset + cluster_num
					cluster_num <- 1L
					c[p[j], i] <- cluster_num + offset
					seeds.index[Q[j]] <- j
					next
				} else if (c[p[j], i] > 0L) { # cluster pre-assigned
					if (invertCenters) {
						c[p[j], i] <- abs(c[c[p[j], i], i])
					} else {
						c[p[j], i] <- c[c[p[j], i], i]
					}
					next
				}
			}
			
			compare <- bL[Q[j]]:bR[Q[j]]
			compare <- as.integer(compare)
			compare <- compare[compare >= 1L & compare <= l]
			if (!ASC && i > 1) { # check bounds
				if (invertCenters) {
					compare <- compare[abs(c[p[seeds.index[compare]], i - 1L]) == abs(c[p[j], i - 1L])]
				} else {
					compare <- compare[c[p[seeds.index[compare]], i - 1L] == c[p[j], i - 1L]]
				}
			}
			compare <- compare[seeds.index[compare] != 0L]
			if (length(compare) == 0L) {
				cluster_num <- cluster_num + 1L
				c[p[j], i] <- cluster_num + offset
				if (invertCenters)
					c[p[j], i] <- -c[p[j], i]
				seeds.index[Q[j]] <- j
				next
			}
			
			compare <- compare[!duplicated(seeds.index[compare])]
			res <- overlap(j, seeds.index[compare])
			m <- mapply(similarity, res, wu[P[j]], wu[P[seeds.index[compare]]])
			w <- head(order(m, decreasing=TRUE), maxAlignments)
			
			for (k in seq_along(w)) {
				if (k > attempts) {
					s <- seq(to=k - 1L, length.out=attempts)
					if (sum(m[w[s]])/attempts + max(0.1, 3*sd(m[w[s]])) < cutoff[i])
						break
				}
				ali <- AlignProfiles(.subset(myXStringSet, p[j]),
					.subset(myXStringSet, p[seeds.index[compare[w[k]]]]),
					anchor=res[[w[k]]])
				m[w[k]] <- .Call("distMatrix",
					ali,
					typeX,
					FALSE, # includeTerminalGaps
					FALSE, # penalizeGapGapMatches
					TRUE, # penalizeGapLetterMatches
					TRUE, # full matrix
					2L, # type = "dist"
					0, # correction
					FALSE, # verbose
					NULL, # progress bar
					1L, # processors
					PACKAGE="DECIPHER")
				if (is.na(m[w[k]])) {
					m[w[k]] <- 0
				} else {
					m[w[k]] <- 1 - m[w[k]]
				}
				if (m[w[k]] >= cutoff[i])
					break
			}
			w <- w[which.max(m[w])]
			
			if (length(w) == 0L ||
				m[w] < cutoff[i]) { # form a new group
				cluster_num <- cluster_num + 1L
				c[p[j], i] <- cluster_num + offset
				if (invertCenters)
					c[p[j], i] <- -c[p[j], i]
				seeds.index[Q[j]] <- j
			} else { # part of an existing group
				if (invertCenters) {
					c[p[j], i] <- -c[p[seeds.index[compare[w]]], i]
				} else {
					c[p[j], i] <- c[p[seeds.index[compare[w]]], i]
				}
				if (!ASC && i < lc) {
					cols <- (i + 1):lc
					cols <- cols[which(m[w] >= cutoff[cols])]
					if (length(cols) > 0L) # assign forward
						c[p[j], cols] <- p[seeds.index[compare[w]]]
				}
			}
		}
		
		c[, i] <- c[x, i]
	}
	myClusters <- as.data.frame(c)
	
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
	
	return(myClusters)
}
