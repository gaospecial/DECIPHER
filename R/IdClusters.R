IdClusters <- function(myXStringSet,
	cutoff=0,
	maxReps=1000,
	maxAttempts=100,
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
	if (!is.numeric(maxAttempts))
		stop("maxAttempts must be a numeric.")
	if (length(maxAttempts) != 1)
		stop("maxAttempts must only be a single number.")
	if (maxAttempts < 1)
		stop("maxAttempts must be at least 1.")
	if (floor(maxAttempts) != maxAttempts)
		stop("maxAttempts must be a whole number.")
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
	denom <- 3 # relative factor of k-mer length differences to allow
	
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
			gsub("^.+\\.", "", cutoff),
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
	wu <- widths[u]
	
	f <- function(x, y) {
		X <- v[[x]]
		OX <- z[[x]]
		size_x <- length(X)
		ly <- length(y)

		res <- vector("list", ly)
		for (i in seq_len(ly)) {
			Y <- v[[y[i]]]
			OY <- z[[y[i]]]
			size_y <- length(Y)
			
			m <- integer(size_x)
			j <- k <- 1L
			while (j <= size_x && k <= size_y) {
				if (is.na(X[OX[j]]) || is.na(Y[OY[k]])) {
					break
				} else if (X[OX[j]] == Y[OY[k]]) {
					m[OX[j]] <- OY[k]
					j <- j + 1L
					k <- k + 1L
				} else if (X[OX[j]] < Y[OY[k]]) {
					j <- j + 1L
				} else {
					k <- k + 1L
				}
			}
			
			# count the number of anchors
			n <- 0L
			new <- TRUE
			for (j in seq_len(size_x)) {
				if (m[j] > 0) {
					if (new) {
						new <- FALSE
						n <- n + 1L
					} else if (m[j] != m[j - 1L] + 1L) {
						n <- n + 1L
					}
				} else {
					new <- TRUE
				}
			}
			
			# record anchor ranges
			t <- matrix(0, nrow=3, ncol=n)
			k <- 0L
			new <- TRUE
			for (j in seq_len(size_x)) {
				if (m[j] > 0) {
					if (new ||
						m[j] != m[j - 1L] + 1L) {
						new <- FALSE
						k <- k + 1L
						t[1L, k] <- j
						t[2L, k] <- m[j]
						t[3L, k] <- wordSize
					} else {
						t[3L, k] <- t[3L, k] + 1L
					}
				} else {
					if (k == n)
						break
					new <- TRUE
				}
			}
			
			# chain anchors
			b <- numeric(n) # traceback
			s <- t[3L,] # cumulative score
			j <- 2L
			p <- 1L
			while (j <= n) {
				for (k in seq_len(j - 1)) {
					if (t[2L, k] < t[2L, j] && # starts within bounds
						((t[2L, k] + t[3L, k] <= t[2L, j] &&
						t[1L, k] + t[3L, k] <= t[1L, j]) ||
						t[1L, j] - t[1L, k] >= t[2L, j] - t[2L, k])) { # ends are compatible (no overlap is ==)
						if (s[k] + t[3L, j] > s[j]) { # extend chain
							s[j] <- s[k] + t[3L, j]
							b[j] <- k
						}
					}
				}
				if (s[j] > s[p]) # higher score
					p <- j
				j <- j + 1L
			}
			
			# rectify anchors
			keep <- logical(n)
			if (n > 0) {
				while (p > 0) { # traceback
					keep[p] <- TRUE
					p <- b[p]
				}
			}
			k <- 0L
			j <- 1L
			while (j <= n) {
				if (keep[j]) {
					if (k > 0L &&
						t[1L, k] + t[3L, k] >= t[1L, j] &&
						t[1L, j] - t[1L, k] == t[2L, j] - t[2L, k]) {
						# merge anchors
						keep[j] <- FALSE
						t[3L, k] <- t[3L, j] + t[1L, j] - t[1L, k]
					} else {
						k <- j
					}
				}
				j <- j + 1L
			}
			
			# record final anchors
			k <- 0L
			for (j in seq_len(n))
				if (keep[j])
					k <- k + 1L
			res[[i]] <- matrix(0, nrow=4, ncol=k)
			k <- 0L
			for (j in seq_len(n)) {
				if (keep[j]) {
					k <- k + 1L
					res[[i]][1L, k] <- t[1L, j]
					res[[i]][2L, k] <- t[1L, j] + t[3L, j] - 1L
					res[[i]][3L, k] <- t[2L, j]
					res[[i]][4L, k] <- t[2L, j] + t[3L, j] - 1L
					# remove overlap
					if (k > 1L) {
						delta <- res[[i]][2L, k - 1L] - res[[i]][1L, k]
						if (delta >= 0) {
							res[[i]][1L, k] <- res[[i]][1L, k] + delta + 1L
							res[[i]][3L, k] <- res[[i]][3L, k] + delta + 1L
						}
						delta <- res[[i]][4L, k - 1L] - res[[i]][3L, k]
						if (delta >= 0) {
							res[[i]][1L, k] <- res[[i]][1L, k] + delta + 1L
							res[[i]][3L, k] <- res[[i]][3L, k] + delta + 1L
						}
					}
				}
			}
		}
		
		res
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
	Z <- z <- lapply(v, order)
	
	# order sequences by approximate similarity
	var <- rep(1, length(v))
	avg_cor <- 0
	for (i in seq_len(maxReps)) {
		repeat {
			r <- sample(length(v), 4L, prob=var)
			res <- f(r[1], r[2:4])
			m <- sapply(res,
				function(x)
					sum(x[2,] - x[1,]))
			w <- which.max(m)
			i1 <- intersect(v[[r[1L]]], v[[r[w + 1L]]])
			w <- (2:4)[-w]
			i2 <- intersect(v[[r[w[1L]]]], v[[r[w[2L]]]])
			k1 <- i1[!(i1 %in% i2)]
			k2 <- i2[!(i2 %in% i1)]
			if (length(k1) > length(i1)/denom &&
				length(k2) > length(i2)/denom &&
				length(k1) > length(k2)/denom &&
				length(k2) > length(k1)/denom)
				break
		}
		
		d <- sapply(v,
			function(x)
				sum(k1 %in% x)/length(k1) - sum(k2 %in% x)/length(k2))
		
		if (i > 1) {
			if (cor(rS, d) < 0) # negatively correlated
				d <- -d
			rS <- rS + d
			
			var <- alpha*abs(d/rS) + (1 - alpha)*var
			w <- which(is.na(var))
			var[w] <- 1
			w <- which(var > 1)
			var[w] <- 1
			avg_cor <- alpha*(1 - sum(var)/length(var)) + (1 - alpha)*avg_cor
		} else {
			rS <- d
		}
		
		if (verbose && interactive()) {
			cat("\riteration = ",
				i,
				" of up to ",
				maxReps,
				" (stability = ",
				round(100*avg_cor),
				"%) ",
				sep="")
		}
		if (avg_cor >= 0.995) # rounds to 100%
			break # reached stasis
	}
	
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
	
	O <- order(rS,
		wu,
		t, # frequency
		decreasing=TRUE,
		method="radix")
	
	cutoff <- 1 - cutoff
	for (i in seq_len(lc)) {
		if (!ASC && i > 1) {
			O <- order(c[u, i - 1],
				rS,
				wu,
				t, # frequency
				decreasing=TRUE,
				method="radix")
		}
		
		v <- V[O]
		z <- Z[O]
		o <- u[O]
		
		cluster_num <- 1L
		offset <- 0L
		c[o[1L], i] <- 1L
		seeds.index <- integer(l)
		seeds.index[1] <- 1L
		
		for (j in seq_along(o)[-1]) {
			if (verbose) {
				value <- round(((i - 1)*l + j)/lc/l, 2)
				if (value > lastValue) {
					lastValue <- value
					setTxtProgressBar(pBar, value)
				}
			}
			
			if (!ASC && i > 1) {
				if (c[o[j], i - 1] != c[o[j - 1], i - 1]) {
					# different clusters in last cutoff
					offset <- offset + cluster_num
					cluster_num <- 1L
					c[o[j], i] <- cluster_num + offset
					seeds.index[1] <- j
					next
				} else if (c[o[j], i] > 0) { # cluster pre-assigned
					c[o[j], i] <- c[c[o[j], i], i]
					next
				}
			}
			
			compare <- seq(to=cluster_num, length.out=min(cluster_num, maxAttempts))
			res <- f(j, seeds.index[compare])
#			m <- .Call("matchListsDual",
#				v[j],
#				v[seeds.index[seq_len(cluster_num)]],
#				FALSE, # verbose
#				NULL, # pBar
#				processors,
#				PACKAGE="DECIPHER")
			m <- sapply(res,
				function(x) {
					n <- ncol(x)
					if (n > 0) {
						diff_s <- x[3L, 1L] - x[1L, 1L] # start difference
						diff_e <- widths[o[j]] - x[4L, n] + x[2L, n] # end difference
						sum(x[2,] - x[1,])/(diff_s + diff_e) # fraction shared
					} else {
						0
					}
				})
			w <- order(m, decreasing=TRUE)
			
			if (m[w[1L]] < cutoff[i]) {
				for (k in seq_along(w)) {
					if (k > attempts) {
						s <- seq(to=k - 1L, length.out=attempts)
						if (sum(m[w[s]])/attempts + 3*sd(m[w[s]]) < cutoff[i])
							break
					}
					ali <- AlignProfiles(.subset(myXStringSet, o[j]),
						.subset(myXStringSet, o[seeds.index[compare[w[k]]]]),
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
			}
			w <- which.max(m)
			
			if (length(w)==0 ||
				m[w] < cutoff[i]) { # form a new group
				cluster_num <- cluster_num + 1L
				c[o[j], i] <- cluster_num
				seeds.index[cluster_num] <- j
			} else { # part of an existing group
				c[o[j], i] <- compare[w] + offset
				if (!ASC && i < lc) {
					cols <- (i + 1):lc
					cols <- cols[which(m[w] >= cutoff[cols])]
					if (length(cols) > 0) # assign forward
						c[o[j], cols] <- o[seeds.index[compare[w]]]
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
