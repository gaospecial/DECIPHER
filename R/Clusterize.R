Clusterize <- function(myXStringSet,
	cutoff=0,
	method="overlap",
	includeTerminalGaps=FALSE,
	penalizeGapLetterMatches=NA,
	minCoverage=0.5,
	maxReps=1000,
	avgComparisons=5000,
	maxAlignments=100,
	invertCenters=FALSE,
	processors=1,
	verbose=TRUE) {
	
	# error checking
	if (!is.numeric(cutoff))
		stop("cutoff must be a numeric.")
	if (is.integer(cutoff))
		cutoff <- as.numeric(cutoff)
	if (any(is.na(cutoff)))
		stop("cutoff must not contain NA values.")
	if (any(cutoff < 0))
		stop("cutoff must be at least zero.")
	if (any(cutoff >= 1))
		stop("cutoff must be less than one.")
	if (any(duplicated(cutoff)))
		stop("cutoff cannot contain duplicated values.")
	METHODS <- c("overlap", "shortest", "longest")
	method <- pmatch(method[1], METHODS)
	if (is.na(method))
		stop("Invalid method.")
	if (method==-1)
		stop("Ambiguous method.")
	if (!is.logical(includeTerminalGaps))
		stop("includeTerminalGaps must be a logical.")
	if (!is.logical(penalizeGapLetterMatches))
		stop("penalizeGapLetterMatches must be a logical.")
	ASC <- TRUE
	if (is.unsorted(cutoff)) {
		if (is.unsorted(rev(cutoff))) {
			stop("cutoff must be sorted.")
		} else {
			ASC <- FALSE
		}
	}
	if (!is.numeric(minCoverage))
		stop("maxCoverage must be a numeric.")
	if (minCoverage < 0)
		stop("minCoverage must be at least zero.")
	if (minCoverage > 1)
		stop("minCoverage can be at most one.")
	if (!is.numeric(maxReps))
		stop("maxReps must be a numeric.")
	if (length(maxReps) != 1)
		stop("maxReps must only be a single number.")
	if (maxReps < 1)
		stop("maxReps must be at least 1.")
	if (floor(maxReps) != maxReps)
		stop("maxReps must be a whole number.")
	if (!is.numeric(avgComparisons))
		stop("avgComparisons must be a numeric.")
	if (length(avgComparisons) != 1)
		stop("avgComparisons must only be a single number.")
	if (avgComparisons < 2)
		stop("avgComparisons must be at least 2.")
	if (floor(avgComparisons) != avgComparisons)
		stop("avgComparisons must be a whole number.")
	if (!is.numeric(maxAlignments))
		stop("maxAlignments must be a numeric.")
	if (length(maxAlignments) != 1)
		stop("maxAlignments must only be a single number.")
	if (maxAlignments > avgComparisons)
		stop("maxAlignments must be less than or equal to avgComparisons.")
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
	alpha <- 0.01 # weight of exponential moving average
	attempts <- 10L # comparison attempts before possibly skipping
	minAttempts <- min(attempts, avgComparisons) # minimum number of comparison attempts
	interval <- 0.999 # prediction interval for k-mer similarity limit (> 0 and < 1)
	minSimilarities <- 1000L # minimum number of similarities to predict limit
	buffer <- 0.01 # offset applied to k-mer simility limit at every cutoff
	minAlignments <- 1L # minimum number of alignments (> 0)
	minSize <- max(10*avgComparisons, 1e3L) # minimum size to continue splitting bins
	maxVariance <- ifelse(typeX == 3L, 0.5, 0.2) # maximum tolerable variance to stop splitting bins
	batchSize <- 50L # number of k-mer similarities to compute per batch (ideally >> processors)
	maxComparisons <- 10*avgComparisons
	totalAlignments <- 10/(1 - interval) # sufficient alignments for predicting similarity
	if (verbose)
		time.1 <- Sys.time()
	
	if (typeX == 3L) { # AAStringSet
		# use rounded PFASUM
		sM <- matrix(c(4L, -1L, -1L, -1L, 0L, -1L, -1L, 0L, -2L, -1L, -1L, -1L, -1L, -2L, -1L, 1L, 0L, -3L, -2L, 0L, -1L, 6L, 0L, -1L, -3L, 2L, 1L, -2L, 1L, -4L, -3L, 3L, -2L, -4L, -1L, 0L, -1L, -2L, -2L, -3L, -1L, 0L, 6L, 2L, -3L, 1L, 1L, 0L, 1L, -4L, -4L, 1L, -3L, -4L, -1L, 1L, 0L, -4L, -2L, -4L, -1L, -1L, 2L, 7L, -4L, 1L, 3L, -1L, 0L, -5L, -5L, 0L, -4L, -5L, -1L, 0L, -1L, -5L, -3L, -5L, 0L, -3L, -3L, -4L, 14L, -3L, -4L, -2L, -2L, -1L, -1L, -4L, -1L, -1L, -4L, 0L, -1L, -2L, -1L, 0L, -1L, 2L, 1L, 1L, -3L, 6L, 2L, -2L, 1L, -3L, -3L, 2L, -1L, -4L, -1L, 0L, 0L, -3L, -2L, -3L, -1L, 1L, 1L, 3L, -4L, 2L, 6L, -2L, 0L, -4L, -4L, 1L, -3L, -5L, -1L, 0L, 0L, -4L, -3L, -3L, 0L, -2L, 0L, -1L, -2L, -2L, -2L, 8L, -2L, -5L, -4L, -2L, -3L, -4L, -1L, 0L, -1L, -4L, -4L, -4L, -2L, 1L, 1L, 0L, -2L, 1L, 0L, -2L, 10L, -3L, -3L, 0L, -2L, -1L, -2L, 0L, -1L, -1L, 2L, -3L, -1L, -4L, -4L, -5L, -1L, -3L, -4L, -5L, -3L, 5L, 3L, -4L, 2L, 1L, -3L, -3L, -1L, -2L, -1L, 3L, -1L, -3L, -4L, -5L, -1L, -3L, -4L, -4L, -3L, 3L, 5L, -3L, 3L, 2L, -3L, -3L, -2L, -1L, -1L, 1L, -1L, 3L, 1L, 0L, -4L, 2L, 1L, -2L, 0L, -4L, -3L, 6L, -2L, -4L, -1L, 0L, 0L, -4L, -2L, -3L, -1L, -2L, -3L, -4L, -1L, -1L, -3L, -3L, -2L, 2L, 3L, -2L, 7L, 1L, -3L, -2L, -1L, -1L, 0L, 1L, -2L, -4L, -4L, -5L, -1L, -4L, -5L, -4L, -1L, 1L, 2L, -4L, 1L, 7L, -4L, -3L, -2L, 3L, 4L, 0L, -1L, -1L, -1L, -1L, -4L, -1L, -1L, -1L, -2L, -3L, -3L, -1L, -3L, -4L, 9L, 0L, -1L, -3L, -3L, -3L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, -3L, -3L, 0L, -2L, -3L, 0L, 4L, 2L, -3L, -2L, -2L, 0L, -1L, 0L, -1L, -1L, 0L, 0L, -1L, -1L, -1L, -2L, 0L, -1L, -2L, -1L, 2L, 5L, -3L, -2L, 0L, -3L, -2L, -4L, -5L, -2L, -3L, -4L, -4L, -1L, -2L, -1L, -4L, -1L, 3L, -3L, -3L, -3L, 14L, 3L, -2L, -2L, -2L, -2L, -3L, -1L, -2L, -3L, -4L, 2L, -1L, -1L, -2L, 0L, 4L, -3L, -2L, -2L, 3L, 9L, -1L, 0L, -3L, -4L, -5L, 0L, -3L, -3L, -4L, -3L, 3L, 1L, -3L, 1L, 0L, -3L, -2L, 0L, -2L, -1L, 5L),
			nrow=20,
			ncol=20)
	} else { # DNAStringSet or RNAStringSet
		# penalize transitions less than transversions
		sM <- matrix(c(3L, -3L, 0L, -3L, -3L, 3L, -3L, 0L, 0L, -3L, 3L, -3L, -3L, 0L, -3L, 3L),
			nrow=4,
			ncol=4)
	}
	
	if (typeX==3L) { # AAStringSet
		wordSize <- ceiling(log(100*quantile(widths, 0.99),
			.Call("alphabetSizeReducedAA",
				myXStringSet,
				0:19,
				PACKAGE="DECIPHER")))
		if (wordSize > 7)
			wordSize <- 7
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
	wu <- widths[u]
	cutoff <- 1 - cutoff
	
	SPKM <- function(x, K, maxVariance, maxIterations=1e3, minSize=0, tol=1e-5) {
		n <- ncol(x) # number of samples
		d <- nrow(x) # number of dimensions
		if (missing(K) && missing(maxVariance))
			stop("Either K or maxVariance must be specified.")
		if (!missing(K) && K > d)
			stop("K can be at most ", d)
		
		# initialize clusters by partitioning variance
		k <- 1L
		y <- rep(1L, n)
		keep <- seq_len(d)
		repeat {
			t <- tabulate(y, k)
			o <- order(t, decreasing=TRUE)
			for (j in seq_along(o)) {
				if (t[o[j]] < minSize)
					break
				w <- which(y == o[j])
				v <- numeric(length(d))
				for (i in seq_len(d))
					v[i] <- var(x[i, w])
				if (missing(maxVariance) || mean(v) >= maxVariance)
					break
			}
			if (t[o[j]] < minSize)
				break
			if (!missing(maxVariance) && mean(v) < maxVariance)
				break
			k <- k + 1L
			y[w[kmeans(x[keep[which.max(v[keep])], w], 2L)$cluster == 2L]] <- k
			keep <- keep[-which.max(v[keep])]
			if ((!missing(K) && k == K) ||
				length(keep) == 0L)
				break
		}
		if (k == 1L)
			return(y)
		if (verbose)
			cat("Binning sequences with spherical ", k, "-means:\n", sep="")
		
		# normalize vectors to unit length
		for (j in seq_len(n))
			x[, j] <- x[, j]/sqrt(sum(x[, j]^2))
		
		.Call("sphericalKmeans",
			x,
			y,
			k,
			tol,
			maxIterations,
			ifelse(verbose, ifelse(interactive(), TRUE, NA), verbose),
			processors,
			PACKAGE="DECIPHER")
	}
	
	overlap <- function(x, y, processors=1L) {
		.Call("matchOverlap",
			x,
			y,
			v,
			wordSize,
			processors,
			PACKAGE="DECIPHER")
	}
	
	similarity <- function(res, w1, w2, processors=1L) {
		.Call("similarities",
			res,
			w1,
			w2,
			includeTerminalGaps,
			penalizeGapLetterMatches,
			minCoverage,
			method,
			processors,
			PACKAGE="DECIPHER")
	}
	
	countHits <- function(x, y, processors=1L) {
		.Call("countOverlap",
			x,
			y,
			v,
			processors,
			PACKAGE="DECIPHER")
	}
	
	dist <- function(ali) {
		.Call("distMatrix",
			ali,
			typeX,
			includeTerminalGaps,
			penalizeGapLetterMatches,
			TRUE, # full matrix
			2L, # type = "dist"
			0, # correction (none)
			minCoverage,
			method,
			FALSE, # verbose
			NULL, # progress bar
			1L, # processors
			PACKAGE="DECIPHER")
	}
	
	align <- function(pattern,
		subject,
		anchor=NULL,
		GO=-10, # gap opening
		GE=-2, # gap extension
		TG=0, # terminal gap
		maxLength=5, # maximum length to skip alignment in equal length regions
		processors=1L) {
		n <- as.integer(length(anchor)/4) # ncol(anchor) but works when anchor is NULL
		start1 <- start2 <- end1 <- end2 <- integer(n + 1L)
		
		l1 <- 0L
		l2 <- 0L
		for (i in seq_len(n)) {
			start1[i] <- l1 + 1L
			end1[i] <- anchor[1L, i] - 1L
			l1 <- anchor[2L, i]
			start2[i] <- l2 + 1L
			end2[i] <- anchor[3L, i] - 1L
			l2 <- anchor[4L, i]
		}
		start1[n + 1L] <- l1 + 1L
		end1[n + 1L] <- width(pattern)
		start2[n + 1L] <- l2 + 1L
		end2[n + 1L] <- width(subject)
		
		results <- .Call("alignPair",
			pattern,
			subject,
			start1,
			end1,
			start2,
			end2,
			GO,
			GE,
			TG,
			maxLength,
			typeX,
			sM,
			processors,
			PACKAGE="DECIPHER")
		
		if (length(results[[3]]) > 0)
			pattern <- .Call("insertGaps",
				pattern,
				as.integer(results[[3]]),
				as.integer(results[[4]]),
				typeX,
				processors,
				PACKAGE="DECIPHER")
		if (length(results[[1]]) > 0)
			subject <- .Call("insertGaps",
				subject,
				as.integer(results[[1]]),
				as.integer(results[[2]]),
				typeX,
				processors,
				PACKAGE="DECIPHER")
		
		.append(pattern, subject)
	}
	
	if (typeX == 3L) { # AAStringSet
		a <- .Call("frequenciesReducedAA",
			.subset(myXStringSet, u),
			3L,
			c(0L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 3L, 3L, 1L, 3L, 2L, 2L, 0L, 0L, 2L, 2L, 3L),
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		a <- .Call("frequencies",
			.subset(myXStringSet, u),
			3L,
			PACKAGE="DECIPHER")
	}
	
	# convert to log-odds
	mode(a) <- "double"
	for (j in seq_len(ncol(a))) {
		w <- which(a[, j] == 0L)
		if (length(w) > 0)
			a[w, j] <- 1/length(w)
		a[, j] <- a[, j]/sum(a[, j])
	}
	a <- log(a/rowMeans(a))
	
	C <- SPKM(a, maxVariance=maxVariance, minSize=minSize)
	
	rm(a)
	gc(FALSE) # necessary for repeatable timing
	G <- vector("list", max(C))
	for (i in seq_along(G))
		G[[i]] <- which(C == i)
	ls <- lengths(G)
	maxReps <- as.integer(min(maxReps, max(ls/2)))
	if (maxReps < minSimilarities/maxAlignments)
		maxReps <- 0L
	
	if (verbose) {
		if (length(ls) > 1L) {
			time.2 <- Sys.time()
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\n")
			time.1 <- time.2
		}
		
		if (maxReps > 0L)
			cat("Ordering sequences by ",
				wordSize,
				"-mer similarity:\n",
				sep="")
		flush.console()
	}
	
	if (typeX==3L) { # AAStringSet
		V <- v <- .Call("enumerateSequenceAA",
			.subset(myXStringSet, u),
			wordSize,
			PACKAGE="DECIPHER")
	} else { # DNAStringSet or RNAStringSet
		V <- v <- .Call("enumerateSequence",
			.subset(myXStringSet, u),
			wordSize,
			TRUE, # mask repeats
			PACKAGE="DECIPHER")
	}
	V <- v <- lapply(V,
		function(x) {
			o <- order(x, na.last=NA)
			list(x[o], seq_along(x)[o])
		})
	sizes <- sapply(v, function(x) length(x[[1L]]))
	
	# order sequences by approximate similarity
	inPlay <- ls >= 2
	ksims <- psims <- vector("list", maxReps)
	batches <- lapply(ls,
		function(l) {
			batches <- seq(0L, l, avgComparisons)
			if (batches[length(batches)] < l)
				batches <- c(batches, l)
			batches
		})
	bins <- seq(0, 1, 0.01) # similarity ranges
	counts <- rep(1L, length(bins))
	var <- ls[C]/2 # doubles to avoid overflow
	keep <- rep(TRUE, l)
	avg_cor <- 0
	resTime1 <- resTime2 <- rep(NA_real_, processors)
	count1 <- count2 <- integer(processors)
	for (i in seq_len(maxReps)) {
		if (verbose && interactive()) {
			cat("\riteration ",
				i,
				" of up to ",
				maxReps,
				" (",
				formatC(100*avg_cor, digits=1, format="f"),
				"% stability) ",
				sep="")
		}
		
		if (i == 1L) {
			optProcessors1 <- processors
		} else if (i == 2L) {
			optProcessors1 <- 1L
		} else {
			w <- which(is.na(resTime1))
			if (length(w) == 0L) {
				optProcessors1 <- sample(processors, 1L, prob=1/resTime1)
			} else if (length(w) == 1L) {
				optProcessors1 <- w
			} else {
				o <- order(resTime1)
				optProcessors1 <- w[which.min(abs(mean(o[1:2]) - w))]
			}
		}
		
		m1 <- m2 <- numeric(l)
		rand <- matrix(NA_integer_, length(ls), 2L)
		time.3 <- Sys.time()
		for (k in which(inPlay)) {
			rand[k, 1L] <- G[[k]][sample(ls[k], 1L, prob=var[G[[k]]]*sizes[G[[k]]]*keep[G[[k]]])]
			keep[rand[k, 1L]] <- FALSE
			
			# process sequences in batches to reduce memory
			ov <- integer(ls[k])
			for (j in seq_len(length(batches[[k]]) - 1L)) {
				b <- (batches[[k]][j] + 1L):batches[[k]][j + 1L]
				s <- G[[k]][b]
				res1 <- overlap(rand[k, 1L], s, optProcessors1)
				ov[b] <- .Call("overlap",
					res1,
					wu[rand[k, 1L]],
					wu[s],
					PACKAGE="DECIPHER")
				m1[s] <- similarity(res1, wu[rand[k, 1L]], wu[s], optProcessors1)
			}
			
			rand[k, 2L] <- G[[k]][sample(ls[k], 1L, prob=(1.001 - m1[G[[k]]])*keep[G[[k]]]/ov)] # maximize distance and overlap to other random sequence
			keep[rand[k, 2L]] <- FALSE
			
			# process sequences in batches to reduce memory
			for (j in seq_len(length(batches[[k]]) - 1L)) {
				s <- G[[k]][(batches[[k]][j] + 1L):batches[[k]][j + 1L]]
				res2 <- overlap(rand[k, 2L], s, optProcessors1)
				m2[s] <- similarity(res2, wu[rand[k, 2L]], wu[s], optProcessors1)
			}
		}
		time.4 <- Sys.time()
		keep[rand] <- FALSE
		if (is.na(resTime1[optProcessors1])) {
			resTime1[optProcessors1] <- difftime(time.4, time.3, units='secs')/sum(ls[inPlay])
		} else {
			resTime1[optProcessors1] <- (count1[optProcessors1]*resTime1[optProcessors1] + difftime(time.4, time.3, units='secs'))/sum(ls[inPlay])/(count1[optProcessors1] + 1L)
		}
		count1[optProcessors1] <- count1[optProcessors1] + 1L
		optProcessors1 <- which.min(resTime1)
		
		if (sum(counts) - length(bins) < totalAlignments) {
			b1 <- .bincode(m1, bins, include.lowest=TRUE)
			w <- m1 > 0 & m1 < 1
			n <- min(sum(w), maxAlignments)
			if (n > 0) {
				s1 <- sample(length(m1),
					n,
					replace=TRUE, # for speed
					prob=w/counts[b1])
				s1 <- s1[!duplicated(b1[s1])]
			} else {
				s1 <- integer()
			}
			b2 <- .bincode(m2, bins, include.lowest=TRUE)
			w <- m2 > 0 & m2 < 1
			n <- min(sum(w), maxAlignments)
			if (n > 0) {
				s2 <- sample(length(m2),
					n,
					replace=TRUE, # for speed
					prob=w/counts[b2])
				s2 <- s2[!duplicated(b2[s2])]
			} else {
				s2 <- integer()
			}
			
			if (i == 1L) {
				optProcessors2 <- processors
			} else if (i == 2L) {
				optProcessors2 <- 1L
			} else {
				w <- which(is.na(resTime2))
				if (length(w) == 0L) {
					optProcessors2 <- sample(processors, 1L, prob=1/resTime2)
				} else if (length(w) == 1L) {
					optProcessors2 <- w
				} else {
					o <- order(resTime2)
					optProcessors2 <- w[which.min(abs(mean(o[1:2]) - w))]
				}
			}
			
			pdist1 <- numeric(length(s1))
			time.5 <- Sys.time()
			for (j in seq_along(s1)) {
				res <- overlap(rand[C[s1[j]], 1L], s1[j], optProcessors1) # recompute
				ali <- align(.subset(myXStringSet, u[rand[C[s1[j]], 1L]]),
					.subset(myXStringSet, u[s1[j]]),
					anchor=res[[1L]],
					processors=optProcessors2)
				pdist1[j] <- dist(ali)
			}
			time.6 <- Sys.time()
			
			pdist2 <- numeric(length(s2))
			time.7 <- Sys.time()
			for (j in seq_along(s2)) {
				res <- overlap(rand[C[s2[j]], 2L], s2[j], optProcessors1) # recompute
				ali <- align(.subset(myXStringSet, u[rand[C[s2[j]], 2L]]),
					.subset(myXStringSet, u[s2[j]]),
					anchor=res[[1L]],
					processors=optProcessors2)
				pdist2[j] <- dist(ali)
			}
			time.8 <- Sys.time()
			if (is.na(resTime2[optProcessors2])) {
				resTime2[optProcessors2] <- difftime(time.8, time.7, units='secs') + difftime(time.6, time.5, units='secs')
			} else {
				resTime2[optProcessors2] <- (count2[optProcessors2]*resTime2[optProcessors2] + difftime(time.8, time.7, units='secs') + difftime(time.6, time.5, units='secs'))/(count2[optProcessors2] + 1L)
			}
			count2[optProcessors2] <- count2[optProcessors2] + 1L
			
			counts <- counts + tabulate(b1[s1], length(counts))
			counts <- counts + tabulate(b2[s2], length(counts))
			ksims[[i]] <- c(m1[s1], m2[s2])
			psims[[i]] <- 1 - c(pdist1, pdist2)
		}
		
		d <- m1 - m2
		if (i > 1L) {
			for (k in which(inPlay)) {
				if (diff(range(d[G[[k]]])) > 0 && cor(rS[G[[k]]], d[G[[k]]]) < 0) { # negatively correlated
					rS[G[[k]]] <- rS[G[[k]]] - d[G[[k]]]
				} else {
					rS[G[[k]]] <- rS[G[[k]]] + d[G[[k]]]
				}
			}
			
			o <- order(C,
				rS,
				wu,
				t, # frequency
				decreasing=TRUE,
				method="radix")
			r <- o
			r[r] <- seq_along(r) # rank
			
			w <- which(d != 0) # only update positions that could have moved
			var[w] <- alpha*abs(R[w] - r[w]) + (1 - alpha)*var[w] # exponential moving average of change in rank order
			avg_cor <- sum(var*keep < avgComparisons/2)/l # exponential moving average of fraction stabilized in rank order
			
			O <- o
			R <- r
			
			for (k in which(inPlay)) {
				if (all(var[G[[k]]] < avgComparisons/2)) {
					inPlay[k] <- FALSE # 100% stability
				} else if (sum(keep[G[[k]]]) < 2) {
					inPlay[k] <- FALSE # no more pairs
				}
			}
		} else {
			rS <- d # relative similaity
			
			# initialize rank order
			O <- order(C,
				rS,
				wu,
				t, # frequency
				decreasing=TRUE,
				method="radix")
			R <- O
			R[R] <- seq_along(R) # rank
		}
		
		if (!any(inPlay) ||
			avg_cor >= 0.9995) # rounds to 100%
			break # reached rank order stability
	}
	if (maxReps > 0L) {
		optProcessors2 <- which.min(resTime2)
	} else {
		optProcessors1 <- optProcessors2 <- processors
		O <- R <- seq_along(u)
		rS <- numeric(length(O))
	}
	
	# fit conversion from approximate to actual similarity
	ksims <- unlist(ksims)
	psims <- unlist(psims)
	w <- which(psims > 0) # potentially insufficiently overlapping
	ksims <- ksims[w]
	psims <- psims[w]
	interval <- ceiling(interval*length(ksims))
	if (interval < minSimilarities) {
		limits <- numeric(length(cutoff)) # no limit
	} else {
		psims <- psims - 1
		ksims <- ksims - 1
		slope <- sum(psims*ksims)/sum(psims^2)
		delta <- 10
		while (delta > 0.0001) {
			temp <- slope + delta
			if (interval > sum(ksims >= temp*(psims - buffer))) {
				slope <- temp
			} else {
				delta <- delta/2
			}
		}
		limits <- slope*(cutoff - 1) + 1 - slope*buffer
	}
	
	var[var > maxComparisons] <- maxComparisons
	var <- var[O]/mean(var) # normalized variability
	bL <- seq_len(l) - as.integer(var*(avgComparisons/2))
	bR <- seq_len(l) + as.integer(var*(avgComparisons/2))
	bL[bL < 1L] <- 1L
	bR[bR > l] <- l
	
	if (verbose) {
		if (maxReps > 0L) {
			cat("\riteration ",
				i,
				" of up to ",
				maxReps,
				" (",
				formatC(100*avg_cor, digits=1, format="f"),
				"% stability) ",
				sep="")
			
			time.2 <- Sys.time()
			cat("\n\n")
			print(round(difftime(time.2,
				time.1,
				units='secs'),
				digits=2))
			cat("\n")
			
			time.1 <- time.2
		}
		lastValue <- 0
		cat("Clustering sequences by similarity:\n",
			sep="")
		flush.console()
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	}
	
	# use slope of relative distance as a tie-breaker
	if (length(rS) >= 4L) {
		rS <- rS[O]
		count <- 0L
		for (i in rev(seq_along(ls))) { # groups reversed because sorted decreasing
			s <- (count + 1L):(count + ls[i])
			S <- smooth.spline(rS[s],
				spar=0.5)
			rS[s] <- predict(S,
				x=seq_along(s),
				deriv=1L)$y
			count <- count + ls[i]
		}
	}
	
	P <- order(sizes[O],
		t[O],
		rS,
		decreasing=TRUE,
		method="radix")
	P <- O[P]
	Q <- R[P]
	
	for (i in seq_len(lc)) {
		if (!ASC && i > 1) {
			P <- order(c[u, i - 1][O],
				sizes[O],
				t[O],
				rS,
				decreasing=TRUE,
				method="radix")
			P <- O[P]
			Q <- R[P]
		}
		
		v <- V[P] # reorder k-mers
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
		recent <- integer(avgComparisons)
		count <- 1L
		recent[count] <- 1L
		
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
					if (invertCenters)
						c[p[j], i] <- -c[p[j], i]
					seeds.index[Q[j]] <- j
					recent <- integer(avgComparisons)
					count <- 1L
					recent[count] <- j
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
			
			compare <- Q[recent] # most recent centers
			compare <- compare[.Call("radixOrder", abs(Q[j] - compare), 1L, PACKAGE="DECIPHER")] # order by proximity
			compare <- c(bL[Q[j]]:bR[Q[j]], # surrounding centers
				compare)
			compare <- compare[seeds.index[compare] != 0L]
			if (!ASC && i > 1) { # check bounds
				if (invertCenters) {
					compare <- compare[abs(c[p[seeds.index[compare]], i - 1L]) == abs(c[p[j], i - 1L])]
				} else {
					compare <- compare[c[p[seeds.index[compare]], i - 1L] == c[p[j], i - 1L]]
				}
			}
			if (length(compare) == 0L) {
				cluster_num <- cluster_num + 1L
				c[p[j], i] <- cluster_num + offset
				if (invertCenters)
					c[p[j], i] <- -c[p[j], i]
				seeds.index[Q[j]] <- j
				count <- count + 1L
				if (count > avgComparisons)
					count <- 1L
				recent[count] <- j
				next
			}
			compare <- head(compare, bR[Q[j]] - bL[Q[j]] + 1L)
			compare <- compare[!duplicated(seeds.index[compare])]
			
			counts <- countHits(j, seeds.index[compare], optProcessors1)
			s <- .Call("radixOrder", counts, 0L, PACKAGE="DECIPHER")
			m <- numeric(length(s))
			res <- vector("list", length(s))
			batches <- seq.int(0L, length(s), batchSize)
			if (batches[length(batches)] < length(s))
				batches <- c(batches, length(s))
			for (k in seq_len(length(batches) - 1L)) {
				b <- (batches[k] + 1L):batches[k + 1L]
				r <- overlap(j, seeds.index[compare[s[b]]], optProcessors1)
				M <- similarity(r, wu[P[j]], wu[P[seeds.index[compare[s[b]]]]], optProcessors1)
				res[s[b]] <- r
				m[s[b]] <- M
				if (any(M >= cutoff[i])) # met similarity threshold
					break
				if (k < length(batches) - 1L) {
					# check for plausibility of reaching limit
					upper <- max(M)
					if (upper < limits[i] &&
						length(M) == batchSize &&
						upper + 2*sd(M) < limits[i])
						break # unlikely to reach limit
					# check for shift to lower similarities
					if (k == 1L) {
						lower <- min(M)
					} else {
						if (lower > upper)
							break # unlikely to exceed previous
					}
				}
			}
			
			b <- seq_len(batches[k + 1L])
			w <- s[sort.list(m[s[b]], method="radix", decreasing=TRUE)]
			if (length(w) > maxAlignments)
				w <- w[seq_len(maxAlignments)]
			W <- which(m[w] >= limits[i])
			if (length(W) >= minAlignments) {
				w <- w[W]
			} else {
				w <- head(w, minAlignments)
			}
			
			if (m[w[1L]] < cutoff[i]) {
				for (k in seq_along(w)) {
					if (k > attempts) {
						s <- seq(to=k - 1L, length.out=attempts)
						if (max(m[w[s]]) + sd(m[w[s]]) < cutoff[i]) {
							w <- w[seq_len(k - 1L)] # prevent matching clusters without alignment
							break
						}
					}
					ali <- align(.subset(myXStringSet, p[j]),
						.subset(myXStringSet, p[seeds.index[compare[w[k]]]]),
						anchor=res[[w[k]]],
						processors=optProcessors2)
					m[w[k]] <- dist(ali)
					if (is.na(m[w[k]])) {
						m[w[k]] <- 0
					} else {
						m[w[k]] <- 1 - m[w[k]]
					}
					if (m[w[k]] >= cutoff[i])
						break
				}
				w <- w[which.max(m[w])]
			} else {
				w <- w[1L] # max similarity
			}
			
			if (length(w) == 0L ||
				m[w] < cutoff[i]) { # form a new group
				cluster_num <- cluster_num + 1L
				c[p[j], i] <- cluster_num + offset
				if (invertCenters)
					c[p[j], i] <- -c[p[j], i]
				seeds.index[Q[j]] <- j
				count <- count + 1L
				if (count > avgComparisons)
					count <- 1L
				recent[count] <- j
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
		setTxtProgressBar(pBar, 1)
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
