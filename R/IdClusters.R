IdClusters <- function(myXStringSet,
	cutoff=0,
	processors=1,
	verbose=TRUE) {
	
	# initialize variables
	time.1 <- Sys.time()
	
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
	a <- vcountPattern("+", myXStringSet)
	if (any(a > 0))
		stop("Mask characters ('+') are not allowed in myXStringSet.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') are not allowed in myXStringSet.")
	if (all(width(myXStringSet)==0L))
		stop("All sequences in myXStringSet are zero width.")
	
	if (verbose) {
		lastValue <- 0
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
	}
	if (typeX==3L) { # AAStringSet
		wordSize <- ceiling(log(500*quantile(width(myXStringSet), 0.99),
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
		wordSize <- ceiling(log(500*quantile(width(myXStringSet), 0.99),
			.Call("alphabetSize",
				myXStringSet,
				PACKAGE="DECIPHER")))
		if (wordSize > 15)
			wordSize <- 15
		if (wordSize < 1)
			wordSize <- 1
		words <- 4^wordSize
	}
	
	l <- length(myXStringSet)
	if (l==0)
		stop("myXStringSet contains no sequences.")
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
	
	cutoff <- cutoff/2 # cluster radius is half the diameter
	threshold <- (1 - cutoff)^wordSize # probability of k-mer match
	for (i in seq_len(lc)) {
		if (!ASC && i > 1) {
			o <- u[order(c[u, i - 1],
				width(myXStringSet)[u],
				t, # frequency
				decreasing=TRUE)]
		} else {
			o <- u[order(width(myXStringSet)[u],
				t, # frequency
				decreasing=TRUE)]
		}
		
		if (typeX==3L) { # AAStringSet
			v <- .Call("enumerateSequenceAA",
				.subset(myXStringSet, o),
				wordSize,
				PACKAGE="DECIPHER")
		} else { # DNAStringSet or RNAStringSet
			v <- .Call("enumerateSequence",
				.subset(myXStringSet, o),
				wordSize,
				FALSE, # mask repeats
				PACKAGE="DECIPHER")
		}
		v <- lapply(v,
			sort,
			na.last=NA)
		
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
			
			m <- .Call("matchListsDual",
				v[j],
				v[seeds.index[seq_len(cluster_num)]],
				FALSE, # verbose
				NULL, # pBar
				processors,
				PACKAGE="DECIPHER")
			w <- which.max(m)
			
			if (length(w)==0 ||
				m[w] < threshold[i]) { # form a new group
				cluster_num <- cluster_num + 1L
				c[o[j], i] <- cluster_num
				seeds.index[cluster_num] <- j
			} else { # part of an existing group
				c[o[j], i] <- w + offset
				if (!ASC && i < lc) {
					cols <- (i + 1):lc
					cols <- cols[which(m[w] >= threshold[cols])]
					if (length(cols) > 0) # assign forward
						c[o[j], cols] <- o[seeds.index[w]]
				}
			}
		}
		
		c[, i] <- c[x, i]
	}
	myClusters <- as.data.frame(c)
	
	if (verbose) {
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
