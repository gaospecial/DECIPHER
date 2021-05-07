DetectRepeats <- function(myXStringSet,
	type="tandem",
	minScore=10,
	allScores=FALSE,
	maxPeriod=10000,
	maxFailures=2,
	maxShifts=5,
	alphabet=AA_REDUCED[[125]],
	processors=1,
	verbose=TRUE,
	...) {
	
	# error checking
	if (length(type)==0)
		stop("No type specified.")
	if (length(type) > 1)
		stop("Only one type may be specified.")
	TYPES <- c("tandem", "interspersed", "both")
	type <- pmatch(type, TYPES)
	if (is.na(type))
		stop("Invalid type.")
	if (type==-1)
		stop("Ambiguous type.")
	if (is(myXStringSet, "DNAStringSet")) {
		xtype <- 1L
	} else if (is(myXStringSet, "RNAStringSet")) {
		xtype <- 2L
	} else if (is(myXStringSet, "AAStringSet")) {
		if (type > 1L)
			stop("type must be 'tandem' when myXStringSet is an AAStringSet.")
		xtype <- 3L
		
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
		n <- as.integer(floor(log(4294967295, sizeAA)))
		alphabet <- alphabet - 1L
	} else {
		stop("myXStringSet must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	}
	if (!is.numeric(minScore) || length(minScore) > 1)
		stop("minScore must be a numeric.")
	if (!is.logical(allScores))
		stop("allScores must be a logical.")
	if (!is.numeric(maxPeriod) || length(maxPeriod) > 1)
		stop("maxPeriod must be a numeric.")
	if (!is.numeric(maxFailures))
		stop("maxFailures must be a numeric.")
	if (maxFailures != floor(maxFailures))
		stop("maxFailures must be a whole number.")
	if (maxFailures < 0)
		stop("maxFailures must be at least zero.")
	if (!is.numeric(maxShifts))
		stop("maxShifts must be a numeric.")
	if (maxShifts != floor(maxShifts))
		stop("maxShifts must be a whole number.")
	if (maxShifts < 0)
		stop("maxShifts must be at least zero.")
	if (maxPeriod < 1)
		stop("maxPeriod must be at least one.")
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	a <- vcountPattern("-", myXStringSet)
	if (any(a > 0))
		stop("Gap characters ('-') are not allowed in myXStringSet.")
	a <- vcountPattern(".", myXStringSet)
	if (any(a > 0))
		stop("Unknown characters ('.') are not allowed in myXStringSet.")
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
	
	if (verbose)
		time.1 <- Sys.time()
	
	# default parameters
	N <- 10 # find k-mers on average N times per maxPeriod
	gapCost <- c(4, 1.5) # cost of opening and extending each gap in the alignment of repeats
	maxExact <- 1e4 # maximum combinations to use exact scoring
	mult <- 2L # multiplier on K for lookahead length
	ratio <- 0.2 # ratio of repeat length to periodicity worth aligning
	transition <- 0.3 # relative similarity for transitions
	
	# multinomial test
	.multinomial <- function(x,
		p,
		gapCost) {
		
		if (xtype == 3L) {
			m <- .Call("consensusProfileAA",
						x,
						rep(1, length(x)),
						NULL,
						PACKAGE="DECIPHER")
			m2 <- rowsum(m[1:20,, drop=FALSE], alphabet)
			m2 <- t(m2)
			gapEx <- (1 - m[24,]*m[27,])*m[27,]
			gapOp <- (m[, 25] > 0) + (m[, 26] > 0)
			m <- m2/(1 - m[24,])*gapEx
		} else {
			m <- .Call("consensusProfile",
						x,
						rep(1, length(x)),
						NULL,
						PACKAGE="DECIPHER")
			m <- t(m)
			gapEx <- (1 - m[, 5]*m[, 8])*m[, 8]
			gapOp <- (m[, 6] > 0) + (m[, 7] > 0)
			m <- m[, 1:4, drop=FALSE]/(1 - m[, 5])*gapEx
		}
		
		m <- m*length(x)
		a <- apply(m,
			1,
			function(x) {
				sx <- sum(x)
				if (sx <= length(ecdfs)) { # use exact
					y <- dmultinom(x, prob=p)
					y <- ecdfs[[ceiling(sx)]][[1L]](y)*y
					y <- y*ecdfs[[ceiling(sx)]][[2L]]
					log(y)
				} else { # approximate
					chisq <- sum(x)*p
					chisq <- sum((x - chisq)^2/chisq)
					pchisq(chisq,
						length(x) - 1,
						lower.tail=FALSE,
						log.p=TRUE)
				}
			})
		a <- a + gapCost[2]*(1 - gapEx)*length(x)
		a <- a + gapCost[1]*gapOp*length(x)
		
		rS <- rowSums(m)
		score <- -sum(log(rS))
		
		score <- score - sum(a, na.rm=TRUE)
		
		similarity <- .Call(functionCall,
			x,
			sM,
			0, # gap opening
			0, # gap extension
			rep(1, length(x)), # weights
			NULL, # structures
			numeric(), # structure matrix
			PACKAGE="DECIPHER")
		rS <- rS*(rS - 1)/2 # normalization constant
		similarity <- similarity/rS
		score*mean(similarity, na.rm=TRUE)
	}
	
	if (xtype == 3L) {
		bg <- alphabetFrequency(myXStringSet,
			as.prob=TRUE,
			collapse=TRUE)[1:20]
		bg <- bg/sum(bg)
		bg <- tapply(bg, alphabet, sum)
		
		sM <- matrix(c(0.83, 0.23, 0.22, 0.21, 0.31, 0.26, 0.25, 0.3, 0.2, 0.23, 0.22, 0.23, 0.26, 0.18, 0.25, 0.33, 0.29, 0.15, 0.17, 0.29, 0, 0.23, 0.95, 0.29, 0.25, 0.12, 0.37, 0.31, 0.18, 0.33, 0.1, 0.13, 0.44, 0.18, 0.09, 0.22, 0.26, 0.26, 0.17, 0.19, 0.13, 0, 0.22, 0.29, 0.95, 0.41, 0.15, 0.33, 0.32, 0.29, 0.34, 0.07, 0.07, 0.33, 0.15, 0.09, 0.23, 0.35, 0.3, 0.09, 0.19, 0.1, 0, 0.21, 0.25, 0.41, 0.97, 0.06, 0.32, 0.42, 0.26, 0.28, 0, 0.02, 0.3, 0.08, 0.01, 0.26, 0.31, 0.26, 0.04, 0.11, 0.05, 0, 0.31, 0.12, 0.15, 0.06, 1.32, 0.11, 0.06, 0.17, 0.19, 0.25, 0.24, 0.09, 0.26, 0.24, 0.1, 0.27, 0.24, 0.19, 0.23, 0.29, 0, 0.26, 0.37, 0.33, 0.32, 0.11, 0.91, 0.4, 0.2, 0.34, 0.11, 0.14, 0.38, 0.23, 0.09, 0.23, 0.3, 0.28, 0.13, 0.19, 0.15, 0, 0.25, 0.31, 0.32, 0.42, 0.06, 0.4, 0.91, 0.2, 0.27, 0.07, 0.08, 0.36, 0.14, 0.03, 0.25, 0.29, 0.27, 0.07, 0.13, 0.11, 0, 0.3, 0.18, 0.29, 0.26, 0.17, 0.2, 0.2, 1.01, 0.19, 0.04, 0.06, 0.2, 0.11, 0.08, 0.21, 0.3, 0.21, 0.09, 0.09, 0.09, 0, 0.2, 0.33, 0.34, 0.28, 0.19, 0.34, 0.27, 0.19, 1.12, 0.11, 0.14, 0.29, 0.19, 0.22, 0.21, 0.26, 0.24, 0.24, 0.37, 0.14, 0, 0.23, 0.1, 0.07, 0, 0.25, 0.11, 0.07, 0.04, 0.11, 0.88, 0.42, 0.1, 0.38, 0.34, 0.11, 0.12, 0.22, 0.21, 0.23, 0.46, 0, 0.22, 0.13, 0.07, 0.02, 0.24, 0.14, 0.08, 0.06, 0.14, 0.42, 0.86, 0.11, 0.42, 0.37, 0.11, 0.12, 0.19, 0.26, 0.25, 0.36, 0, 0.23, 0.44, 0.33, 0.3, 0.09, 0.38, 0.36, 0.2, 0.29, 0.1, 0.11, 0.9, 0.17, 0.06, 0.25, 0.29, 0.28, 0.1, 0.16, 0.13, 0, 0.26, 0.18, 0.15, 0.08, 0.26, 0.23, 0.14, 0.11, 0.19, 0.38, 0.42, 0.17, 0.99, 0.35, 0.13, 0.19, 0.25, 0.26, 0.27, 0.34, 0, 0.18, 0.09, 0.09, 0.01, 0.24, 0.09, 0.03, 0.08, 0.22, 0.34, 0.37, 0.06, 0.35, 1, 0.1, 0.13, 0.17, 0.42, 0.49, 0.3, 0, 0.25, 0.22, 0.23, 0.26, 0.1, 0.23, 0.25, 0.21, 0.21, 0.11, 0.11, 0.25, 0.13, 0.1, 1.09, 0.28, 0.24, 0.11, 0.11, 0.15, 0, 0.33, 0.26, 0.35, 0.31, 0.27, 0.3, 0.29, 0.3, 0.26, 0.12, 0.12, 0.29, 0.19, 0.13, 0.28, 0.84, 0.38, 0.12, 0.17, 0.17, 0, 0.29, 0.26, 0.3, 0.26, 0.24, 0.28, 0.27, 0.21, 0.24, 0.22, 0.19, 0.28, 0.25, 0.17, 0.24, 0.38, 0.87, 0.14, 0.19, 0.27, 0, 0.15, 0.17, 0.09, 0.04, 0.19, 0.13, 0.07, 0.09, 0.24, 0.21, 0.26, 0.1, 0.26, 0.42, 0.11, 0.12, 0.14, 1.33, 0.46, 0.19, 0, 0.17, 0.19, 0.19, 0.11, 0.23, 0.19, 0.13, 0.09, 0.37, 0.23, 0.25, 0.16, 0.27, 0.49, 0.11, 0.17, 0.19, 0.46, 1.07, 0.22, 0, 0.29, 0.13, 0.1, 0.05, 0.29, 0.15, 0.11, 0.09, 0.14, 0.46, 0.36, 0.13, 0.34, 0.3, 0.15, 0.17, 0.27, 0.19, 0.22, 0.86, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.33),
				nrow=21,
				ncol=21,
				dimnames=list(c(AA_STANDARD, "*"), c(AA_STANDARD, "*")))
		
		functionCall <- "colScoresAA"
	} else {
		bg <- oligonucleotideFrequency(myXStringSet,
			width=1,
			as.prob=TRUE,
			simplify.as="collapse")
		
		sM <- diag(length(DNA_BASES))
		diag(sM) <- 1 - log(bg*length(DNA_BASES))
		dimnames(sM) <- list(DNA_BASES, DNA_BASES)
		sM["A", "G"] <- sM["G", "A"] <- transition
		sM["C", "T"] <- sM["T", "C"] <- transition
		
		functionCall <- "colScores"
	}
	
	.combos <- function(n, k) {
		if (n == 1L) {
			k
		} else {
			res <- matrix(rep(0L, n - 1L), 1L)
			for (i in seq_len(k))
				res <- rbind(res,
					.combos(n - 1L, i))
			cbind(res,
				k - rowSums(res))
		}
	}
	sizes <- choose(length(bg) + 0:49, # n + k - 1
		length(bg) - 1) # k - 1
	exact <- which.max(sizes > maxExact) - 1L
	exact <- min(20, exact) # convergence to approximate
	ecdfs <- vector("list", exact)
	for (i in seq_len(exact)) {
		a <- apply(.combos(length(bg), i),
			1L,
			dmultinom,
			prob=bg)
		ecdfs[[i]] <- list(ecdf(a), length(a))
	}
	
	if (type == 1L || type == 3L) { # tandem repeats
		if (verbose) {
			if (type > 1L) {
				cat("Detecting tandem repeats:\n")
				flush.console()
			}
			pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
			totW <- c(0, cumsum(width(myXStringSet)))
		}
		
		if (xtype == 3L) {
			size <- .Call("alphabetSizeReducedAA",
				myXStringSet,
				alphabet,
				PACKAGE="DECIPHER")
		} else {
			size <- .Call("alphabetSize",
				myXStringSet,
				PACKAGE="DECIPHER")
		}
		K <- as.integer(ceiling(log(maxPeriod/N, size)))
		if (xtype == 3L) {
			if (K > n)
				K <- n
			kmers <- .Call("enumerateSequenceReducedAA",
				myXStringSet,
				K,
				alphabet,
				FALSE, # mask repeats
				PACKAGE="DECIPHER")
		} else {
			if (K > 16L)
				K <- 16L
			kmers <- .Call("enumerateSequence",
				myXStringSet,
				K,
				FALSE, # mask repeats
				PACKAGE="DECIPHER")
		}
		
		o <- lapply(kmers, order, method="radix")
		d <- lapply(o, diff)
		d <- mapply(function(a, b) {
					a[a] <- c(b, NA)
					a
				},
			o,
			d,
			SIMPLIFY=FALSE)
		r <- mapply(function(x, y) {
				start <- NA # start of current repeat
				off <- mult*K # lookahead length
				values <- numeric(length(x)) # periodicity
				lengths <- numeric(length(x)) # total evidence
				i <- i2 <- 1L # current and next index
				while (i < length(x)) {
					if (!is.na(x[i]) &&
						!is.na(y[i]) &&
						x[i] > 0) {
						if (is.na(start)) {
							start <- i
							divisor <- x[i]
						} else if (x[i] < divisor) {
							if ((divisor/x[i]) %in% 1:3) {
								divisor <- x[i]
								if (i + 2*divisor <= length(y) &&
									!is.na(y[i + 2*divisor]) &&
									y[i + 2*divisor] == y[i]) {
									i2 <- i + divisor
								} else if (i + off <= length(x) &&
									!is.na(x[i + off]) &&
									x[i + off] == x[i]) {
									i2 <- i + off
								}
							} else {
								if (i2 > i && i2 < i + off) {
									i <- i2
									next
								}
								values[start] <- divisor
								lengths[start] <- i - start
								start <- i
								divisor <- x[i]
							}
						} else if (x[i] > divisor) {
							if ((x[i]/divisor) %in% 1:3) {
								if (i + 2*divisor <= length(y) &&
									!is.na(y[i + 2*divisor]) &&
									y[i + 2*divisor] == y[i]) {
									i2 <- i + divisor
								} else if (i + off <= length(x) &&
									!is.na(x[i + off]) &&
									x[i + off] == x[i]) {
									i2 <- i + off
								}
							} else {
								if (i2 > i && i2 < i + off) {
									i <- i2
									next
								}
								values[start] <- divisor
								lengths[start] <- i - start
								start <- i
								divisor <- x[i]
							}
						} else { # x[i] == divisor
							if (i + 2*divisor <= length(y) &&
								!is.na(y[i + 2*divisor]) &&
								y[i + 2*divisor] == y[i]) {
								i2 <- i + divisor
							} else if (i + off <= length(x) &&
								!is.na(x[i + off]) &&
								x[i + off] == x[i]) {
								i2 <- i + off
							}
						}
					} else if (!is.na(start)) {
						if (i2 > i && i2 < i + off) {
							i <- i2
							next
						}
						values[start] <- divisor
						lengths[start] <- i - start
						start <- NA
					}
					
					i <- i + 1L
				}
				list(values=values, lengths=lengths)
			},
			d,
			kmers,
			SIMPLIFY=FALSE)
		
		shift <- function(x, LnL, posL, posR, maxShifts) {
			# optimize the phase
			tempL <- posL
			tempR <- posR
			POSL <- posL - 1L
			POSR <- posR - 1L
			attempts <- 0L
			while (POSL[1L] >= 1L &&
				attempts <= maxShifts) {
				y <- .Call("replaceGaps",
					x,
					myXString,
					POSL,
					xtype,
					PACKAGE="DECIPHER")
				temp <- .multinomial(y, bg, gapCost)
				if (temp > LnL) {
					x <- y
					LnL <- temp
					posL <- POSL
					posR <- POSR
					attempts <- 0L
				} else {
					attempts <- attempts + 1L
				}
				POSL <- POSL - 1L
				POSR <- POSR - 1L
			}
			
			POSL <- tempL + 1L
			POSR <- tempR + 1L
			attempts <- 0L
			while (POSR[length(POSR)] <= length(myXString) &&
				attempts <= maxShifts) {
				y <- .Call("replaceGaps",
					x,
					myXString,
					POSL,
					xtype,
					PACKAGE="DECIPHER")
				temp <- .multinomial(y, bg, gapCost)
				if (temp > LnL) {
					x <- y
					LnL <- temp
					posL <- POSL
					posR <- POSR
					attempts <- 0L
				} else {
					attempts <- attempts + 1L
				}
				POSL <- POSL + 1L
				POSR <- POSR + 1L
			}
			
			list(x, LnL, posL, posR)
		}
		
		result <- list()
		for (k in seq_along(r)) {
			myXString <- myXStringSet[[k]]
			values <- r[[k]]$values
			lengths <- r[[k]]$lengths
			w <- which(lengths/values > ratio &
				values <= maxPeriod)
			res <- vector("list", length(w))
			for (i in seq_along(w)) {
				posL <- seq(0,
					lengths[w[i]] + values[w[i]],
					values[w[i]])
				posL <- posL + w[i]
				posR <- posL + values[w[i]] - 1L
				keep <- posL <= length(myXString) &
					posR <= length(myXString)
				if (sum(keep) < length(keep)) {
					posL <- posL[keep]
					posR <- posR[keep]
				}
				
				delta <- as.integer(values[w[i]]/2)
				if (length(posR) > 1) { # extend unknown right bound
					posR[length(posR)] <- posR[length(posR)] + delta
					if (posR[length(posR)] > length(myXString))
						posR[length(posR)] <- length(myXString)
				}
				
				if (length(posR) < 2L) {
					res[[i]] <- list(posL, posR, -Inf, k)
					next
				}
				
				# align the repeats
				x <- extractAt(myXString, IRanges(posL, posR))
				ux <- unique(x)
				if (length(ux) > 1) {
					index <- match(x, ux)
					ux <- AlignSeqs(ux,
						terminalGap=c(-1000, -5),
						iterations=0,
						refinements=0,
						anchor=NA,
						processors=processors,
						verbose=FALSE)
					x <- ux[index]
					
					t <- TerminalChar(x)
					off <- min(t[-nrow(t), "trailingChar"])
					x <- subseq(x, end=width(x)[1L] - off)
					posR[length(posR)] <- posR[length(posR)] - off
				}
				
				# calculate the likelihood
				LnL <- .multinomial(x, bg, gapCost)
				out <- shift(x, LnL, posL, posR, min(values[w[i]], maxShifts))
				x <- out[[1L]]
				LnL <- out[[2L]]
				posL <- out[[3L]]
				posR <- out[[4L]]
				if (LnL < minScore/2) {
					res[[i]] <- list(posL, posR, LnL, k)
					next
				}
				
				# try extending left
				start <- posL[1L] - values[w[i]]
				attempts <- 0L
				X <- x
				POSL <- posL
				POSR <- posR
				while (start >= 1L &&
					attempts <= maxFailures) {
					start <- start - delta
					if (start < 1)
						start <- 1L
					y <- AlignProfiles(extractAt(myXString,
							IRanges(start, POSL[1L] - 1L)),
						X,
						anchor=NA,
						processors=processors)
					t <- TerminalChar(y)
					off <- min(t[-1L, "leadingChar"])
					y <- subseq(y, off + 1L)
					temp <- .multinomial(y, bg, gapCost)
					start <- start + off
					POSR <- c(POSL[1L] - 1L, POSR)
					POSL <- c(start, POSL)
					
					out <- shift(y, temp, POSL, POSR, min(values[w[i]], maxShifts))
					y <- out[[1L]]
					temp <- out[[2L]]
					POSL <- out[[3L]]
					POSR <- out[[4L]]
					
					X <- y
					if (temp > LnL) {
						LnL <- temp
						posL <- POSL
						posR <- POSR
						x <- y
						attempts <- 0L
					} else {
						attempts <- attempts + 1L
					}
					start <- POSL[1L] - values[w[i]]
				}
				
				# try extending right
				end <- posR[length(posR)] + values[w[i]]
				attempts <- 0L
				X <- x
				POSL <- posL
				POSR <- posR
				while (end <= length(myXString) &&
					attempts <= maxFailures) {
					end <- end + delta
					if (end > length(myXString))
						end <- length(myXString)
					y <- AlignProfiles(X,
						extractAt(myXString,
							IRanges(POSR[length(POSR)] + 1L, end)),
						anchor=NA,
						processors=processors)
					t <- TerminalChar(y)
					off <- min(t[-nrow(t), "trailingChar"])
					y <- subseq(y, end=width(y)[1L] - off)
					temp <- .multinomial(y, bg, gapCost)
					end <- end - off
					POSL <- c(POSL, POSR[length(POSR)] + 1L)
					POSR <- c(POSR, end)
					
					out <- shift(y, temp, POSL, POSR, min(values[w[i]], maxShifts))
					y <- out[[1L]]
					temp <- out[[2L]]
					POSL <- out[[3L]]
					POSR <- out[[4L]]
					
					X <- y
					if (temp > LnL) {
						LnL <- temp
						posL <- POSL
						posR <- POSR
						x <- y
						attempts <- 0L
					} else {
						attempts <- attempts + 1L
					}
					end <- POSR[length(POSR)] + values[w[i]]
				}
				
				# record the result
				res[[i]] <- list(posL, posR, LnL, k)
				if (verbose)
					setTxtProgressBar(pBar,
						(totW[k] + w[i])/totW[length(totW)])
			}
			
			result <- c(result, res)
		}
		
		result <- unique(result)
		posL <- lapply(result, `[[`, 1L)
		posR <- lapply(result, `[[`, 2L)
		result <- data.frame(Index=sapply(result, `[[`, 4L),
			Begin=sapply(posL, `[`, 1L),
			End=sapply(posR, tail, n=1L),
			Left=I(posL),
			Right=I(posR),
			Score=sapply(result, `[[`, 3L))
		result <- result[result[, "Score"] >= minScore,]
		result <- result[order(result[, "Index"], result[, "Begin"]),]
		
		if (!allScores) {
			# pick highest score when overlapping
			if (nrow(result) > 1) {
				keep <- logical(nrow(result))
				keep [1L] <- TRUE
				i <- 1L
				j <- 2L
				while (j <= nrow(result)) {
					if (result[i, "Index"] == result[j, "Index"] &&
						result[i, "End"] >= result[j, "Begin"]) {
						if (result[i, "Score"] < result[j, "Score"]) {
							keep[i] <- FALSE
							keep[j] <- TRUE
							i <- j
						}
					} else {
						keep[j] <- TRUE
						i <- j
					}
					j <- j + 1L
				}
				result <- result[keep,]
			}
		}
		rownames(result) <- NULL
		
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
			if (type > 1L)
				cat("Detecting interspersed repeats:\n")
			flush.console()
		}
	}
	
	if (type > 1L) { # interspersed repeats
		dbConn <- dbConnect(SQLite(), ":memory:")
		Seqs2DB(myXStringSet,
			"XStringSet",
			dbConn,
			"1",
			verbose=FALSE)
		
		syn <- FindSynteny(dbConn,
			identifier=c("1", "1"),
			verbose=FALSE,
			...)
		
		# remove results within maxPeriod
		w <- which(syn[[2, 1]][, "index1"] == syn[[2, 1]][, "index2"] &
			abs(syn[[2, 1]][, "start1"] - syn[[2, 1]][, "start2"]) > maxPeriod &
			abs(syn[[2, 1]][, "end1"] - syn[[2, 1]][, "end2"]) > maxPeriod |
			syn[[2, 1]][, "index1"] != syn[[2, 1]][, "index2"])
		if (length(w) > 0) {
			hits <- mapply(seq,
				syn[[2, 1]][w, "first_hit"],
				syn[[2, 1]][w, "last_hit"],
				SIMPLIFY=FALSE)
			syn[[2, 1]] <- syn[[2, 1]][w,]
			syn[[1, 2]] <- syn[[1, 2]][unlist(hits),]
			hits <- cumsum(lengths(hits))
			syn[[2, 1]][, "first_hit"] <- c(1L, hits[-length(hits)] + 1L)
			syn[[2, 1]][, "last_hit"] <- hits
			
			ali <- AlignSynteny(syn,
				dbConn,
				verbose=verbose)
			
			ali <- ali[[1L]]
			res2 <- data.frame(syn[[2, 1]][, 1:8])
			if (length(ali) > 0) {
				res2[, "score"] <- sapply(ali,
					.multinomial,
					p=bg,
					gapCost=gapCost)
				res2 <- res2[res2[, "score"] >= minScore,]
				
				if (!allScores) {
					# eliminate overlapping repeats
					keep <- rep(TRUE, nrow(res2))
					o <- order(res2[, "score"], decreasing=TRUE)
					range1 <- IRanges(res2[, "start1"],
						res2[, "end1"])
					range2 <- IRanges(res2[, "start2"],
						res2[, "end2"])
					for (i in seq_along(o)[-1]) {
						w <- o[which(keep[o[seq_len(i - 1)]])]
						int <- intersect(range1[o[i],],
							range1[w[res2[o[i], "index1"] == res2[w, "index1"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
						int <- intersect(range1[o[i],],
							range2[w[res2[o[i], "index1"] == res2[w, "index2"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
						int <- intersect(range2[o[i],],
							range1[w[res2[o[i], "index2"] == res2[w, "index1"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
						int <- intersect(range2[o[i],],
							range2[w[res2[o[i], "index2"] == res2[w, "index2"]],])
						if (length(int) > 0) {
							keep[o[i]] <- FALSE
							next
						}
					}
					res2 <- res2[keep,]
					rownames(res2) <- NULL
				}
			}
		} else {
			res2 <- data.frame(syn[[2, 1]][w, 1:8])
		}
		dbDisconnect(dbConn)
		
		if (type == 2L) { # interspersed
			result <- res2
		} else { # both
			result <- list(result, res2)
		}
	}
	
	return(result)
}
