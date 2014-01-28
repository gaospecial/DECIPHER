.align <- function(t, lp, ls) {
	end.p <- t[1]
	end.s <- t[2]
	start.p <- t[3]
	start.s <- t[4]
	count <- t[5] - 4
	t <- t[-1:-5]
	
	p.inserts <- character()
	s.inserts <- character()
	p.ats <- integer()
	s.ats <- integer()
	
	if (start.s < start.p) {
		s.inserts <- paste(rep("-", start.p - 1), collapse="")
		s.ats <- 1
	} else if (start.p < start.s) {
		p.inserts <- paste(rep("-", start.s - 1), collapse="")
		p.ats <- 1
	}
	
	i <- start.p
	j <- start.s
	
	if (count < length(t)) {
		k <- count + 1
		while (k <= length(t)) {
			l <- length(.Call("multiMatch", t, t[k], as.integer(k), PACKAGE="DECIPHER"))
			if (t[k] == 0) {
				i <- i + l
				j <- j + l
			} else if (t[k] > 0) {
				p.inserts <- c(p.inserts, paste(rep("-", l), collapse=""))
				p.ats <- c(p.ats, i)
				j <- j + l
			} else {
				s.inserts <- c(s.inserts, paste(rep("-", l), collapse=""))
				s.ats <- c(s.ats, j)
				i <- i + l
			}
			k <- k + l
		}
	}
	
	if (ls - end.s > 0) {
		p.inserts <- c(p.inserts, paste(rep("-", ls - end.s), collapse=""))
		p.ats <- c(p.ats, lp + 1)
	}
	
	if (lp - end.p > 0) {
		s.inserts <- c(s.inserts, paste(rep("-", lp - end.p), collapse=""))
		s.ats <- c(s.ats, ls + 1)
	}
	
	return(list(p.ats, p.inserts, s.ats, s.inserts))
}

AlignProfiles <- function(pattern,
	subject,
	p.weight=1,
	s.weight=1,
	perfectMatch=6,
	misMatch=0,
	gapOpening=-9,
	gapExtension=-3,
	terminalGap=-2,
	restrict=-1000,
	anchor=0.7,
	substitutionMatrix="BLOSUM62",
	processors=NULL) {
	
	# error checking
	if (!is(pattern, "DNAStringSet") && !is(pattern, "RNAStringSet") && !is(pattern, "AAStringSet"))
		stop("pattern must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	if (length(pattern) < 1)
		stop("At least one sequence is required in the pattern.")
	if (!is(subject, "DNAStringSet") && !is(subject, "RNAStringSet") && !is(subject, "AAStringSet"))
		stop("subject must be an AAStringSet, DNAStringSet, or RNAStringSet.")
	if (typeof(pattern) != typeof(subject))
		stop("pattern and subject must be of the same type.")
	if (length(subject) < 1)
		stop("At least one sequence is required in the subject.")
	if (!is.numeric(p.weight))
		stop("p.weight must be a numeric.")
	if (length(p.weight)!=1 && length(p.weight)!=length(pattern))
		stop("Length of p.weight must equal one or the length of the pattern.")
	if (length(p.weight)==1)
		p.weight <- rep(1, length(pattern))
	if (!is.numeric(s.weight))
		stop("s.weight must be a numeric.")
	if (length(s.weight)!=1 && length(s.weight)!=length(subject))
		stop("Length of s.weight must equal one or the length of the subject.")
	if (length(s.weight)==1)
		s.weight <- rep(1, length(subject))
	if (!is.numeric(perfectMatch))
		stop("perfectMatch must be a numeric.")
	if (!is.numeric(misMatch))
		stop("misMatch must be a numeric.")
	if (!is.numeric(gapOpening))
		stop("gapOpening must be a numeric.")
	if (!is.numeric(gapExtension))
		stop("gapExtension must be a numeric.")
	if (!all(is.numeric(terminalGap)))
		stop("terminalGap must be a numeric.")
	if (length(terminalGap) > 2 || length(terminalGap) < 1)
		stop("terminalGap must be of length 1 or 2.")
	if (any(is.infinite(terminalGap)))
		stop("terminalGap must be finite.")
	if (length(terminalGap)==1)
		terminalGap[2] <- terminalGap[1]
	if (!all(is.numeric(restrict)))
		stop("restrict must be a numeric.")
	if (restrict >= 0)
		stop("restrict must be less than zero.")
	if (!is.numeric(anchor) && !is.na(anchor))
		stop("anchor must be a numeric.")
	if (is.numeric(anchor) && anchor <= 0)
		stop("anchor must be greater than zero.")
	if (is.numeric(anchor) && anchor > 1)
		stop("anchor must be less than or equal to one.")
	if (!(substitutionMatrix %in% c("BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM100",
		"PAM30", "PAM40", "PAM70", "PAM120", "PAM250")))
		stop("Invalid substitutionMatrix.")
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
	if (length(pattern) > 2147483647)
		stop(paste("Length of pattern (",
			length(pattern),
			") longer than the maximum allowable length (2,147,483,647).",
			sep=""))
	if (length(subject) > 2147483647)
		stop(paste("Length of subject (",
			length(subject),
			") longer than the maximum allowable length (2,147,483,647).",
			sep=""))
	
	if (is(pattern, "AAStringSet")) {
		AAs <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
			"L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*")
		subMatrix <- eval(parse(text=data(list=substitutionMatrix)))
		subMatrix <- subMatrix[AAs, AAs]
	}
	
	f <- function(pattern, subject, tGaps=terminalGap) {
		if (length(pattern)*length(subject) > 2147483647) # maximum when indexing by signed integer
			stop(paste("Alignment larger (",
				length(pattern)*length(subject),
				") than the maximum allowable size (2,147,483,647).",
				sep=""))
		if (is(pattern, "AAStringSet")) {
			p.profile <- .Call("consensusProfileAA",
				pattern,
				p.weight,
				PACKAGE="DECIPHER")
			s.profile <- .Call("consensusProfileAA",
				subject,
				s.weight,
				PACKAGE="DECIPHER")
			
			t <- .Call("alignProfilesAA",
				p.profile,
				s.profile,
				subMatrix,
				perfectMatch,
				misMatch,
				gapOpening,
				gapExtension,
				tGaps[1],
				tGaps[2],
				restrict,
				processors,
				PACKAGE="DECIPHER")
		} else { # DNAStringSet or RNAStringSet
			p.profile <- .Call("consensusProfile",
				pattern,
				p.weight,
				PACKAGE="DECIPHER")
			s.profile <- .Call("consensusProfile",
				subject,
				s.weight,
				PACKAGE="DECIPHER")
			
			t <- .Call("alignProfiles",
				p.profile,
				s.profile,
				perfectMatch,
				misMatch,
				gapOpening,
				gapExtension,
				tGaps[1],
				tGaps[2],
				restrict,
				processors,
				PACKAGE="DECIPHER")
		}
		lp <- dim(p.profile)[2]
		ls <- dim(s.profile)[2]
		
		.align(t, lp, ls)
	}
	
	if (is.na(anchor)) { # don't use anchors
		inserts <- f(pattern, subject)
	} else { # use anchors
		w.p <- width(pattern[1])
		w.s <- width(subject[1])
		
		if (is(pattern, "AAStringSet")) {
				wordSize <- 7
		} else {
				wordSize <- 15
		}
		l <- min(length(pattern), length(subject))
		o.p <- order(p.weight, decreasing=TRUE)
		o.s <- order(s.weight, decreasing=TRUE)
		if (is(pattern, "AAStringSet")) {
			num.p <- .Call("enumerateGappedSequenceAA",
				pattern[o.p[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
			num.s <- .Call("enumerateGappedSequenceAA",
				subject[o.s[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
		} else {
			num.p <- .Call("enumerateGappedSequence",
				pattern[o.p[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
			num.s <- .Call("enumerateGappedSequence",
				subject[o.s[1:l]],
				wordSize,
				PACKAGE="DECIPHER")
		}
		
		anchors <- .Call("matchRanges",
			num.p,
			num.s,
			wordSize,
			w.p,
			anchor,
			PACKAGE="DECIPHER")
		
		numAnchors <- dim(anchors)[2]
		if (numAnchors==0) {
			inserts <- f(pattern, subject)
		} else {
			if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				1L, anchors[1, 1],
				1L, anchors[3, 1],
				PACKAGE="DECIPHER")) {
				inserts <- f(subseq(pattern, 1L, anchors[1, 1]),
					subseq(subject, 1L, anchors[3, 1]),
					tGaps=c(terminalGap[1], -1000))
			} else {
				inserts <- list(integer(), character(),
					integer(), character())
			}
			
			n <- 2L
			while (n <= numAnchors) { # align regions between anchors
				if (!.Call("firstSeqsEqual",
					pattern,
					subject,
					anchors[2, n - 1], anchors[1, n],
					anchors[4, n - 1], anchors[3, n],
					PACKAGE="DECIPHER")) {
					temp <- f(subseq(pattern, anchors[2, n - 1], anchors[1, n]),
						subseq(subject, anchors[4, n - 1], anchors[3, n]),
						tGaps=c(-1000, -1000))
					inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, n - 1] - 1L)
					inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, n - 1] - 1L)
					inserts[[2]] <- c(inserts[[2]], temp[[2]])
					inserts[[4]] <- c(inserts[[4]], temp[[4]])
				}
				n <- n + 1L
			}
			
			n <- 1L
			while (n <= numAnchors) { # align anchor regions
				if (!.Call("firstSeqsEqual",
					pattern,
					subject,
					anchors[1, n], anchors[2, n],
					anchors[3, n], anchors[4, n],
					PACKAGE="DECIPHER")) {
					temp <- f(subseq(pattern, anchors[1, n], anchors[2, n]),
						subseq(subject, anchors[3, n], anchors[4, n]),
						tGaps=c(-1000, -1000))
					inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[1, n] - 1L)
					inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[3, n] - 1L)
					inserts[[2]] <- c(inserts[[2]], temp[[2]])
					inserts[[4]] <- c(inserts[[4]], temp[[4]])
				}
				n <- n + 1L
			}
			
			end.p <- anchors[2, numAnchors] == w.p
			end.s <- anchors[4, numAnchors] == w.s
			if (end.p && !end.s) { # need to add gaps after pattern
				inserts[[1]] <- c(inserts[[1]], w.p + 1L)
				inserts[[2]] <- c(inserts[[2]],
					paste(rep("-", w.s - anchors[4, numAnchors]), collapse=""))
			} else if (end.s && !end.p) { # need to add gaps after subject
				inserts[[3]] <- c(inserts[[3]], w.s + 1L)
				inserts[[4]] <- c(inserts[[4]],
					paste(rep("-", w.p - anchors[2, numAnchors]), collapse=""))
			} else if (!.Call("firstSeqsEqual",
				pattern,
				subject,
				anchors[2, numAnchors], w.p,
				anchors[4, numAnchors], w.s,
				PACKAGE="DECIPHER")) { # need to align
				temp <- f(subseq(pattern, anchors[2, numAnchors], w.p),
					subseq(subject, anchors[4, numAnchors], w.s),
					tGaps=c(-1000, terminalGap[2]))
				inserts[[1]] <- c(inserts[[1]], temp[[1]] + anchors[2, numAnchors] - 1L)
				inserts[[3]] <- c(inserts[[3]], temp[[3]] + anchors[4, numAnchors] - 1L)
				inserts[[2]] <- c(inserts[[2]], temp[[2]])
				inserts[[4]] <- c(inserts[[4]], temp[[4]])
			} # else don't do anything
		}
	}
	
	if (length(inserts[[1]]) > 0)
		pattern <- replaceAt(pattern, at=inserts[1], value=inserts[[2]])
	if (length(inserts[[3]]) > 0)
		subject <- replaceAt(subject, at=inserts[3], value=inserts[[4]])
	
	return(append(pattern, subject))
}