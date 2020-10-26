`[.Genes` <- function(x, i, j, ...) {
	if (class(x)[1] != "Genes")
		stop("x must be an object of class 'Genes'.")
	a <- attributes(x)
	x <- unclass(x)
	x <- x[i, j, drop=FALSE]
	if (missing(j)) {
		class(x) <- "Genes"
		attr(x, "widths") <- a$widths
		attr(x, "geneticCode") <- a$geneticCode
	}
	
	return(x)
}

print.Genes <- function(x, ...) {
	if (!is(x, "Genes"))
		stop("x must be an object of class 'Genes'.")
	N <- nrow(x)
	if (N > 0) {
		w <- which(x[, "Gene"]==1)
		if (length(w) > 0) {
			l <- x[w, "End"] - x[w, "Begin"] + 1L
			n <- length(w)
			s <- suppressWarnings(cor(x[, -c(1:6, 17:18)],
				method="spearman"))
			s <- 100*s[upper.tri(s)]
			cat("Genes object specifying",
				formatC(n,
					big.mark=",",
					format="fg"))
			if (n==1) {
				cat(" gene of",
					formatC(l,
						big.mark=",",
						format="fg"),
					"nucleotides.")
			} else {
				cat(" genes from",
					formatC(min(l),
						big.mark=",",
						format="fg"),
					"to",
					formatC(max(l),
						big.mark=",",
						format="fg"),
					"nucleotides.")
			}
			u <- unique(x[w, "StartScore"])
			u <- u[u != 0]
			cat("\nDensity: ",
				round(100*sum(l)/sum(attr(x, "widths")), 1),
				"%; Initation codons: ",
				length(u),
				sep="")
			if (any(is.na(s))) {
				cat(".\n\n")
			} else {
				cat("; Score correlations: ",
					round(min(s)),
					"% to ",
					round(max(s)),
					"%.\n\n",
					sep="")
			}
		} else {
			l <- x[, "End"] - x[, "Begin"] + 1L
			n <- nrow(x)
			cat("Genes object specifying",
				formatC(n,
					big.mark=",",
					format="fg"))
			if (n==1) {
				cat(" open reading frame of",
					formatC(l,
						big.mark=",",
						format="fg"),
					"nucleotides.\n\n")
			} else {
				cat(" open reading frames from",
					formatC(min(l),
						big.mark=",",
						format="fg"),
					"to",
					formatC(max(l),
						big.mark=",",
						format="fg"),
					"nucleotides.\n\n")
			}
		}
		
		print(cbind(as.data.frame(round(head(x[, 1:5]), 2)),
			`...`="...",
			as.data.frame(round(head(x[, 17:18]), 2))))
		if (N==7) {
			cat("... with 1 more row.\n")
		} else if (N > 6) {
			cat("... with ",
				formatC(N - 6,
					big.mark=",",
					format="fg"),
				" more rows.\n",
				sep="")
		}
	} else {
		cat("Genes object specifying no genes and no open reading frames.\n") 
	}
}

plot.Genes <- function(x,
	xlim=c(1, 1e4),
	ylim=c(-100, 500),
	interact=FALSE,
	...) {
	
	ws <- attr(x, "widths")
	ws <- cumsum(ws) - ws
	x[, "Begin"] <- x[, "Begin"] + ws[x[, "Index"]]
	x[, "End"] <- x[, "End"] + ws[x[, "Index"]]
	ws <- ws[-1L]
	
	pos <- which(x[, "Strand"]==0)
	neg <- which(x[, "Strand"]==1)
	genes <- which(x[, "Gene"]==1)
	pos_genes <- which(x[, "Gene"]==1 & x[, "Strand"]==0)
	neg_genes <- which(x[, "Gene"]==1 & x[, "Strand"]==1)
	
	start <- xlim[1L]
	offset <- diff(xlim) + 1
	
	while(start > 0 &&
		start < max(x[, "End"] + offset)) {
		plot(NA,
			xlim=c(start, start + offset - 1L),
			ylim=ylim,
			xlab="Cumulative genome position",
			ylab="Total score")
		segments(x[pos, "Begin"],
			x[pos, "TotalScore"],
			x[pos, "End"],
			x[pos, "TotalScore"],
			col="darkred")
		segments(x[neg, "Begin"],
			x[neg, "TotalScore"],
			x[neg, "End"],
			x[neg, "TotalScore"],
			col="darkblue")
		segments(x[genes, "Begin"],
			x[genes, "TotalScore"],
			x[genes, "End"],
			x[genes, "TotalScore"],
			col="green4",
			lwd=2)
		abline(v=x[pos_genes, "Begin"])
		abline(v=x[pos_genes, "End"],
			lty=5)
		abline(v=x[neg_genes, "Begin"],
			lty=5)
		abline(v=x[neg_genes, "End"])
		abline(h=0, v=ws, lwd=2, lty=3)
		
		if (!interact)
			break
		
		l <- locator(n=1)
		if (length(l$x) != 1)
			break
		if (l$x > start + offset/2) {
			start <- start + offset - 1L
		} else {
			start <- start - offset + 1L
		}
	}
}
