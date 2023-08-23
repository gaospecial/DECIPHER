#' Create a Matrix Representation of a Microarray
#' 
#' Converts the output of \code{DesignArray} into the sparse matrix format used
#' by \code{NNLS}.
#' 
#' A microarray can be represented by a matrix of hybridization efficiencies,
#' where the rows represent each of the probes and the columns represent each
#' the possible templates.  This matrix is sparse since microarray probes are
#' designed to only target a small subset of the possible templates.
#' 
#' @name Array2Matrix
#' @param probes A set of microarray probes in the format output by
#' \code{DesignArray}.
#' @param verbose Logical indicating whether to display progress.
#' @return A list specifying the hybridization efficiency of each probe to its
#' potential templates.  \item{i}{ Element's row index in the sparse matrix. }
#' \item{j}{ Element's column index in the sparse matrix. } \item{x}{ Non-zero
#' elements' values representing hybridization efficiencies. } \item{dimnames}{
#' A list of two components: the names of each probe, and the names of each
#' template. }
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{DesignArray}}, \code{\link{NNLS}}
#' @references ES Wright et al. (2013) Identification of Bacterial and Archaeal
#' Communities From Source to Tap. Water Research Foundation, Denver, CO.
#' 
#' DR Noguera, et al. (2014). Mathematical tools to optimize the design of
#' oligonucleotide probes and primers. Applied Microbiology and Biotechnology.
#' doi:10.1007/s00253-014-6165-x.
#' @examples
#' 
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' names(dna) <- 1:length(dna)
#' probes <- DesignArray(dna)
#' A <- Array2Matrix(probes)
#' 
#' @export Array2Matrix
Array2Matrix <- function(probes,
	verbose=TRUE) {
	
	if (!is.logical(verbose))
		stop("verbose must be a logical.")
	if (verbose) {
		pBar <- txtProgressBar(style=ifelse(interactive(), 3, 1))
		time.1 <- Sys.time()
	}
	
	MMs <- strsplit(probes$mismatch, ", ", fixed=TRUE)
	
	hyb_effs <- list()
	l <- length(MMs)
	for (i in 1:l) {
		if (length(MMs[[i]])==0)
			next
		
		MMs[[i]] <- unlist(strsplit(MMs[[i]], " ", fixed=TRUE))
		
		num <- length(MMs[[i]])
		
		hyb_effs[[i]] <- MMs[[i]][2*(1:(num/2))]
		hyb_effs[[i]] <- .Call("replaceChar",
			hyb_effs[[i]],
			"%",
			"",
			PACKAGE="DECIPHER")
		hyb_effs[[i]] <- .Call("replaceChar",
			hyb_effs[[i]],
			"(",
			"",
			PACKAGE="DECIPHER")
		hyb_effs[[i]] <- as.numeric(.Call("replaceChar",
			hyb_effs[[i]],
			")",
			"",
			PACKAGE="DECIPHER"))
		
		MMs[[i]] <- MMs[[i]][2*(1:(num/2)) - 1]
		
		if (verbose)
			setTxtProgressBar(pBar, i/l)
	}
	
	index <- lapply(1:l, function(l) {
		rep(l, length(MMs[[l]]))
	})
	
	p_names <- unique(probes$name)
	MMs <- match(unlist(MMs), p_names)
	
	A <- list(i=c(1:l,unlist(index)),
		j=c(match(probes$name, p_names),MMs),
		x=c(as.numeric(probes$hyb_eff),unlist(hyb_effs))/100,
		dimnames=list(1:l, p_names))
	
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
	
	return(A)
}
