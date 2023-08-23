#' Sequential Coordinate-wise Algorithm for the Non-negative Least Squares
#' Problem
#' 
#' Consider the linear system \eqn{\bold{A} x = b} where \eqn{\bold{A} \in
#' R\textsuperscript{m x n}}, \eqn{x \in R\textsuperscript{n}}, and \eqn{b \in
#' R\textsuperscript{m}}.  The technique of least squares proposes to compute
#' \eqn{x} so that the sum of squared residuals is minimized.  \code{NNLS}
#' solves the least squares problem \eqn{\min{||\bold{A} x =
#' b||\textsuperscript{2}}} subject to the constraint \eqn{x \ge 0}.  This
#' implementation of the Sequential Coordinate-wise Algorithm uses a sparse
#' input matrix \eqn{\bold{A}}, which makes it efficient for large sparse
#' problems.
#' 
#' The input \eqn{b} can be either a matrix or a vector of numerics.  If it is
#' a matrix then it is assumed that each column contains a set of observations,
#' and the output \eqn{x} will have the same number of columns.  This allows
#' multiple NNLS problems using the same \eqn{\bold{A}} matrix to be solved
#' simultaneously, and greatly accelerates computation relative to solving each
#' sequentially.
#' 
#' @name NNLS
#' @param A List representing the sparse matrix with integer components i and
#' j, numeric component x.  The fourth component, dimnames, is a list of two
#' components that contains the names for every row (component 1) and column
#' (component 2).
#' @param b Numeric matrix for the set of observed values.  (See details
#' section below.)
#' @param precision The desired accuracy.
#' @param processors The number of processors to use, or \code{NULL} to
#' automatically detect and use all available processors.
#' @param verbose Logical indicating whether to display progress.
#' @return A list of two components: \item{x}{ The matrix of non-negative
#' values that best explains the observed values given by \eqn{b}. }
#' \item{res}{ A matrix of residuals given by \eqn{\bold{A} x - b}. }
#' @seealso \code{\link{Array2Matrix}}, \code{\link{DesignArray}}
#' @references Franc, V., et al. (2005). Sequential coordinate-wise algorithm
#' for the non-negative least squares problem.  Computer Analysis of Images and
#' Patterns, 407-414.
#' @examples
#' 
#' # unconstrained least squares:
#' A <- matrix(c(1, -3, 2, -3, 10, -5, 2, -5, 6), ncol=3)
#' b <- matrix(c(27, -78, 64), ncol=1)
#' x <- solve(crossprod(A), crossprod(A, b))
#' 
#' # Non-negative least squares:
#' w <- which(A > 0, arr.ind=TRUE)
#' A <- list(i=w[,"row"], j=w[,"col"], x=A[w],
#'           dimnames=list(1:dim(A)[1], 1:dim(A)[2]))
#' x_nonneg <- NNLS(A, b)
#' 
#' # compare the unconstrained and constrained solutions:
#' cbind(x, x_nonneg$x)
#' 
#' # the input value "b" can also be a matrix:
#' b2 <- matrix(b, nrow=length(b), ncol=2) # repeat b in two columns
#' x_nonneg <- NNLS(A, b2) # solution is repeated in two output columns
#' 
#' @export NNLS
NNLS <- function(A,
	b,
	precision=sqrt(.Machine$double.eps),
	processors=1,
	verbose=TRUE) {
	
	# error checking:
	if (length(A) != 4)
		stop("A must have four components: A$i, A$j, A$x, and A$dimnames.")
	if (!is.integer(A$i))
		stop("Rows (i) must be a vector of integers.")
	if (!is.integer(A$j))
		stop("Columns (j) must be a vector of integers.")
	if (!is.double(A$x))
		stop("Values (x) must be a vector of doubles.")
	if (length(A$i) != length(A$j))
		stop("The length of columns (j) and rows (i) must be equal.")
	if (length(A$i) != length(A$x))
		stop("The length of rows (i) and values (x) must be equal.")
	if (max(A$j) > length(A$dimnames[[2]]))
		stop("More columns than column names.")
	if (max(A$i) > length(A$dimnames[[1]]))
		stop("More rows than row names.")
	if (!is.numeric(b))
		stop("b must be a numeric vector or matrix.")
	if (!((length(b) %% length(A$dimnames[[1]]))==0))
		stop("The length of b must be a multiple of the number of rows in A.")
	if (is(b, "matrix")) {
		if (nrow(b) != length(A$dimnames[[1]]))
			stop("The number of rows in b must equal the number of rows in A.")
	} else {
		b <- matrix(b, ncol=ncol(x))
	}
	if (!is.numeric(precision))
		stop("precision must be a numeric.")
	if (precision <= 0)
		stop("precision must be a positive number.")
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
	
	if (verbose) {
		time.1 <- Sys.time()
		pBar <- txtProgressBar(max=100, style=ifelse(interactive(), 3, 1))
	} else {
		pBar <- NULL
	}
	
	o <- order(A$i)
	x <- .Call("NNLS",
		A$i[o],
		A$j[o],
		A$x[o],
		length(A$dimnames[[1]]),
		length(A$dimnames[[2]]),
		b,
		precision,
		verbose,
		pBar,
		processors,
		PACKAGE="DECIPHER")
	
	res <- b - .Call("sparseMult",
		A$i,
		A$j,
		A$x,
		length(A$dimnames[[1]]),
		length(A$dimnames[[2]]),
		x,
		PACKAGE="DECIPHER")
	
	rownames(x) <- A$dimnames[[2]]
	colnames(x) <- colnames(b)
	rownames(res) <- A$dimnames[[1]]
	
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
	
	return(list(x=x, residuals=res))
}
