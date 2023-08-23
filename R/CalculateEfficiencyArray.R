#' Predict the Hybridization Efficiency of Probe/Target Sequence Pairs
#' 
#' Calculates the Gibbs free energy and hybridization efficiency of
#' probe/target pairs at varying concentrations of the denaturant formamide.
#' 
#' This function calculates the free energy and hybridization efficiency (HE)
#' for a given formamide concentration ([FA]) using the linear free energy
#' model given by: \deqn{HE = Po*exp[-(dG_0 + m*FA)/RT]/(1+Po*exp[-(dG_0 +
#' m*FA)/RT])}
#' 
#' The \code{probe} and \code{target} input sequences must be aligned in pairs,
#' such that the first probe is aligned to the first target, second-to-second,
#' and so on.  Ambiguity codes in the \code{IUPAC_CODE_MAP} are accepted in
#' probe and target sequences.  Any ambiguities will default to perfect match
#' pairings by inheriting the nucleotide in the same position on the opposite
#' sequence whenever possible.  If the ambiguity results in a mismatch then
#' ``T'', ``G'', ``C'', and ``A'' are substituted, in that order.  For example,
#' if a probe nucleotide is ``S'' (``C'' or ``G'') then it will be considered a
#' ``C'' if the target nucleotide in the same position is a ``C'', otherwise
#' the ambiguity will be interpreted as a ``G''.
#' 
#' If \code{deltaGrules} is NULL then the rules defined in
#' \code{data(deltaGrules)} will be used.  Note that \code{deltaGrules} of the
#' same format may be customized for any application and specified as an input.
#' 
#' @name CalculateEfficiencyArray
#' @param probe A \code{DNAStringSet} object or character vector with
#' pairwise-aligned probe sequences in 5' to 3' orientation.
#' @param target A \code{DNAStringSet} object or character vector with
#' pairwise-aligned target sequences in 5' to 3' orientation.
#' @param FA A vector of one or more formamide concentrations (as percent v/v).
#' @param dGini The initiation free energy.  The default is 1.96 [kcal/mol].
#' @param Po The effective probe concentration.
#' @param m The m-value defining the linear relationship of denaturation in the
#' presence of formamide.
#' @param temp Equilibrium temperature in degrees Celsius.
#' @param deltaGrules Free energy rules for all possible base pairings in
#' quadruplets.  If NULL, defaults to the parameters obtained using NimbleGen
#' microarrays and a Linear Free Energy Model developed by Yilmaz \emph{et al}.
#' @return A \code{matrix} with the predicted Gibbs free energy (dG) and
#' hybridization efficiency (HE) at each concentration of formamide ([FA]).
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{deltaGrules}}
#' @references Yilmaz LS, Loy A, Wright ES, Wagner M, Noguera DR (2012)
#' Modeling Formamide Denaturation of Probe-Target Hybrids for Improved
#' Microarray Probe Design in Microbial Diagnostics. PLoS ONE 7(8): e43862.
#' doi:10.1371/journal.pone.0043862.
#' @examples
#' 
#' probes <- c("AAAAACGGGGAGCGGGGGGATACTG", "AAAAACTCAACCCGAGGAGCGGGGG")
#' targets <- c("CAACCCGGGGAGCGGGGGGATACTG", "TCGGGCTCAACCCGAGGAGCGGGGG")
#' result <- CalculateEfficiencyArray(probes, targets, FA=0:40)
#' dG0 <- result[, "dG_0"]
#' HE0 <- result[, "HybEff_0"]
#' plot(result[1, 1:40], xlab="[FA]", ylab="HE", main="Probe/Target # 1", type="l")
#' 
#' @export CalculateEfficiencyArray
CalculateEfficiencyArray <- function(probe,
	target,
	FA=0,
	dGini=1.96,
	Po=10^-2.0021,
	m=0.1731,
	temp=42,
	deltaGrules=NULL) {
	
	# error checking
	if (is(probe, "DNAStringSet"))
		probe <- strsplit(toString(probe), ", ", fixed=TRUE)[[1]]
	if (is(target, "DNAStringSet"))
		target <- strsplit(toString(target), ", ", fixed=TRUE)[[1]]
	if (!is.character(probe))
		stop("probe must be a character vector.")
	if (!is.character(target))
		stop("target must be a character vector.")
	if (any(nchar(target) != nchar(probe)))
		stop("probe and target must be aligned (equal length).")
	if (!is.numeric(Po))
		stop("Po must be a numeric.")
	if (!(Po > 0))
		stop("Po must be greater than zero.")
	if (!is.numeric(m))
		stop("m must be a numeric.")
	if (!(m > 0))
		stop("m must be greater than zero.")
	if (!is.numeric(dGini))
		stop("dGini must be a numeric.")
	if (!is.numeric(FA))
		stop("FA must be a numeric.")
	if (any(FA < 0))
		stop("FA must be greater than or equal to zero.")
	if (!is.numeric(temp))
		stop("temp must be a numeric.")
	if (temp < -273)
		stop("temp must be greater than or equal to absolute zero.")
	
	if (is.null(deltaGrules)) {
		data("deltaGrules", envir=environment(), package="DECIPHER")
	} else {
		if (!is.numeric(deltaGrules))
			stop("deltaGrules must be numeric.")
		if (length(deltaGrules) != 390625L)
			stop("deltaGrules must be of dimensions 5 x 5 x 5 x 5 x 5 x 5 x 5 x 5.")
	}
	
	l <- length(probe)
	if (l==0)
		stop("No probe specified.")
	if (l!=length(target))
		stop("probe is not the same length as target.")
	
	dG <- dGini + .Call("calculateDeltaG", probe, target, deltaGrules, PACKAGE="DECIPHER")
	
	RT <- .0019871*(273.15 + temp) # [kcal/mol]
	deltaG <- matrix(0, nrow=l, ncol=length(FA), dimnames=list(1:l, paste("dG", FA, sep="_")))
	eff <- matrix(0, nrow=l, ncol=length(FA), dimnames=list(1:l, paste("HybEff", FA, sep="_")))
	for (i in 1:length(FA)) {
		deltaG[, i] <- dG + m*FA[i]
		eff[, i] <- Po*exp(-deltaG[, i]/RT)/(1 + Po*exp(-deltaG[, i]/RT))
	}
	
	return(cbind(eff, deltaG))
}
