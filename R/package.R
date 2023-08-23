

#' Tools for curating, analyzing, and manipulating biological sequences
#' 
#' DECIPHER is a software toolset that can be used for deciphering and managing
#' biological sequences efficiently using the R statistical programming
#' language.  The program is designed to be used with non-destructive workflows
#' for importing, maintaining, analyzing, manipulating, and exporting a massive
#' amount of sequences.
#' 
#' \tabular{ll}{ Package: \tab DECIPHER\cr Type: \tab Package\cr Depends: \tab
#' R (>= 3.5.0), Biostrings (>= 2.59.1), RSQLite (>= 1.1), stats, parallel\cr
#' Imports: \tab methods, DBI, S4Vectors, IRanges, XVector\cr LinkingTo: \tab
#' Biostrings, RSQLite, S4Vectors, IRanges, XVector\cr License: \tab GPL-3\cr
#' LazyLoad: \tab yes\cr }
#' 
#' Index: \preformatted{ AA_REDUCED Reduced amino acid alphabets Add2DB Add
#' Data to a Database AdjustAlignment Improve An Existing Alignment By
#' Adjusting Gap Placements AlignDB Align Two Sets of Aligned Sequences in a
#' Sequence Database AlignProfiles Align Two Sets of Aligned Sequences
#' AlignSeqs Align a Set of Unaligned Sequences AlignSynteny Pairwise Aligns
#' Syntenic Blocks AlignTranslation Align Sequences By Their Amino Acid
#' Translation AmplifyDNA Simulate Amplification of DNA by PCR Array2Matrix
#' Create a Matrix Representation of a Microarray BrowseDB View a Database
#' Table in a Web Browser BrowseSeqs View Sequences in a Web Browser
#' CalculateEfficiencyArray Predict the Hybridization Efficiency of
#' Probe/Target Sequence Pairs CalculateEfficiencyFISH Predict Thermodynamic
#' Parameters of Probe/Target Sequence Pairs CalculateEfficiencyPCR Predict
#' Amplification Efficiency of Primer Sequences Clusterize Cluster Sequences By
#' Distance Codec Compression/Decompression of Character Vectors
#' ConsensusSequence Create a Consensus Sequence Cophenetic Compute cophenetic
#' distances on dendrogram objects CorrectFrameshifts Corrects Frameshift
#' Errors In Protein Coding Sequences CreateChimeras Create Artificial Chimeras
#' DB2Seqs Export Database Sequences to a FASTA or FASTQ File deltaGrules Free
#' Energy of Hybridization of Probe/Target Quadruplets deltaHrules Change in
#' Enthalpy of Hybridization of DNA/DNA Quadruplets in Solution deltaHrulesRNA
#' Change in Enthalpy of Hybridization of RNA/RNA Quadruplets in Solution
#' deltaSrules Change in Entropy of Hybridization of DNA/DNA Quadruplets in
#' Solution deltaSrulesRNA Change in Entropy of Hybridization of RNA/RNA
#' Quadruplets in Solution DesignArray Design a Set of DNA Microarray Probes
#' for Detecting Sequences DesignPrimers Design Primers Targeting a Specific
#' Group of Sequences DesignProbes Design FISH Probes Targeting a Specific
#' Group of Sequences DesignSignatures Design PCR Primers for Amplifying
#' Group-Specific Signatures DetectRepeats Detect Repeats in a Sequence
#' DigestDNA Simulate Restriction Digestion of DNA Disambiguate Expand
#' Ambiguities into All Permutations of a DNAStringSet DistanceMatrix Calculate
#' the Distance Between Sequences ExtractGenes Extract Predicted Genes from a
#' Genome FindChimeras Find Chimeras in a Sequence Database FindGenes Find
#' Genes in a Genome FindNonCoding Find Non-Coding RNAs in a Genome FindSynteny
#' Finds Synteny in a Sequence Database FormGroups Forms Groups By Rank
#' Genes-class Genes objects and accessors HEC_MI Mutual Information for
#' Protein Secondary Structure Prediction IdConsensus Create Consensus
#' Sequences by Groups IdentifyByRank Identify By Taxonomic Rank IdLengths
#' Determine the Number of Characters in Each Sequence of Each Sequence IdTaxa
#' Assign Sequences a Taxonomic Classification LearnNonCoding Learn a
#' Non-Coding RNA Model LearnTaxa Train a Classifier for Assigning Taxonomy
#' MapCharacters Map Changes in Ancestral Character States MaskAlignment Mask
#' Highly Variable Regions of An Alignment MeltDNA Simulate Melting of DNA MIQS
#' MIQS Amino Acid Substitution Matrix MODELS Available Models of Sequence
#' Evolution NNLS Sequential Coordinate-wise Algorithm for the Non-negative
#' Least Squares Problem NonCoding NonCoding Models for Common Non-Coding RNA
#' Families NonCoding-class NonCoding Objects and Methods OrientNucleotides
#' Orient Nucleotide Sequences PFASUM PFASUM Amino Acid Substitution Matrices
#' PredictDBN Predict RNA Secondary Structure in Dot-Bracket Notation
#' PredictHEC Predict Protein Secondary Structure as Helix, Beta-Sheet, or Coil
#' Read Dendrogram Read a Dendrogram from a Newick Formatted File RemoveGaps
#' Remove Gap Characters in Sequences RESTRICTION_ENZYMES Common Restriction
#' Enzyme's Cut Sites ScoreAlignment Score a Multiple Sequence Alignment
#' SearchDB Obtain Specific Sequences from a Database Seqs2DB Add Sequences
#' from Text File to Database StaggerAlignment Produce a Staggered Alignment
#' Synteny-class Synteny blocks and hits Taxa-class Taxa training and testing
#' objects TerminalChar Determine the Number of Terminal Characters TileSeqs
#' Form a Set of Tiles for Each Group of Sequences TrainingSet_16S Training Set
#' for Classification of 16S rRNA Gene Sequences TreeLine Construct a
#' Phylogenetic Tree TrimDNA Trims DNA Sequences to the High Quality Region
#' Between Patterns WriteDendrogram Write a Dendrogram to Newick Format
#' WriteGenes Write Genes to a File }
#' 
#' @name DECIPHER-package
#' @aliases DECIPHER-package DECIPHER
#' @docType package
#' @author Erik Wright
#' 
#' Maintainer: Erik Wright <eswright@@pitt.edu>
#' @keywords package
NULL





#' Free Energy of Hybridization of Probe/Target Quadruplets on Microarrays
#' 
#' An 8D array with four adjacent base pairs of the probe and target sequences
#' at a time.  Each dimension has five elements defining the nucleotide at that
#' position ("A", "C", "G", "T", or "-").  The array contains the standard
#' Gibbs free energy change of probe binding (dG, [kcal/mol]) for every
#' quadruple base pairing.
#' 
#' The first four dimensions correspond to the four probe positions from 5' to
#' 3'.  The fifth to eighth dimensions correspond to the four positions from 5'
#' to 3' of the target sequence.
#' 
#' @name deltaGrules
#' @docType data
#' @format The format is: num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.141 0
#' 0 0 0 ...  - attr(*, "dimnames")=List of 8 ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...
#' @references Yilmaz LS, Loy A, Wright ES, Wagner M, Noguera DR (2012)
#' Modeling Formamide Denaturation of Probe-Target Hybrids for Improved
#' Microarray Probe Design in Microbial Diagnostics. PLoS ONE 7(8): e43862.
#' doi:10.1371/journal.pone.0043862.
#' @source Data obtained using NimbleGen microarrays and a Linear Free Energy
#' Model developed by Yilmaz \emph{et al}.
#' @keywords datasets
#' @examples
#' 
#' data(deltaGrules)
#' # dG of probe = AGCT / target = A-CT pairing
#' deltaGrules["A", "G", "C", "T", "A", "-", "C", "T"]
#' 
NULL





#' Pseudoenergy Parameters for RNA Quadruplets
#' 
#' An 8D array with four adjacent base pairs of the RNA duplex.  Each dimension
#' has five elements defining the nucleotide at that position ("A", "C", "G",
#' "U", or "-").  The array contains the pseudoenergy of duplex formation for
#' every quadruple base pairing.
#' 
#' The first four dimensions correspond to the four top strand positions from
#' 5' to 3'.  The fifth to eighth dimensions correspond to the four bottom
#' strand positions from 5' to 3'.
#' 
#' @name deltaGrulesRNA
#' @docType data
#' @format The format is: num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.776
#' -0.197 -0.197 -0.291 0 ...  - attr(*, "dimnames")=List of 8 ..$ : chr [1:5]
#' "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5]
#' "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5]
#' "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5]
#' "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U" ...
#' @source Psuedoenergy values of ungapped quadruplets are inferred from their
#' log-odds of being found in palindromes within hairpin regions across
#' thousands of non-coding RNA alignments.  Each value represents the log-odds
#' of \emph{in vivo} folding relative to chance.
#' @keywords datasets
#' @examples
#' 
#' data(deltaGrulesRNA)
#' # dG of the duplex AGCU / ACCU pairing (1 mismatch)
#' deltaGrulesRNA["A", "G", "C", "U", "A", "C", "C", "U"]
#' 
NULL





#' Change in Enthalpy of Hybridization of DNA/DNA Quadruplets in Solution
#' 
#' An 8D array with four adjacent base pairs of the DNA duplex.  Each dimension
#' has five elements defining the nucleotide at that position ("A", "C", "G",
#' "T", or "-").  The array contains the standard enthalpy change of probe
#' binding (dH, [kcal/mol]) for every quadruple base pairing.
#' 
#' The first four dimensions correspond to the four top strand positions from
#' 5' to 3'.  The fifth to eighth dimensions correspond to the four bottom
#' strand positions from 5' to 3'.
#' 
#' @name deltaHrules
#' @docType data
#' @format The format is: num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -7.97 0
#' 0 0 0 ...  - attr(*, "dimnames")=List of 8 ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...
#' @references SantaLucia, J., Jr., & Hicks, D. (2004) The Thermodynamics of
#' DNA Structural Motifs. Annual Review of Biophysics and Biomolecular
#' Structure, 33(1), 415-440. doi:10.1146/annurev.biophys.32.110601.141800.
#' @source Data from a variety of publications by SantaLucia \emph{et al}.
#' @keywords datasets
#' @examples
#' 
#' data(deltaHrules)
#' # dH of the duplex AGCT / A-CT pairing
#' deltaHrules["A", "G", "C", "T", "A", "-", "C", "T"]
#' 
NULL





#' Change in Enthalpy of Hybridization of RNA/RNA Quadruplets in Solution
#' 
#' An 8D array with four adjacent base pairs of the RNA duplex.  Each dimension
#' has five elements defining the nucleotide at that position ("A", "C", "G",
#' "U", or "-").  The array contains the standard enthalpy change of probe
#' binding (dH, [kcal/mol]) for every quadruple base pairing.
#' 
#' The first four dimensions correspond to the four top strand positions from
#' 5' to 3'.  The fifth to eighth dimensions correspond to the four bottom
#' strand positions from 5' to 3'.
#' 
#' @name deltaHrulesRNA
#' @docType data
#' @format The format is: num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -6.55 0
#' 0 0 0 ...  - attr(*, "dimnames")=List of 8 ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...
#' @references SantaLucia, J., Jr., & Hicks, D. (2004) The Thermodynamics of
#' DNA Structural Motifs. Annual Review of Biophysics and Biomolecular
#' Structure, 33(1), 415-440. doi:10.1146/annurev.biophys.32.110601.141800.
#' @source Data from a variety of publications by SantaLucia \emph{et al}.
#' @keywords datasets
#' @examples
#' 
#' data(deltaHrulesRNA)
#' # dH of the duplex AGCU / A-CU pairing
#' deltaHrulesRNA["A", "G", "C", "U", "A", "-", "C", "U"]
#' 
NULL





#' Change in Entropy of Hybridization of DNA/DNA Quadruplets in Solution
#' 
#' An 8D array with four adjacent base pairs of the DNA duplex.  Each dimension
#' has five elements defining the nucleotide at that position ("A", "C", "G",
#' "T", or "-").  The array contains the standard entropy change of probe
#' binding (dS, [kcal/mol]) for every quadruple base pairing.
#' 
#' The first four dimensions correspond to the four top strand positions from
#' 5' to 3'.  The fifth to eighth dimensions correspond to the four bottom
#' strand positions from 5' to 3'.
#' 
#' @name deltaSrules
#' @docType data
#' @format The format is: num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.0226
#' 0 0 0 0 ...  - attr(*, "dimnames")=List of 8 ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...  ..$ : chr [1:5] "A" "C" "G" "T"
#' ...  ..$ : chr [1:5] "A" "C" "G" "T" ...
#' @references SantaLucia, J., Jr., & Hicks, D. (2004) The Thermodynamics of
#' DNA Structural Motifs. Annual Review of Biophysics and Biomolecular
#' Structure, 33(1), 415-440. doi:10.1146/annurev.biophys.32.110601.141800.
#' @source Data from a variety of publications by SantaLucia \emph{et al}.
#' @keywords datasets
#' @examples
#' 
#' data(deltaSrules)
#' # dS of the duplex AGCT / A-CT pairing
#' deltaSrules["A", "G", "C", "T", "A", "-", "C", "T"]
#' 
NULL





#' Change in Entropy of Hybridization of RNA/RNA Quadruplets in Solution
#' 
#' An 8D array with four adjacent base pairs of the RNA duplex.  Each dimension
#' has five elements defining the nucleotide at that position ("A", "C", "G",
#' "T", or "-").  The array contains the standard entropy change of probe
#' binding (dS, [kcal/mol]) for every quadruple base pairing.
#' 
#' The first four dimensions correspond to the four top strand positions from
#' 5' to 3'.  The fifth to eighth dimensions correspond to the four bottom
#' strand positions from 5' to 3'.
#' 
#' @name deltaSrulesRNA
#' @docType data
#' @format The format is: num [1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5, 1:5] -0.0182
#' 0 0 0 0 ...  - attr(*, "dimnames")=List of 8 ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...  ..$ : chr [1:5] "A" "C" "G" "U"
#' ...  ..$ : chr [1:5] "A" "C" "G" "U" ...
#' @references SantaLucia, J., Jr., & Hicks, D. (2004) The Thermodynamics of
#' DNA Structural Motifs. Annual Review of Biophysics and Biomolecular
#' Structure, 33(1), 415-440. doi:10.1146/annurev.biophys.32.110601.141800.
#' @source Data from a variety of publications by SantaLucia \emph{et al}.
#' @keywords datasets
#' @examples
#' 
#' data(deltaSrulesRNA)
#' # dS of the duplex AGCU / A-CU pairing
#' deltaSrulesRNA["A", "G", "C", "U", "A", "-", "C", "U"]
#' 
NULL





#' Genes objects and accessors
#' 
#' Gene prediction consist of delimiting the boundaries of regions that
#' function as genes within a genome.  Class \code{Genes} provides objects and
#' functions for storing the boundaries of genes and associated information
#' resulting from gene prediction.
#' 
#' Objects of class \code{Genes} are stored as numeric matrices containing
#' information pertaining to gene predictions.  The matrix columns include the
#' index (\code{"Index"}) of the corresponding sequence in the original genome,
#' the strand (\code{"Strand"}) where the gene is located (either \code{"+"}
#' (\code{0}) or \code{"-"} (\code{1}), the beginning (\code{"Begin"}) and
#' ending (\code{"End"}) positions of the gene, scores acquired during
#' prediction, and whether (\code{!= 0}) or not (\code{0}) the region was
#' predicted to be a gene.  Note that the start of the gene is at the beginning
#' position when the strand is \code{"+"} and end when the strand is
#' \code{"-"}.  By convention, rows with negative values in the \code{"Gene"}
#' column represent non-coding RNAs and rows with positive values represent
#' protein coding genes.
#' 
#' The \code{plot} method will show the total score of each prediction along
#' the genome.  This is most useful when displaying the result of setting
#' \code{allScores} to \code{TRUE} in \code{FindGenes}.  Here, possible genes
#' on either strand will be shown (by default), with the predicted genes
#' highlighted.  The beginning (solid) and ending (dashed) positions are
#' denoted by vertical lines.  Note that the x-axis is cumulative genome
#' position, and changes between genome sequences indices are demarcated by
#' dashed vertical lines.
#' 
#' @name Genes
#' @aliases Genes-class [.Genes print.Genes plot.Genes
#' @param x An object of class \code{Genes}.
#' @param xlim Numeric vector of length 2 specifying the x-axis limits for
#' plotting.
#' @param ylim Numeric vector of length 2 specifying the y-axis limits for
#' plotting.
#' @param interact Logical determining whether the plot is interactive.  If
#' \code{TRUE}, clicking the plot on the right or left side will scroll one
#' frame in that direction.  To end interaction, either right-click, press the
#' escape key, or press the stop button depending on the graphics device in
#' use.
#' @param colorBy Character string indicating the name of the column in
#' \code{x} that should be used for coloring.  Unambiguous abbreviations are
#' also permitted.
#' @param colorRamp A function that will return \code{n} colors when given a
#' number \code{n}.  Examples are \code{rainbow}, \code{heat.colors},
#' \code{terrain.colors}, \code{cm.colors}, or (the default)
#' \code{colorRampPalette}.
#' @param colorGenes Character string specifying the color of genes, or
#' \code{NA} to color genes according to \code{colorBy}.
#' @param i Numeric or character vector of row indices to extract from
#' \code{x}.
#' @param j Numeric or character vector of column indices to extract from
#' \code{x}.  If \code{j} is missing, all columns are included and the returned
#' object will also belong to class \code{Genes}.
#' @param \dots Other optional parameters.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{ExtractGenes}}, \code{\link{FindGenes}},
#' \code{\link{WriteGenes}}
#' @examples
#' 
#' # import a test genome
#' fas <- system.file("extdata",
#' 	"Chlamydia_trachomatis_NC_000117.fas.gz",
#' 	package="DECIPHER")
#' genome <- readDNAStringSet(fas)
#' 
#' x <- FindGenes(genome, allScores=TRUE)
#' x
#' head(unclass(x)) # the underlying structure
#' 
#' plot(x) # default coloring by "Strand"
#' # color by RBS score (blue is weak/low, red is strong/high)
#' plot(x, colorBy="RibosomeBindingSiteScore", colorGenes=NA)
#' # color by fraction of times a gene was chosen
#' plot(x, colorBy="FractionReps", colorGenes=NA)
#' # color by which codon model was selected for each ORF
#' plot(x, colorBy="CodonModel", xlim=c(1, 3e4))
#' 
NULL





#' Genes objects and accessors
#' 
#' Gene prediction consist of delimiting the boundaries of regions that
#' function as genes within a genome.  Class \code{Genes} provides objects and
#' functions for storing the boundaries of genes and associated information
#' resulting from gene prediction.
#' 
#' Objects of class \code{Genes} are stored as numeric matrices containing
#' information pertaining to gene predictions.  The matrix columns include the
#' index (\code{"Index"}) of the corresponding sequence in the original genome,
#' the strand (\code{"Strand"}) where the gene is located (either \code{"+"}
#' (\code{0}) or \code{"-"} (\code{1}), the beginning (\code{"Begin"}) and
#' ending (\code{"End"}) positions of the gene, scores acquired during
#' prediction, and whether (\code{!= 0}) or not (\code{0}) the region was
#' predicted to be a gene.  Note that the start of the gene is at the beginning
#' position when the strand is \code{"+"} and end when the strand is
#' \code{"-"}.  By convention, rows with negative values in the \code{"Gene"}
#' column represent non-coding RNAs and rows with positive values represent
#' protein coding genes.
#' 
#' The \code{plot} method will show the total score of each prediction along
#' the genome.  This is most useful when displaying the result of setting
#' \code{allScores} to \code{TRUE} in \code{FindGenes}.  Here, possible genes
#' on either strand will be shown (by default), with the predicted genes
#' highlighted.  The beginning (solid) and ending (dashed) positions are
#' denoted by vertical lines.  Note that the x-axis is cumulative genome
#' position, and changes between genome sequences indices are demarcated by
#' dashed vertical lines.
#' 
#' @name Genes
#' @aliases Genes Genes-class [.Genes print.Genes plot.Genes
#' @param x An object of class \code{Genes}.
#' @param xlim Numeric vector of length 2 specifying the x-axis limits for
#' plotting.
#' @param ylim Numeric vector of length 2 specifying the y-axis limits for
#' plotting.
#' @param interact Logical determining whether the plot is interactive.  If
#' \code{TRUE}, clicking the plot on the right or left side will scroll one
#' frame in that direction.  To end interaction, either right-click, press the
#' escape key, or press the stop button depending on the graphics device in
#' use.
#' @param colorBy Character string indicating the name of the column in
#' \code{x} that should be used for coloring.  Unambiguous abbreviations are
#' also permitted.
#' @param colorRamp A function that will return \code{n} colors when given a
#' number \code{n}.  Examples are \code{rainbow}, \code{heat.colors},
#' \code{terrain.colors}, \code{cm.colors}, or (the default)
#' \code{colorRampPalette}.
#' @param colorGenes Character string specifying the color of genes, or
#' \code{NA} to color genes according to \code{colorBy}.
#' @param i Numeric or character vector of row indices to extract from
#' \code{x}.
#' @param j Numeric or character vector of column indices to extract from
#' \code{x}.  If \code{j} is missing, all columns are included and the returned
#' object will also belong to class \code{Genes}.
#' @param \dots Other optional parameters.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{ExtractGenes}}, \code{\link{FindGenes}},
#' \code{\link{WriteGenes}}
#' @examples
#' 
#' 
#' # import a test genome
#' fas <- system.file("extdata",
#' 	"Chlamydia_trachomatis_NC_000117.fas.gz",
#' 	package="DECIPHER")
#' genome <- readDNAStringSet(fas)
#' 
#' x <- FindGenes(genome, allScores=TRUE)
#' x
#' head(unclass(x)) # the underlying structure
#' 
#' plot(x) # default coloring by "Strand"
#' # color by RBS score (blue is weak/low, red is strong/high)
#' plot(x, colorBy="RibosomeBindingSiteScore", colorGenes=NA)
#' # color by fraction of times a gene was chosen
#' plot(x, colorBy="FractionReps", colorGenes=NA)
#' # color by which codon model was selected for each ORF
#' plot(x, colorBy="CodonModel", xlim=c(1, 3e4))
#' 
#' 
NULL





#' Mutual Information for Protein Secondary Structure Prediction
#' 
#' Arrays containing values of mutual information for single residues
#' (\code{HEC_MI1}) and pairs of residues (\code{HEC_MI2}) located within 10
#' residues of the position being predicted (position "0").  The arrays have
#' dimensions corresponding to the 20 (standard) amino acids, positions (-10 to
#' 10), and states (helix ("H"), sheet ("E"), or coil ("C")).
#' 
#' The values in each matrix were derived based on a set of 15,201 proteins in
#' the ASTRAL Compendium (Chandonia, 2004).  The 8-states assigned by the
#' Dictionary of Protein Secondary Structure (DSSP) were reduced to 3-states
#' via H = G, H, or I; E = E; and C = B, S, C, or T.
#' 
#' @name HEC_MI
#' @aliases HEC_MI1 HEC_MI2
#' @docType data
#' @format The format of HEC_MI1 is: num [1:20, 1:21, 1:3] 0.04264 -0.00117
#' 0.02641 0.08264 -0.04876 ...  - attr(*, "dimnames")=List of 3 ..$ : chr
#' [1:20] "A" "R" "N" "D" ...  ..$ : chr [1:21] "-10" "-9" "-8" "-7" ...  ..$ :
#' chr [1:3] "H" "E" "C"
#' 
#' The format of HEC_MI2 is: num [1:20, 1:20, 1:21, 1:21, 1:3] 2.56 -Inf -Inf
#' -Inf -Inf ...  - attr(*, "dimnames")=List of 5 ..$ : chr [1:20] "A" "R" "N"
#' "D" ...  ..$ : chr [1:20] "A" "R" "N" "D" ...  ..$ : chr [1:21] "-10" "-9"
#' "-8" "-7" ...  ..$ : chr [1:21] "-10" "-9" "-8" "-7" ...  ..$ : chr [1:3]
#' "H" "E" "C"
#' @references Chandonia, J. M. (2004). The ASTRAL Compendium in 2004.
#' \emph{Nucleic Acids Research}, \bold{32(90001)}, 189D-192.
#' doi:10.1093/nar/gkh034.
#' @keywords datasets
#' @examples
#' 
#' data(HEC_MI1)
#' # the contribution of an arginine ("R")
#' # located 3 residues left of center
#' # to a helical ("H") state at the center
#' HEC_MI1["R", "-3", "H"]
#' 
#' data(HEC_MI2)
#' # the contribution of arginine and lysine ("K")
#' # located at positions -1 and +1, respectively
#' # to a coil ("C") state at the center position
#' HEC_MI2["R", "K", "-1", "1", "C"]
#' 
#' matplot(-10:10, t(HEC_MI1[,, "H"]),
#'        type="l", col=1:8, lty=rep(1:3, each=8),
#'        xlab="Amino Acid Position Relative to Center",
#'        ylab="Log-Odds of Helix at Center Position")
#' legend("bottomleft",
#'        lwd=1, col=1:8, lty=rep(1:3, each=8),
#'        legend=dimnames(HEC_MI1)[[1]], ncol=2)
#' 
NULL





#' MIQS Amino Acid Substitution Matrix
#' 
#' The MIQS amino acid substitution matrix defined by Yamada & Tomii (2014).
#' 
#' Substitution matrix values represent the log-odds of observing an aligned
#' pair of amino acids versus the likelihood of finding the pair by chance.
#' Values in the MIQS matrix are in units of third-bits (\eqn{log(odds\
#' ratio)*3/log(2)}).
#' 
#' @name MIQS
#' @docType data
#' @format The format is: num [1:25, 1:25] 3.2 -1.3 -0.4 -0.4 1.5 -0.2 -0.4 0.4
#' -1.2 -1.3 ...  - attr(*, "dimnames")=List of 2 ..$ : chr [1:25] "A" "R" "N"
#' "D" ...  ..$ : chr [1:25] "A" "R" "N" "D" ...
#' @source Yamada, K., & Tomii, K. (2014). Revisiting amino acid substitution
#' matrices for identifying distantly related proteins. \emph{Bioinformatics},
#' \bold{30(3)}, 317-325.
#' @keywords datasets
#' @examples
#' 
#' data(MIQS)
#' MIQS["A", "R"] # score for A/R pairing
#' 
#' data(BLOSUM62)
#' plot(BLOSUM62[1:20, 1:20], MIQS[1:20, 1:20])
#' abline(a=0, b=1)
#' 
NULL





#' NonCoding Objects and Methods
#' 
#' Non-coding RNAs can be represented by their conserved sequence motifs,
#' secondary structure, and k-mer frequencies.  Class \code{NonCoding} provides
#' objects and functions for representing non-coding RNAs.
#' 
#' Objects of class \code{NonCoding} are stored as lists containing a compact
#' representation of a family of non-coding RNAs.  The first list component is
#' a matrix of sequence motifs that identify the non-coding RNAs, the second is
#' a matrix of hairpin loops that are conserved across the family, the third is
#' a list of k-mer frequencies derived from representative sequences, and the
#' fourth is a vector of log-odds scores for sequence lengths.  An optional
#' fifth list component denotes the log-odds scores for dependencies among
#' patterns.  Patterns are defined by their distance to either end of the
#' non-coding RNA, which helps to identify the boundaries of the non-coding RNA
#' in a genome.
#' 
#' @name NonCoding
#' @aliases NonCoding-class print.NonCoding
#' @param x An object of class \code{NonCoding}.
#' @param \dots Other optional parameters.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{LearnNonCoding}}, \code{\link{FindNonCoding}}
#' @references Wright, E. S. (2021). FindNonCoding: rapid and simple detection
#' of non-coding RNAs in genomes. Bioinformatics.
#' https://doi.org/10.1093/bioinformatics/btab708
#' @examples
#' 
#' data(NonCodingRNA_Bacteria)
#' x <- NonCodingRNA_Bacteria
#' print(x)
#' class(x)
#' attributes(x[[1]])
#' x[[1]] # the first non-coding RNA
#' x[[1]][["motifs"]] # sequence motifs
#' x[[1]][["hairpins"]] # hairpin loops
#' head(x[[1]][["kmers"]]) # k-mer frequencies
#' 
NULL





#' NonCoding Objects and Methods
#' 
#' Non-coding RNAs can be represented by their conserved sequence motifs,
#' secondary structure, and k-mer frequencies.  Class \code{NonCoding} provides
#' objects and functions for representing non-coding RNAs.
#' 
#' Objects of class \code{NonCoding} are stored as lists containing a compact
#' representation of a family of non-coding RNAs.  The first list component is
#' a matrix of sequence motifs that identify the non-coding RNAs, the second is
#' a matrix of hairpin loops that are conserved across the family, the third is
#' a list of k-mer frequencies derived from representative sequences, and the
#' fourth is a vector of log-odds scores for sequence lengths.  An optional
#' fifth list component denotes the log-odds scores for dependencies among
#' patterns.  Patterns are defined by their distance to either end of the
#' non-coding RNA, which helps to identify the boundaries of the non-coding RNA
#' in a genome.
#' 
#' @name NonCoding
#' @aliases NonCoding NonCoding-class print.NonCoding
#' @param x An object of class \code{NonCoding}.
#' @param \dots Other optional parameters.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{LearnNonCoding}}, \code{\link{FindNonCoding}}
#' @references Wright, E. S. (2021). FindNonCoding: rapid and simple detection
#' of non-coding RNAs in genomes. Bioinformatics.
#' https://doi.org/10.1093/bioinformatics/btab708
#' @examples
#' 
#' 
#' data(NonCodingRNA_Bacteria)
#' x <- NonCodingRNA_Bacteria
#' print(x)
#' class(x)
#' attributes(x[[1]])
#' x[[1]] # the first non-coding RNA
#' x[[1]][["motifs"]] # sequence motifs
#' x[[1]][["hairpins"]] # hairpin loops
#' head(x[[1]][["kmers"]]) # k-mer frequencies
#' 
#' 
NULL





#' NonCoding Models for Common Non-Coding RNA Families
#' 
#' Pre-trained with \code{NonCoding} models for common RNA families found in
#' genomes from organisms belonging to each domain of life.
#' 
#' A set of \code{NonCoding} models contained in a named list.  Models were
#' built from up to 1000 representative sequences per non-coding RNA family.
#' 
#' @name NonCodingRNA
#' @aliases NonCodingRNA_Archaea NonCodingRNA_Bacteria NonCodingRNA_Eukarya
#' @docType data
#' @source Models were built from sequences belonging to families in tRNADB-CE
#' (\url{http://trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/index.cgi}) or Rfam
#' (\url{http://rfam.xfam.org}).
#' @keywords datasets
#' @examples
#' 
#' data(NonCodingRNA_Archaea)
#' data(NonCodingRNA_Bacteria)
#' data(NonCodingRNA_Eukarya)
#' names(NonCodingRNA_Bacteria)
#' head(NonCodingRNA_Bacteria)
#' 
NULL





#' PFASUM Amino Acid Substitution Matrices
#' 
#' The PFASUM amino acid substitution matrices defined by Keul, F., et al.
#' (2017).
#' 
#' Substitution matrix values represent the log-odds of observing an aligned
#' pair of amino acids versus the likelihood of finding the pair by chance.
#' The PFASUM substitution matrices are stored as an array named by each
#' sub-matrix's similarity threshold.  (See examples section below.)  In all
#' cases values are in units of third-bits (\eqn{log(odds\ ratio)*3/log(2)}).
#' 
#' @name PFASUM
#' @docType data
#' @format The format is: num [1:25, 1:25, 1:90] 0.9492 -1.7337 0.2764 1.8153
#' 0.0364 ...  - attr(*, "dimnames")=List of 3 ..$ : chr [1:25] "A" "R" "N" "D"
#' ...  ..$ : chr [1:25] "A" "R" "N" "D" ...  ..$ : chr [1:90] "11" "12" "13"
#' "14" ...
#' @source Keul, F., et al. (2017). PFASUM: a substitution matrix from Pfam
#' structural alignments. \emph{BMC Bioinformatics}, \bold{18(1)}, 293.
#' @keywords datasets
#' @examples
#' 
#' data(PFASUM)
#' PFASUM31 <- PFASUM[,, "31"] # the PFASUM31 matrix
#' PFASUM31["A", "R"] # score for A/R pairing
#' 
#' data(BLOSUM62)
#' plot(BLOSUM62[1:20, 1:20], PFASUM31[1:20, 1:20])
#' abline(a=0, b=1)
#' 
NULL





#' Common Restriction Enzyme's Cut Sites
#' 
#' A character vector of common restriction sites named by the restriction
#' enzyme that cuts at each site.  Sequence specificity is listed in 5' to 3'
#' orientation based on the \code{IUPAC_CODE_MAP}.  The cut site is either
#' signified by a ``/'' for palindromic sites, or two numbers giving the
#' position of the top and bottom cut positions relative to the site's 3'-end.
#' 
#' 
#' @name RESTRICTION_ENZYMES
#' @docType data
#' @format The format is: Named chr [1:224] "GACGT/C" "G/GTACC" "GT/MKAC" ...
#' - attr(*, "names")= chr [1:224] "AatII" "Acc65I" "AccI" "AciI" ...
#' @source Restriction enzymes sold by \href{http://www.neb.comNew England
#' BioLabs}.
#' @keywords datasets
#' @examples
#' 
#' data(RESTRICTION_ENZYMES)
#' RESTRICTION_ENZYMES
#' 
NULL





#' Synteny blocks and hits
#' 
#' Syntenic blocks are DNA segments composed of conserved hits occurring in the
#' same order on two sequences.  The two sequences are typically chromosomes of
#' different species that are hypothesized to contain homology.  Class
#' \code{"Synteny"} provides objects and functions for storing and viewing
#' syntenic blocks and hits that are shared between sequences.
#' 
#' Objects of class \code{Synteny} are stored as square matrices of list
#' elements with \code{dimnames} giving the ``identifier'' of the corresponding
#' sequences.  The synteny matrix can be separated into three parts: along,
#' above, and below the diagonal.  Each list element along the diagonal
#' contains an integer vector with the width of the sequence(s) belonging to
#' that ``identifier''.  List elements above the diagonal (column \emph{j} >
#' row \emph{i}) each contain a \code{matrix} with ``hits'' corresponding to
#' matches between sequences \emph{i} and \emph{j}.  List elements below the
#' diagonal each contain a \code{matrix} with ``blocks'' of synteny between
#' sequences \emph{j} and \emph{i}.
#' 
#' The \code{pairs} method creates a scatterplot matrix from a \code{Synteny}
#' object.  Dot plots above the diagonal show hits between identifier \emph{i}
#' and \emph{j}, where forward hits are colored in black, and hits to the
#' reverse strand of identifier \emph{j} are colored in red.  Plots below the
#' diagonal show blocks of synteny colored by their score, from green (highest
#' scoring) to blue to magenta (lowest scoring).
#' 
#' The \code{plot} method displays a bar view of the sequences in the same
#' order as the input object (\code{x}).  The coloring scheme of each bar is
#' determined by the \code{colorBy} argument, and the color palette is set by
#' \code{colorRamp}.  When \code{colorBy} is an index, the sequences are
#' colored according to regions of shared homology with the specified reference
#' sequence (by default \code{1}).  If \code{colorBy} is ``neighbor'' then
#' shared syntenic blocks are connected between neighboring sequences.  If
#' \code{colorBy} is ``frequency'' then positions in each sequence are colored
#' based on the degree of conservation with the other sequences.  In each case,
#' regions that have no correspondence in the other sequence(s) are colored
#' \code{barColor}.
#' 
#' @name Synteny
#' @aliases Synteny-class [.Synteny print.Synteny plot.Synteny pairs.Synteny
#' @param x An object of class \code{Synteny}.
#' @param bounds Logical specifying whether to plot sequence boundaries as
#' horizontal or vertical lines.
#' @param boxBlocks Logical indicating whether to draw a rectangle around hits
#' belonging to the same block of synteny.
#' @param colorBy Numeric giving the index of a reference sequence, or a
#' character string indicating to color by ``neighbor'', ``frequency'', or
#' ``none''.  (See details section below.)
#' @param colorRamp A function that will return \code{n} colors when given a
#' number \code{n}.  Examples are \code{rainbow}, \code{heat.colors},
#' \code{terrain.colors}, \code{cm.colors}, or (the default)
#' \code{colorRampPalette}.
#' @param barColor Character string giving the background color of each bar.
#' @param barSides Logical indicating whether to draw black lines along the
#' long-sides of each bar.
#' @param horizontal Logical indicating whether to plot the sequences
#' horizontally (\code{TRUE}) or vertically (\code{FALSE}).
#' @param labels Character vector providing names corresponding to each
#' ``identifier'' for labels on the diagonal.
#' @param width Numeric giving the fractional width of each bar between zero
#' and one.
#' @param scaleBar Logical controlling whether a scale bar is drawn when
#' \code{colorBy} is ``frequency''.  The scale bar displays the mapping between
#' color and the level of sequence conservation.  Not applicable when
#' \code{colorBy} is a value other than ``frequency''.
#' @param gap Distance between subplots, in margin lines.
#' @param line.main If \code{main} is specified, \code{line.main} provides the
#' \code{line} argument to \code{mtext}.
#' @param cex.labels Magnification of the labels.
#' @param font.labels Font of labels on the diagonal.
#' @param quote Logical indicating whether to print the output surrounded by
#' quotes.
#' @param right Logical specifying whether to right align strings.
#' @param \dots Other graphical parameters for \code{pairs} or \code{plot},
#' including: \code{main}, \code{cex.main}, \code{font.main}, and \code{oma}.
#' Other arguments for \code{print}, including \code{print.gap} and \code{max}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{AlignSynteny}}, \code{\link{FindSynteny}}
#' @examples
#' 
#' # a small example:
#' dbConn <- dbConnect(SQLite(), ":memory:")
#' s1 <- DNAStringSet("ACTAGACCCAGACCGATAAACGGACTGGACAAG")
#' s3 <- reverseComplement(s1)
#' s2 <- c(s1, s3)
#' Seqs2DB(c(c(s1, s2), s3),
#'         "XStringSet",
#'         dbConn,
#'         c("s1", "s2", "s2", "s3"))
#' syn <- FindSynteny(dbConn, minScore=1)
#' syn # Note:  > 100% hits because of sequence reuse across blocks
#' pairs(syn, boxBlocks=TRUE)
#' plot(syn)
#' dbDisconnect(dbConn)
#' 
#' # a larger example:
#' db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
#' synteny <- FindSynteny(db, minScore=50)
#' class(synteny) # 'Synteny'
#' synteny
#' 
#' # accessing parts
#' i <- 1
#' j <- 2
#' synteny[i, i][[1]] # width of sequences in i
#' synteny[j, j][[1]] # width of sequences in j
#' head(synteny[i, j][[1]]) # hits between i & j
#' synteny[j, i][[1]] # blocks between i & j
#' 
#' # plotting
#' pairs(synteny) # dot plots
#' pairs(synteny, boxBlocks=TRUE) # boxes around blocks
#' 
#' plot(synteny) # bar view colored by position in genome 1
#' plot(synteny, barColor="#268FD6") # emphasize missing regions
#' plot(synteny, "frequency") # most regions are shared by all
#' plot(synteny, "frequency", colorRamp=rainbow) # change the colors
#' plot(synteny, "neighbor") # connect neighbors
#' 
NULL





#' Synteny blocks and hits
#' 
#' Syntenic blocks are DNA segments composed of conserved hits occurring in the
#' same order on two sequences.  The two sequences are typically chromosomes of
#' different species that are hypothesized to contain homology.  Class
#' \code{"Synteny"} provides objects and functions for storing and viewing
#' syntenic blocks and hits that are shared between sequences.
#' 
#' Objects of class \code{Synteny} are stored as square matrices of list
#' elements with \code{dimnames} giving the \verb{identifier'' of the
#' corresponding sequences.  The synteny matrix can be separated into three
#' parts: along, above, and below the diagonal.  Each list element along the
#' diagonal contains an integer vector with the width of the sequence(s)
#' belonging to that }identifier''.  List elements above the diagonal (column
#' \emph{j} > row \emph{i}) each contain a \code{matrix} with \verb{hits''
#' corresponding to matches between sequences \emph{i} and \emph{j}.  List
#' elements below the diagonal each contain a \code{matrix} with }blocks'' of
#' synteny between sequences \emph{j} and \emph{i}.
#' 
#' The \code{pairs} method creates a scatterplot matrix from a \code{Synteny}
#' object.  Dot plots above the diagonal show hits between identifier \emph{i}
#' and \emph{j}, where forward hits are colored in black, and hits to the
#' reverse strand of identifier \emph{j} are colored in red.  Plots below the
#' diagonal show blocks of synteny colored by their score, from green (highest
#' scoring) to blue to magenta (lowest scoring).
#' 
#' The \code{plot} method displays a bar view of the sequences in the same
#' order as the input object (\code{x}).  The coloring scheme of each bar is
#' determined by the \code{colorBy} argument, and the color palette is set by
#' \code{colorRamp}.  When \code{colorBy} is an index, the sequences are
#' colored according to regions of shared homology with the specified reference
#' sequence (by default \code{1}).  If \code{colorBy} is \verb{neighbor'' then
#' shared syntenic blocks are connected between neighboring sequences.  If
#' \code{colorBy} is }frequency'' then positions in each sequence are colored
#' based on the degree of conservation with the other sequences.  In each case,
#' regions that have no correspondence in the other sequence(s) are colored
#' \code{barColor}.
#' 
#' @name Synteny
#' @aliases Synteny Synteny-class [.Synteny print.Synteny plot.Synteny
#' pairs.Synteny
#' @param x An object of class \code{Synteny}.
#' @param bounds Logical specifying whether to plot sequence boundaries as
#' horizontal or vertical lines.
#' @param boxBlocks Logical indicating whether to draw a rectangle around hits
#' belonging to the same block of synteny.
#' @param colorBy Numeric giving the index of a reference sequence, or a
#' character string indicating to color by \verb{neighbor'', }frequency'', or
#' ``none''.  (See details section below.)
#' @param colorRamp A function that will return \code{n} colors when given a
#' number \code{n}.  Examples are \code{rainbow}, \code{heat.colors},
#' \code{terrain.colors}, \code{cm.colors}, or (the default)
#' \code{colorRampPalette}.
#' @param barColor Character string giving the background color of each bar.
#' @param barSides Logical indicating whether to draw black lines along the
#' long-sides of each bar.
#' @param horizontal Logical indicating whether to plot the sequences
#' horizontally (\code{TRUE}) or vertically (\code{FALSE}).
#' @param labels Character vector providing names corresponding to each
#' ``identifier'' for labels on the diagonal.
#' @param width Numeric giving the fractional width of each bar between zero
#' and one.
#' @param scaleBar Logical controlling whether a scale bar is drawn when
#' \code{colorBy} is \verb{frequency''.  The scale bar displays the mapping
#' between color and the level of sequence conservation.  Not applicable when
#' \code{colorBy} is a value other than }frequency''.
#' @param gap Distance between subplots, in margin lines.
#' @param line.main If \code{main} is specified, \code{line.main} provides the
#' \code{line} argument to \code{mtext}.
#' @param cex.labels Magnification of the labels.
#' @param font.labels Font of labels on the diagonal.
#' @param quote Logical indicating whether to print the output surrounded by
#' quotes.
#' @param right Logical specifying whether to right align strings.
#' @param \dots Other graphical parameters for \code{pairs} or \code{plot},
#' including: \code{main}, \code{cex.main}, \code{font.main}, and \code{oma}.
#' Other arguments for \code{print}, including \code{print.gap} and \code{max}.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{AlignSynteny}}, \code{\link{FindSynteny}}
#' @examples
#' 
#' 
#' # a small example:
#' dbConn <- dbConnect(SQLite(), ":memory:")
#' s1 <- DNAStringSet("ACTAGACCCAGACCGATAAACGGACTGGACAAG")
#' s3 <- reverseComplement(s1)
#' s2 <- c(s1, s3)
#' Seqs2DB(c(c(s1, s2), s3),
#'         "XStringSet",
#'         dbConn,
#'         c("s1", "s2", "s2", "s3"))
#' syn <- FindSynteny(dbConn, minScore=1)
#' syn # Note:  > 100% hits because of sequence reuse across blocks
#' pairs(syn, boxBlocks=TRUE)
#' plot(syn)
#' dbDisconnect(dbConn)
#' 
#' # a larger example:
#' db <- system.file("extdata", "Influenza.sqlite", package="DECIPHER")
#' synteny <- FindSynteny(db, minScore=50)
#' class(synteny) # 'Synteny'
#' synteny
#' 
#' # accessing parts
#' i <- 1
#' j <- 2
#' synteny[i, i][[1]] # width of sequences in i
#' synteny[j, j][[1]] # width of sequences in j
#' head(synteny[i, j][[1]]) # hits between i & j
#' synteny[j, i][[1]] # blocks between i & j
#' 
#' # plotting
#' pairs(synteny) # dot plots
#' pairs(synteny, boxBlocks=TRUE) # boxes around blocks
#' 
#' plot(synteny) # bar view colored by position in genome 1
#' plot(synteny, barColor="#268FD6") # emphasize missing regions
#' plot(synteny, "frequency") # most regions are shared by all
#' plot(synteny, "frequency", colorRamp=rainbow) # change the colors
#' plot(synteny, "neighbor") # connect neighbors
#' 
#' 
NULL





#' Taxa training and testing objects
#' 
#' Taxonomic classification is the process of assigning an organism a label
#' that is part of a taxonomic hierarchy (e.g., Phylum, Class, Order, Family,
#' Genus).  Here, labels are assigned based on an organism's DNA or RNA
#' sequence at a rank level determined by the classification's confidence.
#' Class \code{Taxa} provides objects and functions for storing and viewing
#' training and testing objects used in taxonomic classification.
#' 
#' Objects of class \code{Taxa} are stored as lists, and can have either
#' subclass \code{Train} or \code{Test}.  The function \code{LearnTaxa} returns
#' an object of subclass \code{Train}, while the function \code{IdTaxa} can
#' return an object of class \code{Test}.
#' 
#' \code{Train}ing objects are built from a set of reference sequences with
#' known taxonomic classifications.  List elements contain information required
#' by \code{IdTaxa} for assigning a classification to test sequences.
#' 
#' \code{Test}ing objects can be generated by \code{IdTaxa} from a
#' \code{Train}ing object and a set of test sequences.  Each list element
#' contains the taxon, confidence, and (optionally) rank name of the taxonomic
#' assignment.
#' 
#' The information stored in \code{Taxa} can be visualized with the \code{plot}
#' function or displayed with \code{print}.  Only objects of subclass
#' \code{Train} can be subsetted without losing their class.
#' 
#' @name Taxa
#' @aliases Taxa-class [.Taxa c.Taxa print.Taxa plot.Taxa
#' @param x An object of class \code{Taxa} with subclass \code{Train} or
#' \code{Test}.
#' @param y An (optional) object of class \code{Taxa} with the opposite
#' subclass as \code{x}.
#' @param showRanks Logical specifying whether to show all rank levels when
#' plotting an object of class \code{Taxa} and subclass \code{Test}.  If
#' \code{TRUE} (the default), then ranks are shown as (colored) concentric
#' rings with radial lines delimiting taxa boundaries.
#' @param n Numeric vector giving the frequency of each classification if
#' \code{x} or \code{y} is an object of subclass \code{Test}, or the default
#' (\code{NULL}) to treat all classifications as occurring once.  Typically,
#' specifying \code{n} is useful when the classifications represent varying
#' numbers of observations, e.g., when only unique sequences were originally
#' classified.
#' @param \dots Other optional parameters.
#' @param i Numeric or character vector of indices to extract from objects of
#' class \code{Taxa} with subclass \code{Test}.
#' @param j Numeric or character vector of rank levels to extract from objects
#' of class \code{Taxa} with subclass \code{Test}.
#' @param threshold Numeric specifying the confidence \code{threshold} at which
#' to truncate the output taxonomic classifications. Note that \code{threshold}
#' must be higher than the original for the classifications to change.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{LearnTaxa}}, \code{\link{IdTaxa}}
#' @examples
#' 
#' data("TrainingSet_16S")
#' plot(TrainingSet_16S)
#' 
#' # import test sequences
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' 
#' # remove any gaps in the sequences
#' dna <- RemoveGaps(dna)
#' 
#' # classify the test sequences
#' ids <- IdTaxa(dna, TrainingSet_16S, strand="top")
#' ids
#' 
#' plot(ids) # plot all rank levels
#' plot(ids[, 1:4]) # plot the first rank levels
#' plot(ids[j=c("rootrank", "class", "genus")]) # plot specific rank levels
#' plot(ids[threshold=70]) # plot high confidence classifications
#' 
NULL





#' Taxa training and testing objects
#' 
#' Taxonomic classification is the process of assigning an organism a label
#' that is part of a taxonomic hierarchy (e.g., Phylum, Class, Order, Family,
#' Genus).  Here, labels are assigned based on an organism's DNA or RNA
#' sequence at a rank level determined by the classification's confidence.
#' Class \code{Taxa} provides objects and functions for storing and viewing
#' training and testing objects used in taxonomic classification.
#' 
#' Objects of class \code{Taxa} are stored as lists, and can have either
#' subclass \code{Train} or \code{Test}.  The function \code{LearnTaxa} returns
#' an object of subclass \code{Train}, while the function \code{IdTaxa} can
#' return an object of class \code{Test}.
#' 
#' \code{Train}ing objects are built from a set of reference sequences with
#' known taxonomic classifications.  List elements contain information required
#' by \code{IdTaxa} for assigning a classification to test sequences.
#' 
#' \code{Test}ing objects can be generated by \code{IdTaxa} from a
#' \code{Train}ing object and a set of test sequences.  Each list element
#' contains the taxon, confidence, and (optionally) rank name of the taxonomic
#' assignment.
#' 
#' The information stored in \code{Taxa} can be visualized with the \code{plot}
#' function or displayed with \code{print}.  Only objects of subclass
#' \code{Train} can be subsetted without losing their class.
#' 
#' @name Taxa
#' @aliases Taxa Taxa-class [.Taxa c.Taxa print.Taxa plot.Taxa
#' @param x An object of class \code{Taxa} with subclass \code{Train} or
#' \code{Test}.
#' @param y An (optional) object of class \code{Taxa} with the opposite
#' subclass as \code{x}.
#' @param showRanks Logical specifying whether to show all rank levels when
#' plotting an object of class \code{Taxa} and subclass \code{Test}.  If
#' \code{TRUE} (the default), then ranks are shown as (colored) concentric
#' rings with radial lines delimiting taxa boundaries.
#' @param n Numeric vector giving the frequency of each classification if
#' \code{x} or \code{y} is an object of subclass \code{Test}, or the default
#' (\code{NULL}) to treat all classifications as occurring once.  Typically,
#' specifying \code{n} is useful when the classifications represent varying
#' numbers of observations, e.g., when only unique sequences were originally
#' classified.
#' @param \dots Other optional parameters.
#' @param i Numeric or character vector of indices to extract from objects of
#' class \code{Taxa} with subclass \code{Test}.
#' @param j Numeric or character vector of rank levels to extract from objects
#' of class \code{Taxa} with subclass \code{Test}.
#' @param threshold Numeric specifying the confidence \code{threshold} at which
#' to truncate the output taxonomic classifications. Note that \code{threshold}
#' must be higher than the original for the classifications to change.
#' @author Erik Wright \email{eswright@@pitt.edu}
#' @seealso \code{\link{LearnTaxa}}, \code{\link{IdTaxa}}
#' @examples
#' 
#' 
#' data("TrainingSet_16S")
#' plot(TrainingSet_16S)
#' 
#' # import test sequences
#' fas <- system.file("extdata", "Bacteria_175seqs.fas", package="DECIPHER")
#' dna <- readDNAStringSet(fas)
#' 
#' # remove any gaps in the sequences
#' dna <- RemoveGaps(dna)
#' 
#' # classify the test sequences
#' ids <- IdTaxa(dna, TrainingSet_16S, strand="top")
#' ids
#' 
#' plot(ids) # plot all rank levels
#' plot(ids[, 1:4]) # plot the first rank levels
#' plot(ids[j=c("rootrank", "class", "genus")]) # plot specific rank levels
#' plot(ids[threshold=70]) # plot high confidence classifications
#' 
#' 
NULL





#' Training Set for Classification of 16S rRNA Gene Sequences
#' 
#' A pre-trained classifier for 16S rRNA gene sequences generated by
#' \code{\link{LearnTaxa}}.
#' 
#' The original training sequences were pruned to a maximum of one sequence per
#' group, as described in the 'Classifying Sequences' vignette.
#' 
#' @name TrainingSet_16S
#' @docType data
#' @format A training set of class 'Taxa' * K-mer size: 8 * Number of rank
#' levels: 10 * Total number of sequences: 2472 * Number of taxonomic groups:
#' 2472 * Number of problem groups: 5 * Number of problem sequences: 8
#' @note This 16S rRNA training set is provided for illustrative purposes only.
#' It is highly recommended to use a more comprehensive training set when
#' classifying real sequences.  Examples of comprehensive training sets can be
#' found at \url{http://DECIPHER.codes/Download.html}.
#' @references Whitman, W.B., Goodfellow, M., Kampfer, P., Busse, H.-J.,
#' Trujillo, M.E., Ludwig, W. & Suzuki, K.-i. (eds., 2012). Bergey's Manual of
#' Systematic Bacteriology, 2nd ed., Springer-Verlag, New York, NY.
#' @source Derived from version 16 of the RDP Training Set
#' (\url{http://rdp.cme.msu.edu}) based on Bergey's Manual.
#' @keywords datasets
#' @examples
#' 
#' data(TrainingSet_16S)
#' TrainingSet_16S
#' plot(TrainingSet_16S)
#' 
NULL



