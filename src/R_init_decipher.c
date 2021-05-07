/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

// for math functions
#include <math.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

/*
 * -- REGISTRATION OF THE .Call ENTRY POINTS ---
 */
static const R_CallMethodDef callMethods[] = { // method call, pointer, num args
	{"consensusSequence", (DL_FUNC) &consensusSequence, 6},
	{"consensusSequenceAA", (DL_FUNC) &consensusSequenceAA, 6},
	{"cluster", (DL_FUNC) &cluster, 7},
	{"reclusterNJ", (DL_FUNC) &reclusterNJ, 2},
	{"reclusterUPGMA", (DL_FUNC) &reclusterUPGMA, 2},
	{"clusterML", (DL_FUNC) &clusterML, 8},
	{"distMatrix", (DL_FUNC) &distMatrix, 10},
	{"gaps", (DL_FUNC) &gaps, 2},
	{"designProbes", (DL_FUNC) &designProbes, 16},
	{"commonGaps", (DL_FUNC) &commonGaps, 1},
	{"multiMatch", (DL_FUNC) &multiMatch, 3},
	{"multiMatchUpper", (DL_FUNC) &multiMatchUpper, 3},
	{"multiMatchCharNotNA", (DL_FUNC) &multiMatchCharNotNA, 1},
	{"replaceChars", (DL_FUNC) &replaceChars, 3},
	{"replaceChar", (DL_FUNC) &replaceChar, 3},
	{"intMatch", (DL_FUNC) &intMatch, 2},
	{"firstMatchUpper", (DL_FUNC) &firstMatchUpper, 3},
	{"terminalMismatch", (DL_FUNC) &terminalMismatch, 5},
	{"NNLS", (DL_FUNC) &NNLS, 10},
	{"sparseMult", (DL_FUNC) &sparseMult, 6},
	{"calculateDeltaG", (DL_FUNC) &calculateDeltaG, 3},
	{"calculateHairpinDeltaG", (DL_FUNC) &calculateHairpinDeltaG, 3},
	{"calculateFISH", (DL_FUNC) &calculateFISH, 2},
	{"alignProfiles", (DL_FUNC) &alignProfiles, 15},
	{"alignProfilesAA", (DL_FUNC) &alignProfilesAA, 12},
	{"consensusProfile", (DL_FUNC) &consensusProfile, 3},
	{"consensusProfileAA", (DL_FUNC) &consensusProfileAA, 3},
	{"matchLists", (DL_FUNC) &matchLists, 4},
	{"adjustHeights", (DL_FUNC) &adjustHeights, 1},
	{"enumerateSequence", (DL_FUNC) &enumerateSequence, 3},
	{"enumerateSequenceAA", (DL_FUNC) &enumerateSequenceAA, 2},
	{"enumerateSequenceReducedAA", (DL_FUNC) &enumerateSequenceReducedAA, 4},
	{"enumerateGappedSequence", (DL_FUNC) &enumerateGappedSequence, 3},
	{"enumerateGappedSequenceAA", (DL_FUNC) &enumerateGappedSequenceAA, 3},
	{"matchOrder", (DL_FUNC) &matchOrder, 4},
	{"matchRanges", (DL_FUNC) &matchRanges, 5},
	{"firstSeqsEqual", (DL_FUNC) &firstSeqsEqual, 8},
	{"matchOrderDual", (DL_FUNC) &matchOrderDual, 3},
	{"boundedMatches", (DL_FUNC) &boundedMatches, 3},
	{"gcContent", (DL_FUNC) &gcContent, 3},
	{"intDist", (DL_FUNC) &intDist, 7},
	{"meltPolymer", (DL_FUNC) &meltPolymer, 4},
	{"insertGaps", (DL_FUNC) &insertGaps, 5},
	{"intMatchOnce", (DL_FUNC) &intMatchOnce, 4},
	{"expandAmbiguities", (DL_FUNC) &expandAmbiguities, 2},
	{"colScoresAA", (DL_FUNC) &colScoresAA, 7},
	{"colScores", (DL_FUNC) &colScores, 7},
	{"shiftGaps", (DL_FUNC) &shiftGaps, 8},
	{"shiftGapsAA", (DL_FUNC) &shiftGapsAA, 8},
	{"removeCommonGaps", (DL_FUNC) &removeCommonGaps, 3},
	{"predictHEC", (DL_FUNC) &predictHEC, 6},
	{"clearIns", (DL_FUNC) &clearIns, 1},
	{"all", (DL_FUNC) &all, 1},
	{"any", (DL_FUNC) &any, 1},
	{"consolidateGaps", (DL_FUNC) &consolidateGaps, 2},
	{"composition", (DL_FUNC) &composition, 1},
	{"matchListsDual", (DL_FUNC) &matchListsDual, 5},
	{"findFrameshifts", (DL_FUNC) &findFrameshifts, 14},
	{"radixOrder", (DL_FUNC) &radixOrder, 2},
	{"fillOverlaps", (DL_FUNC) &fillOverlaps, 2},
	{"indexByContig", (DL_FUNC) &indexByContig, 5},
	{"chainSegments", (DL_FUNC) &chainSegments, 18},
	{"basicTranslate", (DL_FUNC) &basicTranslate, 3},
	{"firstSeqsGapsEqual", (DL_FUNC) &firstSeqsGapsEqual, 9},
	{"positionWeightMatrix", (DL_FUNC) &positionWeightMatrix, 4},
	{"extendSegments", (DL_FUNC) &extendSegments, 13},
	{"firstSeqsPosEqual", (DL_FUNC) &firstSeqsPosEqual, 9},
	{"collapse", (DL_FUNC) &collapse, 3},
	{"nbit", (DL_FUNC) &nbit, 4},
	{"decompress", (DL_FUNC) &decompress, 2},
	{"extractFields", (DL_FUNC) &extractFields, 4},
	{"intDiff", (DL_FUNC) &intDiff, 1},
	{"qbit", (DL_FUNC) &qbit, 3},
	{"movAvg", (DL_FUNC) &movAvg, 7},
	{"getPools", (DL_FUNC) &getPools, 1},
	{"predictDBN", (DL_FUNC) &predictDBN, 14},
	{"informationContent", (DL_FUNC) &informationContent, 4},
	{"informationContentAA", (DL_FUNC) &informationContentAA, 4},
	{"vectorSum", (DL_FUNC) &vectorSum, 4},
	{"parallelMatch", (DL_FUNC) &parallelMatch, 7},
	{"groupMax", (DL_FUNC) &groupMax, 3},
	{"removeGaps", (DL_FUNC) &removeGaps, 3},
	{"alphabetSize", (DL_FUNC) &alphabetSize, 1},
	{"alphabetSizeReducedAA", (DL_FUNC) &alphabetSizeReducedAA, 2},
	{"getORFs", (DL_FUNC) &getORFs, 5},
	{"codonModel", (DL_FUNC) &codonModel, 5},
	{"scoreCodonModel", (DL_FUNC) &scoreCodonModel, 3},
	{"dicodonModel", (DL_FUNC) &dicodonModel, 3},
	{"startCodonModel", (DL_FUNC) &startCodonModel, 4},
	{"scoreStartCodonModel", (DL_FUNC) &scoreStartCodonModel, 3},
	{"initialCodonModel", (DL_FUNC) &initialCodonModel, 4},
	{"scoreInitialCodonModel", (DL_FUNC) &scoreInitialCodonModel, 3},
	{"getRegion", (DL_FUNC) &getRegion, 5},
	{"autocorrelationModel", (DL_FUNC) &autocorrelationModel, 4},
	{"scoreAutocorrelationModel", (DL_FUNC) &scoreAutocorrelationModel, 4},
	{"nucleotideBiasModel", (DL_FUNC) &nucleotideBiasModel, 4},
	{"scoreNucleotideBiasModel", (DL_FUNC) &scoreNucleotideBiasModel, 3},
	{"upstreamMotifModel", (DL_FUNC) &upstreamMotifModel, 6},
	{"scoreUpstreamMotifModel", (DL_FUNC) &scoreUpstreamMotifModel, 6},
	{"runLengthModel", (DL_FUNC) &runLengthModel, 3},
	{"scoreRunLengthModel", (DL_FUNC) &scoreRunLengthModel, 4},
	{"stopCodonModel", (DL_FUNC) &stopCodonModel, 4},
	{"scoreStopCodonModel", (DL_FUNC) &scoreStopCodonModel, 3},
	{"terminationCodonModel", (DL_FUNC) &terminationCodonModel, 4},
	{"scoreTerminationCodonModel", (DL_FUNC) &scoreTerminationCodonModel, 3},
	{"codonFrequencies", (DL_FUNC) &codonFrequencies, 3},
	{"unicodonModel", (DL_FUNC) &unicodonModel, 3},
	{"chainGenes", (DL_FUNC) &chainGenes, 9},
	{"longestORFs", (DL_FUNC) &longestORFs, 1},
	{"getIndex", (DL_FUNC) &getIndex, 4},
	{"inBounds", (DL_FUNC) &inBounds, 5},
	{"getBounds", (DL_FUNC) &getBounds, 12},
	{"addIfElse", (DL_FUNC) &addIfElse, 3},
	{"kmerScores", (DL_FUNC) &kmerScores, 4},
	{"getHits", (DL_FUNC) &getHits, 7},
	{"allZero", (DL_FUNC) &allZero, 6},
	{"couplingModel", (DL_FUNC) &couplingModel, 5},
	{"scoreCouplingModel", (DL_FUNC) &scoreCouplingModel, 4},
	{"maxPerORF", (DL_FUNC) &maxPerORF, 2},
	{"replaceGaps", (DL_FUNC) &replaceGaps, 4},
	{"intMatchSelfOnce", (DL_FUNC) &intMatchSelfOnce, 2},
	{"scorePWM", (DL_FUNC) &scorePWM, 4},
	{NULL, NULL, 0}
};

void R_init_DECIPHER(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
