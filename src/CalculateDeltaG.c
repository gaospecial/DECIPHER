/****************************************************************************
 *               Calculates Delta G of Aligned Probe/Target Pairs           *
 *                           Author: Erik Wright                            *
 ****************************************************************************/

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

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("calculateDeltaG", probe, target, deltaGrules, PACKAGE="DECIPHER")
SEXP calculateDeltaG(SEXP p, SEXP t, SEXP deltaGrules)
{
	int i, j, n, count;
	int s1 = 0;
	int s2 = 0;
	int l = length(p);
	const char *seq1;
	const char *seq2;
	int dG[8];
	double dGov;
	double *rans, *dGrules = REAL(deltaGrules);
	SEXP ans;
	
	PROTECT(ans = allocVector(REALSXP, l));
	rans = REAL(ans);
	
	for (i = 0; i < l; i++) {
		seq1 = CHAR(STRING_ELT(p, i));
		seq2 = CHAR(STRING_ELT(t, i));
		n = length(STRING_ELT(p, i));
		count = 0;
		dGov = 0;
		rans[i] = 0;
		for (j = 0; j < n + 2; j++) {
			if (j < n) {
				switch (seq1[j]) {
					case '-':
						s1 = 4;
						break;
					case 'A':
					case 'a':
						s1 = 0;
						break;
					case 'C':
					case 'c':
						s1 = 1;
						break;
					case 'G':
					case 'g':
						s1 = 2;
						break;
					case 'T':
					case 't':
						s1 = 3;
						break;
					case 'M':
					case 'm':
						// M = A or C
						if (seq2[j]=='A' || seq2[j]=='a') {
							s1 = 0;
						} else { // default to C
							s1 = 1;
						}
						break;
					case 'R':
					case 'r':
						// R = A or G
						if (seq2[j]=='A' || seq2[j]=='a') {
							s1 = 0;
						} else { // default to G
							s1 = 2;
						}
						break;
					case 'S':
					case 's':
						// S = C or G
						if (seq2[j]=='C' || seq2[j]=='c') {
							s1 = 1;
						} else { // default to G
							s1 = 2;
						}
						break;
					case 'W':
					case 'w':
						// W = A or T
						if (seq2[j]=='A' || seq2[j]=='a') {
							s1 = 0;
						} else { // default to T
							s1 = 3;
						}
						break;
					case 'Y':
					case 'y':
						// Y = C or T
						if (seq2[j]=='C' || seq2[j]=='c') {
							s1 = 1;
						} else { // default to T
							s1 = 3;
						}
						break;
					case 'K':
					case 'k':
						// K = G or T
						if (seq2[j]=='G' || seq2[j]=='g') {
							s1 = 2;
						} else { // default to T
							s1 = 3;
						}
						break;
					case 'V':
					case 'v':
						// V = A or C or G
						if (seq2[j]=='A' || seq2[j]=='a') {
							s1 = 0;
						} else if (seq2[j]=='C' || seq2[j]=='c') {
							s1 = 1;
						} else { // default to G
							s1 = 2;
						}
						break;
					case 'H':
					case 'h':
						// H = A or C or T
						if (seq2[j]=='A' || seq2[j]=='a') {
							s1 = 0;
						} else if (seq2[j]=='C' || seq2[j]=='c') {
							s1 = 1;
						} else { // default to T
							s1 = 3;
						}
						break;
					case 'D':
					case 'd':
						// D = A or G or T
						if (seq2[j]=='A' || seq2[j]=='a') {
							s1 = 0;
						} else if (seq2[j]=='G' || seq2[j]=='g') {
							s1 = 2;
						} else { // default to T
							s1 = 3;
						}
						break;
					case 'B':
					case 'b':
						// B = C or G or T
						if (seq2[j]=='C' || seq2[j]=='c') {
							s1 = 1;
						} else if (seq2[j]=='G' || seq2[j]=='g') {
							s1 = 2;
						} else { // default to T
							s1 = 3;
						}
						break;
					case 'N':
					case 'n':
						// N = A or C or G or T
						if (seq2[j]=='A' || seq2[j]=='a') {
							s1 = 0;
						} else if (seq2[j]=='C' || seq2[j]=='c') {
							s1 = 1;
						} else if (seq2[j]=='G' || seq2[j]=='g') {
							s1 = 2;
						} else { // default to T
							s1 = 3;
						}
						break;
					default:
						error("Unrecognized character.");
						break;
				}
				switch (seq2[j]) {
					case '-':
						s2 = 4;
						break;
					case 'A':
					case 'a':
						s2 = 0;
						break;
					case 'C':
					case 'c':
						s2 = 1;
						break;
					case 'G':
					case 'g':
						s2 = 2;
						break;
					case 'T':
					case 't':
						s2 = 3;
						break;
					case 'M':
					case 'm':
						// M = A or C
						if (s1==0) {
							s2 = 0;
						} else { // default to C
							s2 = 1;
						}
						break;
					case 'R':
					case 'r':
						// R = A or G
						if (s1==0) {
							s2 = 0;
						} else { // default to G
							s2 = 2;
						}
						break;
					case 'S':
					case 's':
						// S = C or G
						if (s1==1) {
							s2 = 1;
						} else { // default to G
							s2 = 2;
						}
						break;
					case 'W':
					case 'w':
						// W = A or T
						if (s1==0) {
							s2 = 0;
						} else { // default to T
							s2 = 3;
						}
						break;
					case 'Y':
					case 'y':
						// Y = C or T
						if (s1==1) {
							s2 = 1;
						} else { // default to T
							s2 = 3;
						}
						break;
					case 'K':
					case 'k':
						// K = G or T
						if (s1==2) {
							s2 = 2;
						} else { // default to T
							s2 = 3;
						}
						break;
					case 'V':
					case 'v':
						// V = A or C or G
						if (s1==0) {
							s2 = 0;
						} else if (s1==1) {
							s2 = 1;
						} else { // default to G
							s2 = 2;
						}
						break;
					case 'H':
					case 'h':
						// H = A or C or T
						if (s1==0) {
							s2 = 0;
						} else if (s1==1) {
							s2 = 1;
						} else { // default to T
							s2 = 3;
						}
						break;
					case 'D':
					case 'd':
						// D = A or G or T
						if (s1==0) {
							s2 = 0;
						} else if (s1==2) {
							s2 = 2;
						} else { // default to T
							s2 = 3;
						}
						break;
					case 'B':
					case 'b':
						// B = C or G or T
						if (s1==1) {
							s2 = 1;
						} else if (s1==2) {
							s2 = 2;
						} else { // default to T
							s2 = 3;
						}
						break;
					case 'N':
					case 'n':
						// N = A or C or G or T
						if (s1==0) {
							s2 = 0;
						} else if (s1==1) {
							s2 = 1;
						} else if (s1==2) {
							s2 = 2;
						} else {
							s2 = 3;
						}
						break;
					default:
						error("Unrecognized character.");
						break;
				}
				if (!((s1 == 4) && (s2 == 4))) { // not a common gap
					if (count > 3) { // rotate dG array
						dG[0] = dG[1];
						dG[1] = dG[2];
						dG[2] = dG[3];
						dG[3] = s1;
						dG[4] = dG[5];
						dG[5] = dG[6];
						dG[6] = dG[7];
						dG[7] = s2;
					} else { // fill dG array
						dG[count] = s1;
						dG[count + 4] = s2;
						count++;
					}
				}
			} else {
				// set dG array to gap
				if (count > 3) { // rotate dG array
					dG[0] = dG[1];
					dG[1] = dG[2];
					dG[2] = dG[3];
					dG[3] = 4;
					dG[4] = dG[5];
					dG[5] = dG[6];
					dG[6] = dG[7];
					dG[7] = 4;
				}
			}
			
			if ((!((s1 == 4) && (s2 == 4)) || j >= n) && count > 3) {
				dGov += *(dGrules + dG[0] + dG[1]*5 + dG[2]*25 + dG[3]*125 + dG[4]*625 + dG[5]*3125 + dG[6]*15625 + dG[7]*78125);
				
				if (dGov > 0)
					dGov = 0;
				
				// find minima
				if (dGov < rans[i])
					rans[i] = dGov;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;	
}

SEXP calculateHairpinDeltaG(SEXP x, SEXP arms, SEXP deltaGrules)
{
	int i, j, k, count, s1, s2;
	int *a = INTEGER(arms);
	int dG[8];
	double dGov;
	double *dGrules = REAL(deltaGrules);
	
	// initialize the XStringSet
	XStringSet_holder x_set;
	Chars_holder x_i;
	x_set = hold_XStringSet(x);
	int x_length = get_length_from_XStringSet_holder(&x_set);
	
	SEXP ans;
	PROTECT(ans = allocVector(REALSXP, x_length));
	double *rans = REAL(ans);
	
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		count = 0;
		dGov = 0;
		rans[i] = 0;
		k = x_i.length - 1; // other end of sequence
		for (j = 0; j < a[i] + 2; j++) {
			if (j < a[i]) { // within hairpin
				switch (x_i.ptr[j]) {
					case 1: // A
						s1 = 0;
						break;
					case 2: // C
						s1 = 1;
						break;
					case 4: // G
						s1 = 2;
						break;
					case 8: // T/U
						s1 = 3;
						break;
					default: // treat others as gap
						s1 = 4;
						break;
				}
				switch (x_i.ptr[k--]) {
					case 1: // A
						s2 = 3;
						break;
					case 2: // C
						s2 = 2;
						break;
					case 4: // G
						s2 = 1;
						break;
					case 8: // T/U
						s2 = 0;
						break;
					default: // treat others as gap
						s2 = 4;
						break;
				}
				if (!((s1 == 4) && (s2 == 4))) { // not a common gap
					if (count > 3) { // rotate dG array
						dG[0] = dG[1];
						dG[1] = dG[2];
						dG[2] = dG[3];
						dG[3] = s1;
						dG[4] = dG[5];
						dG[5] = dG[6];
						dG[6] = dG[7];
						dG[7] = s2;
					} else { // fill dG array
						dG[count] = s1;
						dG[count + 4] = s2;
						count++;
					}
				}
			} else {
				// set dG array to gap
				if (count > 3) { // rotate dG array
					dG[0] = dG[1];
					dG[1] = dG[2];
					dG[2] = dG[3];
					dG[3] = 4;
					dG[4] = dG[5];
					dG[5] = dG[6];
					dG[6] = dG[7];
					dG[7] = 4;
				}
			}
			
			if ((!((s1 == 4) && (s2 == 4)) || j >= a[i]) && count > 3) {
				dGov += *(dGrules + dG[0] + dG[1]*5 + dG[2]*25 + dG[3]*125 + dG[4]*625 + dG[5]*3125 + dG[6]*15625 + dG[7]*78125);
				
				if (dGov > 0)
					dGov = 0;
				
				// find minima
				if (dGov < rans[i])
					rans[i] = dGov;
			}
		}
	}
	
	UNPROTECT(1);
	
	return ans;	
}
