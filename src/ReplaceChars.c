/****************************************************************************
 *                        Quickly Replace Characters                        *
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

// strcpy
#include <string.h>

/* for Calloc/Free */
#include <R_ext/RS.h>

// for math functions
#include <math.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"
#include "XVector_interface.h"
#include "S4Vectors_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

//ans_start <- .Call("replaceChars", sequences, replaceChar, type, PACKAGE="DECIPHER")
SEXP replaceChars(SEXP x, SEXP r, SEXP t)
{
	int i, j, l, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	const char *repChar = CHAR(STRING_ELT(r, 0));
	int fail = (STRING_ELT(r, 0)==NA_STRING);
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char *s = Calloc(longest + 1, char); // each sequence
	
	// write new character vector
	if (asInteger(t)==1) {
		for (i = 0; i < n; i++) {
			l = length(STRING_ELT(x, i));
			seq = CHAR(STRING_ELT(x, i));
			count = 0;
			for (j = 0; j < l; j++) {
				if (seq[j]!='U' && seq[j]!='u') {
					switch (seq[j]) {
						case '-':
						case 'A':
						case 'a':
						case 'C':
						case 'c':
						case 'G':
						case 'g':
						case 'T':
						case 't':
						case 'N':
						case 'n':
						case 'M':
						case 'm':
						case 'R':
						case 'r':
						case 'W':
						case 'w':
						case 'S':
						case 's':
						case 'Y':
						case 'y':
						case 'K':
						case 'k':
						case 'V':
						case 'v':
						case 'H':
						case 'h':
						case 'D':
						case 'd':
						case 'B':
						case 'b':
						case '+':
						case '.':
							s[count] = seq[j];
							count++;
							break;
						default:
							if (fail) {
								error("Incompatible character ('%c') found in DNAStringSet when replaceChar = NA.", seq[j]);
							} else if (repChar[0] != '\0') {
								s[count] = repChar[0];
								count++;
							}
							break;
					}
				} else {
					s[count] = 'T';
					count++;
				}
			}
			s[count] = '\0'; // null-terminate
			SET_STRING_ELT(seqs, i, mkChar(s));
		}
	} else if (asInteger(t)==2) {
		for (i = 0; i < n; i++) {
			l = length(STRING_ELT(x, i));
			seq = CHAR(STRING_ELT(x, i));
			count = 0;
			for (j = 0; j < l; j++) {
				if (seq[j]!='T' && seq[j]!='t') {
					switch (seq[j]) {
						case '-':
						case 'A':
						case 'a':
						case 'C':
						case 'c':
						case 'G':
						case 'g':
						case 'U':
						case 'u':
						case 'N':
						case 'n':
						case 'M':
						case 'm':
						case 'R':
						case 'r':
						case 'W':
						case 'w':
						case 'S':
						case 's':
						case 'Y':
						case 'y':
						case 'K':
						case 'k':
						case 'V':
						case 'v':
						case 'H':
						case 'h':
						case 'D':
						case 'd':
						case 'B':
						case 'b':
						case '+':
						case '.':
							s[count] = seq[j];
							count++;
							break;
						default:
							if (fail) {
								error("Incompatible character ('%c') in RNAStringSet found when replaceChar = NA.", seq[j]);
							} else if (repChar[0] != '\0') {
								s[count] = repChar[0];
								count++;
							}
							break;
					}
				} else {
					s[count] = 'U';
					count++;
				}
			}
			s[count] = '\0'; // null-terminate
			SET_STRING_ELT(seqs, i, mkChar(s));
		}
	} else {
		for (i = 0; i < n; i++) {
			l = length(STRING_ELT(x, i));
			seq = CHAR(STRING_ELT(x, i));
			count = 0;
			for (j = 0; j < l; j++) {
				switch (seq[j]) {
					case '-':
					case 'A':
					case 'R':
					case 'N':
					case 'D':
					case 'C':
					case 'Q':
					case 'E':
					case 'G':
					case 'H':
					case 'I':
					case 'L':
					case 'K':
					case 'M':
					case 'F':
					case 'P':
					case 'S':
					case 'T':
					case 'W':
					case 'Y':
					case 'V':
					case 'U':
					case 'O':
					case 'B':
					case 'Z':
					case 'X':
					case '*':
					case '+':
					case '.':
						s[count] = seq[j];
						count++;
						break;
					case 'a':
						s[count] = 'A';
						count++;
						break;
					case 'r':
						s[count] = 'R';
						count++;
						break;
					case 'n':
						s[count] = 'N';
						count++;
						break;
					case 'd':
						s[count] = 'D';
						count++;
						break;
					case 'c':
						s[count] = 'C';
						count++;
						break;
					case 'q':
						s[count] = 'Q';
						count++;
						break;
					case 'e':
						s[count] = 'E';
						count++;
						break;
					case 'g':
						s[count] = 'G';
						count++;
						break;
					case 'h':
						s[count] = 'H';
						count++;
						break;
					case 'i':
						s[count] = 'I';
						count++;
						break;
					case 'l':
						s[count] = 'L';
						count++;
						break;
					case 'k':
						s[count] = 'K';
						count++;
						break;
					case 'm':
						s[count] = 'M';
						count++;
						break;
					case 'f':
						s[count] = 'F';
						count++;
						break;
					case 'p':
						s[count] = 'P';
						count++;
						break;
					case 's':
						s[count] = 'S';
						count++;
						break;
					case 't':
						s[count] = 'T';
						count++;
						break;
					case 'w':
						s[count] = 'W';
						count++;
						break;
					case 'y':
						s[count] = 'Y';
						count++;
						break;
					case 'v':
						s[count] = 'V';
						count++;
						break;
					case 'u':
						s[count] = 'U';
						count++;
						break;
					case 'o':
						s[count] = 'O';
						count++;
						break;
					case 'b':
						s[count] = 'B';
						count++;
						break;
					case 'j':
						s[count] = 'J';
						count++;
						break;
					case 'z':
						s[count] = 'Z';
						count++;
						break;
					case 'x':
						s[count] = 'X';
						count++;
						break;
					default:
						if (fail) {
							error("Incompatible character ('%c') in AAStringSet found when replaceChar = NA.", seq[j]);
						} else if (repChar[0] != '\0') {
							s[count] = repChar[0];
							count++;
						}
						break;
				}
			}
			s[count] = '\0'; // null-terminate
			SET_STRING_ELT(seqs, i, mkChar(s));
		}
	}
	
	Free(s);
	
	UNPROTECT(1);
	
	return seqs;
}

//ans_start <- .Call("replaceChar", sequences, charReplace, replaceChar, PACKAGE="DECIPHER")
SEXP replaceChar(SEXP x, SEXP c, SEXP r)
{
	int i, j, l, count;
	int n = length(x);
	int longest = 0;
	const char *seq;
	const char *repChar = CHAR(STRING_ELT(r, 0));
	const char *charRep = CHAR(STRING_ELT(c, 0));
	
	// find longest string
	for (i = 0; i < n; i++)
		if (length(STRING_ELT(x, i)) > longest)
			longest = length(STRING_ELT(x, i));
	
	SEXP seqs;
	PROTECT(seqs = allocVector(STRSXP, n));
	char *s = Calloc(longest + 1, char); // each sequence
	
	// write new character vector
	for (i = 0; i < n; i++) {
		l = length(STRING_ELT(x, i));
		seq = CHAR(STRING_ELT(x, i));
		count = 0;
		for (j = 0; j < l; j++) {
			if (seq[j]==charRep[0]) {
				if (repChar[0] != '\0') {
					s[count] = repChar[0];
					count++;
				}
			} else {
				s[count] = seq[j];
				count++;
			}
		}
		s[count] = '\0'; // null-terminate
		SET_STRING_ELT(seqs, i, mkChar(s));
	}
	
	Free(s);
	
	UNPROTECT(1);
	
	return seqs;
}

//ans_start <- .Call("replaceGaps", sequences, newsequence, start, type, PACKAGE="DECIPHER")
SEXP replaceGaps(SEXP x, SEXP y, SEXP start, SEXP type)
{
	int i, j, x_length;
	SEXP ans_width, ans;
	int s = asInteger(start) - 1; // starting position in y
	int t = asInteger(type);
	
	// determine the element type of the XStringSet
	const char *ans_element_type;
	ans_element_type = get_List_elementType(x);
	
	// determine the length of the XStringSet
	Chars_holder x_i, y_holder, ans_elt_holder;
	XStringSet_holder x_set, ans_holder;
	y_holder = hold_XRaw(y);
	x_set = hold_XStringSet(x);
	x_length = get_length_from_XStringSet_holder(&x_set);
	
	// count the sequence lengths
	PROTECT(ans_width = NEW_INTEGER(x_length));
	int *width = INTEGER(ans_width);
	for (i = 0; i < x_length; i++) {
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		width[i] = x_i.length;
	}
	
	// set the class of the XStringSet
	char ans_classname[40];
	if (t==1) {
		strcpy(ans_classname, "DNAStringSet");
	} else if (t==2) {
		strcpy(ans_classname, "RNAStringSet");
	} else { // t==3
		strcpy(ans_classname, "AAStringSet");
	}
	
	// initialize a new XStringSet
	PROTECT(ans = alloc_XRawList(ans_classname, ans_element_type, ans_width));
	ans_holder = hold_XVectorList(ans);
	
	for (i = 0; i < x_length; i++) {
		ans_elt_holder = get_elt_from_XStringSet_holder(&ans_holder, i);
		x_i = get_elt_from_XStringSet_holder(&x_set, i);
		if (t==3) { // AAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (!(x_i.ptr[j] ^ 0x2D) || !(x_i.ptr[j] ^ 0x2E)) { // position is a gap
					memcpy((char *) ans_elt_holder.ptr + j, x_i.ptr + j, 1 * sizeof(char));
				} else {
					memcpy((char *) ans_elt_holder.ptr + j, y_holder.ptr + s, 1 * sizeof(char));
					s++;
				}
			}
		} else { // DNAStringSet or RNAStringSet
			for (j = 0; j < x_i.length; j++) {
				if (x_i.ptr[j] & 0x10 || x_i.ptr[j] & 0x40) { // position is a gap
					memcpy((char *) ans_elt_holder.ptr + j, x_i.ptr + j, 1 * sizeof(char));
				} else {
					memcpy((char *) ans_elt_holder.ptr + j, y_holder.ptr + s, 1 * sizeof(char));
					s++;
				}
			}
		}
	}
	
	UNPROTECT(2);
	return ans;
}
