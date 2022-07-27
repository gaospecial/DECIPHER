/****************************************************************************
 *                        Cluster Maximum Likelihood                        *
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

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

/* for Calloc/Free */
#include <R_ext/RS.h>

// for math functions
#include <math.h>

// for OpenMP parallel processing
#ifdef SUPPORT_OPENMP
#include <omp.h>
#endif

// for calloc/free
#include <stdlib.h>

// for floating point limits
#include <float.h>

/*
 * Biostrings_interface.h is needed for the DNAencode(), get_XString_asRoSeq(),
 * init_match_reporting(), report_match() and reported_matches_asSEXP()
 * protoypes, and for the COUNT_MRMODE and START_MRMODE constant symbols.
 */
#include "Biostrings_interface.h"

// DECIPHER header file
#include "DECIPHER.h"

// for SIMD
//#if defined(__AVX2__) && __AVX2__
//#define _L_unknown_5 L_unknown_5_SIMD
//#define _L_unknown L_unknown_SIMD
//#define _L_unknown_AA_5 L_unknown_AA_5
//#define _L_unknown_AA L_unknown_AA
//#else
#define _L_unknown_5 L_unknown_5
#define _L_unknown L_unknown
#define _L_unknown_AA_5 L_unknown_AA_5
#define _L_unknown_AA L_unknown_AA
//#endif

static void L_known(const char *p, double *Ls, int gaps)
{
	switch (*p) {
		case 1: // A
			*(Ls) = 1;
			break;
		case 2: // C
			*(Ls + 1) = 1;
			break;
		case 3: // M
			*(Ls) = 1/(double)2; *(Ls + 1) = 1/(double)2; // AC
			break;
		case 4: // G
			*(Ls + 2) = 1;
			break;
		case 5: // R
			*(Ls) = 1/(double)2; *(Ls + 2) = 1/(double)2; // AG
			break;
		case 6: // S
			*(Ls + 1) = 1/(double)2; *(Ls + 2) = 1/(double)2; // CG
			break;
		case 7: // V
			*(Ls) = 1/(double)3; *(Ls + 1) = 1/(double)3; *(Ls + 2) = 1/(double)3; // ACG
			break;
		case 8: // T
			*(Ls + 3) = 1;
			break;
		case 9: // W
			*(Ls) = 1/(double)2; *(Ls + 3) = 1/(double)2; // AT
			break;
		case 10: // Y
			*(Ls + 1) = 1/(double)2; *(Ls + 3) = 1/(double)2; // CT
			break;
		case 11: // H
			*(Ls) = 1/(double)3; *(Ls + 1) = 1/(double)3; *(Ls + 3) = 1/(double)3; // ACT
			break;
		case 12: // K
			*(Ls + 2) = 1/(double)2; *(Ls + 3) = 1/(double)2; // GT
			break;
		case 13: // D
			*(Ls) = 1/(double)3; *(Ls + 2) = 1/(double)3; *(Ls + 3) = 1/(double)3; // AGT
			break;
		case 14: // B
			*(Ls + 1) = 1/(double)3; *(Ls + 2) = 1/(double)3; *(Ls + 3) = 1/(double)3; // CGT
			break;
		case 15: // N
			// ACGT (unknown)
			break;
		case 16: // -
			if (gaps)
				*(Ls + 4) = 1;
			break;
		case 32: // +
			// mask
			break;
		case 64: // .
			if (gaps)
				*(Ls + 4) = 1;
			break;
		default:
			error("not nucleotides!");
			break;
	}
}

static void L_known_AA(const char *p, double *Ls, int gaps)
{
	switch (*p) {
		case 65: // A
			*(Ls) = 1;
			break;
		case 82: // R
			*(Ls + 1) = 1;
			break;
		case 78: // N
			*(Ls + 2) = 1;
			break;
		case 68: // D
			*(Ls + 3) = 1;
			break;
		case 67: // C
			*(Ls + 4) = 1;
			break;
		case 81: // Q
			*(Ls + 5) = 1;
			break;
		case 69: // E
			*(Ls + 6) = 1;
			break;
		case 71: // G
			*(Ls + 7) = 1;
			break;
		case 72: // H
			*(Ls + 8) = 1;
			break;
		case 73: // I
			*(Ls + 9) = 1;
			break;
		case 76: // L
			*(Ls + 10) = 1;
			break;
		case 75: // K
			*(Ls + 11) = 1;
			break;
		case 77: // M
			*(Ls + 12) = 1;
			break;
		case 70: // F
			*(Ls + 13) = 1;
			break;
		case 80: // P
			*(Ls + 14) = 1;
			break;
		case 83: // S
			*(Ls + 15) = 1;
			break;
		case 84: // T
			*(Ls + 16) = 1;
			break;
		case 87: // W
			*(Ls + 17) = 1;
			break;
		case 89: // Y
			*(Ls + 18) = 1;
			break;
		case 86: // V
			*(Ls + 19) = 1;
			break;
		case 85: // U
			// ignore
			break;
		case 79: // O
			// ignore
			break;
		case 66: // B = N or D
			*(Ls + 2) = 1/(double)2; *(Ls + 3) = 1/(double)2;
			break;
		case 90: // Z = Q or E
			*(Ls + 5) = 1/(double)2; *(Ls + 6) = 1/(double)2;
			break;
		case 74: // J = I or L
			*(Ls + 9) = 1/(double)2; *(Ls + 10) = 1/(double)2;
			break;
		case 88: // X = any letter
			// ignore
			break;
		case 42: // * (stop)
			// ignore
			break;
		case 45: // -
			if (gaps)
				*(Ls + 20) = 1;
			break;
		case 43: // +
			// mask
			break;
		case 46: // .
			if (gaps)
				*(Ls + 20) = 1;
			break;
		default:
			error("not AA!");
			break;
	}
}

static void L_unknown(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	double L1, L2;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lc1 = *(Ls1 + 1);
	const double Lg1 = *(Ls1 + 2);
	const double Lt1 = *(Ls1 + 3);
	const double La2 = *(Ls2);
	const double Lc2 = *(Ls2 + 1);
	const double Lg2 = *(Ls2 + 2);
	const double Lt2 = *(Ls2 + 3);
	
	if (root == 0 &&
		(La1 != 0 ||
		Lc1 != 0 ||
		Lg1 != 0 ||
		Lt1 != 0)) {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0) {
			// neither branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 5)*(Lc1); // Lac
			L1 += *(P1 + 10)*(Lg1); // Lag
			L1 += *(P1 + 15)*(Lt1); // Lat
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 5)*(Lc2); // Lac
			L2 += *(P2 + 10)*(Lg2); // Lag
			L2 += *(P2 + 15)*(Lt2); // Lat
			*(Ls3) = L1*L2;
			
			// L(C)
			L1 = *(P1 + 1)*(La1); // Lca
			L1 += *(P1 + 6)*(Lc1); // Lcc
			L1 += *(P1 + 11)*(Lg1); // Lcg
			L1 += *(P1 + 16)*(Lt1); // Lct
			L2 = *(P2 + 1)*(La2); // Lca
			L2 += *(P2 + 6)*(Lc2); // Lcc
			L2 += *(P2 + 11)*(Lg2); // Lcg
			L2 += *(P2 + 16)*(Lt2); // Lct
			*(Ls3 + 1) = L1*L2;
			
			// L(G)
			L1 = *(P1 + 2)*(La1); // Lga
			L1 += *(P1 + 7)*(Lc1); // Lgc
			L1 += *(P1 + 12)*(Lg1); // Lgg
			L1 += *(P1 + 17)*(Lt1); // Lgt
			L2 = *(P2 + 2)*(La2); // Lga
			L2 += *(P2 + 7)*(Lc2); // Lgc
			L2 += *(P2 + 12)*(Lg2); // Lgg
			L2 += *(P2 + 17)*(Lt2); // Lgt
			*(Ls3 + 2) = L1*L2;
			
			// L(T)
			L1 = *(P1 + 3)*(La1); // Lta
			L1 += *(P1 + 8)*(Lc1); // Ltc
			L1 += *(P1 + 13)*(Lg1); // Ltg
			L1 += *(P1 + 18)*(Lt1); // Ltt
			L2 = *(P2 + 3)*(La2); // Lta
			L2 += *(P2 + 8)*(Lc2); // Ltc
			L2 += *(P2 + 13)*(Lg2); // Ltg
			L2 += *(P2 + 18)*(Lt2); // Ltt
			*(Ls3 + 3) = L1*L2;
			
			*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 5) += 1;
			}
		} else {
			// second branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 5)*(Lc1); // Lac
			L1 += *(P1 + 10)*(Lg1); // Lag
			L1 += *(P1 + 15)*(Lt1); // Lat
			*(Ls3) = L1;
			
			// L(C)
			L1 = *(P1 + 1)*(La1); // Lca
			L1 += *(P1 + 6)*(Lc1); // Lcc
			L1 += *(P1 + 11)*(Lg1); // Lcg
			L1 += *(P1 + 16)*(Lt1); // Lct
			*(Ls3 + 1) = L1;
			
			// L(G)
			L1 = *(P1 + 2)*(La1); // Lga
			L1 += *(P1 + 7)*(Lc1); // Lgc
			L1 += *(P1 + 12)*(Lg1); // Lgg
			L1 += *(P1 + 17)*(Lt1); // Lgt
			*(Ls3 + 2) = L1;
			
			// L(T)
			L1 = *(P1 + 3)*(La1); // Lta
			L1 += *(P1 + 8)*(Lc1); // Ltc
			L1 += *(P1 + 13)*(Lg1); // Ltg
			L1 += *(P1 + 18)*(Lt1); // Ltt
			*(Ls3 + 3) = L1;
			
			*(Ls3 + 5) = *(Ls1 + 5);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 5) += 1;
			}
		}
	} else {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0) {
			// first branch can be disregarded
			
			// L(A)
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 5)*(Lc2); // Lac
			L2 += *(P2 + 10)*(Lg2); // Lag
			L2 += *(P2 + 15)*(Lt2); // Lat
			*(Ls3) = L2;
			
			// L(C)
			L2 = *(P2 + 1)*(La2); // Lca
			L2 += *(P2 + 6)*(Lc2); // Lcc
			L2 += *(P2 + 11)*(Lg2); // Lcg
			L2 += *(P2 + 16)*(Lt2); // Lct
			*(Ls3 + 1) = L2;
			
			// L(G)
			L2 = *(P2 + 2)*(La2); // Lga
			L2 += *(P2 + 7)*(Lc2); // Lgc
			L2 += *(P2 + 12)*(Lg2); // Lgg
			L2 += *(P2 + 17)*(Lt2); // Lgt
			*(Ls3 + 2) = L2;
			
			// L(T)
			L2 = *(P2 + 3)*(La2); // Lta
			L2 += *(P2 + 8)*(Lc2); // Ltc
			L2 += *(P2 + 13)*(Lg2); // Ltg
			L2 += *(P2 + 18)*(Lt2); // Ltt
			*(Ls3 + 3) = L2;
			
			if (root &&
				(La1 != 0 ||
				Lc1 != 0 ||
				Lg1 != 0 ||
				Lt1 != 0)) {
				*(Ls3) *= La1;
				*(Ls3 + 1) *= Lc1;
				*(Ls3 + 2) *= Lg1;
				*(Ls3 + 3) *= Lt1;
				*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
			} else {
				*(Ls3 + 5) = *(Ls2 + 5);
			}
			
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 5) += 1;
			}
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lc1;
			*(Ls3 + 2) = Lg1;
			*(Ls3 + 3) = Lt1;
			*(Ls3 + 5) = *(Ls1 + 5);
		}
	}
	*(Ls3 + 4) = 0;
}

static void L_unknown_AA(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	double L1, L2;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lr1 = *(Ls1 + 1);
	const double Ln1 = *(Ls1 + 2);
	const double Ld1 = *(Ls1 + 3);
	const double Lc1 = *(Ls1 + 4);
	const double Lq1 = *(Ls1 + 5);
	const double Le1 = *(Ls1 + 6);
	const double Lg1 = *(Ls1 + 7);
	const double Lh1 = *(Ls1 + 8);
	const double Li1 = *(Ls1 + 9);
	const double Ll1 = *(Ls1 + 10);
	const double Lk1 = *(Ls1 + 11);
	const double Lm1 = *(Ls1 + 12);
	const double Lf1 = *(Ls1 + 13);
	const double Lp1 = *(Ls1 + 14);
	const double Lz1 = *(Ls1 + 15);
	const double Lt1 = *(Ls1 + 16);
	const double Lw1 = *(Ls1 + 17);
	const double Ly1 = *(Ls1 + 18);
	const double Lv1 = *(Ls1 + 19);
	const double La2 = *(Ls2);
	const double Lr2 = *(Ls2 + 1);
	const double Ln2 = *(Ls2 + 2);
	const double Ld2 = *(Ls2 + 3);
	const double Lc2 = *(Ls2 + 4);
	const double Lq2 = *(Ls2 + 5);
	const double Le2 = *(Ls2 + 6);
	const double Lg2 = *(Ls2 + 7);
	const double Lh2 = *(Ls2 + 8);
	const double Li2 = *(Ls2 + 9);
	const double Ll2 = *(Ls2 + 10);
	const double Lk2 = *(Ls2 + 11);
	const double Lm2 = *(Ls2 + 12);
	const double Lf2 = *(Ls2 + 13);
	const double Lp2 = *(Ls2 + 14);
	const double Lz2 = *(Ls2 + 15);
	const double Lt2 = *(Ls2 + 16);
	const double Lw2 = *(Ls2 + 17);
	const double Ly2 = *(Ls2 + 18);
	const double Lv2 = *(Ls2 + 19);
	
	if (root == 0 &&
		(La1 != 0 ||
		Lr1 != 0 ||
		Ln1 != 0 ||
		Ld1 != 0 ||
		Lc1 != 0 ||
		Lq1 != 0 ||
		Le1 != 0 ||
		Lg1 != 0 ||
		Lh1 != 0 ||
		Li1 != 0 ||
		Ll1 != 0 ||
		Lk1 != 0 ||
		Lm1 != 0 ||
		Lf1 != 0 ||
		Lp1 != 0 ||
		Lz1 != 0 ||
		Lt1 != 0 ||
		Lw1 != 0 ||
		Ly1 != 0 ||
		Lv1 != 0)) {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0) {
			// neither branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 21)*(Lr1); // Lar
			L1 += *(P1 + 42)*(Ln1); // Lan
			L1 += *(P1 + 63)*(Ld1); // Lad
			L1 += *(P1 + 84)*(Lc1); // Lac
			L1 += *(P1 + 105)*(Lq1); // Laq
			L1 += *(P1 + 126)*(Le1); // Lae
			L1 += *(P1 + 147)*(Lg1); // Lag
			L1 += *(P1 + 168)*(Lh1); // Lah
			L1 += *(P1 + 189)*(Li1); // Lai
			L1 += *(P1 + 210)*(Ll1); // Lal
			L1 += *(P1 + 231)*(Lk1); // Lak
			L1 += *(P1 + 252)*(Lm1); // Lam
			L1 += *(P1 + 273)*(Lf1); // Laf
			L1 += *(P1 + 294)*(Lp1); // Lap
			L1 += *(P1 + 315)*(Lz1); // Las
			L1 += *(P1 + 336)*(Lt1); // Lat
			L1 += *(P1 + 357)*(Lw1); // Law
			L1 += *(P1 + 378)*(Ly1); // Lay
			L1 += *(P1 + 399)*(Lv1); // Lav
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 21)*(Lr2); // Lar
			L2 += *(P2 + 42)*(Ln2); // Lan
			L2 += *(P2 + 63)*(Ld2); // Lad
			L2 += *(P2 + 84)*(Lc2); // Lac
			L2 += *(P2 + 105)*(Lq2); // Laq
			L2 += *(P2 + 126)*(Le2); // Lae
			L2 += *(P2 + 147)*(Lg2); // Lag
			L2 += *(P2 + 168)*(Lh2); // Lah
			L2 += *(P2 + 189)*(Li2); // Lai
			L2 += *(P2 + 210)*(Ll2); // Lal
			L2 += *(P2 + 231)*(Lk2); // Lak
			L2 += *(P2 + 252)*(Lm2); // Lam
			L2 += *(P2 + 273)*(Lf2); // Laf
			L2 += *(P2 + 294)*(Lp2); // Lap
			L2 += *(P2 + 315)*(Lz2); // Las
			L2 += *(P2 + 336)*(Lt2); // Lat
			L2 += *(P2 + 357)*(Lw2); // Law
			L2 += *(P2 + 378)*(Ly2); // Lay
			L2 += *(P2 + 399)*(Lv2); // Lav
			*(Ls3 + 0) = L1*L2;
			
			// L(R)
			L1 = *(P1 + 1)*(La1); // Lra
			L1 += *(P1 + 22)*(Lr1); // Lrr
			L1 += *(P1 + 43)*(Ln1); // Lrn
			L1 += *(P1 + 64)*(Ld1); // Lrd
			L1 += *(P1 + 85)*(Lc1); // Lrc
			L1 += *(P1 + 106)*(Lq1); // Lrq
			L1 += *(P1 + 127)*(Le1); // Lre
			L1 += *(P1 + 148)*(Lg1); // Lrg
			L1 += *(P1 + 169)*(Lh1); // Lrh
			L1 += *(P1 + 190)*(Li1); // Lri
			L1 += *(P1 + 211)*(Ll1); // Lrl
			L1 += *(P1 + 232)*(Lk1); // Lrk
			L1 += *(P1 + 253)*(Lm1); // Lrm
			L1 += *(P1 + 274)*(Lf1); // Lrf
			L1 += *(P1 + 295)*(Lp1); // Lrp
			L1 += *(P1 + 316)*(Lz1); // Lrs
			L1 += *(P1 + 337)*(Lt1); // Lrt
			L1 += *(P1 + 358)*(Lw1); // Lrw
			L1 += *(P1 + 379)*(Ly1); // Lry
			L1 += *(P1 + 400)*(Lv1); // Lrv
			L2 = *(P2 + 1)*(La2); // Lra
			L2 += *(P2 + 22)*(Lr2); // Lrr
			L2 += *(P2 + 43)*(Ln2); // Lrn
			L2 += *(P2 + 64)*(Ld2); // Lrd
			L2 += *(P2 + 85)*(Lc2); // Lrc
			L2 += *(P2 + 106)*(Lq2); // Lrq
			L2 += *(P2 + 127)*(Le2); // Lre
			L2 += *(P2 + 148)*(Lg2); // Lrg
			L2 += *(P2 + 169)*(Lh2); // Lrh
			L2 += *(P2 + 190)*(Li2); // Lri
			L2 += *(P2 + 211)*(Ll2); // Lrl
			L2 += *(P2 + 232)*(Lk2); // Lrk
			L2 += *(P2 + 253)*(Lm2); // Lrm
			L2 += *(P2 + 274)*(Lf2); // Lrf
			L2 += *(P2 + 295)*(Lp2); // Lrp
			L2 += *(P2 + 316)*(Lz2); // Lrs
			L2 += *(P2 + 337)*(Lt2); // Lrt
			L2 += *(P2 + 358)*(Lw2); // Lrw
			L2 += *(P2 + 379)*(Ly2); // Lry
			L2 += *(P2 + 400)*(Lv2); // Lrv
			*(Ls3 + 1) = L1*L2;
			
			// L(N)
			L1 = *(P1 + 2)*(La1); // Lna
			L1 += *(P1 + 23)*(Lr1); // Lnr
			L1 += *(P1 + 44)*(Ln1); // Lnn
			L1 += *(P1 + 65)*(Ld1); // Lnd
			L1 += *(P1 + 86)*(Lc1); // Lnc
			L1 += *(P1 + 107)*(Lq1); // Lnq
			L1 += *(P1 + 128)*(Le1); // Lne
			L1 += *(P1 + 149)*(Lg1); // Lng
			L1 += *(P1 + 170)*(Lh1); // Lnh
			L1 += *(P1 + 191)*(Li1); // Lni
			L1 += *(P1 + 212)*(Ll1); // Lnl
			L1 += *(P1 + 233)*(Lk1); // Lnk
			L1 += *(P1 + 254)*(Lm1); // Lnm
			L1 += *(P1 + 275)*(Lf1); // Lnf
			L1 += *(P1 + 296)*(Lp1); // Lnp
			L1 += *(P1 + 317)*(Lz1); // Lns
			L1 += *(P1 + 338)*(Lt1); // Lnt
			L1 += *(P1 + 359)*(Lw1); // Lnw
			L1 += *(P1 + 380)*(Ly1); // Lny
			L1 += *(P1 + 401)*(Lv1); // Lnv
			L2 = *(P2 + 2)*(La2); // Lna
			L2 += *(P2 + 23)*(Lr2); // Lnr
			L2 += *(P2 + 44)*(Ln2); // Lnn
			L2 += *(P2 + 65)*(Ld2); // Lnd
			L2 += *(P2 + 86)*(Lc2); // Lnc
			L2 += *(P2 + 107)*(Lq2); // Lnq
			L2 += *(P2 + 128)*(Le2); // Lne
			L2 += *(P2 + 149)*(Lg2); // Lng
			L2 += *(P2 + 170)*(Lh2); // Lnh
			L2 += *(P2 + 191)*(Li2); // Lni
			L2 += *(P2 + 212)*(Ll2); // Lnl
			L2 += *(P2 + 233)*(Lk2); // Lnk
			L2 += *(P2 + 254)*(Lm2); // Lnm
			L2 += *(P2 + 275)*(Lf2); // Lnf
			L2 += *(P2 + 296)*(Lp2); // Lnp
			L2 += *(P2 + 317)*(Lz2); // Lns
			L2 += *(P2 + 338)*(Lt2); // Lnt
			L2 += *(P2 + 359)*(Lw2); // Lnw
			L2 += *(P2 + 380)*(Ly2); // Lny
			L2 += *(P2 + 401)*(Lv2); // Lnv
			*(Ls3 + 2) = L1*L2;
			
			// L(D)
			L1 = *(P1 + 3)*(La1); // Lda
			L1 += *(P1 + 24)*(Lr1); // Ldr
			L1 += *(P1 + 45)*(Ln1); // Ldn
			L1 += *(P1 + 66)*(Ld1); // Ldd
			L1 += *(P1 + 87)*(Lc1); // Ldc
			L1 += *(P1 + 108)*(Lq1); // Ldq
			L1 += *(P1 + 129)*(Le1); // Lde
			L1 += *(P1 + 150)*(Lg1); // Ldg
			L1 += *(P1 + 171)*(Lh1); // Ldh
			L1 += *(P1 + 192)*(Li1); // Ldi
			L1 += *(P1 + 213)*(Ll1); // Ldl
			L1 += *(P1 + 234)*(Lk1); // Ldk
			L1 += *(P1 + 255)*(Lm1); // Ldm
			L1 += *(P1 + 276)*(Lf1); // Ldf
			L1 += *(P1 + 297)*(Lp1); // Ldp
			L1 += *(P1 + 318)*(Lz1); // Lds
			L1 += *(P1 + 339)*(Lt1); // Ldt
			L1 += *(P1 + 360)*(Lw1); // Ldw
			L1 += *(P1 + 381)*(Ly1); // Ldy
			L1 += *(P1 + 402)*(Lv1); // Ldv
			L2 = *(P2 + 3)*(La2); // Lda
			L2 += *(P2 + 24)*(Lr2); // Ldr
			L2 += *(P2 + 45)*(Ln2); // Ldn
			L2 += *(P2 + 66)*(Ld2); // Ldd
			L2 += *(P2 + 87)*(Lc2); // Ldc
			L2 += *(P2 + 108)*(Lq2); // Ldq
			L2 += *(P2 + 129)*(Le2); // Lde
			L2 += *(P2 + 150)*(Lg2); // Ldg
			L2 += *(P2 + 171)*(Lh2); // Ldh
			L2 += *(P2 + 192)*(Li2); // Ldi
			L2 += *(P2 + 213)*(Ll2); // Ldl
			L2 += *(P2 + 234)*(Lk2); // Ldk
			L2 += *(P2 + 255)*(Lm2); // Ldm
			L2 += *(P2 + 276)*(Lf2); // Ldf
			L2 += *(P2 + 297)*(Lp2); // Ldp
			L2 += *(P2 + 318)*(Lz2); // Lds
			L2 += *(P2 + 339)*(Lt2); // Ldt
			L2 += *(P2 + 360)*(Lw2); // Ldw
			L2 += *(P2 + 381)*(Ly2); // Ldy
			L2 += *(P2 + 402)*(Lv2); // Ldv
			*(Ls3 + 3) = L1*L2;
			
			// L(C)
			L1 = *(P1 + 4)*(La1); // Lca
			L1 += *(P1 + 25)*(Lr1); // Lcr
			L1 += *(P1 + 46)*(Ln1); // Lcn
			L1 += *(P1 + 67)*(Ld1); // Lcd
			L1 += *(P1 + 88)*(Lc1); // Lcc
			L1 += *(P1 + 109)*(Lq1); // Lcq
			L1 += *(P1 + 130)*(Le1); // Lce
			L1 += *(P1 + 151)*(Lg1); // Lcg
			L1 += *(P1 + 172)*(Lh1); // Lch
			L1 += *(P1 + 193)*(Li1); // Lci
			L1 += *(P1 + 214)*(Ll1); // Lcl
			L1 += *(P1 + 235)*(Lk1); // Lck
			L1 += *(P1 + 256)*(Lm1); // Lcm
			L1 += *(P1 + 277)*(Lf1); // Lcf
			L1 += *(P1 + 298)*(Lp1); // Lcp
			L1 += *(P1 + 319)*(Lz1); // Lcs
			L1 += *(P1 + 340)*(Lt1); // Lct
			L1 += *(P1 + 361)*(Lw1); // Lcw
			L1 += *(P1 + 382)*(Ly1); // Lcy
			L1 += *(P1 + 403)*(Lv1); // Lcv
			L2 = *(P2 + 4)*(La2); // Lca
			L2 += *(P2 + 25)*(Lr2); // Lcr
			L2 += *(P2 + 46)*(Ln2); // Lcn
			L2 += *(P2 + 67)*(Ld2); // Lcd
			L2 += *(P2 + 88)*(Lc2); // Lcc
			L2 += *(P2 + 109)*(Lq2); // Lcq
			L2 += *(P2 + 130)*(Le2); // Lce
			L2 += *(P2 + 151)*(Lg2); // Lcg
			L2 += *(P2 + 172)*(Lh2); // Lch
			L2 += *(P2 + 193)*(Li2); // Lci
			L2 += *(P2 + 214)*(Ll2); // Lcl
			L2 += *(P2 + 235)*(Lk2); // Lck
			L2 += *(P2 + 256)*(Lm2); // Lcm
			L2 += *(P2 + 277)*(Lf2); // Lcf
			L2 += *(P2 + 298)*(Lp2); // Lcp
			L2 += *(P2 + 319)*(Lz2); // Lcs
			L2 += *(P2 + 340)*(Lt2); // Lct
			L2 += *(P2 + 361)*(Lw2); // Lcw
			L2 += *(P2 + 382)*(Ly2); // Lcy
			L2 += *(P2 + 403)*(Lv2); // Lcv
			*(Ls3 + 4) = L1*L2;
			
			// L(Q)
			L1 = *(P1 + 5)*(La1); // Lqa
			L1 += *(P1 + 26)*(Lr1); // Lqr
			L1 += *(P1 + 47)*(Ln1); // Lqn
			L1 += *(P1 + 68)*(Ld1); // Lqd
			L1 += *(P1 + 89)*(Lc1); // Lqc
			L1 += *(P1 + 110)*(Lq1); // Lqq
			L1 += *(P1 + 131)*(Le1); // Lqe
			L1 += *(P1 + 152)*(Lg1); // Lqg
			L1 += *(P1 + 173)*(Lh1); // Lqh
			L1 += *(P1 + 194)*(Li1); // Lqi
			L1 += *(P1 + 215)*(Ll1); // Lql
			L1 += *(P1 + 236)*(Lk1); // Lqk
			L1 += *(P1 + 257)*(Lm1); // Lqm
			L1 += *(P1 + 278)*(Lf1); // Lqf
			L1 += *(P1 + 299)*(Lp1); // Lqp
			L1 += *(P1 + 320)*(Lz1); // Lqs
			L1 += *(P1 + 341)*(Lt1); // Lqt
			L1 += *(P1 + 362)*(Lw1); // Lqw
			L1 += *(P1 + 383)*(Ly1); // Lqy
			L1 += *(P1 + 404)*(Lv1); // Lqv
			L2 = *(P2 + 5)*(La2); // Lqa
			L2 += *(P2 + 26)*(Lr2); // Lqr
			L2 += *(P2 + 47)*(Ln2); // Lqn
			L2 += *(P2 + 68)*(Ld2); // Lqd
			L2 += *(P2 + 89)*(Lc2); // Lqc
			L2 += *(P2 + 110)*(Lq2); // Lqq
			L2 += *(P2 + 131)*(Le2); // Lqe
			L2 += *(P2 + 152)*(Lg2); // Lqg
			L2 += *(P2 + 173)*(Lh2); // Lqh
			L2 += *(P2 + 194)*(Li2); // Lqi
			L2 += *(P2 + 215)*(Ll2); // Lql
			L2 += *(P2 + 236)*(Lk2); // Lqk
			L2 += *(P2 + 257)*(Lm2); // Lqm
			L2 += *(P2 + 278)*(Lf2); // Lqf
			L2 += *(P2 + 299)*(Lp2); // Lqp
			L2 += *(P2 + 320)*(Lz2); // Lqs
			L2 += *(P2 + 341)*(Lt2); // Lqt
			L2 += *(P2 + 362)*(Lw2); // Lqw
			L2 += *(P2 + 383)*(Ly2); // Lqy
			L2 += *(P2 + 404)*(Lv2); // Lqv
			*(Ls3 + 5) = L1*L2;
			
			// L(E)
			L1 = *(P1 + 6)*(La1); // Lea
			L1 += *(P1 + 27)*(Lr1); // Ler
			L1 += *(P1 + 48)*(Ln1); // Len
			L1 += *(P1 + 69)*(Ld1); // Led
			L1 += *(P1 + 90)*(Lc1); // Lec
			L1 += *(P1 + 111)*(Lq1); // Leq
			L1 += *(P1 + 132)*(Le1); // Lee
			L1 += *(P1 + 153)*(Lg1); // Leg
			L1 += *(P1 + 174)*(Lh1); // Leh
			L1 += *(P1 + 195)*(Li1); // Lei
			L1 += *(P1 + 216)*(Ll1); // Lel
			L1 += *(P1 + 237)*(Lk1); // Lek
			L1 += *(P1 + 258)*(Lm1); // Lem
			L1 += *(P1 + 279)*(Lf1); // Lef
			L1 += *(P1 + 300)*(Lp1); // Lep
			L1 += *(P1 + 321)*(Lz1); // Les
			L1 += *(P1 + 342)*(Lt1); // Let
			L1 += *(P1 + 363)*(Lw1); // Lew
			L1 += *(P1 + 384)*(Ly1); // Ley
			L1 += *(P1 + 405)*(Lv1); // Lev
			L2 = *(P2 + 6)*(La2); // Lea
			L2 += *(P2 + 27)*(Lr2); // Ler
			L2 += *(P2 + 48)*(Ln2); // Len
			L2 += *(P2 + 69)*(Ld2); // Led
			L2 += *(P2 + 90)*(Lc2); // Lec
			L2 += *(P2 + 111)*(Lq2); // Leq
			L2 += *(P2 + 132)*(Le2); // Lee
			L2 += *(P2 + 153)*(Lg2); // Leg
			L2 += *(P2 + 174)*(Lh2); // Leh
			L2 += *(P2 + 195)*(Li2); // Lei
			L2 += *(P2 + 216)*(Ll2); // Lel
			L2 += *(P2 + 237)*(Lk2); // Lek
			L2 += *(P2 + 258)*(Lm2); // Lem
			L2 += *(P2 + 279)*(Lf2); // Lef
			L2 += *(P2 + 300)*(Lp2); // Lep
			L2 += *(P2 + 321)*(Lz2); // Les
			L2 += *(P2 + 342)*(Lt2); // Let
			L2 += *(P2 + 363)*(Lw2); // Lew
			L2 += *(P2 + 384)*(Ly2); // Ley
			L2 += *(P2 + 405)*(Lv2); // Lev
			*(Ls3 + 6) = L1*L2;
			
			// L(G)
			L1 = *(P1 + 7)*(La1); // Lga
			L1 += *(P1 + 28)*(Lr1); // Lgr
			L1 += *(P1 + 49)*(Ln1); // Lgn
			L1 += *(P1 + 70)*(Ld1); // Lgd
			L1 += *(P1 + 91)*(Lc1); // Lgc
			L1 += *(P1 + 112)*(Lq1); // Lgq
			L1 += *(P1 + 133)*(Le1); // Lge
			L1 += *(P1 + 154)*(Lg1); // Lgg
			L1 += *(P1 + 175)*(Lh1); // Lgh
			L1 += *(P1 + 196)*(Li1); // Lgi
			L1 += *(P1 + 217)*(Ll1); // Lgl
			L1 += *(P1 + 238)*(Lk1); // Lgk
			L1 += *(P1 + 259)*(Lm1); // Lgm
			L1 += *(P1 + 280)*(Lf1); // Lgf
			L1 += *(P1 + 301)*(Lp1); // Lgp
			L1 += *(P1 + 322)*(Lz1); // Lgs
			L1 += *(P1 + 343)*(Lt1); // Lgt
			L1 += *(P1 + 364)*(Lw1); // Lgw
			L1 += *(P1 + 385)*(Ly1); // Lgy
			L1 += *(P1 + 406)*(Lv1); // Lgv
			L2 = *(P2 + 7)*(La2); // Lga
			L2 += *(P2 + 28)*(Lr2); // Lgr
			L2 += *(P2 + 49)*(Ln2); // Lgn
			L2 += *(P2 + 70)*(Ld2); // Lgd
			L2 += *(P2 + 91)*(Lc2); // Lgc
			L2 += *(P2 + 112)*(Lq2); // Lgq
			L2 += *(P2 + 133)*(Le2); // Lge
			L2 += *(P2 + 154)*(Lg2); // Lgg
			L2 += *(P2 + 175)*(Lh2); // Lgh
			L2 += *(P2 + 196)*(Li2); // Lgi
			L2 += *(P2 + 217)*(Ll2); // Lgl
			L2 += *(P2 + 238)*(Lk2); // Lgk
			L2 += *(P2 + 259)*(Lm2); // Lgm
			L2 += *(P2 + 280)*(Lf2); // Lgf
			L2 += *(P2 + 301)*(Lp2); // Lgp
			L2 += *(P2 + 322)*(Lz2); // Lgs
			L2 += *(P2 + 343)*(Lt2); // Lgt
			L2 += *(P2 + 364)*(Lw2); // Lgw
			L2 += *(P2 + 385)*(Ly2); // Lgy
			L2 += *(P2 + 406)*(Lv2); // Lgv
			*(Ls3 + 7) = L1*L2;
			
			// L(H)
			L1 = *(P1 + 8)*(La1); // Lha
			L1 += *(P1 + 29)*(Lr1); // Lhr
			L1 += *(P1 + 50)*(Ln1); // Lhn
			L1 += *(P1 + 71)*(Ld1); // Lhd
			L1 += *(P1 + 92)*(Lc1); // Lhc
			L1 += *(P1 + 113)*(Lq1); // Lhq
			L1 += *(P1 + 134)*(Le1); // Lhe
			L1 += *(P1 + 155)*(Lg1); // Lhg
			L1 += *(P1 + 176)*(Lh1); // Lhh
			L1 += *(P1 + 197)*(Li1); // Lhi
			L1 += *(P1 + 218)*(Ll1); // Lhl
			L1 += *(P1 + 239)*(Lk1); // Lhk
			L1 += *(P1 + 260)*(Lm1); // Lhm
			L1 += *(P1 + 281)*(Lf1); // Lhf
			L1 += *(P1 + 302)*(Lp1); // Lhp
			L1 += *(P1 + 323)*(Lz1); // Lhs
			L1 += *(P1 + 344)*(Lt1); // Lht
			L1 += *(P1 + 365)*(Lw1); // Lhw
			L1 += *(P1 + 386)*(Ly1); // Lhy
			L1 += *(P1 + 407)*(Lv1); // Lhv
			L2 = *(P2 + 8)*(La2); // Lha
			L2 += *(P2 + 29)*(Lr2); // Lhr
			L2 += *(P2 + 50)*(Ln2); // Lhn
			L2 += *(P2 + 71)*(Ld2); // Lhd
			L2 += *(P2 + 92)*(Lc2); // Lhc
			L2 += *(P2 + 113)*(Lq2); // Lhq
			L2 += *(P2 + 134)*(Le2); // Lhe
			L2 += *(P2 + 155)*(Lg2); // Lhg
			L2 += *(P2 + 176)*(Lh2); // Lhh
			L2 += *(P2 + 197)*(Li2); // Lhi
			L2 += *(P2 + 218)*(Ll2); // Lhl
			L2 += *(P2 + 239)*(Lk2); // Lhk
			L2 += *(P2 + 260)*(Lm2); // Lhm
			L2 += *(P2 + 281)*(Lf2); // Lhf
			L2 += *(P2 + 302)*(Lp2); // Lhp
			L2 += *(P2 + 323)*(Lz2); // Lhs
			L2 += *(P2 + 344)*(Lt2); // Lht
			L2 += *(P2 + 365)*(Lw2); // Lhw
			L2 += *(P2 + 386)*(Ly2); // Lhy
			L2 += *(P2 + 407)*(Lv2); // Lhv
			*(Ls3 + 8) = L1*L2;
			
			// L(I)
			L1 = *(P1 + 9)*(La1); // Lia
			L1 += *(P1 + 30)*(Lr1); // Lir
			L1 += *(P1 + 51)*(Ln1); // Lin
			L1 += *(P1 + 72)*(Ld1); // Lid
			L1 += *(P1 + 93)*(Lc1); // Lic
			L1 += *(P1 + 114)*(Lq1); // Liq
			L1 += *(P1 + 135)*(Le1); // Lie
			L1 += *(P1 + 156)*(Lg1); // Lig
			L1 += *(P1 + 177)*(Lh1); // Lih
			L1 += *(P1 + 198)*(Li1); // Lii
			L1 += *(P1 + 219)*(Ll1); // Lil
			L1 += *(P1 + 240)*(Lk1); // Lik
			L1 += *(P1 + 261)*(Lm1); // Lim
			L1 += *(P1 + 282)*(Lf1); // Lif
			L1 += *(P1 + 303)*(Lp1); // Lip
			L1 += *(P1 + 324)*(Lz1); // Lis
			L1 += *(P1 + 345)*(Lt1); // Lit
			L1 += *(P1 + 366)*(Lw1); // Liw
			L1 += *(P1 + 387)*(Ly1); // Liy
			L1 += *(P1 + 408)*(Lv1); // Liv
			L2 = *(P2 + 9)*(La2); // Lia
			L2 += *(P2 + 30)*(Lr2); // Lir
			L2 += *(P2 + 51)*(Ln2); // Lin
			L2 += *(P2 + 72)*(Ld2); // Lid
			L2 += *(P2 + 93)*(Lc2); // Lic
			L2 += *(P2 + 114)*(Lq2); // Liq
			L2 += *(P2 + 135)*(Le2); // Lie
			L2 += *(P2 + 156)*(Lg2); // Lig
			L2 += *(P2 + 177)*(Lh2); // Lih
			L2 += *(P2 + 198)*(Li2); // Lii
			L2 += *(P2 + 219)*(Ll2); // Lil
			L2 += *(P2 + 240)*(Lk2); // Lik
			L2 += *(P2 + 261)*(Lm2); // Lim
			L2 += *(P2 + 282)*(Lf2); // Lif
			L2 += *(P2 + 303)*(Lp2); // Lip
			L2 += *(P2 + 324)*(Lz2); // Lis
			L2 += *(P2 + 345)*(Lt2); // Lit
			L2 += *(P2 + 366)*(Lw2); // Liw
			L2 += *(P2 + 387)*(Ly2); // Liy
			L2 += *(P2 + 408)*(Lv2); // Liv
			*(Ls3 + 9) = L1*L2;
			
			// L(L)
			L1 = *(P1 + 10)*(La1); // Lla
			L1 += *(P1 + 31)*(Lr1); // Llr
			L1 += *(P1 + 52)*(Ln1); // Lln
			L1 += *(P1 + 73)*(Ld1); // Lld
			L1 += *(P1 + 94)*(Lc1); // Llc
			L1 += *(P1 + 115)*(Lq1); // Llq
			L1 += *(P1 + 136)*(Le1); // Lle
			L1 += *(P1 + 157)*(Lg1); // Llg
			L1 += *(P1 + 178)*(Lh1); // Llh
			L1 += *(P1 + 199)*(Li1); // Lli
			L1 += *(P1 + 220)*(Ll1); // Lll
			L1 += *(P1 + 241)*(Lk1); // Llk
			L1 += *(P1 + 262)*(Lm1); // Llm
			L1 += *(P1 + 283)*(Lf1); // Llf
			L1 += *(P1 + 304)*(Lp1); // Llp
			L1 += *(P1 + 325)*(Lz1); // Lls
			L1 += *(P1 + 346)*(Lt1); // Llt
			L1 += *(P1 + 367)*(Lw1); // Llw
			L1 += *(P1 + 388)*(Ly1); // Lly
			L1 += *(P1 + 409)*(Lv1); // Llv
			L2 = *(P2 + 10)*(La2); // Lla
			L2 += *(P2 + 31)*(Lr2); // Llr
			L2 += *(P2 + 52)*(Ln2); // Lln
			L2 += *(P2 + 73)*(Ld2); // Lld
			L2 += *(P2 + 94)*(Lc2); // Llc
			L2 += *(P2 + 115)*(Lq2); // Llq
			L2 += *(P2 + 136)*(Le2); // Lle
			L2 += *(P2 + 157)*(Lg2); // Llg
			L2 += *(P2 + 178)*(Lh2); // Llh
			L2 += *(P2 + 199)*(Li2); // Lli
			L2 += *(P2 + 220)*(Ll2); // Lll
			L2 += *(P2 + 241)*(Lk2); // Llk
			L2 += *(P2 + 262)*(Lm2); // Llm
			L2 += *(P2 + 283)*(Lf2); // Llf
			L2 += *(P2 + 304)*(Lp2); // Llp
			L2 += *(P2 + 325)*(Lz2); // Lls
			L2 += *(P2 + 346)*(Lt2); // Llt
			L2 += *(P2 + 367)*(Lw2); // Llw
			L2 += *(P2 + 388)*(Ly2); // Lly
			L2 += *(P2 + 409)*(Lv2); // Llv
			*(Ls3 + 10) = L1*L2;
			
			// L(K)
			L1 = *(P1 + 11)*(La1); // Lka
			L1 += *(P1 + 32)*(Lr1); // Lkr
			L1 += *(P1 + 53)*(Ln1); // Lkn
			L1 += *(P1 + 74)*(Ld1); // Lkd
			L1 += *(P1 + 95)*(Lc1); // Lkc
			L1 += *(P1 + 116)*(Lq1); // Lkq
			L1 += *(P1 + 137)*(Le1); // Lke
			L1 += *(P1 + 158)*(Lg1); // Lkg
			L1 += *(P1 + 179)*(Lh1); // Lkh
			L1 += *(P1 + 200)*(Li1); // Lki
			L1 += *(P1 + 221)*(Ll1); // Lkl
			L1 += *(P1 + 242)*(Lk1); // Lkk
			L1 += *(P1 + 263)*(Lm1); // Lkm
			L1 += *(P1 + 284)*(Lf1); // Lkf
			L1 += *(P1 + 305)*(Lp1); // Lkp
			L1 += *(P1 + 326)*(Lz1); // Lks
			L1 += *(P1 + 347)*(Lt1); // Lkt
			L1 += *(P1 + 368)*(Lw1); // Lkw
			L1 += *(P1 + 389)*(Ly1); // Lky
			L1 += *(P1 + 410)*(Lv1); // Lkv
			L2 = *(P2 + 11)*(La2); // Lka
			L2 += *(P2 + 32)*(Lr2); // Lkr
			L2 += *(P2 + 53)*(Ln2); // Lkn
			L2 += *(P2 + 74)*(Ld2); // Lkd
			L2 += *(P2 + 95)*(Lc2); // Lkc
			L2 += *(P2 + 116)*(Lq2); // Lkq
			L2 += *(P2 + 137)*(Le2); // Lke
			L2 += *(P2 + 158)*(Lg2); // Lkg
			L2 += *(P2 + 179)*(Lh2); // Lkh
			L2 += *(P2 + 200)*(Li2); // Lki
			L2 += *(P2 + 221)*(Ll2); // Lkl
			L2 += *(P2 + 242)*(Lk2); // Lkk
			L2 += *(P2 + 263)*(Lm2); // Lkm
			L2 += *(P2 + 284)*(Lf2); // Lkf
			L2 += *(P2 + 305)*(Lp2); // Lkp
			L2 += *(P2 + 326)*(Lz2); // Lks
			L2 += *(P2 + 347)*(Lt2); // Lkt
			L2 += *(P2 + 368)*(Lw2); // Lkw
			L2 += *(P2 + 389)*(Ly2); // Lky
			L2 += *(P2 + 410)*(Lv2); // Lkv
			*(Ls3 + 11) = L1*L2;
			
			// L(M)
			L1 = *(P1 + 12)*(La1); // Lma
			L1 += *(P1 + 33)*(Lr1); // Lmr
			L1 += *(P1 + 54)*(Ln1); // Lmn
			L1 += *(P1 + 75)*(Ld1); // Lmd
			L1 += *(P1 + 96)*(Lc1); // Lmc
			L1 += *(P1 + 117)*(Lq1); // Lmq
			L1 += *(P1 + 138)*(Le1); // Lme
			L1 += *(P1 + 159)*(Lg1); // Lmg
			L1 += *(P1 + 180)*(Lh1); // Lmh
			L1 += *(P1 + 201)*(Li1); // Lmi
			L1 += *(P1 + 222)*(Ll1); // Lml
			L1 += *(P1 + 243)*(Lk1); // Lmk
			L1 += *(P1 + 264)*(Lm1); // Lmm
			L1 += *(P1 + 285)*(Lf1); // Lmf
			L1 += *(P1 + 306)*(Lp1); // Lmp
			L1 += *(P1 + 327)*(Lz1); // Lms
			L1 += *(P1 + 348)*(Lt1); // Lmt
			L1 += *(P1 + 369)*(Lw1); // Lmw
			L1 += *(P1 + 390)*(Ly1); // Lmy
			L1 += *(P1 + 411)*(Lv1); // Lmv
			L2 = *(P2 + 12)*(La2); // Lma
			L2 += *(P2 + 33)*(Lr2); // Lmr
			L2 += *(P2 + 54)*(Ln2); // Lmn
			L2 += *(P2 + 75)*(Ld2); // Lmd
			L2 += *(P2 + 96)*(Lc2); // Lmc
			L2 += *(P2 + 117)*(Lq2); // Lmq
			L2 += *(P2 + 138)*(Le2); // Lme
			L2 += *(P2 + 159)*(Lg2); // Lmg
			L2 += *(P2 + 180)*(Lh2); // Lmh
			L2 += *(P2 + 201)*(Li2); // Lmi
			L2 += *(P2 + 222)*(Ll2); // Lml
			L2 += *(P2 + 243)*(Lk2); // Lmk
			L2 += *(P2 + 264)*(Lm2); // Lmm
			L2 += *(P2 + 285)*(Lf2); // Lmf
			L2 += *(P2 + 306)*(Lp2); // Lmp
			L2 += *(P2 + 327)*(Lz2); // Lms
			L2 += *(P2 + 348)*(Lt2); // Lmt
			L2 += *(P2 + 369)*(Lw2); // Lmw
			L2 += *(P2 + 390)*(Ly2); // Lmy
			L2 += *(P2 + 411)*(Lv2); // Lmv
			*(Ls3 + 12) = L1*L2;
			
			// L(F)
			L1 = *(P1 + 13)*(La1); // Lfa
			L1 += *(P1 + 34)*(Lr1); // Lfr
			L1 += *(P1 + 55)*(Ln1); // Lfn
			L1 += *(P1 + 76)*(Ld1); // Lfd
			L1 += *(P1 + 97)*(Lc1); // Lfc
			L1 += *(P1 + 118)*(Lq1); // Lfq
			L1 += *(P1 + 139)*(Le1); // Lfe
			L1 += *(P1 + 160)*(Lg1); // Lfg
			L1 += *(P1 + 181)*(Lh1); // Lfh
			L1 += *(P1 + 202)*(Li1); // Lfi
			L1 += *(P1 + 223)*(Ll1); // Lfl
			L1 += *(P1 + 244)*(Lk1); // Lfk
			L1 += *(P1 + 265)*(Lm1); // Lfm
			L1 += *(P1 + 286)*(Lf1); // Lff
			L1 += *(P1 + 307)*(Lp1); // Lfp
			L1 += *(P1 + 328)*(Lz1); // Lfs
			L1 += *(P1 + 349)*(Lt1); // Lft
			L1 += *(P1 + 370)*(Lw1); // Lfw
			L1 += *(P1 + 391)*(Ly1); // Lfy
			L1 += *(P1 + 412)*(Lv1); // Lfv
			L2 = *(P2 + 13)*(La2); // Lfa
			L2 += *(P2 + 34)*(Lr2); // Lfr
			L2 += *(P2 + 55)*(Ln2); // Lfn
			L2 += *(P2 + 76)*(Ld2); // Lfd
			L2 += *(P2 + 97)*(Lc2); // Lfc
			L2 += *(P2 + 118)*(Lq2); // Lfq
			L2 += *(P2 + 139)*(Le2); // Lfe
			L2 += *(P2 + 160)*(Lg2); // Lfg
			L2 += *(P2 + 181)*(Lh2); // Lfh
			L2 += *(P2 + 202)*(Li2); // Lfi
			L2 += *(P2 + 223)*(Ll2); // Lfl
			L2 += *(P2 + 244)*(Lk2); // Lfk
			L2 += *(P2 + 265)*(Lm2); // Lfm
			L2 += *(P2 + 286)*(Lf2); // Lff
			L2 += *(P2 + 307)*(Lp2); // Lfp
			L2 += *(P2 + 328)*(Lz2); // Lfs
			L2 += *(P2 + 349)*(Lt2); // Lft
			L2 += *(P2 + 370)*(Lw2); // Lfw
			L2 += *(P2 + 391)*(Ly2); // Lfy
			L2 += *(P2 + 412)*(Lv2); // Lfv
			*(Ls3 + 13) = L1*L2;
			
			// L(P)
			L1 = *(P1 + 14)*(La1); // Lpa
			L1 += *(P1 + 35)*(Lr1); // Lpr
			L1 += *(P1 + 56)*(Ln1); // Lpn
			L1 += *(P1 + 77)*(Ld1); // Lpd
			L1 += *(P1 + 98)*(Lc1); // Lpc
			L1 += *(P1 + 119)*(Lq1); // Lpq
			L1 += *(P1 + 140)*(Le1); // Lpe
			L1 += *(P1 + 161)*(Lg1); // Lpg
			L1 += *(P1 + 182)*(Lh1); // Lph
			L1 += *(P1 + 203)*(Li1); // Lpi
			L1 += *(P1 + 224)*(Ll1); // Lpl
			L1 += *(P1 + 245)*(Lk1); // Lpk
			L1 += *(P1 + 266)*(Lm1); // Lpm
			L1 += *(P1 + 287)*(Lf1); // Lpf
			L1 += *(P1 + 308)*(Lp1); // Lpp
			L1 += *(P1 + 329)*(Lz1); // Lps
			L1 += *(P1 + 350)*(Lt1); // Lpt
			L1 += *(P1 + 371)*(Lw1); // Lpw
			L1 += *(P1 + 392)*(Ly1); // Lpy
			L1 += *(P1 + 413)*(Lv1); // Lpv
			L2 = *(P2 + 14)*(La2); // Lpa
			L2 += *(P2 + 35)*(Lr2); // Lpr
			L2 += *(P2 + 56)*(Ln2); // Lpn
			L2 += *(P2 + 77)*(Ld2); // Lpd
			L2 += *(P2 + 98)*(Lc2); // Lpc
			L2 += *(P2 + 119)*(Lq2); // Lpq
			L2 += *(P2 + 140)*(Le2); // Lpe
			L2 += *(P2 + 161)*(Lg2); // Lpg
			L2 += *(P2 + 182)*(Lh2); // Lph
			L2 += *(P2 + 203)*(Li2); // Lpi
			L2 += *(P2 + 224)*(Ll2); // Lpl
			L2 += *(P2 + 245)*(Lk2); // Lpk
			L2 += *(P2 + 266)*(Lm2); // Lpm
			L2 += *(P2 + 287)*(Lf2); // Lpf
			L2 += *(P2 + 308)*(Lp2); // Lpp
			L2 += *(P2 + 329)*(Lz2); // Lps
			L2 += *(P2 + 350)*(Lt2); // Lpt
			L2 += *(P2 + 371)*(Lw2); // Lpw
			L2 += *(P2 + 392)*(Ly2); // Lpy
			L2 += *(P2 + 413)*(Lv2); // Lpv
			*(Ls3 + 14) = L1*L2;
			
			// L(S)
			L1 = *(P1 + 15)*(La1); // Lsa
			L1 += *(P1 + 36)*(Lr1); // Lsr
			L1 += *(P1 + 57)*(Ln1); // Lsn
			L1 += *(P1 + 78)*(Ld1); // Lsd
			L1 += *(P1 + 99)*(Lc1); // Lsc
			L1 += *(P1 + 120)*(Lq1); // Lsq
			L1 += *(P1 + 141)*(Le1); // Lse
			L1 += *(P1 + 162)*(Lg1); // Lsg
			L1 += *(P1 + 183)*(Lh1); // Lsh
			L1 += *(P1 + 204)*(Li1); // Lsi
			L1 += *(P1 + 225)*(Ll1); // Lsl
			L1 += *(P1 + 246)*(Lk1); // Lsk
			L1 += *(P1 + 267)*(Lm1); // Lsm
			L1 += *(P1 + 288)*(Lf1); // Lsf
			L1 += *(P1 + 309)*(Lp1); // Lsp
			L1 += *(P1 + 330)*(Lz1); // Lss
			L1 += *(P1 + 351)*(Lt1); // Lst
			L1 += *(P1 + 372)*(Lw1); // Lsw
			L1 += *(P1 + 393)*(Ly1); // Lsy
			L1 += *(P1 + 414)*(Lv1); // Lsv
			L2 = *(P2 + 15)*(La2); // Lsa
			L2 += *(P2 + 36)*(Lr2); // Lsr
			L2 += *(P2 + 57)*(Ln2); // Lsn
			L2 += *(P2 + 78)*(Ld2); // Lsd
			L2 += *(P2 + 99)*(Lc2); // Lsc
			L2 += *(P2 + 120)*(Lq2); // Lsq
			L2 += *(P2 + 141)*(Le2); // Lse
			L2 += *(P2 + 162)*(Lg2); // Lsg
			L2 += *(P2 + 183)*(Lh2); // Lsh
			L2 += *(P2 + 204)*(Li2); // Lsi
			L2 += *(P2 + 225)*(Ll2); // Lsl
			L2 += *(P2 + 246)*(Lk2); // Lsk
			L2 += *(P2 + 267)*(Lm2); // Lsm
			L2 += *(P2 + 288)*(Lf2); // Lsf
			L2 += *(P2 + 309)*(Lp2); // Lsp
			L2 += *(P2 + 330)*(Lz2); // Lss
			L2 += *(P2 + 351)*(Lt2); // Lst
			L2 += *(P2 + 372)*(Lw2); // Lsw
			L2 += *(P2 + 393)*(Ly2); // Lsy
			L2 += *(P2 + 414)*(Lv2); // Lsv
			*(Ls3 + 15) = L1*L2;
			
			// L(T)
			L1 = *(P1 + 16)*(La1); // Lta
			L1 += *(P1 + 37)*(Lr1); // Ltr
			L1 += *(P1 + 58)*(Ln1); // Ltn
			L1 += *(P1 + 79)*(Ld1); // Ltd
			L1 += *(P1 + 100)*(Lc1); // Ltc
			L1 += *(P1 + 121)*(Lq1); // Ltq
			L1 += *(P1 + 142)*(Le1); // Lte
			L1 += *(P1 + 163)*(Lg1); // Ltg
			L1 += *(P1 + 184)*(Lh1); // Lth
			L1 += *(P1 + 205)*(Li1); // Lti
			L1 += *(P1 + 226)*(Ll1); // Ltl
			L1 += *(P1 + 247)*(Lk1); // Ltk
			L1 += *(P1 + 268)*(Lm1); // Ltm
			L1 += *(P1 + 289)*(Lf1); // Ltf
			L1 += *(P1 + 310)*(Lp1); // Ltp
			L1 += *(P1 + 331)*(Lz1); // Lts
			L1 += *(P1 + 352)*(Lt1); // Ltt
			L1 += *(P1 + 373)*(Lw1); // Ltw
			L1 += *(P1 + 394)*(Ly1); // Lty
			L1 += *(P1 + 415)*(Lv1); // Ltv
			L2 = *(P2 + 16)*(La2); // Lta
			L2 += *(P2 + 37)*(Lr2); // Ltr
			L2 += *(P2 + 58)*(Ln2); // Ltn
			L2 += *(P2 + 79)*(Ld2); // Ltd
			L2 += *(P2 + 100)*(Lc2); // Ltc
			L2 += *(P2 + 121)*(Lq2); // Ltq
			L2 += *(P2 + 142)*(Le2); // Lte
			L2 += *(P2 + 163)*(Lg2); // Ltg
			L2 += *(P2 + 184)*(Lh2); // Lth
			L2 += *(P2 + 205)*(Li2); // Lti
			L2 += *(P2 + 226)*(Ll2); // Ltl
			L2 += *(P2 + 247)*(Lk2); // Ltk
			L2 += *(P2 + 268)*(Lm2); // Ltm
			L2 += *(P2 + 289)*(Lf2); // Ltf
			L2 += *(P2 + 310)*(Lp2); // Ltp
			L2 += *(P2 + 331)*(Lz2); // Lts
			L2 += *(P2 + 352)*(Lt2); // Ltt
			L2 += *(P2 + 373)*(Lw2); // Ltw
			L2 += *(P2 + 394)*(Ly2); // Lty
			L2 += *(P2 + 415)*(Lv2); // Ltv
			*(Ls3 + 16) = L1*L2;
			
			// L(W)
			L1 = *(P1 + 17)*(La1); // Lwa
			L1 += *(P1 + 38)*(Lr1); // Lwr
			L1 += *(P1 + 59)*(Ln1); // Lwn
			L1 += *(P1 + 80)*(Ld1); // Lwd
			L1 += *(P1 + 101)*(Lc1); // Lwc
			L1 += *(P1 + 122)*(Lq1); // Lwq
			L1 += *(P1 + 143)*(Le1); // Lwe
			L1 += *(P1 + 164)*(Lg1); // Lwg
			L1 += *(P1 + 185)*(Lh1); // Lwh
			L1 += *(P1 + 206)*(Li1); // Lwi
			L1 += *(P1 + 227)*(Ll1); // Lwl
			L1 += *(P1 + 248)*(Lk1); // Lwk
			L1 += *(P1 + 269)*(Lm1); // Lwm
			L1 += *(P1 + 290)*(Lf1); // Lwf
			L1 += *(P1 + 311)*(Lp1); // Lwp
			L1 += *(P1 + 332)*(Lz1); // Lws
			L1 += *(P1 + 353)*(Lt1); // Lwt
			L1 += *(P1 + 374)*(Lw1); // Lww
			L1 += *(P1 + 395)*(Ly1); // Lwy
			L1 += *(P1 + 416)*(Lv1); // Lwv
			L2 = *(P2 + 17)*(La2); // Lwa
			L2 += *(P2 + 38)*(Lr2); // Lwr
			L2 += *(P2 + 59)*(Ln2); // Lwn
			L2 += *(P2 + 80)*(Ld2); // Lwd
			L2 += *(P2 + 101)*(Lc2); // Lwc
			L2 += *(P2 + 122)*(Lq2); // Lwq
			L2 += *(P2 + 143)*(Le2); // Lwe
			L2 += *(P2 + 164)*(Lg2); // Lwg
			L2 += *(P2 + 185)*(Lh2); // Lwh
			L2 += *(P2 + 206)*(Li2); // Lwi
			L2 += *(P2 + 227)*(Ll2); // Lwl
			L2 += *(P2 + 248)*(Lk2); // Lwk
			L2 += *(P2 + 269)*(Lm2); // Lwm
			L2 += *(P2 + 290)*(Lf2); // Lwf
			L2 += *(P2 + 311)*(Lp2); // Lwp
			L2 += *(P2 + 332)*(Lz2); // Lws
			L2 += *(P2 + 353)*(Lt2); // Lwt
			L2 += *(P2 + 374)*(Lw2); // Lww
			L2 += *(P2 + 395)*(Ly2); // Lwy
			L2 += *(P2 + 416)*(Lv2); // Lwv
			*(Ls3 + 17) = L1*L2;
			
			// L(Y)
			L1 = *(P1 + 18)*(La1); // Lya
			L1 += *(P1 + 39)*(Lr1); // Lyr
			L1 += *(P1 + 60)*(Ln1); // Lyn
			L1 += *(P1 + 81)*(Ld1); // Lyd
			L1 += *(P1 + 102)*(Lc1); // Lyc
			L1 += *(P1 + 123)*(Lq1); // Lyq
			L1 += *(P1 + 144)*(Le1); // Lye
			L1 += *(P1 + 165)*(Lg1); // Lyg
			L1 += *(P1 + 186)*(Lh1); // Lyh
			L1 += *(P1 + 207)*(Li1); // Lyi
			L1 += *(P1 + 228)*(Ll1); // Lyl
			L1 += *(P1 + 249)*(Lk1); // Lyk
			L1 += *(P1 + 270)*(Lm1); // Lym
			L1 += *(P1 + 291)*(Lf1); // Lyf
			L1 += *(P1 + 312)*(Lp1); // Lyp
			L1 += *(P1 + 333)*(Lz1); // Lys
			L1 += *(P1 + 354)*(Lt1); // Lyt
			L1 += *(P1 + 375)*(Lw1); // Lyw
			L1 += *(P1 + 396)*(Ly1); // Lyy
			L1 += *(P1 + 417)*(Lv1); // Lyv
			L2 = *(P2 + 18)*(La2); // Lya
			L2 += *(P2 + 39)*(Lr2); // Lyr
			L2 += *(P2 + 60)*(Ln2); // Lyn
			L2 += *(P2 + 81)*(Ld2); // Lyd
			L2 += *(P2 + 102)*(Lc2); // Lyc
			L2 += *(P2 + 123)*(Lq2); // Lyq
			L2 += *(P2 + 144)*(Le2); // Lye
			L2 += *(P2 + 165)*(Lg2); // Lyg
			L2 += *(P2 + 186)*(Lh2); // Lyh
			L2 += *(P2 + 207)*(Li2); // Lyi
			L2 += *(P2 + 228)*(Ll2); // Lyl
			L2 += *(P2 + 249)*(Lk2); // Lyk
			L2 += *(P2 + 270)*(Lm2); // Lym
			L2 += *(P2 + 291)*(Lf2); // Lyf
			L2 += *(P2 + 312)*(Lp2); // Lyp
			L2 += *(P2 + 333)*(Lz2); // Lys
			L2 += *(P2 + 354)*(Lt2); // Lyt
			L2 += *(P2 + 375)*(Lw2); // Lyw
			L2 += *(P2 + 396)*(Ly2); // Lyy
			L2 += *(P2 + 417)*(Lv2); // Lyv
			*(Ls3 + 18) = L1*L2;
			
			// L(V)
			L1 = *(P1 + 19)*(La1); // Lva
			L1 += *(P1 + 40)*(Lr1); // Lvr
			L1 += *(P1 + 61)*(Ln1); // Lvn
			L1 += *(P1 + 82)*(Ld1); // Lvd
			L1 += *(P1 + 103)*(Lc1); // Lvc
			L1 += *(P1 + 124)*(Lq1); // Lvq
			L1 += *(P1 + 145)*(Le1); // Lve
			L1 += *(P1 + 166)*(Lg1); // Lvg
			L1 += *(P1 + 187)*(Lh1); // Lvh
			L1 += *(P1 + 208)*(Li1); // Lvi
			L1 += *(P1 + 229)*(Ll1); // Lvl
			L1 += *(P1 + 250)*(Lk1); // Lvk
			L1 += *(P1 + 271)*(Lm1); // Lvm
			L1 += *(P1 + 292)*(Lf1); // Lvf
			L1 += *(P1 + 313)*(Lp1); // Lvp
			L1 += *(P1 + 334)*(Lz1); // Lvs
			L1 += *(P1 + 355)*(Lt1); // Lvt
			L1 += *(P1 + 376)*(Lw1); // Lvw
			L1 += *(P1 + 397)*(Ly1); // Lvy
			L1 += *(P1 + 418)*(Lv1); // Lvv
			L2 = *(P2 + 19)*(La2); // Lva
			L2 += *(P2 + 40)*(Lr2); // Lvr
			L2 += *(P2 + 61)*(Ln2); // Lvn
			L2 += *(P2 + 82)*(Ld2); // Lvd
			L2 += *(P2 + 103)*(Lc2); // Lvc
			L2 += *(P2 + 124)*(Lq2); // Lvq
			L2 += *(P2 + 145)*(Le2); // Lve
			L2 += *(P2 + 166)*(Lg2); // Lvg
			L2 += *(P2 + 187)*(Lh2); // Lvh
			L2 += *(P2 + 208)*(Li2); // Lvi
			L2 += *(P2 + 229)*(Ll2); // Lvl
			L2 += *(P2 + 250)*(Lk2); // Lvk
			L2 += *(P2 + 271)*(Lm2); // Lvm
			L2 += *(P2 + 292)*(Lf2); // Lvf
			L2 += *(P2 + 313)*(Lp2); // Lvp
			L2 += *(P2 + 334)*(Lz2); // Lvs
			L2 += *(P2 + 355)*(Lt2); // Lvt
			L2 += *(P2 + 376)*(Lw2); // Lvw
			L2 += *(P2 + 397)*(Ly2); // Lvy
			L2 += *(P2 + 418)*(Lv2); // Lvv
			*(Ls3 + 19) = L1*L2;
			
			*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon) ||
				(*(Ls3 + 5) > 0 && *(Ls3 + 5) < inv_epsilon) ||
				(*(Ls3 + 6) > 0 && *(Ls3 + 6) < inv_epsilon) ||
				(*(Ls3 + 7) > 0 && *(Ls3 + 7) < inv_epsilon) ||
				(*(Ls3 + 8) > 0 && *(Ls3 + 8) < inv_epsilon) ||
				(*(Ls3 + 9) > 0 && *(Ls3 + 9) < inv_epsilon) ||
				(*(Ls3 + 10) > 0 && *(Ls3 + 10) < inv_epsilon) ||
				(*(Ls3 + 11) > 0 && *(Ls3 + 11) < inv_epsilon) ||
				(*(Ls3 + 12) > 0 && *(Ls3 + 12) < inv_epsilon) ||
				(*(Ls3 + 13) > 0 && *(Ls3 + 13) < inv_epsilon) ||
				(*(Ls3 + 14) > 0 && *(Ls3 + 14) < inv_epsilon) ||
				(*(Ls3 + 15) > 0 && *(Ls3 + 15) < inv_epsilon) ||
				(*(Ls3 + 16) > 0 && *(Ls3 + 16) < inv_epsilon) ||
				(*(Ls3 + 17) > 0 && *(Ls3 + 17) < inv_epsilon) ||
				(*(Ls3 + 18) > 0 && *(Ls3 + 18) < inv_epsilon) ||
				(*(Ls3 + 19) > 0 && *(Ls3 + 19) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) *= epsilon;
				*(Ls3 + 6) *= epsilon;
				*(Ls3 + 7) *= epsilon;
				*(Ls3 + 8) *= epsilon;
				*(Ls3 + 9) *= epsilon;
				*(Ls3 + 10) *= epsilon;
				*(Ls3 + 11) *= epsilon;
				*(Ls3 + 12) *= epsilon;
				*(Ls3 + 13) *= epsilon;
				*(Ls3 + 14) *= epsilon;
				*(Ls3 + 15) *= epsilon;
				*(Ls3 + 16) *= epsilon;
				*(Ls3 + 17) *= epsilon;
				*(Ls3 + 18) *= epsilon;
				*(Ls3 + 19) *= epsilon;
				*(Ls3 + 21) += 1;
			}
		} else {
			// second branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 21)*(Lr1); // Lar
			L1 += *(P1 + 42)*(Ln1); // Lan
			L1 += *(P1 + 63)*(Ld1); // Lad
			L1 += *(P1 + 84)*(Lc1); // Lac
			L1 += *(P1 + 105)*(Lq1); // Laq
			L1 += *(P1 + 126)*(Le1); // Lae
			L1 += *(P1 + 147)*(Lg1); // Lag
			L1 += *(P1 + 168)*(Lh1); // Lah
			L1 += *(P1 + 189)*(Li1); // Lai
			L1 += *(P1 + 210)*(Ll1); // Lal
			L1 += *(P1 + 231)*(Lk1); // Lak
			L1 += *(P1 + 252)*(Lm1); // Lam
			L1 += *(P1 + 273)*(Lf1); // Laf
			L1 += *(P1 + 294)*(Lp1); // Lap
			L1 += *(P1 + 315)*(Lz1); // Las
			L1 += *(P1 + 336)*(Lt1); // Lat
			L1 += *(P1 + 357)*(Lw1); // Law
			L1 += *(P1 + 378)*(Ly1); // Lay
			L1 += *(P1 + 399)*(Lv1); // Lav
			*(Ls3 + 0) = L1;
			
			// L(R)
			L1 = *(P1 + 1)*(La1); // Lra
			L1 += *(P1 + 22)*(Lr1); // Lrr
			L1 += *(P1 + 43)*(Ln1); // Lrn
			L1 += *(P1 + 64)*(Ld1); // Lrd
			L1 += *(P1 + 85)*(Lc1); // Lrc
			L1 += *(P1 + 106)*(Lq1); // Lrq
			L1 += *(P1 + 127)*(Le1); // Lre
			L1 += *(P1 + 148)*(Lg1); // Lrg
			L1 += *(P1 + 169)*(Lh1); // Lrh
			L1 += *(P1 + 190)*(Li1); // Lri
			L1 += *(P1 + 211)*(Ll1); // Lrl
			L1 += *(P1 + 232)*(Lk1); // Lrk
			L1 += *(P1 + 253)*(Lm1); // Lrm
			L1 += *(P1 + 274)*(Lf1); // Lrf
			L1 += *(P1 + 295)*(Lp1); // Lrp
			L1 += *(P1 + 316)*(Lz1); // Lrs
			L1 += *(P1 + 337)*(Lt1); // Lrt
			L1 += *(P1 + 358)*(Lw1); // Lrw
			L1 += *(P1 + 379)*(Ly1); // Lry
			L1 += *(P1 + 400)*(Lv1); // Lrv
			*(Ls3 + 1) = L1;
			
			// L(N)
			L1 = *(P1 + 2)*(La1); // Lna
			L1 += *(P1 + 23)*(Lr1); // Lnr
			L1 += *(P1 + 44)*(Ln1); // Lnn
			L1 += *(P1 + 65)*(Ld1); // Lnd
			L1 += *(P1 + 86)*(Lc1); // Lnc
			L1 += *(P1 + 107)*(Lq1); // Lnq
			L1 += *(P1 + 128)*(Le1); // Lne
			L1 += *(P1 + 149)*(Lg1); // Lng
			L1 += *(P1 + 170)*(Lh1); // Lnh
			L1 += *(P1 + 191)*(Li1); // Lni
			L1 += *(P1 + 212)*(Ll1); // Lnl
			L1 += *(P1 + 233)*(Lk1); // Lnk
			L1 += *(P1 + 254)*(Lm1); // Lnm
			L1 += *(P1 + 275)*(Lf1); // Lnf
			L1 += *(P1 + 296)*(Lp1); // Lnp
			L1 += *(P1 + 317)*(Lz1); // Lns
			L1 += *(P1 + 338)*(Lt1); // Lnt
			L1 += *(P1 + 359)*(Lw1); // Lnw
			L1 += *(P1 + 380)*(Ly1); // Lny
			L1 += *(P1 + 401)*(Lv1); // Lnv
			*(Ls3 + 2) = L1;
			
			// L(D)
			L1 = *(P1 + 3)*(La1); // Lda
			L1 += *(P1 + 24)*(Lr1); // Ldr
			L1 += *(P1 + 45)*(Ln1); // Ldn
			L1 += *(P1 + 66)*(Ld1); // Ldd
			L1 += *(P1 + 87)*(Lc1); // Ldc
			L1 += *(P1 + 108)*(Lq1); // Ldq
			L1 += *(P1 + 129)*(Le1); // Lde
			L1 += *(P1 + 150)*(Lg1); // Ldg
			L1 += *(P1 + 171)*(Lh1); // Ldh
			L1 += *(P1 + 192)*(Li1); // Ldi
			L1 += *(P1 + 213)*(Ll1); // Ldl
			L1 += *(P1 + 234)*(Lk1); // Ldk
			L1 += *(P1 + 255)*(Lm1); // Ldm
			L1 += *(P1 + 276)*(Lf1); // Ldf
			L1 += *(P1 + 297)*(Lp1); // Ldp
			L1 += *(P1 + 318)*(Lz1); // Lds
			L1 += *(P1 + 339)*(Lt1); // Ldt
			L1 += *(P1 + 360)*(Lw1); // Ldw
			L1 += *(P1 + 381)*(Ly1); // Ldy
			L1 += *(P1 + 402)*(Lv1); // Ldv
			*(Ls3 + 3) = L1;
			
			// L(C)
			L1 = *(P1 + 4)*(La1); // Lca
			L1 += *(P1 + 25)*(Lr1); // Lcr
			L1 += *(P1 + 46)*(Ln1); // Lcn
			L1 += *(P1 + 67)*(Ld1); // Lcd
			L1 += *(P1 + 88)*(Lc1); // Lcc
			L1 += *(P1 + 109)*(Lq1); // Lcq
			L1 += *(P1 + 130)*(Le1); // Lce
			L1 += *(P1 + 151)*(Lg1); // Lcg
			L1 += *(P1 + 172)*(Lh1); // Lch
			L1 += *(P1 + 193)*(Li1); // Lci
			L1 += *(P1 + 214)*(Ll1); // Lcl
			L1 += *(P1 + 235)*(Lk1); // Lck
			L1 += *(P1 + 256)*(Lm1); // Lcm
			L1 += *(P1 + 277)*(Lf1); // Lcf
			L1 += *(P1 + 298)*(Lp1); // Lcp
			L1 += *(P1 + 319)*(Lz1); // Lcs
			L1 += *(P1 + 340)*(Lt1); // Lct
			L1 += *(P1 + 361)*(Lw1); // Lcw
			L1 += *(P1 + 382)*(Ly1); // Lcy
			L1 += *(P1 + 403)*(Lv1); // Lcv
			*(Ls3 + 4) = L1;
			
			// L(Q)
			L1 = *(P1 + 5)*(La1); // Lqa
			L1 += *(P1 + 26)*(Lr1); // Lqr
			L1 += *(P1 + 47)*(Ln1); // Lqn
			L1 += *(P1 + 68)*(Ld1); // Lqd
			L1 += *(P1 + 89)*(Lc1); // Lqc
			L1 += *(P1 + 110)*(Lq1); // Lqq
			L1 += *(P1 + 131)*(Le1); // Lqe
			L1 += *(P1 + 152)*(Lg1); // Lqg
			L1 += *(P1 + 173)*(Lh1); // Lqh
			L1 += *(P1 + 194)*(Li1); // Lqi
			L1 += *(P1 + 215)*(Ll1); // Lql
			L1 += *(P1 + 236)*(Lk1); // Lqk
			L1 += *(P1 + 257)*(Lm1); // Lqm
			L1 += *(P1 + 278)*(Lf1); // Lqf
			L1 += *(P1 + 299)*(Lp1); // Lqp
			L1 += *(P1 + 320)*(Lz1); // Lqs
			L1 += *(P1 + 341)*(Lt1); // Lqt
			L1 += *(P1 + 362)*(Lw1); // Lqw
			L1 += *(P1 + 383)*(Ly1); // Lqy
			L1 += *(P1 + 404)*(Lv1); // Lqv
			*(Ls3 + 5) = L1;
			
			// L(E)
			L1 = *(P1 + 6)*(La1); // Lea
			L1 += *(P1 + 27)*(Lr1); // Ler
			L1 += *(P1 + 48)*(Ln1); // Len
			L1 += *(P1 + 69)*(Ld1); // Led
			L1 += *(P1 + 90)*(Lc1); // Lec
			L1 += *(P1 + 111)*(Lq1); // Leq
			L1 += *(P1 + 132)*(Le1); // Lee
			L1 += *(P1 + 153)*(Lg1); // Leg
			L1 += *(P1 + 174)*(Lh1); // Leh
			L1 += *(P1 + 195)*(Li1); // Lei
			L1 += *(P1 + 216)*(Ll1); // Lel
			L1 += *(P1 + 237)*(Lk1); // Lek
			L1 += *(P1 + 258)*(Lm1); // Lem
			L1 += *(P1 + 279)*(Lf1); // Lef
			L1 += *(P1 + 300)*(Lp1); // Lep
			L1 += *(P1 + 321)*(Lz1); // Les
			L1 += *(P1 + 342)*(Lt1); // Let
			L1 += *(P1 + 363)*(Lw1); // Lew
			L1 += *(P1 + 384)*(Ly1); // Ley
			L1 += *(P1 + 405)*(Lv1); // Lev
			*(Ls3 + 6) = L1;
			
			// L(G)
			L1 = *(P1 + 7)*(La1); // Lga
			L1 += *(P1 + 28)*(Lr1); // Lgr
			L1 += *(P1 + 49)*(Ln1); // Lgn
			L1 += *(P1 + 70)*(Ld1); // Lgd
			L1 += *(P1 + 91)*(Lc1); // Lgc
			L1 += *(P1 + 112)*(Lq1); // Lgq
			L1 += *(P1 + 133)*(Le1); // Lge
			L1 += *(P1 + 154)*(Lg1); // Lgg
			L1 += *(P1 + 175)*(Lh1); // Lgh
			L1 += *(P1 + 196)*(Li1); // Lgi
			L1 += *(P1 + 217)*(Ll1); // Lgl
			L1 += *(P1 + 238)*(Lk1); // Lgk
			L1 += *(P1 + 259)*(Lm1); // Lgm
			L1 += *(P1 + 280)*(Lf1); // Lgf
			L1 += *(P1 + 301)*(Lp1); // Lgp
			L1 += *(P1 + 322)*(Lz1); // Lgs
			L1 += *(P1 + 343)*(Lt1); // Lgt
			L1 += *(P1 + 364)*(Lw1); // Lgw
			L1 += *(P1 + 385)*(Ly1); // Lgy
			L1 += *(P1 + 406)*(Lv1); // Lgv
			*(Ls3 + 7) = L1;
			
			// L(H)
			L1 = *(P1 + 8)*(La1); // Lha
			L1 += *(P1 + 29)*(Lr1); // Lhr
			L1 += *(P1 + 50)*(Ln1); // Lhn
			L1 += *(P1 + 71)*(Ld1); // Lhd
			L1 += *(P1 + 92)*(Lc1); // Lhc
			L1 += *(P1 + 113)*(Lq1); // Lhq
			L1 += *(P1 + 134)*(Le1); // Lhe
			L1 += *(P1 + 155)*(Lg1); // Lhg
			L1 += *(P1 + 176)*(Lh1); // Lhh
			L1 += *(P1 + 197)*(Li1); // Lhi
			L1 += *(P1 + 218)*(Ll1); // Lhl
			L1 += *(P1 + 239)*(Lk1); // Lhk
			L1 += *(P1 + 260)*(Lm1); // Lhm
			L1 += *(P1 + 281)*(Lf1); // Lhf
			L1 += *(P1 + 302)*(Lp1); // Lhp
			L1 += *(P1 + 323)*(Lz1); // Lhs
			L1 += *(P1 + 344)*(Lt1); // Lht
			L1 += *(P1 + 365)*(Lw1); // Lhw
			L1 += *(P1 + 386)*(Ly1); // Lhy
			L1 += *(P1 + 407)*(Lv1); // Lhv
			*(Ls3 + 8) = L1;
			
			// L(I)
			L1 = *(P1 + 9)*(La1); // Lia
			L1 += *(P1 + 30)*(Lr1); // Lir
			L1 += *(P1 + 51)*(Ln1); // Lin
			L1 += *(P1 + 72)*(Ld1); // Lid
			L1 += *(P1 + 93)*(Lc1); // Lic
			L1 += *(P1 + 114)*(Lq1); // Liq
			L1 += *(P1 + 135)*(Le1); // Lie
			L1 += *(P1 + 156)*(Lg1); // Lig
			L1 += *(P1 + 177)*(Lh1); // Lih
			L1 += *(P1 + 198)*(Li1); // Lii
			L1 += *(P1 + 219)*(Ll1); // Lil
			L1 += *(P1 + 240)*(Lk1); // Lik
			L1 += *(P1 + 261)*(Lm1); // Lim
			L1 += *(P1 + 282)*(Lf1); // Lif
			L1 += *(P1 + 303)*(Lp1); // Lip
			L1 += *(P1 + 324)*(Lz1); // Lis
			L1 += *(P1 + 345)*(Lt1); // Lit
			L1 += *(P1 + 366)*(Lw1); // Liw
			L1 += *(P1 + 387)*(Ly1); // Liy
			L1 += *(P1 + 408)*(Lv1); // Liv
			*(Ls3 + 9) = L1;
			
			// L(L)
			L1 = *(P1 + 10)*(La1); // Lla
			L1 += *(P1 + 31)*(Lr1); // Llr
			L1 += *(P1 + 52)*(Ln1); // Lln
			L1 += *(P1 + 73)*(Ld1); // Lld
			L1 += *(P1 + 94)*(Lc1); // Llc
			L1 += *(P1 + 115)*(Lq1); // Llq
			L1 += *(P1 + 136)*(Le1); // Lle
			L1 += *(P1 + 157)*(Lg1); // Llg
			L1 += *(P1 + 178)*(Lh1); // Llh
			L1 += *(P1 + 199)*(Li1); // Lli
			L1 += *(P1 + 220)*(Ll1); // Lll
			L1 += *(P1 + 241)*(Lk1); // Llk
			L1 += *(P1 + 262)*(Lm1); // Llm
			L1 += *(P1 + 283)*(Lf1); // Llf
			L1 += *(P1 + 304)*(Lp1); // Llp
			L1 += *(P1 + 325)*(Lz1); // Lls
			L1 += *(P1 + 346)*(Lt1); // Llt
			L1 += *(P1 + 367)*(Lw1); // Llw
			L1 += *(P1 + 388)*(Ly1); // Lly
			L1 += *(P1 + 409)*(Lv1); // Llv
			*(Ls3 + 10) = L1;
			
			// L(K)
			L1 = *(P1 + 11)*(La1); // Lka
			L1 += *(P1 + 32)*(Lr1); // Lkr
			L1 += *(P1 + 53)*(Ln1); // Lkn
			L1 += *(P1 + 74)*(Ld1); // Lkd
			L1 += *(P1 + 95)*(Lc1); // Lkc
			L1 += *(P1 + 116)*(Lq1); // Lkq
			L1 += *(P1 + 137)*(Le1); // Lke
			L1 += *(P1 + 158)*(Lg1); // Lkg
			L1 += *(P1 + 179)*(Lh1); // Lkh
			L1 += *(P1 + 200)*(Li1); // Lki
			L1 += *(P1 + 221)*(Ll1); // Lkl
			L1 += *(P1 + 242)*(Lk1); // Lkk
			L1 += *(P1 + 263)*(Lm1); // Lkm
			L1 += *(P1 + 284)*(Lf1); // Lkf
			L1 += *(P1 + 305)*(Lp1); // Lkp
			L1 += *(P1 + 326)*(Lz1); // Lks
			L1 += *(P1 + 347)*(Lt1); // Lkt
			L1 += *(P1 + 368)*(Lw1); // Lkw
			L1 += *(P1 + 389)*(Ly1); // Lky
			L1 += *(P1 + 410)*(Lv1); // Lkv
			*(Ls3 + 11) = L1;
			
			// L(M)
			L1 = *(P1 + 12)*(La1); // Lma
			L1 += *(P1 + 33)*(Lr1); // Lmr
			L1 += *(P1 + 54)*(Ln1); // Lmn
			L1 += *(P1 + 75)*(Ld1); // Lmd
			L1 += *(P1 + 96)*(Lc1); // Lmc
			L1 += *(P1 + 117)*(Lq1); // Lmq
			L1 += *(P1 + 138)*(Le1); // Lme
			L1 += *(P1 + 159)*(Lg1); // Lmg
			L1 += *(P1 + 180)*(Lh1); // Lmh
			L1 += *(P1 + 201)*(Li1); // Lmi
			L1 += *(P1 + 222)*(Ll1); // Lml
			L1 += *(P1 + 243)*(Lk1); // Lmk
			L1 += *(P1 + 264)*(Lm1); // Lmm
			L1 += *(P1 + 285)*(Lf1); // Lmf
			L1 += *(P1 + 306)*(Lp1); // Lmp
			L1 += *(P1 + 327)*(Lz1); // Lms
			L1 += *(P1 + 348)*(Lt1); // Lmt
			L1 += *(P1 + 369)*(Lw1); // Lmw
			L1 += *(P1 + 390)*(Ly1); // Lmy
			L1 += *(P1 + 411)*(Lv1); // Lmv
			*(Ls3 + 12) = L1;
			
			// L(F)
			L1 = *(P1 + 13)*(La1); // Lfa
			L1 += *(P1 + 34)*(Lr1); // Lfr
			L1 += *(P1 + 55)*(Ln1); // Lfn
			L1 += *(P1 + 76)*(Ld1); // Lfd
			L1 += *(P1 + 97)*(Lc1); // Lfc
			L1 += *(P1 + 118)*(Lq1); // Lfq
			L1 += *(P1 + 139)*(Le1); // Lfe
			L1 += *(P1 + 160)*(Lg1); // Lfg
			L1 += *(P1 + 181)*(Lh1); // Lfh
			L1 += *(P1 + 202)*(Li1); // Lfi
			L1 += *(P1 + 223)*(Ll1); // Lfl
			L1 += *(P1 + 244)*(Lk1); // Lfk
			L1 += *(P1 + 265)*(Lm1); // Lfm
			L1 += *(P1 + 286)*(Lf1); // Lff
			L1 += *(P1 + 307)*(Lp1); // Lfp
			L1 += *(P1 + 328)*(Lz1); // Lfs
			L1 += *(P1 + 349)*(Lt1); // Lft
			L1 += *(P1 + 370)*(Lw1); // Lfw
			L1 += *(P1 + 391)*(Ly1); // Lfy
			L1 += *(P1 + 412)*(Lv1); // Lfv
			*(Ls3 + 13) = L1;
			
			// L(P)
			L1 = *(P1 + 14)*(La1); // Lpa
			L1 += *(P1 + 35)*(Lr1); // Lpr
			L1 += *(P1 + 56)*(Ln1); // Lpn
			L1 += *(P1 + 77)*(Ld1); // Lpd
			L1 += *(P1 + 98)*(Lc1); // Lpc
			L1 += *(P1 + 119)*(Lq1); // Lpq
			L1 += *(P1 + 140)*(Le1); // Lpe
			L1 += *(P1 + 161)*(Lg1); // Lpg
			L1 += *(P1 + 182)*(Lh1); // Lph
			L1 += *(P1 + 203)*(Li1); // Lpi
			L1 += *(P1 + 224)*(Ll1); // Lpl
			L1 += *(P1 + 245)*(Lk1); // Lpk
			L1 += *(P1 + 266)*(Lm1); // Lpm
			L1 += *(P1 + 287)*(Lf1); // Lpf
			L1 += *(P1 + 308)*(Lp1); // Lpp
			L1 += *(P1 + 329)*(Lz1); // Lps
			L1 += *(P1 + 350)*(Lt1); // Lpt
			L1 += *(P1 + 371)*(Lw1); // Lpw
			L1 += *(P1 + 392)*(Ly1); // Lpy
			L1 += *(P1 + 413)*(Lv1); // Lpv
			*(Ls3 + 14) = L1;
			
			// L(S)
			L1 = *(P1 + 15)*(La1); // Lsa
			L1 += *(P1 + 36)*(Lr1); // Lsr
			L1 += *(P1 + 57)*(Ln1); // Lsn
			L1 += *(P1 + 78)*(Ld1); // Lsd
			L1 += *(P1 + 99)*(Lc1); // Lsc
			L1 += *(P1 + 120)*(Lq1); // Lsq
			L1 += *(P1 + 141)*(Le1); // Lse
			L1 += *(P1 + 162)*(Lg1); // Lsg
			L1 += *(P1 + 183)*(Lh1); // Lsh
			L1 += *(P1 + 204)*(Li1); // Lsi
			L1 += *(P1 + 225)*(Ll1); // Lsl
			L1 += *(P1 + 246)*(Lk1); // Lsk
			L1 += *(P1 + 267)*(Lm1); // Lsm
			L1 += *(P1 + 288)*(Lf1); // Lsf
			L1 += *(P1 + 309)*(Lp1); // Lsp
			L1 += *(P1 + 330)*(Lz1); // Lss
			L1 += *(P1 + 351)*(Lt1); // Lst
			L1 += *(P1 + 372)*(Lw1); // Lsw
			L1 += *(P1 + 393)*(Ly1); // Lsy
			L1 += *(P1 + 414)*(Lv1); // Lsv
			*(Ls3 + 15) = L1;
			
			// L(T)
			L1 = *(P1 + 16)*(La1); // Lta
			L1 += *(P1 + 37)*(Lr1); // Ltr
			L1 += *(P1 + 58)*(Ln1); // Ltn
			L1 += *(P1 + 79)*(Ld1); // Ltd
			L1 += *(P1 + 100)*(Lc1); // Ltc
			L1 += *(P1 + 121)*(Lq1); // Ltq
			L1 += *(P1 + 142)*(Le1); // Lte
			L1 += *(P1 + 163)*(Lg1); // Ltg
			L1 += *(P1 + 184)*(Lh1); // Lth
			L1 += *(P1 + 205)*(Li1); // Lti
			L1 += *(P1 + 226)*(Ll1); // Ltl
			L1 += *(P1 + 247)*(Lk1); // Ltk
			L1 += *(P1 + 268)*(Lm1); // Ltm
			L1 += *(P1 + 289)*(Lf1); // Ltf
			L1 += *(P1 + 310)*(Lp1); // Ltp
			L1 += *(P1 + 331)*(Lz1); // Lts
			L1 += *(P1 + 352)*(Lt1); // Ltt
			L1 += *(P1 + 373)*(Lw1); // Ltw
			L1 += *(P1 + 394)*(Ly1); // Lty
			L1 += *(P1 + 415)*(Lv1); // Ltv
			*(Ls3 + 16) = L1;
			
			// L(W)
			L1 = *(P1 + 17)*(La1); // Lwa
			L1 += *(P1 + 38)*(Lr1); // Lwr
			L1 += *(P1 + 59)*(Ln1); // Lwn
			L1 += *(P1 + 80)*(Ld1); // Lwd
			L1 += *(P1 + 101)*(Lc1); // Lwc
			L1 += *(P1 + 122)*(Lq1); // Lwq
			L1 += *(P1 + 143)*(Le1); // Lwe
			L1 += *(P1 + 164)*(Lg1); // Lwg
			L1 += *(P1 + 185)*(Lh1); // Lwh
			L1 += *(P1 + 206)*(Li1); // Lwi
			L1 += *(P1 + 227)*(Ll1); // Lwl
			L1 += *(P1 + 248)*(Lk1); // Lwk
			L1 += *(P1 + 269)*(Lm1); // Lwm
			L1 += *(P1 + 290)*(Lf1); // Lwf
			L1 += *(P1 + 311)*(Lp1); // Lwp
			L1 += *(P1 + 332)*(Lz1); // Lws
			L1 += *(P1 + 353)*(Lt1); // Lwt
			L1 += *(P1 + 374)*(Lw1); // Lww
			L1 += *(P1 + 395)*(Ly1); // Lwy
			L1 += *(P1 + 416)*(Lv1); // Lwv
			*(Ls3 + 17) = L1;
			
			// L(Y)
			L1 = *(P1 + 18)*(La1); // Lya
			L1 += *(P1 + 39)*(Lr1); // Lyr
			L1 += *(P1 + 60)*(Ln1); // Lyn
			L1 += *(P1 + 81)*(Ld1); // Lyd
			L1 += *(P1 + 102)*(Lc1); // Lyc
			L1 += *(P1 + 123)*(Lq1); // Lyq
			L1 += *(P1 + 144)*(Le1); // Lye
			L1 += *(P1 + 165)*(Lg1); // Lyg
			L1 += *(P1 + 186)*(Lh1); // Lyh
			L1 += *(P1 + 207)*(Li1); // Lyi
			L1 += *(P1 + 228)*(Ll1); // Lyl
			L1 += *(P1 + 249)*(Lk1); // Lyk
			L1 += *(P1 + 270)*(Lm1); // Lym
			L1 += *(P1 + 291)*(Lf1); // Lyf
			L1 += *(P1 + 312)*(Lp1); // Lyp
			L1 += *(P1 + 333)*(Lz1); // Lys
			L1 += *(P1 + 354)*(Lt1); // Lyt
			L1 += *(P1 + 375)*(Lw1); // Lyw
			L1 += *(P1 + 396)*(Ly1); // Lyy
			L1 += *(P1 + 417)*(Lv1); // Lyv
			*(Ls3 + 18) = L1;
			
			// L(V)
			L1 = *(P1 + 19)*(La1); // Lva
			L1 += *(P1 + 40)*(Lr1); // Lvr
			L1 += *(P1 + 61)*(Ln1); // Lvn
			L1 += *(P1 + 82)*(Ld1); // Lvd
			L1 += *(P1 + 103)*(Lc1); // Lvc
			L1 += *(P1 + 124)*(Lq1); // Lvq
			L1 += *(P1 + 145)*(Le1); // Lve
			L1 += *(P1 + 166)*(Lg1); // Lvg
			L1 += *(P1 + 187)*(Lh1); // Lvh
			L1 += *(P1 + 208)*(Li1); // Lvi
			L1 += *(P1 + 229)*(Ll1); // Lvl
			L1 += *(P1 + 250)*(Lk1); // Lvk
			L1 += *(P1 + 271)*(Lm1); // Lvm
			L1 += *(P1 + 292)*(Lf1); // Lvf
			L1 += *(P1 + 313)*(Lp1); // Lvp
			L1 += *(P1 + 334)*(Lz1); // Lvs
			L1 += *(P1 + 355)*(Lt1); // Lvt
			L1 += *(P1 + 376)*(Lw1); // Lvw
			L1 += *(P1 + 397)*(Ly1); // Lvy
			L1 += *(P1 + 418)*(Lv1); // Lvv
			*(Ls3 + 19) = L1;
			
			*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon) ||
				(*(Ls3 + 5) > 0 && *(Ls3 + 5) < inv_epsilon) ||
				(*(Ls3 + 6) > 0 && *(Ls3 + 6) < inv_epsilon) ||
				(*(Ls3 + 7) > 0 && *(Ls3 + 7) < inv_epsilon) ||
				(*(Ls3 + 8) > 0 && *(Ls3 + 8) < inv_epsilon) ||
				(*(Ls3 + 9) > 0 && *(Ls3 + 9) < inv_epsilon) ||
				(*(Ls3 + 10) > 0 && *(Ls3 + 10) < inv_epsilon) ||
				(*(Ls3 + 11) > 0 && *(Ls3 + 11) < inv_epsilon) ||
				(*(Ls3 + 12) > 0 && *(Ls3 + 12) < inv_epsilon) ||
				(*(Ls3 + 13) > 0 && *(Ls3 + 13) < inv_epsilon) ||
				(*(Ls3 + 14) > 0 && *(Ls3 + 14) < inv_epsilon) ||
				(*(Ls3 + 15) > 0 && *(Ls3 + 15) < inv_epsilon) ||
				(*(Ls3 + 16) > 0 && *(Ls3 + 16) < inv_epsilon) ||
				(*(Ls3 + 17) > 0 && *(Ls3 + 17) < inv_epsilon) ||
				(*(Ls3 + 18) > 0 && *(Ls3 + 18) < inv_epsilon) ||
				(*(Ls3 + 19) > 0 && *(Ls3 + 19) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) *= epsilon;
				*(Ls3 + 6) *= epsilon;
				*(Ls3 + 7) *= epsilon;
				*(Ls3 + 8) *= epsilon;
				*(Ls3 + 9) *= epsilon;
				*(Ls3 + 10) *= epsilon;
				*(Ls3 + 11) *= epsilon;
				*(Ls3 + 12) *= epsilon;
				*(Ls3 + 13) *= epsilon;
				*(Ls3 + 14) *= epsilon;
				*(Ls3 + 15) *= epsilon;
				*(Ls3 + 16) *= epsilon;
				*(Ls3 + 17) *= epsilon;
				*(Ls3 + 18) *= epsilon;
				*(Ls3 + 19) *= epsilon;
				*(Ls3 + 21) += 1;
			}
		}
	} else {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0) {
			// first branch can be disregarded
			
			// L(A)
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 21)*(Lr2); // Lar
			L2 += *(P2 + 42)*(Ln2); // Lan
			L2 += *(P2 + 63)*(Ld2); // Lad
			L2 += *(P2 + 84)*(Lc2); // Lac
			L2 += *(P2 + 105)*(Lq2); // Laq
			L2 += *(P2 + 126)*(Le2); // Lae
			L2 += *(P2 + 147)*(Lg2); // Lag
			L2 += *(P2 + 168)*(Lh2); // Lah
			L2 += *(P2 + 189)*(Li2); // Lai
			L2 += *(P2 + 210)*(Ll2); // Lal
			L2 += *(P2 + 231)*(Lk2); // Lak
			L2 += *(P2 + 252)*(Lm2); // Lam
			L2 += *(P2 + 273)*(Lf2); // Laf
			L2 += *(P2 + 294)*(Lp2); // Lap
			L2 += *(P2 + 315)*(Lz2); // Las
			L2 += *(P2 + 336)*(Lt2); // Lat
			L2 += *(P2 + 357)*(Lw2); // Law
			L2 += *(P2 + 378)*(Ly2); // Lay
			L2 += *(P2 + 399)*(Lv2); // Lav
			*(Ls3 + 0) = L2;
			
			// L(R)
			L2 = *(P2 + 1)*(La2); // Lra
			L2 += *(P2 + 22)*(Lr2); // Lrr
			L2 += *(P2 + 43)*(Ln2); // Lrn
			L2 += *(P2 + 64)*(Ld2); // Lrd
			L2 += *(P2 + 85)*(Lc2); // Lrc
			L2 += *(P2 + 106)*(Lq2); // Lrq
			L2 += *(P2 + 127)*(Le2); // Lre
			L2 += *(P2 + 148)*(Lg2); // Lrg
			L2 += *(P2 + 169)*(Lh2); // Lrh
			L2 += *(P2 + 190)*(Li2); // Lri
			L2 += *(P2 + 211)*(Ll2); // Lrl
			L2 += *(P2 + 232)*(Lk2); // Lrk
			L2 += *(P2 + 253)*(Lm2); // Lrm
			L2 += *(P2 + 274)*(Lf2); // Lrf
			L2 += *(P2 + 295)*(Lp2); // Lrp
			L2 += *(P2 + 316)*(Lz2); // Lrs
			L2 += *(P2 + 337)*(Lt2); // Lrt
			L2 += *(P2 + 358)*(Lw2); // Lrw
			L2 += *(P2 + 379)*(Ly2); // Lry
			L2 += *(P2 + 400)*(Lv2); // Lrv
			*(Ls3 + 1) = L2;
			
			// L(N)
			L2 = *(P2 + 2)*(La2); // Lna
			L2 += *(P2 + 23)*(Lr2); // Lnr
			L2 += *(P2 + 44)*(Ln2); // Lnn
			L2 += *(P2 + 65)*(Ld2); // Lnd
			L2 += *(P2 + 86)*(Lc2); // Lnc
			L2 += *(P2 + 107)*(Lq2); // Lnq
			L2 += *(P2 + 128)*(Le2); // Lne
			L2 += *(P2 + 149)*(Lg2); // Lng
			L2 += *(P2 + 170)*(Lh2); // Lnh
			L2 += *(P2 + 191)*(Li2); // Lni
			L2 += *(P2 + 212)*(Ll2); // Lnl
			L2 += *(P2 + 233)*(Lk2); // Lnk
			L2 += *(P2 + 254)*(Lm2); // Lnm
			L2 += *(P2 + 275)*(Lf2); // Lnf
			L2 += *(P2 + 296)*(Lp2); // Lnp
			L2 += *(P2 + 317)*(Lz2); // Lns
			L2 += *(P2 + 338)*(Lt2); // Lnt
			L2 += *(P2 + 359)*(Lw2); // Lnw
			L2 += *(P2 + 380)*(Ly2); // Lny
			L2 += *(P2 + 401)*(Lv2); // Lnv
			*(Ls3 + 2) = L2;
			
			// L(D)
			L2 = *(P2 + 3)*(La2); // Lda
			L2 += *(P2 + 24)*(Lr2); // Ldr
			L2 += *(P2 + 45)*(Ln2); // Ldn
			L2 += *(P2 + 66)*(Ld2); // Ldd
			L2 += *(P2 + 87)*(Lc2); // Ldc
			L2 += *(P2 + 108)*(Lq2); // Ldq
			L2 += *(P2 + 129)*(Le2); // Lde
			L2 += *(P2 + 150)*(Lg2); // Ldg
			L2 += *(P2 + 171)*(Lh2); // Ldh
			L2 += *(P2 + 192)*(Li2); // Ldi
			L2 += *(P2 + 213)*(Ll2); // Ldl
			L2 += *(P2 + 234)*(Lk2); // Ldk
			L2 += *(P2 + 255)*(Lm2); // Ldm
			L2 += *(P2 + 276)*(Lf2); // Ldf
			L2 += *(P2 + 297)*(Lp2); // Ldp
			L2 += *(P2 + 318)*(Lz2); // Lds
			L2 += *(P2 + 339)*(Lt2); // Ldt
			L2 += *(P2 + 360)*(Lw2); // Ldw
			L2 += *(P2 + 381)*(Ly2); // Ldy
			L2 += *(P2 + 402)*(Lv2); // Ldv
			*(Ls3 + 3) = L2;
			
			// L(C)
			L2 = *(P2 + 4)*(La2); // Lca
			L2 += *(P2 + 25)*(Lr2); // Lcr
			L2 += *(P2 + 46)*(Ln2); // Lcn
			L2 += *(P2 + 67)*(Ld2); // Lcd
			L2 += *(P2 + 88)*(Lc2); // Lcc
			L2 += *(P2 + 109)*(Lq2); // Lcq
			L2 += *(P2 + 130)*(Le2); // Lce
			L2 += *(P2 + 151)*(Lg2); // Lcg
			L2 += *(P2 + 172)*(Lh2); // Lch
			L2 += *(P2 + 193)*(Li2); // Lci
			L2 += *(P2 + 214)*(Ll2); // Lcl
			L2 += *(P2 + 235)*(Lk2); // Lck
			L2 += *(P2 + 256)*(Lm2); // Lcm
			L2 += *(P2 + 277)*(Lf2); // Lcf
			L2 += *(P2 + 298)*(Lp2); // Lcp
			L2 += *(P2 + 319)*(Lz2); // Lcs
			L2 += *(P2 + 340)*(Lt2); // Lct
			L2 += *(P2 + 361)*(Lw2); // Lcw
			L2 += *(P2 + 382)*(Ly2); // Lcy
			L2 += *(P2 + 403)*(Lv2); // Lcv
			*(Ls3 + 4) = L2;
			
			// L(Q)
			L2 = *(P2 + 5)*(La2); // Lqa
			L2 += *(P2 + 26)*(Lr2); // Lqr
			L2 += *(P2 + 47)*(Ln2); // Lqn
			L2 += *(P2 + 68)*(Ld2); // Lqd
			L2 += *(P2 + 89)*(Lc2); // Lqc
			L2 += *(P2 + 110)*(Lq2); // Lqq
			L2 += *(P2 + 131)*(Le2); // Lqe
			L2 += *(P2 + 152)*(Lg2); // Lqg
			L2 += *(P2 + 173)*(Lh2); // Lqh
			L2 += *(P2 + 194)*(Li2); // Lqi
			L2 += *(P2 + 215)*(Ll2); // Lql
			L2 += *(P2 + 236)*(Lk2); // Lqk
			L2 += *(P2 + 257)*(Lm2); // Lqm
			L2 += *(P2 + 278)*(Lf2); // Lqf
			L2 += *(P2 + 299)*(Lp2); // Lqp
			L2 += *(P2 + 320)*(Lz2); // Lqs
			L2 += *(P2 + 341)*(Lt2); // Lqt
			L2 += *(P2 + 362)*(Lw2); // Lqw
			L2 += *(P2 + 383)*(Ly2); // Lqy
			L2 += *(P2 + 404)*(Lv2); // Lqv
			*(Ls3 + 5) = L2;
			
			// L(E)
			L2 = *(P2 + 6)*(La2); // Lea
			L2 += *(P2 + 27)*(Lr2); // Ler
			L2 += *(P2 + 48)*(Ln2); // Len
			L2 += *(P2 + 69)*(Ld2); // Led
			L2 += *(P2 + 90)*(Lc2); // Lec
			L2 += *(P2 + 111)*(Lq2); // Leq
			L2 += *(P2 + 132)*(Le2); // Lee
			L2 += *(P2 + 153)*(Lg2); // Leg
			L2 += *(P2 + 174)*(Lh2); // Leh
			L2 += *(P2 + 195)*(Li2); // Lei
			L2 += *(P2 + 216)*(Ll2); // Lel
			L2 += *(P2 + 237)*(Lk2); // Lek
			L2 += *(P2 + 258)*(Lm2); // Lem
			L2 += *(P2 + 279)*(Lf2); // Lef
			L2 += *(P2 + 300)*(Lp2); // Lep
			L2 += *(P2 + 321)*(Lz2); // Les
			L2 += *(P2 + 342)*(Lt2); // Let
			L2 += *(P2 + 363)*(Lw2); // Lew
			L2 += *(P2 + 384)*(Ly2); // Ley
			L2 += *(P2 + 405)*(Lv2); // Lev
			*(Ls3 + 6) = L2;
			
			// L(G)
			L2 = *(P2 + 7)*(La2); // Lga
			L2 += *(P2 + 28)*(Lr2); // Lgr
			L2 += *(P2 + 49)*(Ln2); // Lgn
			L2 += *(P2 + 70)*(Ld2); // Lgd
			L2 += *(P2 + 91)*(Lc2); // Lgc
			L2 += *(P2 + 112)*(Lq2); // Lgq
			L2 += *(P2 + 133)*(Le2); // Lge
			L2 += *(P2 + 154)*(Lg2); // Lgg
			L2 += *(P2 + 175)*(Lh2); // Lgh
			L2 += *(P2 + 196)*(Li2); // Lgi
			L2 += *(P2 + 217)*(Ll2); // Lgl
			L2 += *(P2 + 238)*(Lk2); // Lgk
			L2 += *(P2 + 259)*(Lm2); // Lgm
			L2 += *(P2 + 280)*(Lf2); // Lgf
			L2 += *(P2 + 301)*(Lp2); // Lgp
			L2 += *(P2 + 322)*(Lz2); // Lgs
			L2 += *(P2 + 343)*(Lt2); // Lgt
			L2 += *(P2 + 364)*(Lw2); // Lgw
			L2 += *(P2 + 385)*(Ly2); // Lgy
			L2 += *(P2 + 406)*(Lv2); // Lgv
			*(Ls3 + 7) = L2;
			
			// L(H)
			L2 = *(P2 + 8)*(La2); // Lha
			L2 += *(P2 + 29)*(Lr2); // Lhr
			L2 += *(P2 + 50)*(Ln2); // Lhn
			L2 += *(P2 + 71)*(Ld2); // Lhd
			L2 += *(P2 + 92)*(Lc2); // Lhc
			L2 += *(P2 + 113)*(Lq2); // Lhq
			L2 += *(P2 + 134)*(Le2); // Lhe
			L2 += *(P2 + 155)*(Lg2); // Lhg
			L2 += *(P2 + 176)*(Lh2); // Lhh
			L2 += *(P2 + 197)*(Li2); // Lhi
			L2 += *(P2 + 218)*(Ll2); // Lhl
			L2 += *(P2 + 239)*(Lk2); // Lhk
			L2 += *(P2 + 260)*(Lm2); // Lhm
			L2 += *(P2 + 281)*(Lf2); // Lhf
			L2 += *(P2 + 302)*(Lp2); // Lhp
			L2 += *(P2 + 323)*(Lz2); // Lhs
			L2 += *(P2 + 344)*(Lt2); // Lht
			L2 += *(P2 + 365)*(Lw2); // Lhw
			L2 += *(P2 + 386)*(Ly2); // Lhy
			L2 += *(P2 + 407)*(Lv2); // Lhv
			*(Ls3 + 8) = L2;
			
			// L(I)
			L2 = *(P2 + 9)*(La2); // Lia
			L2 += *(P2 + 30)*(Lr2); // Lir
			L2 += *(P2 + 51)*(Ln2); // Lin
			L2 += *(P2 + 72)*(Ld2); // Lid
			L2 += *(P2 + 93)*(Lc2); // Lic
			L2 += *(P2 + 114)*(Lq2); // Liq
			L2 += *(P2 + 135)*(Le2); // Lie
			L2 += *(P2 + 156)*(Lg2); // Lig
			L2 += *(P2 + 177)*(Lh2); // Lih
			L2 += *(P2 + 198)*(Li2); // Lii
			L2 += *(P2 + 219)*(Ll2); // Lil
			L2 += *(P2 + 240)*(Lk2); // Lik
			L2 += *(P2 + 261)*(Lm2); // Lim
			L2 += *(P2 + 282)*(Lf2); // Lif
			L2 += *(P2 + 303)*(Lp2); // Lip
			L2 += *(P2 + 324)*(Lz2); // Lis
			L2 += *(P2 + 345)*(Lt2); // Lit
			L2 += *(P2 + 366)*(Lw2); // Liw
			L2 += *(P2 + 387)*(Ly2); // Liy
			L2 += *(P2 + 408)*(Lv2); // Liv
			*(Ls3 + 9) = L2;
			
			// L(L)
			L2 = *(P2 + 10)*(La2); // Lla
			L2 += *(P2 + 31)*(Lr2); // Llr
			L2 += *(P2 + 52)*(Ln2); // Lln
			L2 += *(P2 + 73)*(Ld2); // Lld
			L2 += *(P2 + 94)*(Lc2); // Llc
			L2 += *(P2 + 115)*(Lq2); // Llq
			L2 += *(P2 + 136)*(Le2); // Lle
			L2 += *(P2 + 157)*(Lg2); // Llg
			L2 += *(P2 + 178)*(Lh2); // Llh
			L2 += *(P2 + 199)*(Li2); // Lli
			L2 += *(P2 + 220)*(Ll2); // Lll
			L2 += *(P2 + 241)*(Lk2); // Llk
			L2 += *(P2 + 262)*(Lm2); // Llm
			L2 += *(P2 + 283)*(Lf2); // Llf
			L2 += *(P2 + 304)*(Lp2); // Llp
			L2 += *(P2 + 325)*(Lz2); // Lls
			L2 += *(P2 + 346)*(Lt2); // Llt
			L2 += *(P2 + 367)*(Lw2); // Llw
			L2 += *(P2 + 388)*(Ly2); // Lly
			L2 += *(P2 + 409)*(Lv2); // Llv
			*(Ls3 + 10) = L2;
			
			// L(K)
			L2 = *(P2 + 11)*(La2); // Lka
			L2 += *(P2 + 32)*(Lr2); // Lkr
			L2 += *(P2 + 53)*(Ln2); // Lkn
			L2 += *(P2 + 74)*(Ld2); // Lkd
			L2 += *(P2 + 95)*(Lc2); // Lkc
			L2 += *(P2 + 116)*(Lq2); // Lkq
			L2 += *(P2 + 137)*(Le2); // Lke
			L2 += *(P2 + 158)*(Lg2); // Lkg
			L2 += *(P2 + 179)*(Lh2); // Lkh
			L2 += *(P2 + 200)*(Li2); // Lki
			L2 += *(P2 + 221)*(Ll2); // Lkl
			L2 += *(P2 + 242)*(Lk2); // Lkk
			L2 += *(P2 + 263)*(Lm2); // Lkm
			L2 += *(P2 + 284)*(Lf2); // Lkf
			L2 += *(P2 + 305)*(Lp2); // Lkp
			L2 += *(P2 + 326)*(Lz2); // Lks
			L2 += *(P2 + 347)*(Lt2); // Lkt
			L2 += *(P2 + 368)*(Lw2); // Lkw
			L2 += *(P2 + 389)*(Ly2); // Lky
			L2 += *(P2 + 410)*(Lv2); // Lkv
			*(Ls3 + 11) = L2;
			
			// L(M)
			L2 = *(P2 + 12)*(La2); // Lma
			L2 += *(P2 + 33)*(Lr2); // Lmr
			L2 += *(P2 + 54)*(Ln2); // Lmn
			L2 += *(P2 + 75)*(Ld2); // Lmd
			L2 += *(P2 + 96)*(Lc2); // Lmc
			L2 += *(P2 + 117)*(Lq2); // Lmq
			L2 += *(P2 + 138)*(Le2); // Lme
			L2 += *(P2 + 159)*(Lg2); // Lmg
			L2 += *(P2 + 180)*(Lh2); // Lmh
			L2 += *(P2 + 201)*(Li2); // Lmi
			L2 += *(P2 + 222)*(Ll2); // Lml
			L2 += *(P2 + 243)*(Lk2); // Lmk
			L2 += *(P2 + 264)*(Lm2); // Lmm
			L2 += *(P2 + 285)*(Lf2); // Lmf
			L2 += *(P2 + 306)*(Lp2); // Lmp
			L2 += *(P2 + 327)*(Lz2); // Lms
			L2 += *(P2 + 348)*(Lt2); // Lmt
			L2 += *(P2 + 369)*(Lw2); // Lmw
			L2 += *(P2 + 390)*(Ly2); // Lmy
			L2 += *(P2 + 411)*(Lv2); // Lmv
			*(Ls3 + 12) = L2;
			
			// L(F)
			L2 = *(P2 + 13)*(La2); // Lfa
			L2 += *(P2 + 34)*(Lr2); // Lfr
			L2 += *(P2 + 55)*(Ln2); // Lfn
			L2 += *(P2 + 76)*(Ld2); // Lfd
			L2 += *(P2 + 97)*(Lc2); // Lfc
			L2 += *(P2 + 118)*(Lq2); // Lfq
			L2 += *(P2 + 139)*(Le2); // Lfe
			L2 += *(P2 + 160)*(Lg2); // Lfg
			L2 += *(P2 + 181)*(Lh2); // Lfh
			L2 += *(P2 + 202)*(Li2); // Lfi
			L2 += *(P2 + 223)*(Ll2); // Lfl
			L2 += *(P2 + 244)*(Lk2); // Lfk
			L2 += *(P2 + 265)*(Lm2); // Lfm
			L2 += *(P2 + 286)*(Lf2); // Lff
			L2 += *(P2 + 307)*(Lp2); // Lfp
			L2 += *(P2 + 328)*(Lz2); // Lfs
			L2 += *(P2 + 349)*(Lt2); // Lft
			L2 += *(P2 + 370)*(Lw2); // Lfw
			L2 += *(P2 + 391)*(Ly2); // Lfy
			L2 += *(P2 + 412)*(Lv2); // Lfv
			*(Ls3 + 13) = L2;
			
			// L(P)
			L2 = *(P2 + 14)*(La2); // Lpa
			L2 += *(P2 + 35)*(Lr2); // Lpr
			L2 += *(P2 + 56)*(Ln2); // Lpn
			L2 += *(P2 + 77)*(Ld2); // Lpd
			L2 += *(P2 + 98)*(Lc2); // Lpc
			L2 += *(P2 + 119)*(Lq2); // Lpq
			L2 += *(P2 + 140)*(Le2); // Lpe
			L2 += *(P2 + 161)*(Lg2); // Lpg
			L2 += *(P2 + 182)*(Lh2); // Lph
			L2 += *(P2 + 203)*(Li2); // Lpi
			L2 += *(P2 + 224)*(Ll2); // Lpl
			L2 += *(P2 + 245)*(Lk2); // Lpk
			L2 += *(P2 + 266)*(Lm2); // Lpm
			L2 += *(P2 + 287)*(Lf2); // Lpf
			L2 += *(P2 + 308)*(Lp2); // Lpp
			L2 += *(P2 + 329)*(Lz2); // Lps
			L2 += *(P2 + 350)*(Lt2); // Lpt
			L2 += *(P2 + 371)*(Lw2); // Lpw
			L2 += *(P2 + 392)*(Ly2); // Lpy
			L2 += *(P2 + 413)*(Lv2); // Lpv
			*(Ls3 + 14) = L2;
			
			// L(S)
			L2 = *(P2 + 15)*(La2); // Lsa
			L2 += *(P2 + 36)*(Lr2); // Lsr
			L2 += *(P2 + 57)*(Ln2); // Lsn
			L2 += *(P2 + 78)*(Ld2); // Lsd
			L2 += *(P2 + 99)*(Lc2); // Lsc
			L2 += *(P2 + 120)*(Lq2); // Lsq
			L2 += *(P2 + 141)*(Le2); // Lse
			L2 += *(P2 + 162)*(Lg2); // Lsg
			L2 += *(P2 + 183)*(Lh2); // Lsh
			L2 += *(P2 + 204)*(Li2); // Lsi
			L2 += *(P2 + 225)*(Ll2); // Lsl
			L2 += *(P2 + 246)*(Lk2); // Lsk
			L2 += *(P2 + 267)*(Lm2); // Lsm
			L2 += *(P2 + 288)*(Lf2); // Lsf
			L2 += *(P2 + 309)*(Lp2); // Lsp
			L2 += *(P2 + 330)*(Lz2); // Lss
			L2 += *(P2 + 351)*(Lt2); // Lst
			L2 += *(P2 + 372)*(Lw2); // Lsw
			L2 += *(P2 + 393)*(Ly2); // Lsy
			L2 += *(P2 + 414)*(Lv2); // Lsv
			*(Ls3 + 15) = L2;
			
			// L(T)
			L2 = *(P2 + 16)*(La2); // Lta
			L2 += *(P2 + 37)*(Lr2); // Ltr
			L2 += *(P2 + 58)*(Ln2); // Ltn
			L2 += *(P2 + 79)*(Ld2); // Ltd
			L2 += *(P2 + 100)*(Lc2); // Ltc
			L2 += *(P2 + 121)*(Lq2); // Ltq
			L2 += *(P2 + 142)*(Le2); // Lte
			L2 += *(P2 + 163)*(Lg2); // Ltg
			L2 += *(P2 + 184)*(Lh2); // Lth
			L2 += *(P2 + 205)*(Li2); // Lti
			L2 += *(P2 + 226)*(Ll2); // Ltl
			L2 += *(P2 + 247)*(Lk2); // Ltk
			L2 += *(P2 + 268)*(Lm2); // Ltm
			L2 += *(P2 + 289)*(Lf2); // Ltf
			L2 += *(P2 + 310)*(Lp2); // Ltp
			L2 += *(P2 + 331)*(Lz2); // Lts
			L2 += *(P2 + 352)*(Lt2); // Ltt
			L2 += *(P2 + 373)*(Lw2); // Ltw
			L2 += *(P2 + 394)*(Ly2); // Lty
			L2 += *(P2 + 415)*(Lv2); // Ltv
			*(Ls3 + 16) = L2;
			
			// L(W)
			L2 = *(P2 + 17)*(La2); // Lwa
			L2 += *(P2 + 38)*(Lr2); // Lwr
			L2 += *(P2 + 59)*(Ln2); // Lwn
			L2 += *(P2 + 80)*(Ld2); // Lwd
			L2 += *(P2 + 101)*(Lc2); // Lwc
			L2 += *(P2 + 122)*(Lq2); // Lwq
			L2 += *(P2 + 143)*(Le2); // Lwe
			L2 += *(P2 + 164)*(Lg2); // Lwg
			L2 += *(P2 + 185)*(Lh2); // Lwh
			L2 += *(P2 + 206)*(Li2); // Lwi
			L2 += *(P2 + 227)*(Ll2); // Lwl
			L2 += *(P2 + 248)*(Lk2); // Lwk
			L2 += *(P2 + 269)*(Lm2); // Lwm
			L2 += *(P2 + 290)*(Lf2); // Lwf
			L2 += *(P2 + 311)*(Lp2); // Lwp
			L2 += *(P2 + 332)*(Lz2); // Lws
			L2 += *(P2 + 353)*(Lt2); // Lwt
			L2 += *(P2 + 374)*(Lw2); // Lww
			L2 += *(P2 + 395)*(Ly2); // Lwy
			L2 += *(P2 + 416)*(Lv2); // Lwv
			*(Ls3 + 17) = L2;
			
			// L(Y)
			L2 = *(P2 + 18)*(La2); // Lya
			L2 += *(P2 + 39)*(Lr2); // Lyr
			L2 += *(P2 + 60)*(Ln2); // Lyn
			L2 += *(P2 + 81)*(Ld2); // Lyd
			L2 += *(P2 + 102)*(Lc2); // Lyc
			L2 += *(P2 + 123)*(Lq2); // Lyq
			L2 += *(P2 + 144)*(Le2); // Lye
			L2 += *(P2 + 165)*(Lg2); // Lyg
			L2 += *(P2 + 186)*(Lh2); // Lyh
			L2 += *(P2 + 207)*(Li2); // Lyi
			L2 += *(P2 + 228)*(Ll2); // Lyl
			L2 += *(P2 + 249)*(Lk2); // Lyk
			L2 += *(P2 + 270)*(Lm2); // Lym
			L2 += *(P2 + 291)*(Lf2); // Lyf
			L2 += *(P2 + 312)*(Lp2); // Lyp
			L2 += *(P2 + 333)*(Lz2); // Lys
			L2 += *(P2 + 354)*(Lt2); // Lyt
			L2 += *(P2 + 375)*(Lw2); // Lyw
			L2 += *(P2 + 396)*(Ly2); // Lyy
			L2 += *(P2 + 417)*(Lv2); // Lyv
			*(Ls3 + 18) = L2;
			
			// L(V)
			L2 = *(P2 + 19)*(La2); // Lva
			L2 += *(P2 + 40)*(Lr2); // Lvr
			L2 += *(P2 + 61)*(Ln2); // Lvn
			L2 += *(P2 + 82)*(Ld2); // Lvd
			L2 += *(P2 + 103)*(Lc2); // Lvc
			L2 += *(P2 + 124)*(Lq2); // Lvq
			L2 += *(P2 + 145)*(Le2); // Lve
			L2 += *(P2 + 166)*(Lg2); // Lvg
			L2 += *(P2 + 187)*(Lh2); // Lvh
			L2 += *(P2 + 208)*(Li2); // Lvi
			L2 += *(P2 + 229)*(Ll2); // Lvl
			L2 += *(P2 + 250)*(Lk2); // Lvk
			L2 += *(P2 + 271)*(Lm2); // Lvm
			L2 += *(P2 + 292)*(Lf2); // Lvf
			L2 += *(P2 + 313)*(Lp2); // Lvp
			L2 += *(P2 + 334)*(Lz2); // Lvs
			L2 += *(P2 + 355)*(Lt2); // Lvt
			L2 += *(P2 + 376)*(Lw2); // Lvw
			L2 += *(P2 + 397)*(Ly2); // Lvy
			L2 += *(P2 + 418)*(Lv2); // Lvv
			*(Ls3 + 19) = L2;
			
			if (root &&
				(La1 != 0 ||
				Lr1 != 0 ||
				Ln1 != 0 ||
				Ld1 != 0 ||
				Lc1 != 0 ||
				Lq1 != 0 ||
				Le1 != 0 ||
				Lg1 != 0 ||
				Lh1 != 0 ||
				Li1 != 0 ||
				Ll1 != 0 ||
				Lk1 != 0 ||
				Lm1 != 0 ||
				Lf1 != 0 ||
				Lp1 != 0 ||
				Lz1 != 0 ||
				Lt1 != 0 ||
				Lw1 != 0 ||
				Ly1 != 0 ||
				Lv1 != 0)) {
				*(Ls3) *= La1;
				*(Ls3 + 1) *= Lr1;
				*(Ls3 + 2) *= Ln1;
				*(Ls3 + 3) *= Ld1;
				*(Ls3 + 4) *= Lc1;
				*(Ls3 + 5) *= Lq1;
				*(Ls3 + 6) *= Le1;
				*(Ls3 + 7) *= Lg1;
				*(Ls3 + 8) *= Lh1;
				*(Ls3 + 9) *= Li1;
				*(Ls3 + 10) *= Ll1;
				*(Ls3 + 11) *= Lk1;
				*(Ls3 + 12) *= Lm1;
				*(Ls3 + 13) *= Lf1;
				*(Ls3 + 14) *= Lp1;
				*(Ls3 + 15) *= Lz1;
				*(Ls3 + 16) *= Lt1;
				*(Ls3 + 17) *= Lw1;
				*(Ls3 + 18) *= Ly1;
				*(Ls3 + 19) *= Lv1;
				*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			} else {
				*(Ls3 + 21) = *(Ls2 + 21);
			}
			
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon) ||
				(*(Ls3 + 5) > 0 && *(Ls3 + 5) < inv_epsilon) ||
				(*(Ls3 + 6) > 0 && *(Ls3 + 6) < inv_epsilon) ||
				(*(Ls3 + 7) > 0 && *(Ls3 + 7) < inv_epsilon) ||
				(*(Ls3 + 8) > 0 && *(Ls3 + 8) < inv_epsilon) ||
				(*(Ls3 + 9) > 0 && *(Ls3 + 9) < inv_epsilon) ||
				(*(Ls3 + 10) > 0 && *(Ls3 + 10) < inv_epsilon) ||
				(*(Ls3 + 11) > 0 && *(Ls3 + 11) < inv_epsilon) ||
				(*(Ls3 + 12) > 0 && *(Ls3 + 12) < inv_epsilon) ||
				(*(Ls3 + 13) > 0 && *(Ls3 + 13) < inv_epsilon) ||
				(*(Ls3 + 14) > 0 && *(Ls3 + 14) < inv_epsilon) ||
				(*(Ls3 + 15) > 0 && *(Ls3 + 15) < inv_epsilon) ||
				(*(Ls3 + 16) > 0 && *(Ls3 + 16) < inv_epsilon) ||
				(*(Ls3 + 17) > 0 && *(Ls3 + 17) < inv_epsilon) ||
				(*(Ls3 + 18) > 0 && *(Ls3 + 18) < inv_epsilon) ||
				(*(Ls3 + 19) > 0 && *(Ls3 + 19) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) *= epsilon;
				*(Ls3 + 6) *= epsilon;
				*(Ls3 + 7) *= epsilon;
				*(Ls3 + 8) *= epsilon;
				*(Ls3 + 9) *= epsilon;
				*(Ls3 + 10) *= epsilon;
				*(Ls3 + 11) *= epsilon;
				*(Ls3 + 12) *= epsilon;
				*(Ls3 + 13) *= epsilon;
				*(Ls3 + 14) *= epsilon;
				*(Ls3 + 15) *= epsilon;
				*(Ls3 + 16) *= epsilon;
				*(Ls3 + 17) *= epsilon;
				*(Ls3 + 18) *= epsilon;
				*(Ls3 + 19) *= epsilon;
				*(Ls3 + 21) += 1;
			}
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lr1;
			*(Ls3 + 2) = Ln1;
			*(Ls3 + 3) = Ld1;
			*(Ls3 + 4) = Lc1;
			*(Ls3 + 5) = Lq1;
			*(Ls3 + 6) = Le1;
			*(Ls3 + 7) = Lg1;
			*(Ls3 + 8) = Lh1;
			*(Ls3 + 9) = Li1;
			*(Ls3 + 10) = Ll1;
			*(Ls3 + 11) = Lk1;
			*(Ls3 + 12) = Lm1;
			*(Ls3 + 13) = Lf1;
			*(Ls3 + 14) = Lp1;
			*(Ls3 + 15) = Lz1;
			*(Ls3 + 16) = Lt1;
			*(Ls3 + 17) = Lw1;
			*(Ls3 + 18) = Ly1;
			*(Ls3 + 19) = Lv1;
			*(Ls3 + 21) = *(Ls1 + 21);
		}
	}
	*(Ls3 + 20) = 0;
}

static void L_unknown_5(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	double L1, L2;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lc1 = *(Ls1 + 1);
	const double Lg1 = *(Ls1 + 2);
	const double Lt1 = *(Ls1 + 3);
	const double Li1 = *(Ls1 + 4);
	const double La2 = *(Ls2);
	const double Lc2 = *(Ls2 + 1);
	const double Lg2 = *(Ls2 + 2);
	const double Lt2 = *(Ls2 + 3);
	const double Li2 = *(Ls2 + 4);
	
	if (root == 0 &&
		(La1 != 0 ||
		Lc1 != 0 ||
		Lg1 != 0 ||
		Lt1 != 0 ||
		Li1 != 0)) {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0 ||
			Li2 != 0) {
			// neither branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 5)*(Lc1); // Lac
			L1 += *(P1 + 10)*(Lg1); // Lag
			L1 += *(P1 + 15)*(Lt1); // Lat
			L1 += *(P1 + 20)*(Li1); // La-
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 5)*(Lc2); // Lac
			L2 += *(P2 + 10)*(Lg2); // Lag
			L2 += *(P2 + 15)*(Lt2); // Lat
			L2 += *(P2 + 20)*(Li2); // La-
			*(Ls3) = L1*L2;
			
			// L(C)
			L1 = *(P1 + 1)*(La1); // Lca
			L1 += *(P1 + 6)*(Lc1); // Lcc
			L1 += *(P1 + 11)*(Lg1); // Lcg
			L1 += *(P1 + 16)*(Lt1); // Lct
			L1 += *(P1 + 21)*(Li1); // Lc-
			L2 = *(P2 + 1)*(La2); // Lca
			L2 += *(P2 + 6)*(Lc2); // Lcc
			L2 += *(P2 + 11)*(Lg2); // Lcg
			L2 += *(P2 + 16)*(Lt2); // Lct
			L2 += *(P2 + 21)*(Li2); // Lc-
			*(Ls3 + 1) = L1*L2;
			
			// L(G)
			L1 = *(P1 + 2)*(La1); // Lga
			L1 += *(P1 + 7)*(Lc1); // Lgc
			L1 += *(P1 + 12)*(Lg1); // Lgg
			L1 += *(P1 + 17)*(Lt1); // Lgt
			L1 += *(P1 + 22)*(Li1); // Lg-
			L2 = *(P2 + 2)*(La2); // Lga
			L2 += *(P2 + 7)*(Lc2); // Lgc
			L2 += *(P2 + 12)*(Lg2); // Lgg
			L2 += *(P2 + 17)*(Lt2); // Lgt
			L2 += *(P2 + 22)*(Li2); // Lg-
			*(Ls3 + 2) = L1*L2;
			
			// L(T)
			L1 = *(P1 + 3)*(La1); // Lta
			L1 += *(P1 + 8)*(Lc1); // Ltc
			L1 += *(P1 + 13)*(Lg1); // Ltg
			L1 += *(P1 + 18)*(Lt1); // Ltt
			L1 += *(P1 + 23)*(Li1); // Lt-
			L2 = *(P2 + 3)*(La2); // Lta
			L2 += *(P2 + 8)*(Lc2); // Ltc
			L2 += *(P2 + 13)*(Lg2); // Ltg
			L2 += *(P2 + 18)*(Lt2); // Ltt
			L2 += *(P2 + 23)*(Li2); // Lt-
			*(Ls3 + 3) = L1*L2;
			
			// L(-)
			L1 = *(P1 + 4)*(La1); // L-a
			L1 += *(P1 + 9)*(Lc1); // L-c
			L1 += *(P1 + 14)*(Lg1); // L-g
			L1 += *(P1 + 19)*(Lt1); // L-t
			L1 += *(P1 + 24)*(Li1); // L--
			L2 = *(P2 + 4)*(La2); // L-a
			L2 += *(P2 + 9)*(Lc2); // L-c
			L2 += *(P2 + 14)*(Lg2); // L-g
			L2 += *(P2 + 19)*(Lt2); // L-t
			L2 += *(P2 + 24)*(Li2); // L--
			*(Ls3 + 4) = L1*L2;
			
			*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) += 1;
			}
		} else {
			// second branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 5)*(Lc1); // Lac
			L1 += *(P1 + 10)*(Lg1); // Lag
			L1 += *(P1 + 15)*(Lt1); // Lat
			L1 += *(P1 + 20)*(Li1); // La-
			*(Ls3) = L1;
			
			// L(C)
			L1 = *(P1 + 1)*(La1); // Lca
			L1 += *(P1 + 6)*(Lc1); // Lcc
			L1 += *(P1 + 11)*(Lg1); // Lcg
			L1 += *(P1 + 16)*(Lt1); // Lct
			L1 += *(P1 + 21)*(Li1); // Lc-
			*(Ls3 + 1) = L1;
			
			// L(G)
			L1 = *(P1 + 2)*(La1); // Lga
			L1 += *(P1 + 7)*(Lc1); // Lgc
			L1 += *(P1 + 12)*(Lg1); // Lgg
			L1 += *(P1 + 17)*(Lt1); // Lgt
			L1 += *(P1 + 22)*(Li1); // Lg-
			*(Ls3 + 2) = L1;
			
			// L(T)
			L1 = *(P1 + 3)*(La1); // Lta
			L1 += *(P1 + 8)*(Lc1); // Ltc
			L1 += *(P1 + 13)*(Lg1); // Ltg
			L1 += *(P1 + 18)*(Lt1); // Ltt
			L1 += *(P1 + 23)*(Li1); // Lt-
			*(Ls3 + 3) = L1;
			
			// L(-)
			L1 = *(P1 + 4)*(La1); // L-a
			L1 += *(P1 + 9)*(Lc1); // L-c
			L1 += *(P1 + 14)*(Lg1); // L-g
			L1 += *(P1 + 19)*(Lt1); // L-t
			L1 += *(P1 + 24)*(Li1); // L--
			*(Ls3 + 4) = L1;
			
			*(Ls3 + 5) = *(Ls1 + 5);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) += 1;
			}
		}
	} else {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0 ||
			Li2 != 0) {
			// first branch can be disregarded
			
			// L(A)
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 5)*(Lc2); // Lac
			L2 += *(P2 + 10)*(Lg2); // Lag
			L2 += *(P2 + 15)*(Lt2); // Lat
			L2 += *(P2 + 20)*(Li2); // La-
			*(Ls3) = L2;
			
			// L(C)
			L2 = *(P2 + 1)*(La2); // Lca
			L2 += *(P2 + 6)*(Lc2); // Lcc
			L2 += *(P2 + 11)*(Lg2); // Lcg
			L2 += *(P2 + 16)*(Lt2); // Lct
			L2 += *(P2 + 21)*(Li2); // Lc-
			*(Ls3 + 1) = L2;
			
			// L(G)
			L2 = *(P2 + 2)*(La2); // Lga
			L2 += *(P2 + 7)*(Lc2); // Lgc
			L2 += *(P2 + 12)*(Lg2); // Lgg
			L2 += *(P2 + 17)*(Lt2); // Lgt
			L2 += *(P2 + 22)*(Li2); // Lg-
			*(Ls3 + 2) = L2;
			
			// L(T)
			L2 = *(P2 + 3)*(La2); // Lta
			L2 += *(P2 + 8)*(Lc2); // Ltc
			L2 += *(P2 + 13)*(Lg2); // Ltg
			L2 += *(P2 + 18)*(Lt2); // Ltt
			L2 += *(P2 + 23)*(Li2); // Lt-
			*(Ls3 + 3) = L2;
			
			// L(-)
			L2 = *(P2 + 4)*(La2); // L-a
			L2 += *(P2 + 9)*(Lc2); // L-c
			L2 += *(P2 + 14)*(Lg2); // L-g
			L2 += *(P2 + 19)*(Lt2); // L-t
			L2 += *(P2 + 24)*(Li2); // L--
			*(Ls3 + 4) = L2;
			
			if (root &&
				(La1 != 0 ||
				Lc1 != 0 ||
				Lg1 != 0 ||
				Lt1 != 0 ||
				Li1 != 0)) {
				*(Ls3) *= La1;
				*(Ls3 + 1) *= Lc1;
				*(Ls3 + 2) *= Lg1;
				*(Ls3 + 3) *= Lt1;
				*(Ls3 + 4) *= Li1;
				*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
			} else {
				*(Ls3 + 5) = *(Ls2 + 5);
			}
			
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) += 1;
			}
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lc1;
			*(Ls3 + 2) = Lg1;
			*(Ls3 + 3) = Lt1;
			*(Ls3 + 4) = Li1;
			*(Ls3 + 5) = *(Ls1 + 5);
		}
	}
}

static void L_unknown_AA_5(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	double L1, L2;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lr1 = *(Ls1 + 1);
	const double Ln1 = *(Ls1 + 2);
	const double Ld1 = *(Ls1 + 3);
	const double Lc1 = *(Ls1 + 4);
	const double Lq1 = *(Ls1 + 5);
	const double Le1 = *(Ls1 + 6);
	const double Lg1 = *(Ls1 + 7);
	const double Lh1 = *(Ls1 + 8);
	const double Li1 = *(Ls1 + 9);
	const double Ll1 = *(Ls1 + 10);
	const double Lk1 = *(Ls1 + 11);
	const double Lm1 = *(Ls1 + 12);
	const double Lf1 = *(Ls1 + 13);
	const double Lp1 = *(Ls1 + 14);
	const double Lz1 = *(Ls1 + 15);
	const double Lt1 = *(Ls1 + 16);
	const double Lw1 = *(Ls1 + 17);
	const double Ly1 = *(Ls1 + 18);
	const double Lv1 = *(Ls1 + 19);
	const double Lb1 = *(Ls1 + 20);
	const double La2 = *(Ls2);
	const double Lr2 = *(Ls2 + 1);
	const double Ln2 = *(Ls2 + 2);
	const double Ld2 = *(Ls2 + 3);
	const double Lc2 = *(Ls2 + 4);
	const double Lq2 = *(Ls2 + 5);
	const double Le2 = *(Ls2 + 6);
	const double Lg2 = *(Ls2 + 7);
	const double Lh2 = *(Ls2 + 8);
	const double Li2 = *(Ls2 + 9);
	const double Ll2 = *(Ls2 + 10);
	const double Lk2 = *(Ls2 + 11);
	const double Lm2 = *(Ls2 + 12);
	const double Lf2 = *(Ls2 + 13);
	const double Lp2 = *(Ls2 + 14);
	const double Lz2 = *(Ls2 + 15);
	const double Lt2 = *(Ls2 + 16);
	const double Lw2 = *(Ls2 + 17);
	const double Ly2 = *(Ls2 + 18);
	const double Lv2 = *(Ls2 + 19);
	const double Lb2 = *(Ls2 + 20);
	
	if (root == 0 &&
		(La1 != 0 ||
		Lr1 != 0 ||
		Ln1 != 0 ||
		Ld1 != 0 ||
		Lc1 != 0 ||
		Lq1 != 0 ||
		Le1 != 0 ||
		Lg1 != 0 ||
		Lh1 != 0 ||
		Li1 != 0 ||
		Ll1 != 0 ||
		Lk1 != 0 ||
		Lm1 != 0 ||
		Lf1 != 0 ||
		Lp1 != 0 ||
		Lz1 != 0 ||
		Lt1 != 0 ||
		Lw1 != 0 ||
		Ly1 != 0 ||
		Lv1 != 0 ||
		Lb1 != 0)) {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0 ||
			Lb2 != 0) {
			// neither branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 21)*(Lr1); // Lar
			L1 += *(P1 + 42)*(Ln1); // Lan
			L1 += *(P1 + 63)*(Ld1); // Lad
			L1 += *(P1 + 84)*(Lc1); // Lac
			L1 += *(P1 + 105)*(Lq1); // Laq
			L1 += *(P1 + 126)*(Le1); // Lae
			L1 += *(P1 + 147)*(Lg1); // Lag
			L1 += *(P1 + 168)*(Lh1); // Lah
			L1 += *(P1 + 189)*(Li1); // Lai
			L1 += *(P1 + 210)*(Ll1); // Lal
			L1 += *(P1 + 231)*(Lk1); // Lak
			L1 += *(P1 + 252)*(Lm1); // Lam
			L1 += *(P1 + 273)*(Lf1); // Laf
			L1 += *(P1 + 294)*(Lp1); // Lap
			L1 += *(P1 + 315)*(Lz1); // Las
			L1 += *(P1 + 336)*(Lt1); // Lat
			L1 += *(P1 + 357)*(Lw1); // Law
			L1 += *(P1 + 378)*(Ly1); // Lay
			L1 += *(P1 + 399)*(Lv1); // Lav
			L1 += *(P1 + 420)*(Lb1); // Lab
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 21)*(Lr2); // Lar
			L2 += *(P2 + 42)*(Ln2); // Lan
			L2 += *(P2 + 63)*(Ld2); // Lad
			L2 += *(P2 + 84)*(Lc2); // Lac
			L2 += *(P2 + 105)*(Lq2); // Laq
			L2 += *(P2 + 126)*(Le2); // Lae
			L2 += *(P2 + 147)*(Lg2); // Lag
			L2 += *(P2 + 168)*(Lh2); // Lah
			L2 += *(P2 + 189)*(Li2); // Lai
			L2 += *(P2 + 210)*(Ll2); // Lal
			L2 += *(P2 + 231)*(Lk2); // Lak
			L2 += *(P2 + 252)*(Lm2); // Lam
			L2 += *(P2 + 273)*(Lf2); // Laf
			L2 += *(P2 + 294)*(Lp2); // Lap
			L2 += *(P2 + 315)*(Lz2); // Las
			L2 += *(P2 + 336)*(Lt2); // Lat
			L2 += *(P2 + 357)*(Lw2); // Law
			L2 += *(P2 + 378)*(Ly2); // Lay
			L2 += *(P2 + 399)*(Lv2); // Lav
			L2 += *(P2 + 420)*(Lb2); // Lab
			*(Ls3 + 0) = L1*L2;
			
			// L(R)
			L1 = *(P1 + 1)*(La1); // Lra
			L1 += *(P1 + 22)*(Lr1); // Lrr
			L1 += *(P1 + 43)*(Ln1); // Lrn
			L1 += *(P1 + 64)*(Ld1); // Lrd
			L1 += *(P1 + 85)*(Lc1); // Lrc
			L1 += *(P1 + 106)*(Lq1); // Lrq
			L1 += *(P1 + 127)*(Le1); // Lre
			L1 += *(P1 + 148)*(Lg1); // Lrg
			L1 += *(P1 + 169)*(Lh1); // Lrh
			L1 += *(P1 + 190)*(Li1); // Lri
			L1 += *(P1 + 211)*(Ll1); // Lrl
			L1 += *(P1 + 232)*(Lk1); // Lrk
			L1 += *(P1 + 253)*(Lm1); // Lrm
			L1 += *(P1 + 274)*(Lf1); // Lrf
			L1 += *(P1 + 295)*(Lp1); // Lrp
			L1 += *(P1 + 316)*(Lz1); // Lrs
			L1 += *(P1 + 337)*(Lt1); // Lrt
			L1 += *(P1 + 358)*(Lw1); // Lrw
			L1 += *(P1 + 379)*(Ly1); // Lry
			L1 += *(P1 + 400)*(Lv1); // Lrv
			L1 += *(P1 + 421)*(Lb1); // Lrb
			L2 = *(P2 + 1)*(La2); // Lra
			L2 += *(P2 + 22)*(Lr2); // Lrr
			L2 += *(P2 + 43)*(Ln2); // Lrn
			L2 += *(P2 + 64)*(Ld2); // Lrd
			L2 += *(P2 + 85)*(Lc2); // Lrc
			L2 += *(P2 + 106)*(Lq2); // Lrq
			L2 += *(P2 + 127)*(Le2); // Lre
			L2 += *(P2 + 148)*(Lg2); // Lrg
			L2 += *(P2 + 169)*(Lh2); // Lrh
			L2 += *(P2 + 190)*(Li2); // Lri
			L2 += *(P2 + 211)*(Ll2); // Lrl
			L2 += *(P2 + 232)*(Lk2); // Lrk
			L2 += *(P2 + 253)*(Lm2); // Lrm
			L2 += *(P2 + 274)*(Lf2); // Lrf
			L2 += *(P2 + 295)*(Lp2); // Lrp
			L2 += *(P2 + 316)*(Lz2); // Lrs
			L2 += *(P2 + 337)*(Lt2); // Lrt
			L2 += *(P2 + 358)*(Lw2); // Lrw
			L2 += *(P2 + 379)*(Ly2); // Lry
			L2 += *(P2 + 400)*(Lv2); // Lrv
			L2 += *(P2 + 421)*(Lb2); // Lrb
			*(Ls3 + 1) = L1*L2;
			
			// L(N)
			L1 = *(P1 + 2)*(La1); // Lna
			L1 += *(P1 + 23)*(Lr1); // Lnr
			L1 += *(P1 + 44)*(Ln1); // Lnn
			L1 += *(P1 + 65)*(Ld1); // Lnd
			L1 += *(P1 + 86)*(Lc1); // Lnc
			L1 += *(P1 + 107)*(Lq1); // Lnq
			L1 += *(P1 + 128)*(Le1); // Lne
			L1 += *(P1 + 149)*(Lg1); // Lng
			L1 += *(P1 + 170)*(Lh1); // Lnh
			L1 += *(P1 + 191)*(Li1); // Lni
			L1 += *(P1 + 212)*(Ll1); // Lnl
			L1 += *(P1 + 233)*(Lk1); // Lnk
			L1 += *(P1 + 254)*(Lm1); // Lnm
			L1 += *(P1 + 275)*(Lf1); // Lnf
			L1 += *(P1 + 296)*(Lp1); // Lnp
			L1 += *(P1 + 317)*(Lz1); // Lns
			L1 += *(P1 + 338)*(Lt1); // Lnt
			L1 += *(P1 + 359)*(Lw1); // Lnw
			L1 += *(P1 + 380)*(Ly1); // Lny
			L1 += *(P1 + 401)*(Lv1); // Lnv
			L1 += *(P1 + 422)*(Lb1); // Lnb
			L2 = *(P2 + 2)*(La2); // Lna
			L2 += *(P2 + 23)*(Lr2); // Lnr
			L2 += *(P2 + 44)*(Ln2); // Lnn
			L2 += *(P2 + 65)*(Ld2); // Lnd
			L2 += *(P2 + 86)*(Lc2); // Lnc
			L2 += *(P2 + 107)*(Lq2); // Lnq
			L2 += *(P2 + 128)*(Le2); // Lne
			L2 += *(P2 + 149)*(Lg2); // Lng
			L2 += *(P2 + 170)*(Lh2); // Lnh
			L2 += *(P2 + 191)*(Li2); // Lni
			L2 += *(P2 + 212)*(Ll2); // Lnl
			L2 += *(P2 + 233)*(Lk2); // Lnk
			L2 += *(P2 + 254)*(Lm2); // Lnm
			L2 += *(P2 + 275)*(Lf2); // Lnf
			L2 += *(P2 + 296)*(Lp2); // Lnp
			L2 += *(P2 + 317)*(Lz2); // Lns
			L2 += *(P2 + 338)*(Lt2); // Lnt
			L2 += *(P2 + 359)*(Lw2); // Lnw
			L2 += *(P2 + 380)*(Ly2); // Lny
			L2 += *(P2 + 401)*(Lv2); // Lnv
			L2 += *(P2 + 422)*(Lb2); // Lnb
			*(Ls3 + 2) = L1*L2;
			
			// L(D)
			L1 = *(P1 + 3)*(La1); // Lda
			L1 += *(P1 + 24)*(Lr1); // Ldr
			L1 += *(P1 + 45)*(Ln1); // Ldn
			L1 += *(P1 + 66)*(Ld1); // Ldd
			L1 += *(P1 + 87)*(Lc1); // Ldc
			L1 += *(P1 + 108)*(Lq1); // Ldq
			L1 += *(P1 + 129)*(Le1); // Lde
			L1 += *(P1 + 150)*(Lg1); // Ldg
			L1 += *(P1 + 171)*(Lh1); // Ldh
			L1 += *(P1 + 192)*(Li1); // Ldi
			L1 += *(P1 + 213)*(Ll1); // Ldl
			L1 += *(P1 + 234)*(Lk1); // Ldk
			L1 += *(P1 + 255)*(Lm1); // Ldm
			L1 += *(P1 + 276)*(Lf1); // Ldf
			L1 += *(P1 + 297)*(Lp1); // Ldp
			L1 += *(P1 + 318)*(Lz1); // Lds
			L1 += *(P1 + 339)*(Lt1); // Ldt
			L1 += *(P1 + 360)*(Lw1); // Ldw
			L1 += *(P1 + 381)*(Ly1); // Ldy
			L1 += *(P1 + 402)*(Lv1); // Ldv
			L1 += *(P1 + 423)*(Lb1); // Ldb
			L2 = *(P2 + 3)*(La2); // Lda
			L2 += *(P2 + 24)*(Lr2); // Ldr
			L2 += *(P2 + 45)*(Ln2); // Ldn
			L2 += *(P2 + 66)*(Ld2); // Ldd
			L2 += *(P2 + 87)*(Lc2); // Ldc
			L2 += *(P2 + 108)*(Lq2); // Ldq
			L2 += *(P2 + 129)*(Le2); // Lde
			L2 += *(P2 + 150)*(Lg2); // Ldg
			L2 += *(P2 + 171)*(Lh2); // Ldh
			L2 += *(P2 + 192)*(Li2); // Ldi
			L2 += *(P2 + 213)*(Ll2); // Ldl
			L2 += *(P2 + 234)*(Lk2); // Ldk
			L2 += *(P2 + 255)*(Lm2); // Ldm
			L2 += *(P2 + 276)*(Lf2); // Ldf
			L2 += *(P2 + 297)*(Lp2); // Ldp
			L2 += *(P2 + 318)*(Lz2); // Lds
			L2 += *(P2 + 339)*(Lt2); // Ldt
			L2 += *(P2 + 360)*(Lw2); // Ldw
			L2 += *(P2 + 381)*(Ly2); // Ldy
			L2 += *(P2 + 402)*(Lv2); // Ldv
			L2 += *(P2 + 423)*(Lb2); // Ldb
			*(Ls3 + 3) = L1*L2;
			
			// L(C)
			L1 = *(P1 + 4)*(La1); // Lca
			L1 += *(P1 + 25)*(Lr1); // Lcr
			L1 += *(P1 + 46)*(Ln1); // Lcn
			L1 += *(P1 + 67)*(Ld1); // Lcd
			L1 += *(P1 + 88)*(Lc1); // Lcc
			L1 += *(P1 + 109)*(Lq1); // Lcq
			L1 += *(P1 + 130)*(Le1); // Lce
			L1 += *(P1 + 151)*(Lg1); // Lcg
			L1 += *(P1 + 172)*(Lh1); // Lch
			L1 += *(P1 + 193)*(Li1); // Lci
			L1 += *(P1 + 214)*(Ll1); // Lcl
			L1 += *(P1 + 235)*(Lk1); // Lck
			L1 += *(P1 + 256)*(Lm1); // Lcm
			L1 += *(P1 + 277)*(Lf1); // Lcf
			L1 += *(P1 + 298)*(Lp1); // Lcp
			L1 += *(P1 + 319)*(Lz1); // Lcs
			L1 += *(P1 + 340)*(Lt1); // Lct
			L1 += *(P1 + 361)*(Lw1); // Lcw
			L1 += *(P1 + 382)*(Ly1); // Lcy
			L1 += *(P1 + 403)*(Lv1); // Lcv
			L1 += *(P1 + 424)*(Lb1); // Lcb
			L2 = *(P2 + 4)*(La2); // Lca
			L2 += *(P2 + 25)*(Lr2); // Lcr
			L2 += *(P2 + 46)*(Ln2); // Lcn
			L2 += *(P2 + 67)*(Ld2); // Lcd
			L2 += *(P2 + 88)*(Lc2); // Lcc
			L2 += *(P2 + 109)*(Lq2); // Lcq
			L2 += *(P2 + 130)*(Le2); // Lce
			L2 += *(P2 + 151)*(Lg2); // Lcg
			L2 += *(P2 + 172)*(Lh2); // Lch
			L2 += *(P2 + 193)*(Li2); // Lci
			L2 += *(P2 + 214)*(Ll2); // Lcl
			L2 += *(P2 + 235)*(Lk2); // Lck
			L2 += *(P2 + 256)*(Lm2); // Lcm
			L2 += *(P2 + 277)*(Lf2); // Lcf
			L2 += *(P2 + 298)*(Lp2); // Lcp
			L2 += *(P2 + 319)*(Lz2); // Lcs
			L2 += *(P2 + 340)*(Lt2); // Lct
			L2 += *(P2 + 361)*(Lw2); // Lcw
			L2 += *(P2 + 382)*(Ly2); // Lcy
			L2 += *(P2 + 403)*(Lv2); // Lcv
			L2 += *(P2 + 424)*(Lb2); // Lcb
			*(Ls3 + 4) = L1*L2;
			
			// L(Q)
			L1 = *(P1 + 5)*(La1); // Lqa
			L1 += *(P1 + 26)*(Lr1); // Lqr
			L1 += *(P1 + 47)*(Ln1); // Lqn
			L1 += *(P1 + 68)*(Ld1); // Lqd
			L1 += *(P1 + 89)*(Lc1); // Lqc
			L1 += *(P1 + 110)*(Lq1); // Lqq
			L1 += *(P1 + 131)*(Le1); // Lqe
			L1 += *(P1 + 152)*(Lg1); // Lqg
			L1 += *(P1 + 173)*(Lh1); // Lqh
			L1 += *(P1 + 194)*(Li1); // Lqi
			L1 += *(P1 + 215)*(Ll1); // Lql
			L1 += *(P1 + 236)*(Lk1); // Lqk
			L1 += *(P1 + 257)*(Lm1); // Lqm
			L1 += *(P1 + 278)*(Lf1); // Lqf
			L1 += *(P1 + 299)*(Lp1); // Lqp
			L1 += *(P1 + 320)*(Lz1); // Lqs
			L1 += *(P1 + 341)*(Lt1); // Lqt
			L1 += *(P1 + 362)*(Lw1); // Lqw
			L1 += *(P1 + 383)*(Ly1); // Lqy
			L1 += *(P1 + 404)*(Lv1); // Lqv
			L1 += *(P1 + 425)*(Lb1); // Lqb
			L2 = *(P2 + 5)*(La2); // Lqa
			L2 += *(P2 + 26)*(Lr2); // Lqr
			L2 += *(P2 + 47)*(Ln2); // Lqn
			L2 += *(P2 + 68)*(Ld2); // Lqd
			L2 += *(P2 + 89)*(Lc2); // Lqc
			L2 += *(P2 + 110)*(Lq2); // Lqq
			L2 += *(P2 + 131)*(Le2); // Lqe
			L2 += *(P2 + 152)*(Lg2); // Lqg
			L2 += *(P2 + 173)*(Lh2); // Lqh
			L2 += *(P2 + 194)*(Li2); // Lqi
			L2 += *(P2 + 215)*(Ll2); // Lql
			L2 += *(P2 + 236)*(Lk2); // Lqk
			L2 += *(P2 + 257)*(Lm2); // Lqm
			L2 += *(P2 + 278)*(Lf2); // Lqf
			L2 += *(P2 + 299)*(Lp2); // Lqp
			L2 += *(P2 + 320)*(Lz2); // Lqs
			L2 += *(P2 + 341)*(Lt2); // Lqt
			L2 += *(P2 + 362)*(Lw2); // Lqw
			L2 += *(P2 + 383)*(Ly2); // Lqy
			L2 += *(P2 + 404)*(Lv2); // Lqv
			L2 += *(P2 + 425)*(Lb2); // Lqb
			*(Ls3 + 5) = L1*L2;
			
			// L(E)
			L1 = *(P1 + 6)*(La1); // Lea
			L1 += *(P1 + 27)*(Lr1); // Ler
			L1 += *(P1 + 48)*(Ln1); // Len
			L1 += *(P1 + 69)*(Ld1); // Led
			L1 += *(P1 + 90)*(Lc1); // Lec
			L1 += *(P1 + 111)*(Lq1); // Leq
			L1 += *(P1 + 132)*(Le1); // Lee
			L1 += *(P1 + 153)*(Lg1); // Leg
			L1 += *(P1 + 174)*(Lh1); // Leh
			L1 += *(P1 + 195)*(Li1); // Lei
			L1 += *(P1 + 216)*(Ll1); // Lel
			L1 += *(P1 + 237)*(Lk1); // Lek
			L1 += *(P1 + 258)*(Lm1); // Lem
			L1 += *(P1 + 279)*(Lf1); // Lef
			L1 += *(P1 + 300)*(Lp1); // Lep
			L1 += *(P1 + 321)*(Lz1); // Les
			L1 += *(P1 + 342)*(Lt1); // Let
			L1 += *(P1 + 363)*(Lw1); // Lew
			L1 += *(P1 + 384)*(Ly1); // Ley
			L1 += *(P1 + 405)*(Lv1); // Lev
			L1 += *(P1 + 426)*(Lb1); // Leb
			L2 = *(P2 + 6)*(La2); // Lea
			L2 += *(P2 + 27)*(Lr2); // Ler
			L2 += *(P2 + 48)*(Ln2); // Len
			L2 += *(P2 + 69)*(Ld2); // Led
			L2 += *(P2 + 90)*(Lc2); // Lec
			L2 += *(P2 + 111)*(Lq2); // Leq
			L2 += *(P2 + 132)*(Le2); // Lee
			L2 += *(P2 + 153)*(Lg2); // Leg
			L2 += *(P2 + 174)*(Lh2); // Leh
			L2 += *(P2 + 195)*(Li2); // Lei
			L2 += *(P2 + 216)*(Ll2); // Lel
			L2 += *(P2 + 237)*(Lk2); // Lek
			L2 += *(P2 + 258)*(Lm2); // Lem
			L2 += *(P2 + 279)*(Lf2); // Lef
			L2 += *(P2 + 300)*(Lp2); // Lep
			L2 += *(P2 + 321)*(Lz2); // Les
			L2 += *(P2 + 342)*(Lt2); // Let
			L2 += *(P2 + 363)*(Lw2); // Lew
			L2 += *(P2 + 384)*(Ly2); // Ley
			L2 += *(P2 + 405)*(Lv2); // Lev
			L2 += *(P2 + 426)*(Lb2); // Leb
			*(Ls3 + 6) = L1*L2;
			
			// L(G)
			L1 = *(P1 + 7)*(La1); // Lga
			L1 += *(P1 + 28)*(Lr1); // Lgr
			L1 += *(P1 + 49)*(Ln1); // Lgn
			L1 += *(P1 + 70)*(Ld1); // Lgd
			L1 += *(P1 + 91)*(Lc1); // Lgc
			L1 += *(P1 + 112)*(Lq1); // Lgq
			L1 += *(P1 + 133)*(Le1); // Lge
			L1 += *(P1 + 154)*(Lg1); // Lgg
			L1 += *(P1 + 175)*(Lh1); // Lgh
			L1 += *(P1 + 196)*(Li1); // Lgi
			L1 += *(P1 + 217)*(Ll1); // Lgl
			L1 += *(P1 + 238)*(Lk1); // Lgk
			L1 += *(P1 + 259)*(Lm1); // Lgm
			L1 += *(P1 + 280)*(Lf1); // Lgf
			L1 += *(P1 + 301)*(Lp1); // Lgp
			L1 += *(P1 + 322)*(Lz1); // Lgs
			L1 += *(P1 + 343)*(Lt1); // Lgt
			L1 += *(P1 + 364)*(Lw1); // Lgw
			L1 += *(P1 + 385)*(Ly1); // Lgy
			L1 += *(P1 + 406)*(Lv1); // Lgv
			L1 += *(P1 + 427)*(Lb1); // Lgb
			L2 = *(P2 + 7)*(La2); // Lga
			L2 += *(P2 + 28)*(Lr2); // Lgr
			L2 += *(P2 + 49)*(Ln2); // Lgn
			L2 += *(P2 + 70)*(Ld2); // Lgd
			L2 += *(P2 + 91)*(Lc2); // Lgc
			L2 += *(P2 + 112)*(Lq2); // Lgq
			L2 += *(P2 + 133)*(Le2); // Lge
			L2 += *(P2 + 154)*(Lg2); // Lgg
			L2 += *(P2 + 175)*(Lh2); // Lgh
			L2 += *(P2 + 196)*(Li2); // Lgi
			L2 += *(P2 + 217)*(Ll2); // Lgl
			L2 += *(P2 + 238)*(Lk2); // Lgk
			L2 += *(P2 + 259)*(Lm2); // Lgm
			L2 += *(P2 + 280)*(Lf2); // Lgf
			L2 += *(P2 + 301)*(Lp2); // Lgp
			L2 += *(P2 + 322)*(Lz2); // Lgs
			L2 += *(P2 + 343)*(Lt2); // Lgt
			L2 += *(P2 + 364)*(Lw2); // Lgw
			L2 += *(P2 + 385)*(Ly2); // Lgy
			L2 += *(P2 + 406)*(Lv2); // Lgv
			L2 += *(P2 + 427)*(Lb2); // Lgb
			*(Ls3 + 7) = L1*L2;
			
			// L(H)
			L1 = *(P1 + 8)*(La1); // Lha
			L1 += *(P1 + 29)*(Lr1); // Lhr
			L1 += *(P1 + 50)*(Ln1); // Lhn
			L1 += *(P1 + 71)*(Ld1); // Lhd
			L1 += *(P1 + 92)*(Lc1); // Lhc
			L1 += *(P1 + 113)*(Lq1); // Lhq
			L1 += *(P1 + 134)*(Le1); // Lhe
			L1 += *(P1 + 155)*(Lg1); // Lhg
			L1 += *(P1 + 176)*(Lh1); // Lhh
			L1 += *(P1 + 197)*(Li1); // Lhi
			L1 += *(P1 + 218)*(Ll1); // Lhl
			L1 += *(P1 + 239)*(Lk1); // Lhk
			L1 += *(P1 + 260)*(Lm1); // Lhm
			L1 += *(P1 + 281)*(Lf1); // Lhf
			L1 += *(P1 + 302)*(Lp1); // Lhp
			L1 += *(P1 + 323)*(Lz1); // Lhs
			L1 += *(P1 + 344)*(Lt1); // Lht
			L1 += *(P1 + 365)*(Lw1); // Lhw
			L1 += *(P1 + 386)*(Ly1); // Lhy
			L1 += *(P1 + 407)*(Lv1); // Lhv
			L1 += *(P1 + 428)*(Lb1); // Lhb
			L2 = *(P2 + 8)*(La2); // Lha
			L2 += *(P2 + 29)*(Lr2); // Lhr
			L2 += *(P2 + 50)*(Ln2); // Lhn
			L2 += *(P2 + 71)*(Ld2); // Lhd
			L2 += *(P2 + 92)*(Lc2); // Lhc
			L2 += *(P2 + 113)*(Lq2); // Lhq
			L2 += *(P2 + 134)*(Le2); // Lhe
			L2 += *(P2 + 155)*(Lg2); // Lhg
			L2 += *(P2 + 176)*(Lh2); // Lhh
			L2 += *(P2 + 197)*(Li2); // Lhi
			L2 += *(P2 + 218)*(Ll2); // Lhl
			L2 += *(P2 + 239)*(Lk2); // Lhk
			L2 += *(P2 + 260)*(Lm2); // Lhm
			L2 += *(P2 + 281)*(Lf2); // Lhf
			L2 += *(P2 + 302)*(Lp2); // Lhp
			L2 += *(P2 + 323)*(Lz2); // Lhs
			L2 += *(P2 + 344)*(Lt2); // Lht
			L2 += *(P2 + 365)*(Lw2); // Lhw
			L2 += *(P2 + 386)*(Ly2); // Lhy
			L2 += *(P2 + 407)*(Lv2); // Lhv
			L2 += *(P2 + 428)*(Lb2); // Lhb
			*(Ls3 + 8) = L1*L2;
			
			// L(I)
			L1 = *(P1 + 9)*(La1); // Lia
			L1 += *(P1 + 30)*(Lr1); // Lir
			L1 += *(P1 + 51)*(Ln1); // Lin
			L1 += *(P1 + 72)*(Ld1); // Lid
			L1 += *(P1 + 93)*(Lc1); // Lic
			L1 += *(P1 + 114)*(Lq1); // Liq
			L1 += *(P1 + 135)*(Le1); // Lie
			L1 += *(P1 + 156)*(Lg1); // Lig
			L1 += *(P1 + 177)*(Lh1); // Lih
			L1 += *(P1 + 198)*(Li1); // Lii
			L1 += *(P1 + 219)*(Ll1); // Lil
			L1 += *(P1 + 240)*(Lk1); // Lik
			L1 += *(P1 + 261)*(Lm1); // Lim
			L1 += *(P1 + 282)*(Lf1); // Lif
			L1 += *(P1 + 303)*(Lp1); // Lip
			L1 += *(P1 + 324)*(Lz1); // Lis
			L1 += *(P1 + 345)*(Lt1); // Lit
			L1 += *(P1 + 366)*(Lw1); // Liw
			L1 += *(P1 + 387)*(Ly1); // Liy
			L1 += *(P1 + 408)*(Lv1); // Liv
			L1 += *(P1 + 429)*(Lb1); // Lib
			L2 = *(P2 + 9)*(La2); // Lia
			L2 += *(P2 + 30)*(Lr2); // Lir
			L2 += *(P2 + 51)*(Ln2); // Lin
			L2 += *(P2 + 72)*(Ld2); // Lid
			L2 += *(P2 + 93)*(Lc2); // Lic
			L2 += *(P2 + 114)*(Lq2); // Liq
			L2 += *(P2 + 135)*(Le2); // Lie
			L2 += *(P2 + 156)*(Lg2); // Lig
			L2 += *(P2 + 177)*(Lh2); // Lih
			L2 += *(P2 + 198)*(Li2); // Lii
			L2 += *(P2 + 219)*(Ll2); // Lil
			L2 += *(P2 + 240)*(Lk2); // Lik
			L2 += *(P2 + 261)*(Lm2); // Lim
			L2 += *(P2 + 282)*(Lf2); // Lif
			L2 += *(P2 + 303)*(Lp2); // Lip
			L2 += *(P2 + 324)*(Lz2); // Lis
			L2 += *(P2 + 345)*(Lt2); // Lit
			L2 += *(P2 + 366)*(Lw2); // Liw
			L2 += *(P2 + 387)*(Ly2); // Liy
			L2 += *(P2 + 408)*(Lv2); // Liv
			L2 += *(P2 + 429)*(Lb2); // Lib
			*(Ls3 + 9) = L1*L2;
			
			// L(L)
			L1 = *(P1 + 10)*(La1); // Lla
			L1 += *(P1 + 31)*(Lr1); // Llr
			L1 += *(P1 + 52)*(Ln1); // Lln
			L1 += *(P1 + 73)*(Ld1); // Lld
			L1 += *(P1 + 94)*(Lc1); // Llc
			L1 += *(P1 + 115)*(Lq1); // Llq
			L1 += *(P1 + 136)*(Le1); // Lle
			L1 += *(P1 + 157)*(Lg1); // Llg
			L1 += *(P1 + 178)*(Lh1); // Llh
			L1 += *(P1 + 199)*(Li1); // Lli
			L1 += *(P1 + 220)*(Ll1); // Lll
			L1 += *(P1 + 241)*(Lk1); // Llk
			L1 += *(P1 + 262)*(Lm1); // Llm
			L1 += *(P1 + 283)*(Lf1); // Llf
			L1 += *(P1 + 304)*(Lp1); // Llp
			L1 += *(P1 + 325)*(Lz1); // Lls
			L1 += *(P1 + 346)*(Lt1); // Llt
			L1 += *(P1 + 367)*(Lw1); // Llw
			L1 += *(P1 + 388)*(Ly1); // Lly
			L1 += *(P1 + 409)*(Lv1); // Llv
			L1 += *(P1 + 430)*(Lb1); // Llb
			L2 = *(P2 + 10)*(La2); // Lla
			L2 += *(P2 + 31)*(Lr2); // Llr
			L2 += *(P2 + 52)*(Ln2); // Lln
			L2 += *(P2 + 73)*(Ld2); // Lld
			L2 += *(P2 + 94)*(Lc2); // Llc
			L2 += *(P2 + 115)*(Lq2); // Llq
			L2 += *(P2 + 136)*(Le2); // Lle
			L2 += *(P2 + 157)*(Lg2); // Llg
			L2 += *(P2 + 178)*(Lh2); // Llh
			L2 += *(P2 + 199)*(Li2); // Lli
			L2 += *(P2 + 220)*(Ll2); // Lll
			L2 += *(P2 + 241)*(Lk2); // Llk
			L2 += *(P2 + 262)*(Lm2); // Llm
			L2 += *(P2 + 283)*(Lf2); // Llf
			L2 += *(P2 + 304)*(Lp2); // Llp
			L2 += *(P2 + 325)*(Lz2); // Lls
			L2 += *(P2 + 346)*(Lt2); // Llt
			L2 += *(P2 + 367)*(Lw2); // Llw
			L2 += *(P2 + 388)*(Ly2); // Lly
			L2 += *(P2 + 409)*(Lv2); // Llv
			L2 += *(P2 + 430)*(Lb2); // Llb
			*(Ls3 + 10) = L1*L2;
			
			// L(K)
			L1 = *(P1 + 11)*(La1); // Lka
			L1 += *(P1 + 32)*(Lr1); // Lkr
			L1 += *(P1 + 53)*(Ln1); // Lkn
			L1 += *(P1 + 74)*(Ld1); // Lkd
			L1 += *(P1 + 95)*(Lc1); // Lkc
			L1 += *(P1 + 116)*(Lq1); // Lkq
			L1 += *(P1 + 137)*(Le1); // Lke
			L1 += *(P1 + 158)*(Lg1); // Lkg
			L1 += *(P1 + 179)*(Lh1); // Lkh
			L1 += *(P1 + 200)*(Li1); // Lki
			L1 += *(P1 + 221)*(Ll1); // Lkl
			L1 += *(P1 + 242)*(Lk1); // Lkk
			L1 += *(P1 + 263)*(Lm1); // Lkm
			L1 += *(P1 + 284)*(Lf1); // Lkf
			L1 += *(P1 + 305)*(Lp1); // Lkp
			L1 += *(P1 + 326)*(Lz1); // Lks
			L1 += *(P1 + 347)*(Lt1); // Lkt
			L1 += *(P1 + 368)*(Lw1); // Lkw
			L1 += *(P1 + 389)*(Ly1); // Lky
			L1 += *(P1 + 410)*(Lv1); // Lkv
			L1 += *(P1 + 431)*(Lb1); // Lkb
			L2 = *(P2 + 11)*(La2); // Lka
			L2 += *(P2 + 32)*(Lr2); // Lkr
			L2 += *(P2 + 53)*(Ln2); // Lkn
			L2 += *(P2 + 74)*(Ld2); // Lkd
			L2 += *(P2 + 95)*(Lc2); // Lkc
			L2 += *(P2 + 116)*(Lq2); // Lkq
			L2 += *(P2 + 137)*(Le2); // Lke
			L2 += *(P2 + 158)*(Lg2); // Lkg
			L2 += *(P2 + 179)*(Lh2); // Lkh
			L2 += *(P2 + 200)*(Li2); // Lki
			L2 += *(P2 + 221)*(Ll2); // Lkl
			L2 += *(P2 + 242)*(Lk2); // Lkk
			L2 += *(P2 + 263)*(Lm2); // Lkm
			L2 += *(P2 + 284)*(Lf2); // Lkf
			L2 += *(P2 + 305)*(Lp2); // Lkp
			L2 += *(P2 + 326)*(Lz2); // Lks
			L2 += *(P2 + 347)*(Lt2); // Lkt
			L2 += *(P2 + 368)*(Lw2); // Lkw
			L2 += *(P2 + 389)*(Ly2); // Lky
			L2 += *(P2 + 410)*(Lv2); // Lkv
			L2 += *(P2 + 431)*(Lb2); // Lkb
			*(Ls3 + 11) = L1*L2;
			
			// L(M)
			L1 = *(P1 + 12)*(La1); // Lma
			L1 += *(P1 + 33)*(Lr1); // Lmr
			L1 += *(P1 + 54)*(Ln1); // Lmn
			L1 += *(P1 + 75)*(Ld1); // Lmd
			L1 += *(P1 + 96)*(Lc1); // Lmc
			L1 += *(P1 + 117)*(Lq1); // Lmq
			L1 += *(P1 + 138)*(Le1); // Lme
			L1 += *(P1 + 159)*(Lg1); // Lmg
			L1 += *(P1 + 180)*(Lh1); // Lmh
			L1 += *(P1 + 201)*(Li1); // Lmi
			L1 += *(P1 + 222)*(Ll1); // Lml
			L1 += *(P1 + 243)*(Lk1); // Lmk
			L1 += *(P1 + 264)*(Lm1); // Lmm
			L1 += *(P1 + 285)*(Lf1); // Lmf
			L1 += *(P1 + 306)*(Lp1); // Lmp
			L1 += *(P1 + 327)*(Lz1); // Lms
			L1 += *(P1 + 348)*(Lt1); // Lmt
			L1 += *(P1 + 369)*(Lw1); // Lmw
			L1 += *(P1 + 390)*(Ly1); // Lmy
			L1 += *(P1 + 411)*(Lv1); // Lmv
			L1 += *(P1 + 432)*(Lb1); // Lmb
			L2 = *(P2 + 12)*(La2); // Lma
			L2 += *(P2 + 33)*(Lr2); // Lmr
			L2 += *(P2 + 54)*(Ln2); // Lmn
			L2 += *(P2 + 75)*(Ld2); // Lmd
			L2 += *(P2 + 96)*(Lc2); // Lmc
			L2 += *(P2 + 117)*(Lq2); // Lmq
			L2 += *(P2 + 138)*(Le2); // Lme
			L2 += *(P2 + 159)*(Lg2); // Lmg
			L2 += *(P2 + 180)*(Lh2); // Lmh
			L2 += *(P2 + 201)*(Li2); // Lmi
			L2 += *(P2 + 222)*(Ll2); // Lml
			L2 += *(P2 + 243)*(Lk2); // Lmk
			L2 += *(P2 + 264)*(Lm2); // Lmm
			L2 += *(P2 + 285)*(Lf2); // Lmf
			L2 += *(P2 + 306)*(Lp2); // Lmp
			L2 += *(P2 + 327)*(Lz2); // Lms
			L2 += *(P2 + 348)*(Lt2); // Lmt
			L2 += *(P2 + 369)*(Lw2); // Lmw
			L2 += *(P2 + 390)*(Ly2); // Lmy
			L2 += *(P2 + 411)*(Lv2); // LmvLz
			L2 += *(P2 + 432)*(Lb2); // Lmb
			*(Ls3 + 12) = L1*L2;
			
			// L(F)
			L1 = *(P1 + 13)*(La1); // Lfa
			L1 += *(P1 + 34)*(Lr1); // Lfr
			L1 += *(P1 + 55)*(Ln1); // Lfn
			L1 += *(P1 + 76)*(Ld1); // Lfd
			L1 += *(P1 + 97)*(Lc1); // Lfc
			L1 += *(P1 + 118)*(Lq1); // Lfq
			L1 += *(P1 + 139)*(Le1); // Lfe
			L1 += *(P1 + 160)*(Lg1); // Lfg
			L1 += *(P1 + 181)*(Lh1); // Lfh
			L1 += *(P1 + 202)*(Li1); // Lfi
			L1 += *(P1 + 223)*(Ll1); // Lfl
			L1 += *(P1 + 244)*(Lk1); // Lfk
			L1 += *(P1 + 265)*(Lm1); // Lfm
			L1 += *(P1 + 286)*(Lf1); // Lff
			L1 += *(P1 + 307)*(Lp1); // Lfp
			L1 += *(P1 + 328)*(Lz1); // Lfs
			L1 += *(P1 + 349)*(Lt1); // Lft
			L1 += *(P1 + 370)*(Lw1); // Lfw
			L1 += *(P1 + 391)*(Ly1); // Lfy
			L1 += *(P1 + 412)*(Lv1); // Lfv
			L1 += *(P1 + 433)*(Lb1); // Lfb
			L2 = *(P2 + 13)*(La2); // Lfa
			L2 += *(P2 + 34)*(Lr2); // Lfr
			L2 += *(P2 + 55)*(Ln2); // Lfn
			L2 += *(P2 + 76)*(Ld2); // Lfd
			L2 += *(P2 + 97)*(Lc2); // Lfc
			L2 += *(P2 + 118)*(Lq2); // Lfq
			L2 += *(P2 + 139)*(Le2); // Lfe
			L2 += *(P2 + 160)*(Lg2); // Lfg
			L2 += *(P2 + 181)*(Lh2); // Lfh
			L2 += *(P2 + 202)*(Li2); // Lfi
			L2 += *(P2 + 223)*(Ll2); // Lfl
			L2 += *(P2 + 244)*(Lk2); // Lfk
			L2 += *(P2 + 265)*(Lm2); // Lfm
			L2 += *(P2 + 286)*(Lf2); // Lff
			L2 += *(P2 + 307)*(Lp2); // Lfp
			L2 += *(P2 + 328)*(Lz2); // Lfs
			L2 += *(P2 + 349)*(Lt2); // Lft
			L2 += *(P2 + 370)*(Lw2); // Lfw
			L2 += *(P2 + 391)*(Ly2); // Lfy
			L2 += *(P2 + 412)*(Lv2); // Lfv
			L2 += *(P2 + 433)*(Lb2); // Lfb
			*(Ls3 + 13) = L1*L2;
			
			// L(P)
			L1 = *(P1 + 14)*(La1); // Lpa
			L1 += *(P1 + 35)*(Lr1); // Lpr
			L1 += *(P1 + 56)*(Ln1); // Lpn
			L1 += *(P1 + 77)*(Ld1); // Lpd
			L1 += *(P1 + 98)*(Lc1); // Lpc
			L1 += *(P1 + 119)*(Lq1); // Lpq
			L1 += *(P1 + 140)*(Le1); // Lpe
			L1 += *(P1 + 161)*(Lg1); // Lpg
			L1 += *(P1 + 182)*(Lh1); // Lph
			L1 += *(P1 + 203)*(Li1); // Lpi
			L1 += *(P1 + 224)*(Ll1); // Lpl
			L1 += *(P1 + 245)*(Lk1); // Lpk
			L1 += *(P1 + 266)*(Lm1); // Lpm
			L1 += *(P1 + 287)*(Lf1); // Lpf
			L1 += *(P1 + 308)*(Lp1); // Lpp
			L1 += *(P1 + 329)*(Lz1); // Lps
			L1 += *(P1 + 350)*(Lt1); // Lpt
			L1 += *(P1 + 371)*(Lw1); // Lpw
			L1 += *(P1 + 392)*(Ly1); // Lpy
			L1 += *(P1 + 413)*(Lv1); // Lpv
			L1 += *(P1 + 434)*(Lb1); // Lpb
			L2 = *(P2 + 14)*(La2); // Lpa
			L2 += *(P2 + 35)*(Lr2); // Lpr
			L2 += *(P2 + 56)*(Ln2); // Lpn
			L2 += *(P2 + 77)*(Ld2); // Lpd
			L2 += *(P2 + 98)*(Lc2); // Lpc
			L2 += *(P2 + 119)*(Lq2); // Lpq
			L2 += *(P2 + 140)*(Le2); // Lpe
			L2 += *(P2 + 161)*(Lg2); // Lpg
			L2 += *(P2 + 182)*(Lh2); // Lph
			L2 += *(P2 + 203)*(Li2); // Lpi
			L2 += *(P2 + 224)*(Ll2); // Lpl
			L2 += *(P2 + 245)*(Lk2); // Lpk
			L2 += *(P2 + 266)*(Lm2); // Lpm
			L2 += *(P2 + 287)*(Lf2); // Lpf
			L2 += *(P2 + 308)*(Lp2); // Lpp
			L2 += *(P2 + 329)*(Lz2); // Lps
			L2 += *(P2 + 350)*(Lt2); // Lpt
			L2 += *(P2 + 371)*(Lw2); // Lpw
			L2 += *(P2 + 392)*(Ly2); // Lpy
			L2 += *(P2 + 413)*(Lv2); // Lpv
			L2 += *(P2 + 434)*(Lb2); // Lpb
			*(Ls3 + 14) = L1*L2;
			
			// L(S)
			L1 = *(P1 + 15)*(La1); // Lsa
			L1 += *(P1 + 36)*(Lr1); // Lsr
			L1 += *(P1 + 57)*(Ln1); // Lsn
			L1 += *(P1 + 78)*(Ld1); // Lsd
			L1 += *(P1 + 99)*(Lc1); // Lsc
			L1 += *(P1 + 120)*(Lq1); // Lsq
			L1 += *(P1 + 141)*(Le1); // Lse
			L1 += *(P1 + 162)*(Lg1); // Lsg
			L1 += *(P1 + 183)*(Lh1); // Lsh
			L1 += *(P1 + 204)*(Li1); // Lsi
			L1 += *(P1 + 225)*(Ll1); // Lsl
			L1 += *(P1 + 246)*(Lk1); // Lsk
			L1 += *(P1 + 267)*(Lm1); // Lsm
			L1 += *(P1 + 288)*(Lf1); // Lsf
			L1 += *(P1 + 309)*(Lp1); // Lsp
			L1 += *(P1 + 330)*(Lz1); // Lss
			L1 += *(P1 + 351)*(Lt1); // Lst
			L1 += *(P1 + 372)*(Lw1); // Lsw
			L1 += *(P1 + 393)*(Ly1); // Lsy
			L1 += *(P1 + 414)*(Lv1); // Lsv
			L1 += *(P1 + 435)*(Lb1); // Lsb
			L2 = *(P2 + 15)*(La2); // Lsa
			L2 += *(P2 + 36)*(Lr2); // Lsr
			L2 += *(P2 + 57)*(Ln2); // Lsn
			L2 += *(P2 + 78)*(Ld2); // Lsd
			L2 += *(P2 + 99)*(Lc2); // Lsc
			L2 += *(P2 + 120)*(Lq2); // Lsq
			L2 += *(P2 + 141)*(Le2); // Lse
			L2 += *(P2 + 162)*(Lg2); // Lsg
			L2 += *(P2 + 183)*(Lh2); // Lsh
			L2 += *(P2 + 204)*(Li2); // Lsi
			L2 += *(P2 + 225)*(Ll2); // Lsl
			L2 += *(P2 + 246)*(Lk2); // Lsk
			L2 += *(P2 + 267)*(Lm2); // Lsm
			L2 += *(P2 + 288)*(Lf2); // Lsf
			L2 += *(P2 + 309)*(Lp2); // Lsp
			L2 += *(P2 + 330)*(Lz2); // Lss
			L2 += *(P2 + 351)*(Lt2); // Lst
			L2 += *(P2 + 372)*(Lw2); // Lsw
			L2 += *(P2 + 393)*(Ly2); // Lsy
			L2 += *(P2 + 414)*(Lv2); // Lsv
			L2 += *(P2 + 435)*(Lb2); // Lsb
			*(Ls3 + 15) = L1*L2;
			
			// L(T)
			L1 = *(P1 + 16)*(La1); // Lta
			L1 += *(P1 + 37)*(Lr1); // Ltr
			L1 += *(P1 + 58)*(Ln1); // Ltn
			L1 += *(P1 + 79)*(Ld1); // Ltd
			L1 += *(P1 + 100)*(Lc1); // Ltc
			L1 += *(P1 + 121)*(Lq1); // Ltq
			L1 += *(P1 + 142)*(Le1); // Lte
			L1 += *(P1 + 163)*(Lg1); // Ltg
			L1 += *(P1 + 184)*(Lh1); // Lth
			L1 += *(P1 + 205)*(Li1); // Lti
			L1 += *(P1 + 226)*(Ll1); // Ltl
			L1 += *(P1 + 247)*(Lk1); // Ltk
			L1 += *(P1 + 268)*(Lm1); // Ltm
			L1 += *(P1 + 289)*(Lf1); // Ltf
			L1 += *(P1 + 310)*(Lp1); // Ltp
			L1 += *(P1 + 331)*(Lz1); // Lts
			L1 += *(P1 + 352)*(Lt1); // Ltt
			L1 += *(P1 + 373)*(Lw1); // Ltw
			L1 += *(P1 + 394)*(Ly1); // Lty
			L1 += *(P1 + 415)*(Lv1); // Ltv
			L1 += *(P1 + 436)*(Lb1); // Ltb
			L2 = *(P2 + 16)*(La2); // Lta
			L2 += *(P2 + 37)*(Lr2); // Ltr
			L2 += *(P2 + 58)*(Ln2); // Ltn
			L2 += *(P2 + 79)*(Ld2); // Ltd
			L2 += *(P2 + 100)*(Lc2); // Ltc
			L2 += *(P2 + 121)*(Lq2); // Ltq
			L2 += *(P2 + 142)*(Le2); // Lte
			L2 += *(P2 + 163)*(Lg2); // Ltg
			L2 += *(P2 + 184)*(Lh2); // Lth
			L2 += *(P2 + 205)*(Li2); // Lti
			L2 += *(P2 + 226)*(Ll2); // Ltl
			L2 += *(P2 + 247)*(Lk2); // Ltk
			L2 += *(P2 + 268)*(Lm2); // Ltm
			L2 += *(P2 + 289)*(Lf2); // Ltf
			L2 += *(P2 + 310)*(Lp2); // Ltp
			L2 += *(P2 + 331)*(Lz2); // Lts
			L2 += *(P2 + 352)*(Lt2); // Ltt
			L2 += *(P2 + 373)*(Lw2); // Ltw
			L2 += *(P2 + 394)*(Ly2); // Lty
			L2 += *(P2 + 415)*(Lv2); // Ltv
			L2 += *(P2 + 436)*(Lb2); // Ltb
			*(Ls3 + 16) = L1*L2;
			
			// L(W)
			L1 = *(P1 + 17)*(La1); // Lwa
			L1 += *(P1 + 38)*(Lr1); // Lwr
			L1 += *(P1 + 59)*(Ln1); // Lwn
			L1 += *(P1 + 80)*(Ld1); // Lwd
			L1 += *(P1 + 101)*(Lc1); // Lwc
			L1 += *(P1 + 122)*(Lq1); // Lwq
			L1 += *(P1 + 143)*(Le1); // Lwe
			L1 += *(P1 + 164)*(Lg1); // Lwg
			L1 += *(P1 + 185)*(Lh1); // Lwh
			L1 += *(P1 + 206)*(Li1); // Lwi
			L1 += *(P1 + 227)*(Ll1); // Lwl
			L1 += *(P1 + 248)*(Lk1); // Lwk
			L1 += *(P1 + 269)*(Lm1); // Lwm
			L1 += *(P1 + 290)*(Lf1); // Lwf
			L1 += *(P1 + 311)*(Lp1); // Lwp
			L1 += *(P1 + 332)*(Lz1); // Lws
			L1 += *(P1 + 353)*(Lt1); // Lwt
			L1 += *(P1 + 374)*(Lw1); // Lww
			L1 += *(P1 + 395)*(Ly1); // Lwy
			L1 += *(P1 + 416)*(Lv1); // Lwv
			L1 += *(P1 + 437)*(Lb1); // Lwb
			L2 = *(P2 + 17)*(La2); // Lwa
			L2 += *(P2 + 38)*(Lr2); // Lwr
			L2 += *(P2 + 59)*(Ln2); // Lwn
			L2 += *(P2 + 80)*(Ld2); // Lwd
			L2 += *(P2 + 101)*(Lc2); // Lwc
			L2 += *(P2 + 122)*(Lq2); // Lwq
			L2 += *(P2 + 143)*(Le2); // Lwe
			L2 += *(P2 + 164)*(Lg2); // Lwg
			L2 += *(P2 + 185)*(Lh2); // Lwh
			L2 += *(P2 + 206)*(Li2); // Lwi
			L2 += *(P2 + 227)*(Ll2); // Lwl
			L2 += *(P2 + 248)*(Lk2); // Lwk
			L2 += *(P2 + 269)*(Lm2); // Lwm
			L2 += *(P2 + 290)*(Lf2); // Lwf
			L2 += *(P2 + 311)*(Lp2); // Lwp
			L2 += *(P2 + 332)*(Lz2); // Lws
			L2 += *(P2 + 353)*(Lt2); // Lwt
			L2 += *(P2 + 374)*(Lw2); // Lww
			L2 += *(P2 + 395)*(Ly2); // Lwy
			L2 += *(P2 + 416)*(Lv2); // Lwv
			L2 += *(P2 + 437)*(Lb2); // Lwb
			*(Ls3 + 17) = L1*L2;
			
			// L(Y)
			L1 = *(P1 + 18)*(La1); // Lya
			L1 += *(P1 + 39)*(Lr1); // Lyr
			L1 += *(P1 + 60)*(Ln1); // Lyn
			L1 += *(P1 + 81)*(Ld1); // Lyd
			L1 += *(P1 + 102)*(Lc1); // Lyc
			L1 += *(P1 + 123)*(Lq1); // Lyq
			L1 += *(P1 + 144)*(Le1); // Lye
			L1 += *(P1 + 165)*(Lg1); // Lyg
			L1 += *(P1 + 186)*(Lh1); // Lyh
			L1 += *(P1 + 207)*(Li1); // Lyi
			L1 += *(P1 + 228)*(Ll1); // Lyl
			L1 += *(P1 + 249)*(Lk1); // Lyk
			L1 += *(P1 + 270)*(Lm1); // Lym
			L1 += *(P1 + 291)*(Lf1); // Lyf
			L1 += *(P1 + 312)*(Lp1); // Lyp
			L1 += *(P1 + 333)*(Lz1); // Lys
			L1 += *(P1 + 354)*(Lt1); // Lyt
			L1 += *(P1 + 375)*(Lw1); // Lyw
			L1 += *(P1 + 396)*(Ly1); // Lyy
			L1 += *(P1 + 417)*(Lv1); // Lyv
			L1 += *(P1 + 438)*(Lb1); // Lyb
			L2 = *(P2 + 18)*(La2); // Lya
			L2 += *(P2 + 39)*(Lr2); // Lyr
			L2 += *(P2 + 60)*(Ln2); // Lyn
			L2 += *(P2 + 81)*(Ld2); // Lyd
			L2 += *(P2 + 102)*(Lc2); // Lyc
			L2 += *(P2 + 123)*(Lq2); // Lyq
			L2 += *(P2 + 144)*(Le2); // Lye
			L2 += *(P2 + 165)*(Lg2); // Lyg
			L2 += *(P2 + 186)*(Lh2); // Lyh
			L2 += *(P2 + 207)*(Li2); // Lyi
			L2 += *(P2 + 228)*(Ll2); // Lyl
			L2 += *(P2 + 249)*(Lk2); // Lyk
			L2 += *(P2 + 270)*(Lm2); // Lym
			L2 += *(P2 + 291)*(Lf2); // Lyf
			L2 += *(P2 + 312)*(Lp2); // Lyp
			L2 += *(P2 + 333)*(Lz2); // Lys
			L2 += *(P2 + 354)*(Lt2); // Lyt
			L2 += *(P2 + 375)*(Lw2); // Lyw
			L2 += *(P2 + 396)*(Ly2); // Lyy
			L2 += *(P2 + 417)*(Lv2); // Lyv
			L2 += *(P2 + 438)*(Lb2); // Lyb
			*(Ls3 + 18) = L1*L2;
			
			// L(V)
			L1 = *(P1 + 19)*(La1); // Lva
			L1 += *(P1 + 40)*(Lr1); // Lvr
			L1 += *(P1 + 61)*(Ln1); // Lvn
			L1 += *(P1 + 82)*(Ld1); // Lvd
			L1 += *(P1 + 103)*(Lc1); // Lvc
			L1 += *(P1 + 124)*(Lq1); // Lvq
			L1 += *(P1 + 145)*(Le1); // Lve
			L1 += *(P1 + 166)*(Lg1); // Lvg
			L1 += *(P1 + 187)*(Lh1); // Lvh
			L1 += *(P1 + 208)*(Li1); // Lvi
			L1 += *(P1 + 229)*(Ll1); // Lvl
			L1 += *(P1 + 250)*(Lk1); // Lvk
			L1 += *(P1 + 271)*(Lm1); // Lvm
			L1 += *(P1 + 292)*(Lf1); // Lvf
			L1 += *(P1 + 313)*(Lp1); // Lvp
			L1 += *(P1 + 334)*(Lz1); // Lvs
			L1 += *(P1 + 355)*(Lt1); // Lvt
			L1 += *(P1 + 376)*(Lw1); // Lvw
			L1 += *(P1 + 397)*(Ly1); // Lvy
			L1 += *(P1 + 418)*(Lv1); // Lvv
			L1 += *(P1 + 439)*(Lb1); // Lvb
			L2 = *(P2 + 19)*(La2); // Lva
			L2 += *(P2 + 40)*(Lr2); // Lvr
			L2 += *(P2 + 61)*(Ln2); // Lvn
			L2 += *(P2 + 82)*(Ld2); // Lvd
			L2 += *(P2 + 103)*(Lc2); // Lvc
			L2 += *(P2 + 124)*(Lq2); // Lvq
			L2 += *(P2 + 145)*(Le2); // Lve
			L2 += *(P2 + 166)*(Lg2); // Lvg
			L2 += *(P2 + 187)*(Lh2); // Lvh
			L2 += *(P2 + 208)*(Li2); // Lvi
			L2 += *(P2 + 229)*(Ll2); // Lvl
			L2 += *(P2 + 250)*(Lk2); // Lvk
			L2 += *(P2 + 271)*(Lm2); // Lvm
			L2 += *(P2 + 292)*(Lf2); // Lvf
			L2 += *(P2 + 313)*(Lp2); // Lvp
			L2 += *(P2 + 334)*(Lz2); // Lvs
			L2 += *(P2 + 355)*(Lt2); // Lvt
			L2 += *(P2 + 376)*(Lw2); // Lvw
			L2 += *(P2 + 397)*(Ly2); // Lvy
			L2 += *(P2 + 418)*(Lv2); // Lvv
			L2 += *(P2 + 439)*(Lb2); // Lvb
			*(Ls3 + 19) = L1*L2;
			
			// L(Indels)
			L1 = *(P1 + 20)*(La1); // Lba
			L1 += *(P1 + 41)*(Lr1); // Lbr
			L1 += *(P1 + 62)*(Ln1); // Lbn
			L1 += *(P1 + 83)*(Ld1); // Lbd
			L1 += *(P1 + 104)*(Lc1); // Lbc
			L1 += *(P1 + 125)*(Lq1); // Lbq
			L1 += *(P1 + 146)*(Le1); // Lbe
			L1 += *(P1 + 167)*(Lg1); // Lbg
			L1 += *(P1 + 188)*(Lh1); // Lbh
			L1 += *(P1 + 209)*(Li1); // Lbi
			L1 += *(P1 + 230)*(Ll1); // Lbl
			L1 += *(P1 + 251)*(Lk1); // Lbk
			L1 += *(P1 + 272)*(Lm1); // Lbm
			L1 += *(P1 + 293)*(Lf1); // Lbf
			L1 += *(P1 + 314)*(Lp1); // Lbp
			L1 += *(P1 + 335)*(Lz1); // Lbs
			L1 += *(P1 + 356)*(Lt1); // Lbt
			L1 += *(P1 + 377)*(Lw1); // Lbw
			L1 += *(P1 + 398)*(Ly1); // Lby
			L1 += *(P1 + 419)*(Lv1); // Lbv
			L1 += *(P1 + 440)*(Lb1); // Lbb
			L2 = *(P2 + 20)*(La2); // Lba
			L2 += *(P2 + 41)*(Lr2); // Lbr
			L2 += *(P2 + 62)*(Ln2); // Lbn
			L2 += *(P2 + 83)*(Ld2); // Lbd
			L2 += *(P2 + 104)*(Lc2); // Lbc
			L2 += *(P2 + 125)*(Lq2); // Lbq
			L2 += *(P2 + 146)*(Le2); // Lbe
			L2 += *(P2 + 167)*(Lg2); // Lbg
			L2 += *(P2 + 188)*(Lh2); // Lbh
			L2 += *(P2 + 209)*(Li2); // Lbi
			L2 += *(P2 + 230)*(Ll2); // Lbl
			L2 += *(P2 + 251)*(Lk2); // Lbk
			L2 += *(P2 + 272)*(Lm2); // Lbm
			L2 += *(P2 + 293)*(Lf2); // Lbf
			L2 += *(P2 + 314)*(Lp2); // Lbp
			L2 += *(P2 + 335)*(Lz2); // Lbs
			L2 += *(P2 + 356)*(Lt2); // Lbt
			L2 += *(P2 + 377)*(Lw2); // Lbw
			L2 += *(P2 + 398)*(Ly2); // Lby
			L2 += *(P2 + 419)*(Lv2); // Lbv
			L2 += *(P2 + 440)*(Lb2); // Lbb
			*(Ls3 + 20) = L1*L2;
			
			*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon) ||
				(*(Ls3 + 5) > 0 && *(Ls3 + 5) < inv_epsilon) ||
				(*(Ls3 + 6) > 0 && *(Ls3 + 6) < inv_epsilon) ||
				(*(Ls3 + 7) > 0 && *(Ls3 + 7) < inv_epsilon) ||
				(*(Ls3 + 8) > 0 && *(Ls3 + 8) < inv_epsilon) ||
				(*(Ls3 + 9) > 0 && *(Ls3 + 9) < inv_epsilon) ||
				(*(Ls3 + 10) > 0 && *(Ls3 + 10) < inv_epsilon) ||
				(*(Ls3 + 11) > 0 && *(Ls3 + 11) < inv_epsilon) ||
				(*(Ls3 + 12) > 0 && *(Ls3 + 12) < inv_epsilon) ||
				(*(Ls3 + 13) > 0 && *(Ls3 + 13) < inv_epsilon) ||
				(*(Ls3 + 14) > 0 && *(Ls3 + 14) < inv_epsilon) ||
				(*(Ls3 + 15) > 0 && *(Ls3 + 15) < inv_epsilon) ||
				(*(Ls3 + 16) > 0 && *(Ls3 + 16) < inv_epsilon) ||
				(*(Ls3 + 17) > 0 && *(Ls3 + 17) < inv_epsilon) ||
				(*(Ls3 + 18) > 0 && *(Ls3 + 18) < inv_epsilon) ||
				(*(Ls3 + 19) > 0 && *(Ls3 + 19) < inv_epsilon) ||
				(*(Ls3 + 20) > 0 && *(Ls3 + 20) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) *= epsilon;
				*(Ls3 + 6) *= epsilon;
				*(Ls3 + 7) *= epsilon;
				*(Ls3 + 8) *= epsilon;
				*(Ls3 + 9) *= epsilon;
				*(Ls3 + 10) *= epsilon;
				*(Ls3 + 11) *= epsilon;
				*(Ls3 + 12) *= epsilon;
				*(Ls3 + 13) *= epsilon;
				*(Ls3 + 14) *= epsilon;
				*(Ls3 + 15) *= epsilon;
				*(Ls3 + 16) *= epsilon;
				*(Ls3 + 17) *= epsilon;
				*(Ls3 + 18) *= epsilon;
				*(Ls3 + 19) *= epsilon;
				*(Ls3 + 20) *= epsilon;
				*(Ls3 + 21) += 1;
			}
		} else {
			// second branch can be disregarded
			
			// L(A)
			L1 = *(P1 + 0)*(La1); // Laa
			L1 += *(P1 + 21)*(Lr1); // Lar
			L1 += *(P1 + 42)*(Ln1); // Lan
			L1 += *(P1 + 63)*(Ld1); // Lad
			L1 += *(P1 + 84)*(Lc1); // Lac
			L1 += *(P1 + 105)*(Lq1); // Laq
			L1 += *(P1 + 126)*(Le1); // Lae
			L1 += *(P1 + 147)*(Lg1); // Lag
			L1 += *(P1 + 168)*(Lh1); // Lah
			L1 += *(P1 + 189)*(Li1); // Lai
			L1 += *(P1 + 210)*(Ll1); // Lal
			L1 += *(P1 + 231)*(Lk1); // Lak
			L1 += *(P1 + 252)*(Lm1); // Lam
			L1 += *(P1 + 273)*(Lf1); // Laf
			L1 += *(P1 + 294)*(Lp1); // Lap
			L1 += *(P1 + 315)*(Lz1); // Las
			L1 += *(P1 + 336)*(Lt1); // Lat
			L1 += *(P1 + 357)*(Lw1); // Law
			L1 += *(P1 + 378)*(Ly1); // Lay
			L1 += *(P1 + 399)*(Lv1); // Lav
			L1 += *(P1 + 420)*(Lb1); // Lab
			*(Ls3 + 0) = L1;
			
			// L(R)
			L1 = *(P1 + 1)*(La1); // Lra
			L1 += *(P1 + 22)*(Lr1); // Lrr
			L1 += *(P1 + 43)*(Ln1); // Lrn
			L1 += *(P1 + 64)*(Ld1); // Lrd
			L1 += *(P1 + 85)*(Lc1); // Lrc
			L1 += *(P1 + 106)*(Lq1); // Lrq
			L1 += *(P1 + 127)*(Le1); // Lre
			L1 += *(P1 + 148)*(Lg1); // Lrg
			L1 += *(P1 + 169)*(Lh1); // Lrh
			L1 += *(P1 + 190)*(Li1); // Lri
			L1 += *(P1 + 211)*(Ll1); // Lrl
			L1 += *(P1 + 232)*(Lk1); // Lrk
			L1 += *(P1 + 253)*(Lm1); // Lrm
			L1 += *(P1 + 274)*(Lf1); // Lrf
			L1 += *(P1 + 295)*(Lp1); // Lrp
			L1 += *(P1 + 316)*(Lz1); // Lrs
			L1 += *(P1 + 337)*(Lt1); // Lrt
			L1 += *(P1 + 358)*(Lw1); // Lrw
			L1 += *(P1 + 379)*(Ly1); // Lry
			L1 += *(P1 + 400)*(Lv1); // Lrv
			L1 += *(P1 + 421)*(Lb1); // Lrb
			*(Ls3 + 1) = L1;
			
			// L(N)
			L1 = *(P1 + 2)*(La1); // Lna
			L1 += *(P1 + 23)*(Lr1); // Lnr
			L1 += *(P1 + 44)*(Ln1); // Lnn
			L1 += *(P1 + 65)*(Ld1); // Lnd
			L1 += *(P1 + 86)*(Lc1); // Lnc
			L1 += *(P1 + 107)*(Lq1); // Lnq
			L1 += *(P1 + 128)*(Le1); // Lne
			L1 += *(P1 + 149)*(Lg1); // Lng
			L1 += *(P1 + 170)*(Lh1); // Lnh
			L1 += *(P1 + 191)*(Li1); // Lni
			L1 += *(P1 + 212)*(Ll1); // Lnl
			L1 += *(P1 + 233)*(Lk1); // Lnk
			L1 += *(P1 + 254)*(Lm1); // Lnm
			L1 += *(P1 + 275)*(Lf1); // Lnf
			L1 += *(P1 + 296)*(Lp1); // Lnp
			L1 += *(P1 + 317)*(Lz1); // Lns
			L1 += *(P1 + 338)*(Lt1); // Lnt
			L1 += *(P1 + 359)*(Lw1); // Lnw
			L1 += *(P1 + 380)*(Ly1); // Lny
			L1 += *(P1 + 401)*(Lv1); // Lnv
			L1 += *(P1 + 422)*(Lb1); // Lnb
			*(Ls3 + 2) = L1;
			
			// L(D)
			L1 = *(P1 + 3)*(La1); // Lda
			L1 += *(P1 + 24)*(Lr1); // Ldr
			L1 += *(P1 + 45)*(Ln1); // Ldn
			L1 += *(P1 + 66)*(Ld1); // Ldd
			L1 += *(P1 + 87)*(Lc1); // Ldc
			L1 += *(P1 + 108)*(Lq1); // Ldq
			L1 += *(P1 + 129)*(Le1); // Lde
			L1 += *(P1 + 150)*(Lg1); // Ldg
			L1 += *(P1 + 171)*(Lh1); // Ldh
			L1 += *(P1 + 192)*(Li1); // Ldi
			L1 += *(P1 + 213)*(Ll1); // Ldl
			L1 += *(P1 + 234)*(Lk1); // Ldk
			L1 += *(P1 + 255)*(Lm1); // Ldm
			L1 += *(P1 + 276)*(Lf1); // Ldf
			L1 += *(P1 + 297)*(Lp1); // Ldp
			L1 += *(P1 + 318)*(Lz1); // Lds
			L1 += *(P1 + 339)*(Lt1); // Ldt
			L1 += *(P1 + 360)*(Lw1); // Ldw
			L1 += *(P1 + 381)*(Ly1); // Ldy
			L1 += *(P1 + 402)*(Lv1); // Ldv
			L1 += *(P1 + 423)*(Lb1); // Ldb
			*(Ls3 + 3) = L1;
			
			// L(C)
			L1 = *(P1 + 4)*(La1); // Lca
			L1 += *(P1 + 25)*(Lr1); // Lcr
			L1 += *(P1 + 46)*(Ln1); // Lcn
			L1 += *(P1 + 67)*(Ld1); // Lcd
			L1 += *(P1 + 88)*(Lc1); // Lcc
			L1 += *(P1 + 109)*(Lq1); // Lcq
			L1 += *(P1 + 130)*(Le1); // Lce
			L1 += *(P1 + 151)*(Lg1); // Lcg
			L1 += *(P1 + 172)*(Lh1); // Lch
			L1 += *(P1 + 193)*(Li1); // Lci
			L1 += *(P1 + 214)*(Ll1); // Lcl
			L1 += *(P1 + 235)*(Lk1); // Lck
			L1 += *(P1 + 256)*(Lm1); // Lcm
			L1 += *(P1 + 277)*(Lf1); // Lcf
			L1 += *(P1 + 298)*(Lp1); // Lcp
			L1 += *(P1 + 319)*(Lz1); // Lcs
			L1 += *(P1 + 340)*(Lt1); // Lct
			L1 += *(P1 + 361)*(Lw1); // Lcw
			L1 += *(P1 + 382)*(Ly1); // Lcy
			L1 += *(P1 + 403)*(Lv1); // Lcv
			L1 += *(P1 + 424)*(Lb1); // Lcb
			*(Ls3 + 4) = L1;
			
			// L(Q)
			L1 = *(P1 + 5)*(La1); // Lqa
			L1 += *(P1 + 26)*(Lr1); // Lqr
			L1 += *(P1 + 47)*(Ln1); // Lqn
			L1 += *(P1 + 68)*(Ld1); // Lqd
			L1 += *(P1 + 89)*(Lc1); // Lqc
			L1 += *(P1 + 110)*(Lq1); // Lqq
			L1 += *(P1 + 131)*(Le1); // Lqe
			L1 += *(P1 + 152)*(Lg1); // Lqg
			L1 += *(P1 + 173)*(Lh1); // Lqh
			L1 += *(P1 + 194)*(Li1); // Lqi
			L1 += *(P1 + 215)*(Ll1); // Lql
			L1 += *(P1 + 236)*(Lk1); // Lqk
			L1 += *(P1 + 257)*(Lm1); // Lqm
			L1 += *(P1 + 278)*(Lf1); // Lqf
			L1 += *(P1 + 299)*(Lp1); // Lqp
			L1 += *(P1 + 320)*(Lz1); // Lqs
			L1 += *(P1 + 341)*(Lt1); // Lqt
			L1 += *(P1 + 362)*(Lw1); // Lqw
			L1 += *(P1 + 383)*(Ly1); // Lqy
			L1 += *(P1 + 404)*(Lv1); // Lqv
			L1 += *(P1 + 425)*(Lb1); // Lqb
			*(Ls3 + 5) = L1;
			
			// L(E)
			L1 = *(P1 + 6)*(La1); // Lea
			L1 += *(P1 + 27)*(Lr1); // Ler
			L1 += *(P1 + 48)*(Ln1); // Len
			L1 += *(P1 + 69)*(Ld1); // Led
			L1 += *(P1 + 90)*(Lc1); // Lec
			L1 += *(P1 + 111)*(Lq1); // Leq
			L1 += *(P1 + 132)*(Le1); // Lee
			L1 += *(P1 + 153)*(Lg1); // Leg
			L1 += *(P1 + 174)*(Lh1); // Leh
			L1 += *(P1 + 195)*(Li1); // Lei
			L1 += *(P1 + 216)*(Ll1); // Lel
			L1 += *(P1 + 237)*(Lk1); // Lek
			L1 += *(P1 + 258)*(Lm1); // Lem
			L1 += *(P1 + 279)*(Lf1); // Lef
			L1 += *(P1 + 300)*(Lp1); // Lep
			L1 += *(P1 + 321)*(Lz1); // Les
			L1 += *(P1 + 342)*(Lt1); // Let
			L1 += *(P1 + 363)*(Lw1); // Lew
			L1 += *(P1 + 384)*(Ly1); // Ley
			L1 += *(P1 + 405)*(Lv1); // Lev
			L1 += *(P1 + 426)*(Lb1); // Leb
			*(Ls3 + 6) = L1;
			
			// L(G)
			L1 = *(P1 + 7)*(La1); // Lga
			L1 += *(P1 + 28)*(Lr1); // Lgr
			L1 += *(P1 + 49)*(Ln1); // Lgn
			L1 += *(P1 + 70)*(Ld1); // Lgd
			L1 += *(P1 + 91)*(Lc1); // Lgc
			L1 += *(P1 + 112)*(Lq1); // Lgq
			L1 += *(P1 + 133)*(Le1); // Lge
			L1 += *(P1 + 154)*(Lg1); // Lgg
			L1 += *(P1 + 175)*(Lh1); // Lgh
			L1 += *(P1 + 196)*(Li1); // Lgi
			L1 += *(P1 + 217)*(Ll1); // Lgl
			L1 += *(P1 + 238)*(Lk1); // Lgk
			L1 += *(P1 + 259)*(Lm1); // Lgm
			L1 += *(P1 + 280)*(Lf1); // Lgf
			L1 += *(P1 + 301)*(Lp1); // Lgp
			L1 += *(P1 + 322)*(Lz1); // Lgs
			L1 += *(P1 + 343)*(Lt1); // Lgt
			L1 += *(P1 + 364)*(Lw1); // Lgw
			L1 += *(P1 + 385)*(Ly1); // Lgy
			L1 += *(P1 + 406)*(Lv1); // Lgv
			L1 += *(P1 + 427)*(Lb1); // Lgb
			*(Ls3 + 7) = L1;
			
			// L(H)
			L1 = *(P1 + 8)*(La1); // Lha
			L1 += *(P1 + 29)*(Lr1); // Lhr
			L1 += *(P1 + 50)*(Ln1); // Lhn
			L1 += *(P1 + 71)*(Ld1); // Lhd
			L1 += *(P1 + 92)*(Lc1); // Lhc
			L1 += *(P1 + 113)*(Lq1); // Lhq
			L1 += *(P1 + 134)*(Le1); // Lhe
			L1 += *(P1 + 155)*(Lg1); // Lhg
			L1 += *(P1 + 176)*(Lh1); // Lhh
			L1 += *(P1 + 197)*(Li1); // Lhi
			L1 += *(P1 + 218)*(Ll1); // Lhl
			L1 += *(P1 + 239)*(Lk1); // Lhk
			L1 += *(P1 + 260)*(Lm1); // Lhm
			L1 += *(P1 + 281)*(Lf1); // Lhf
			L1 += *(P1 + 302)*(Lp1); // Lhp
			L1 += *(P1 + 323)*(Lz1); // Lhs
			L1 += *(P1 + 344)*(Lt1); // Lht
			L1 += *(P1 + 365)*(Lw1); // Lhw
			L1 += *(P1 + 386)*(Ly1); // Lhy
			L1 += *(P1 + 407)*(Lv1); // Lhv
			L1 += *(P1 + 428)*(Lb1); // Lhb
			*(Ls3 + 8) = L1;
			
			// L(I)
			L1 = *(P1 + 9)*(La1); // Lia
			L1 += *(P1 + 30)*(Lr1); // Lir
			L1 += *(P1 + 51)*(Ln1); // Lin
			L1 += *(P1 + 72)*(Ld1); // Lid
			L1 += *(P1 + 93)*(Lc1); // Lic
			L1 += *(P1 + 114)*(Lq1); // Liq
			L1 += *(P1 + 135)*(Le1); // Lie
			L1 += *(P1 + 156)*(Lg1); // Lig
			L1 += *(P1 + 177)*(Lh1); // Lih
			L1 += *(P1 + 198)*(Li1); // Lii
			L1 += *(P1 + 219)*(Ll1); // Lil
			L1 += *(P1 + 240)*(Lk1); // Lik
			L1 += *(P1 + 261)*(Lm1); // Lim
			L1 += *(P1 + 282)*(Lf1); // Lif
			L1 += *(P1 + 303)*(Lp1); // Lip
			L1 += *(P1 + 324)*(Lz1); // Lis
			L1 += *(P1 + 345)*(Lt1); // Lit
			L1 += *(P1 + 366)*(Lw1); // Liw
			L1 += *(P1 + 387)*(Ly1); // Liy
			L1 += *(P1 + 408)*(Lv1); // Liv
			L1 += *(P1 + 429)*(Lb1); // Lib
			*(Ls3 + 9) = L1;
			
			// L(L)
			L1 = *(P1 + 10)*(La1); // Lla
			L1 += *(P1 + 31)*(Lr1); // Llr
			L1 += *(P1 + 52)*(Ln1); // Lln
			L1 += *(P1 + 73)*(Ld1); // Lld
			L1 += *(P1 + 94)*(Lc1); // Llc
			L1 += *(P1 + 115)*(Lq1); // Llq
			L1 += *(P1 + 136)*(Le1); // Lle
			L1 += *(P1 + 157)*(Lg1); // Llg
			L1 += *(P1 + 178)*(Lh1); // Llh
			L1 += *(P1 + 199)*(Li1); // Lli
			L1 += *(P1 + 220)*(Ll1); // Lll
			L1 += *(P1 + 241)*(Lk1); // Llk
			L1 += *(P1 + 262)*(Lm1); // Llm
			L1 += *(P1 + 283)*(Lf1); // Llf
			L1 += *(P1 + 304)*(Lp1); // Llp
			L1 += *(P1 + 325)*(Lz1); // Lls
			L1 += *(P1 + 346)*(Lt1); // Llt
			L1 += *(P1 + 367)*(Lw1); // Llw
			L1 += *(P1 + 388)*(Ly1); // Lly
			L1 += *(P1 + 409)*(Lv1); // Llv
			L1 += *(P1 + 430)*(Lb1); // Llb
			*(Ls3 + 10) = L1;
			
			// L(K)
			L1 = *(P1 + 11)*(La1); // Lka
			L1 += *(P1 + 32)*(Lr1); // Lkr
			L1 += *(P1 + 53)*(Ln1); // Lkn
			L1 += *(P1 + 74)*(Ld1); // Lkd
			L1 += *(P1 + 95)*(Lc1); // Lkc
			L1 += *(P1 + 116)*(Lq1); // Lkq
			L1 += *(P1 + 137)*(Le1); // Lke
			L1 += *(P1 + 158)*(Lg1); // Lkg
			L1 += *(P1 + 179)*(Lh1); // Lkh
			L1 += *(P1 + 200)*(Li1); // Lki
			L1 += *(P1 + 221)*(Ll1); // Lkl
			L1 += *(P1 + 242)*(Lk1); // Lkk
			L1 += *(P1 + 263)*(Lm1); // Lkm
			L1 += *(P1 + 284)*(Lf1); // Lkf
			L1 += *(P1 + 305)*(Lp1); // Lkp
			L1 += *(P1 + 326)*(Lz1); // Lks
			L1 += *(P1 + 347)*(Lt1); // Lkt
			L1 += *(P1 + 368)*(Lw1); // Lkw
			L1 += *(P1 + 389)*(Ly1); // Lky
			L1 += *(P1 + 410)*(Lv1); // Lkv
			L1 += *(P1 + 431)*(Lb1); // Lkb
			*(Ls3 + 11) = L1;
			
			// L(M)
			L1 = *(P1 + 12)*(La1); // Lma
			L1 += *(P1 + 33)*(Lr1); // Lmr
			L1 += *(P1 + 54)*(Ln1); // Lmn
			L1 += *(P1 + 75)*(Ld1); // Lmd
			L1 += *(P1 + 96)*(Lc1); // Lmc
			L1 += *(P1 + 117)*(Lq1); // Lmq
			L1 += *(P1 + 138)*(Le1); // Lme
			L1 += *(P1 + 159)*(Lg1); // Lmg
			L1 += *(P1 + 180)*(Lh1); // Lmh
			L1 += *(P1 + 201)*(Li1); // Lmi
			L1 += *(P1 + 222)*(Ll1); // Lml
			L1 += *(P1 + 243)*(Lk1); // Lmk
			L1 += *(P1 + 264)*(Lm1); // Lmm
			L1 += *(P1 + 285)*(Lf1); // Lmf
			L1 += *(P1 + 306)*(Lp1); // Lmp
			L1 += *(P1 + 327)*(Lz1); // Lms
			L1 += *(P1 + 348)*(Lt1); // Lmt
			L1 += *(P1 + 369)*(Lw1); // Lmw
			L1 += *(P1 + 390)*(Ly1); // Lmy
			L1 += *(P1 + 411)*(Lv1); // Lmv
			L1 += *(P1 + 432)*(Lb1); // Lmb
			*(Ls3 + 12) = L1;
			
			// L(F)
			L1 = *(P1 + 13)*(La1); // Lfa
			L1 += *(P1 + 34)*(Lr1); // Lfr
			L1 += *(P1 + 55)*(Ln1); // Lfn
			L1 += *(P1 + 76)*(Ld1); // Lfd
			L1 += *(P1 + 97)*(Lc1); // Lfc
			L1 += *(P1 + 118)*(Lq1); // Lfq
			L1 += *(P1 + 139)*(Le1); // Lfe
			L1 += *(P1 + 160)*(Lg1); // Lfg
			L1 += *(P1 + 181)*(Lh1); // Lfh
			L1 += *(P1 + 202)*(Li1); // Lfi
			L1 += *(P1 + 223)*(Ll1); // Lfl
			L1 += *(P1 + 244)*(Lk1); // Lfk
			L1 += *(P1 + 265)*(Lm1); // Lfm
			L1 += *(P1 + 286)*(Lf1); // Lff
			L1 += *(P1 + 307)*(Lp1); // Lfp
			L1 += *(P1 + 328)*(Lz1); // Lfs
			L1 += *(P1 + 349)*(Lt1); // Lft
			L1 += *(P1 + 370)*(Lw1); // Lfw
			L1 += *(P1 + 391)*(Ly1); // Lfy
			L1 += *(P1 + 412)*(Lv1); // Lfv
			L1 += *(P1 + 433)*(Lb1); // Lfb
			*(Ls3 + 13) = L1;
			
			// L(P)
			L1 = *(P1 + 14)*(La1); // Lpa
			L1 += *(P1 + 35)*(Lr1); // Lpr
			L1 += *(P1 + 56)*(Ln1); // Lpn
			L1 += *(P1 + 77)*(Ld1); // Lpd
			L1 += *(P1 + 98)*(Lc1); // Lpc
			L1 += *(P1 + 119)*(Lq1); // Lpq
			L1 += *(P1 + 140)*(Le1); // Lpe
			L1 += *(P1 + 161)*(Lg1); // Lpg
			L1 += *(P1 + 182)*(Lh1); // Lph
			L1 += *(P1 + 203)*(Li1); // Lpi
			L1 += *(P1 + 224)*(Ll1); // Lpl
			L1 += *(P1 + 245)*(Lk1); // Lpk
			L1 += *(P1 + 266)*(Lm1); // Lpm
			L1 += *(P1 + 287)*(Lf1); // Lpf
			L1 += *(P1 + 308)*(Lp1); // Lpp
			L1 += *(P1 + 329)*(Lz1); // Lps
			L1 += *(P1 + 350)*(Lt1); // Lpt
			L1 += *(P1 + 371)*(Lw1); // Lpw
			L1 += *(P1 + 392)*(Ly1); // Lpy
			L1 += *(P1 + 413)*(Lv1); // Lpv
			L1 += *(P1 + 434)*(Lb1); // Lpb
			*(Ls3 + 14) = L1;
			
			// L(S)
			L1 = *(P1 + 15)*(La1); // Lsa
			L1 += *(P1 + 36)*(Lr1); // Lsr
			L1 += *(P1 + 57)*(Ln1); // Lsn
			L1 += *(P1 + 78)*(Ld1); // Lsd
			L1 += *(P1 + 99)*(Lc1); // Lsc
			L1 += *(P1 + 120)*(Lq1); // Lsq
			L1 += *(P1 + 141)*(Le1); // Lse
			L1 += *(P1 + 162)*(Lg1); // Lsg
			L1 += *(P1 + 183)*(Lh1); // Lsh
			L1 += *(P1 + 204)*(Li1); // Lsi
			L1 += *(P1 + 225)*(Ll1); // Lsl
			L1 += *(P1 + 246)*(Lk1); // Lsk
			L1 += *(P1 + 267)*(Lm1); // Lsm
			L1 += *(P1 + 288)*(Lf1); // Lsf
			L1 += *(P1 + 309)*(Lp1); // Lsp
			L1 += *(P1 + 330)*(Lz1); // Lss
			L1 += *(P1 + 351)*(Lt1); // Lst
			L1 += *(P1 + 372)*(Lw1); // Lsw
			L1 += *(P1 + 393)*(Ly1); // Lsy
			L1 += *(P1 + 414)*(Lv1); // Lsv
			L1 += *(P1 + 435)*(Lb1); // Lsb
			*(Ls3 + 15) = L1;
			
			// L(T)
			L1 = *(P1 + 16)*(La1); // Lta
			L1 += *(P1 + 37)*(Lr1); // Ltr
			L1 += *(P1 + 58)*(Ln1); // Ltn
			L1 += *(P1 + 79)*(Ld1); // Ltd
			L1 += *(P1 + 100)*(Lc1); // Ltc
			L1 += *(P1 + 121)*(Lq1); // Ltq
			L1 += *(P1 + 142)*(Le1); // Lte
			L1 += *(P1 + 163)*(Lg1); // Ltg
			L1 += *(P1 + 184)*(Lh1); // Lth
			L1 += *(P1 + 205)*(Li1); // Lti
			L1 += *(P1 + 226)*(Ll1); // Ltl
			L1 += *(P1 + 247)*(Lk1); // Ltk
			L1 += *(P1 + 268)*(Lm1); // Ltm
			L1 += *(P1 + 289)*(Lf1); // Ltf
			L1 += *(P1 + 310)*(Lp1); // Ltp
			L1 += *(P1 + 331)*(Lz1); // Lts
			L1 += *(P1 + 352)*(Lt1); // Ltt
			L1 += *(P1 + 373)*(Lw1); // Ltw
			L1 += *(P1 + 394)*(Ly1); // Lty
			L1 += *(P1 + 415)*(Lv1); // Ltv
			L1 += *(P1 + 436)*(Lb1); // Ltb
			*(Ls3 + 16) = L1;
			
			// L(W)
			L1 = *(P1 + 17)*(La1); // Lwa
			L1 += *(P1 + 38)*(Lr1); // Lwr
			L1 += *(P1 + 59)*(Ln1); // Lwn
			L1 += *(P1 + 80)*(Ld1); // Lwd
			L1 += *(P1 + 101)*(Lc1); // Lwc
			L1 += *(P1 + 122)*(Lq1); // Lwq
			L1 += *(P1 + 143)*(Le1); // Lwe
			L1 += *(P1 + 164)*(Lg1); // Lwg
			L1 += *(P1 + 185)*(Lh1); // Lwh
			L1 += *(P1 + 206)*(Li1); // Lwi
			L1 += *(P1 + 227)*(Ll1); // Lwl
			L1 += *(P1 + 248)*(Lk1); // Lwk
			L1 += *(P1 + 269)*(Lm1); // Lwm
			L1 += *(P1 + 290)*(Lf1); // Lwf
			L1 += *(P1 + 311)*(Lp1); // Lwp
			L1 += *(P1 + 332)*(Lz1); // Lws
			L1 += *(P1 + 353)*(Lt1); // Lwt
			L1 += *(P1 + 374)*(Lw1); // Lww
			L1 += *(P1 + 395)*(Ly1); // Lwy
			L1 += *(P1 + 416)*(Lv1); // Lwv
			L1 += *(P1 + 437)*(Lb1); // Lwb
			*(Ls3 + 17) = L1;
			
			// L(Y)
			L1 = *(P1 + 18)*(La1); // Lya
			L1 += *(P1 + 39)*(Lr1); // Lyr
			L1 += *(P1 + 60)*(Ln1); // Lyn
			L1 += *(P1 + 81)*(Ld1); // Lyd
			L1 += *(P1 + 102)*(Lc1); // Lyc
			L1 += *(P1 + 123)*(Lq1); // Lyq
			L1 += *(P1 + 144)*(Le1); // Lye
			L1 += *(P1 + 165)*(Lg1); // Lyg
			L1 += *(P1 + 186)*(Lh1); // Lyh
			L1 += *(P1 + 207)*(Li1); // Lyi
			L1 += *(P1 + 228)*(Ll1); // Lyl
			L1 += *(P1 + 249)*(Lk1); // Lyk
			L1 += *(P1 + 270)*(Lm1); // Lym
			L1 += *(P1 + 291)*(Lf1); // Lyf
			L1 += *(P1 + 312)*(Lp1); // Lyp
			L1 += *(P1 + 333)*(Lz1); // Lys
			L1 += *(P1 + 354)*(Lt1); // Lyt
			L1 += *(P1 + 375)*(Lw1); // Lyw
			L1 += *(P1 + 396)*(Ly1); // Lyy
			L1 += *(P1 + 417)*(Lv1); // Lyv
			L1 += *(P1 + 438)*(Lb1); // Lyb
			*(Ls3 + 18) = L1;
			
			// L(V)
			L1 = *(P1 + 19)*(La1); // Lva
			L1 += *(P1 + 40)*(Lr1); // Lvr
			L1 += *(P1 + 61)*(Ln1); // Lvn
			L1 += *(P1 + 82)*(Ld1); // Lvd
			L1 += *(P1 + 103)*(Lc1); // Lvc
			L1 += *(P1 + 124)*(Lq1); // Lvq
			L1 += *(P1 + 145)*(Le1); // Lve
			L1 += *(P1 + 166)*(Lg1); // Lvg
			L1 += *(P1 + 187)*(Lh1); // Lvh
			L1 += *(P1 + 208)*(Li1); // Lvi
			L1 += *(P1 + 229)*(Ll1); // Lvl
			L1 += *(P1 + 250)*(Lk1); // Lvk
			L1 += *(P1 + 271)*(Lm1); // Lvm
			L1 += *(P1 + 292)*(Lf1); // Lvf
			L1 += *(P1 + 313)*(Lp1); // Lvp
			L1 += *(P1 + 334)*(Lz1); // Lvs
			L1 += *(P1 + 355)*(Lt1); // Lvt
			L1 += *(P1 + 376)*(Lw1); // Lvw
			L1 += *(P1 + 397)*(Ly1); // Lvy
			L1 += *(P1 + 418)*(Lv1); // Lvv
			L1 += *(P1 + 439)*(Lb1); // Lvb
			*(Ls3 + 19) = L1;
			
			// L(Indels)
			L1 = *(P1 + 20)*(La1); // Lba
			L1 += *(P1 + 41)*(Lr1); // Lbr
			L1 += *(P1 + 62)*(Ln1); // Lbn
			L1 += *(P1 + 83)*(Ld1); // Lbd
			L1 += *(P1 + 104)*(Lc1); // Lbc
			L1 += *(P1 + 125)*(Lq1); // Lbq
			L1 += *(P1 + 146)*(Le1); // Lbe
			L1 += *(P1 + 167)*(Lg1); // Lbg
			L1 += *(P1 + 188)*(Lh1); // Lbh
			L1 += *(P1 + 209)*(Li1); // Lbi
			L1 += *(P1 + 230)*(Ll1); // Lbl
			L1 += *(P1 + 251)*(Lk1); // Lbk
			L1 += *(P1 + 272)*(Lm1); // Lbm
			L1 += *(P1 + 293)*(Lf1); // Lbf
			L1 += *(P1 + 314)*(Lp1); // Lbp
			L1 += *(P1 + 335)*(Lz1); // Lbs
			L1 += *(P1 + 356)*(Lt1); // Lbt
			L1 += *(P1 + 377)*(Lw1); // Lbw
			L1 += *(P1 + 398)*(Ly1); // Lby
			L1 += *(P1 + 419)*(Lv1); // Lbv
			L1 += *(P1 + 440)*(Lb1); // Lbb
			*(Ls3 + 20) = L1;
			
			*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon) ||
				(*(Ls3 + 5) > 0 && *(Ls3 + 5) < inv_epsilon) ||
				(*(Ls3 + 6) > 0 && *(Ls3 + 6) < inv_epsilon) ||
				(*(Ls3 + 7) > 0 && *(Ls3 + 7) < inv_epsilon) ||
				(*(Ls3 + 8) > 0 && *(Ls3 + 8) < inv_epsilon) ||
				(*(Ls3 + 9) > 0 && *(Ls3 + 9) < inv_epsilon) ||
				(*(Ls3 + 10) > 0 && *(Ls3 + 10) < inv_epsilon) ||
				(*(Ls3 + 11) > 0 && *(Ls3 + 11) < inv_epsilon) ||
				(*(Ls3 + 12) > 0 && *(Ls3 + 12) < inv_epsilon) ||
				(*(Ls3 + 13) > 0 && *(Ls3 + 13) < inv_epsilon) ||
				(*(Ls3 + 14) > 0 && *(Ls3 + 14) < inv_epsilon) ||
				(*(Ls3 + 15) > 0 && *(Ls3 + 15) < inv_epsilon) ||
				(*(Ls3 + 16) > 0 && *(Ls3 + 16) < inv_epsilon) ||
				(*(Ls3 + 17) > 0 && *(Ls3 + 17) < inv_epsilon) ||
				(*(Ls3 + 18) > 0 && *(Ls3 + 18) < inv_epsilon) ||
				(*(Ls3 + 19) > 0 && *(Ls3 + 19) < inv_epsilon) ||
				(*(Ls3 + 20) > 0 && *(Ls3 + 20) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) *= epsilon;
				*(Ls3 + 6) *= epsilon;
				*(Ls3 + 7) *= epsilon;
				*(Ls3 + 8) *= epsilon;
				*(Ls3 + 9) *= epsilon;
				*(Ls3 + 10) *= epsilon;
				*(Ls3 + 11) *= epsilon;
				*(Ls3 + 12) *= epsilon;
				*(Ls3 + 13) *= epsilon;
				*(Ls3 + 14) *= epsilon;
				*(Ls3 + 15) *= epsilon;
				*(Ls3 + 16) *= epsilon;
				*(Ls3 + 17) *= epsilon;
				*(Ls3 + 18) *= epsilon;
				*(Ls3 + 19) *= epsilon;
				*(Ls3 + 20) *= epsilon;
				*(Ls3 + 21) += 1;
			}
		}
	} else {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0 ||
			Lb2 != 0) {
			// first branch can be disregarded
			
			// L(A)
			L2 = *(P2 + 0)*(La2); // Laa
			L2 += *(P2 + 21)*(Lr2); // Lar
			L2 += *(P2 + 42)*(Ln2); // Lan
			L2 += *(P2 + 63)*(Ld2); // Lad
			L2 += *(P2 + 84)*(Lc2); // Lac
			L2 += *(P2 + 105)*(Lq2); // Laq
			L2 += *(P2 + 126)*(Le2); // Lae
			L2 += *(P2 + 147)*(Lg2); // Lag
			L2 += *(P2 + 168)*(Lh2); // Lah
			L2 += *(P2 + 189)*(Li2); // Lai
			L2 += *(P2 + 210)*(Ll2); // Lal
			L2 += *(P2 + 231)*(Lk2); // Lak
			L2 += *(P2 + 252)*(Lm2); // Lam
			L2 += *(P2 + 273)*(Lf2); // Laf
			L2 += *(P2 + 294)*(Lp2); // Lap
			L2 += *(P2 + 315)*(Lz2); // Las
			L2 += *(P2 + 336)*(Lt2); // Lat
			L2 += *(P2 + 357)*(Lw2); // Law
			L2 += *(P2 + 378)*(Ly2); // Lay
			L2 += *(P2 + 399)*(Lv2); // Lav
			L2 += *(P2 + 420)*(Lb2); // Lab
			*(Ls3 + 0) = L2;
			
			// L(R)
			L2 = *(P2 + 1)*(La2); // Lra
			L2 += *(P2 + 22)*(Lr2); // Lrr
			L2 += *(P2 + 43)*(Ln2); // Lrn
			L2 += *(P2 + 64)*(Ld2); // Lrd
			L2 += *(P2 + 85)*(Lc2); // Lrc
			L2 += *(P2 + 106)*(Lq2); // Lrq
			L2 += *(P2 + 127)*(Le2); // Lre
			L2 += *(P2 + 148)*(Lg2); // Lrg
			L2 += *(P2 + 169)*(Lh2); // Lrh
			L2 += *(P2 + 190)*(Li2); // Lri
			L2 += *(P2 + 211)*(Ll2); // Lrl
			L2 += *(P2 + 232)*(Lk2); // Lrk
			L2 += *(P2 + 253)*(Lm2); // Lrm
			L2 += *(P2 + 274)*(Lf2); // Lrf
			L2 += *(P2 + 295)*(Lp2); // Lrp
			L2 += *(P2 + 316)*(Lz2); // Lrs
			L2 += *(P2 + 337)*(Lt2); // Lrt
			L2 += *(P2 + 358)*(Lw2); // Lrw
			L2 += *(P2 + 379)*(Ly2); // Lry
			L2 += *(P2 + 400)*(Lv2); // Lrv
			L2 += *(P2 + 421)*(Lb2); // Lrb
			*(Ls3 + 1) = L2;
			
			// L(N)
			L2 = *(P2 + 2)*(La2); // Lna
			L2 += *(P2 + 23)*(Lr2); // Lnr
			L2 += *(P2 + 44)*(Ln2); // Lnn
			L2 += *(P2 + 65)*(Ld2); // Lnd
			L2 += *(P2 + 86)*(Lc2); // Lnc
			L2 += *(P2 + 107)*(Lq2); // Lnq
			L2 += *(P2 + 128)*(Le2); // Lne
			L2 += *(P2 + 149)*(Lg2); // Lng
			L2 += *(P2 + 170)*(Lh2); // Lnh
			L2 += *(P2 + 191)*(Li2); // Lni
			L2 += *(P2 + 212)*(Ll2); // Lnl
			L2 += *(P2 + 233)*(Lk2); // Lnk
			L2 += *(P2 + 254)*(Lm2); // Lnm
			L2 += *(P2 + 275)*(Lf2); // Lnf
			L2 += *(P2 + 296)*(Lp2); // Lnp
			L2 += *(P2 + 317)*(Lz2); // Lns
			L2 += *(P2 + 338)*(Lt2); // Lnt
			L2 += *(P2 + 359)*(Lw2); // Lnw
			L2 += *(P2 + 380)*(Ly2); // Lny
			L2 += *(P2 + 401)*(Lv2); // Lnv
			L2 += *(P2 + 422)*(Lb2); // Lnb
			*(Ls3 + 2) = L2;
			
			// L(D)
			L2 = *(P2 + 3)*(La2); // Lda
			L2 += *(P2 + 24)*(Lr2); // Ldr
			L2 += *(P2 + 45)*(Ln2); // Ldn
			L2 += *(P2 + 66)*(Ld2); // Ldd
			L2 += *(P2 + 87)*(Lc2); // Ldc
			L2 += *(P2 + 108)*(Lq2); // Ldq
			L2 += *(P2 + 129)*(Le2); // Lde
			L2 += *(P2 + 150)*(Lg2); // Ldg
			L2 += *(P2 + 171)*(Lh2); // Ldh
			L2 += *(P2 + 192)*(Li2); // Ldi
			L2 += *(P2 + 213)*(Ll2); // Ldl
			L2 += *(P2 + 234)*(Lk2); // Ldk
			L2 += *(P2 + 255)*(Lm2); // Ldm
			L2 += *(P2 + 276)*(Lf2); // Ldf
			L2 += *(P2 + 297)*(Lp2); // Ldp
			L2 += *(P2 + 318)*(Lz2); // Lds
			L2 += *(P2 + 339)*(Lt2); // Ldt
			L2 += *(P2 + 360)*(Lw2); // Ldw
			L2 += *(P2 + 381)*(Ly2); // Ldy
			L2 += *(P2 + 402)*(Lv2); // Ldv
			L2 += *(P2 + 423)*(Lb2); // Ldb
			*(Ls3 + 3) = L2;
			
			// L(C)
			L2 = *(P2 + 4)*(La2); // Lca
			L2 += *(P2 + 25)*(Lr2); // Lcr
			L2 += *(P2 + 46)*(Ln2); // Lcn
			L2 += *(P2 + 67)*(Ld2); // Lcd
			L2 += *(P2 + 88)*(Lc2); // Lcc
			L2 += *(P2 + 109)*(Lq2); // Lcq
			L2 += *(P2 + 130)*(Le2); // Lce
			L2 += *(P2 + 151)*(Lg2); // Lcg
			L2 += *(P2 + 172)*(Lh2); // Lch
			L2 += *(P2 + 193)*(Li2); // Lci
			L2 += *(P2 + 214)*(Ll2); // Lcl
			L2 += *(P2 + 235)*(Lk2); // Lck
			L2 += *(P2 + 256)*(Lm2); // Lcm
			L2 += *(P2 + 277)*(Lf2); // Lcf
			L2 += *(P2 + 298)*(Lp2); // Lcp
			L2 += *(P2 + 319)*(Lz2); // Lcs
			L2 += *(P2 + 340)*(Lt2); // Lct
			L2 += *(P2 + 361)*(Lw2); // Lcw
			L2 += *(P2 + 382)*(Ly2); // Lcy
			L2 += *(P2 + 403)*(Lv2); // Lcv
			L2 += *(P2 + 424)*(Lb2); // Lcb
			*(Ls3 + 4) = L2;
			
			// L(Q)
			L2 = *(P2 + 5)*(La2); // Lqa
			L2 += *(P2 + 26)*(Lr2); // Lqr
			L2 += *(P2 + 47)*(Ln2); // Lqn
			L2 += *(P2 + 68)*(Ld2); // Lqd
			L2 += *(P2 + 89)*(Lc2); // Lqc
			L2 += *(P2 + 110)*(Lq2); // Lqq
			L2 += *(P2 + 131)*(Le2); // Lqe
			L2 += *(P2 + 152)*(Lg2); // Lqg
			L2 += *(P2 + 173)*(Lh2); // Lqh
			L2 += *(P2 + 194)*(Li2); // Lqi
			L2 += *(P2 + 215)*(Ll2); // Lql
			L2 += *(P2 + 236)*(Lk2); // Lqk
			L2 += *(P2 + 257)*(Lm2); // Lqm
			L2 += *(P2 + 278)*(Lf2); // Lqf
			L2 += *(P2 + 299)*(Lp2); // Lqp
			L2 += *(P2 + 320)*(Lz2); // Lqs
			L2 += *(P2 + 341)*(Lt2); // Lqt
			L2 += *(P2 + 362)*(Lw2); // Lqw
			L2 += *(P2 + 383)*(Ly2); // Lqy
			L2 += *(P2 + 404)*(Lv2); // Lqv
			L2 += *(P2 + 425)*(Lb2); // Lqb
			*(Ls3 + 5) = L2;
			
			// L(E)
			L2 = *(P2 + 6)*(La2); // Lea
			L2 += *(P2 + 27)*(Lr2); // Ler
			L2 += *(P2 + 48)*(Ln2); // Len
			L2 += *(P2 + 69)*(Ld2); // Led
			L2 += *(P2 + 90)*(Lc2); // Lec
			L2 += *(P2 + 111)*(Lq2); // Leq
			L2 += *(P2 + 132)*(Le2); // Lee
			L2 += *(P2 + 153)*(Lg2); // Leg
			L2 += *(P2 + 174)*(Lh2); // Leh
			L2 += *(P2 + 195)*(Li2); // Lei
			L2 += *(P2 + 216)*(Ll2); // Lel
			L2 += *(P2 + 237)*(Lk2); // Lek
			L2 += *(P2 + 258)*(Lm2); // Lem
			L2 += *(P2 + 279)*(Lf2); // Lef
			L2 += *(P2 + 300)*(Lp2); // Lep
			L2 += *(P2 + 321)*(Lz2); // Les
			L2 += *(P2 + 342)*(Lt2); // Let
			L2 += *(P2 + 363)*(Lw2); // Lew
			L2 += *(P2 + 384)*(Ly2); // Ley
			L2 += *(P2 + 405)*(Lv2); // Lev
			L2 += *(P2 + 426)*(Lb2); // Leb
			*(Ls3 + 6) = L2;
			
			// L(G)
			L2 = *(P2 + 7)*(La2); // Lga
			L2 += *(P2 + 28)*(Lr2); // Lgr
			L2 += *(P2 + 49)*(Ln2); // Lgn
			L2 += *(P2 + 70)*(Ld2); // Lgd
			L2 += *(P2 + 91)*(Lc2); // Lgc
			L2 += *(P2 + 112)*(Lq2); // Lgq
			L2 += *(P2 + 133)*(Le2); // Lge
			L2 += *(P2 + 154)*(Lg2); // Lgg
			L2 += *(P2 + 175)*(Lh2); // Lgh
			L2 += *(P2 + 196)*(Li2); // Lgi
			L2 += *(P2 + 217)*(Ll2); // Lgl
			L2 += *(P2 + 238)*(Lk2); // Lgk
			L2 += *(P2 + 259)*(Lm2); // Lgm
			L2 += *(P2 + 280)*(Lf2); // Lgf
			L2 += *(P2 + 301)*(Lp2); // Lgp
			L2 += *(P2 + 322)*(Lz2); // Lgs
			L2 += *(P2 + 343)*(Lt2); // Lgt
			L2 += *(P2 + 364)*(Lw2); // Lgw
			L2 += *(P2 + 385)*(Ly2); // Lgy
			L2 += *(P2 + 406)*(Lv2); // Lgv
			L2 += *(P2 + 427)*(Lb2); // Lgb
			*(Ls3 + 7) = L2;
			
			// L(H)
			L2 = *(P2 + 8)*(La2); // Lha
			L2 += *(P2 + 29)*(Lr2); // Lhr
			L2 += *(P2 + 50)*(Ln2); // Lhn
			L2 += *(P2 + 71)*(Ld2); // Lhd
			L2 += *(P2 + 92)*(Lc2); // Lhc
			L2 += *(P2 + 113)*(Lq2); // Lhq
			L2 += *(P2 + 134)*(Le2); // Lhe
			L2 += *(P2 + 155)*(Lg2); // Lhg
			L2 += *(P2 + 176)*(Lh2); // Lhh
			L2 += *(P2 + 197)*(Li2); // Lhi
			L2 += *(P2 + 218)*(Ll2); // Lhl
			L2 += *(P2 + 239)*(Lk2); // Lhk
			L2 += *(P2 + 260)*(Lm2); // Lhm
			L2 += *(P2 + 281)*(Lf2); // Lhf
			L2 += *(P2 + 302)*(Lp2); // Lhp
			L2 += *(P2 + 323)*(Lz2); // Lhs
			L2 += *(P2 + 344)*(Lt2); // Lht
			L2 += *(P2 + 365)*(Lw2); // Lhw
			L2 += *(P2 + 386)*(Ly2); // Lhy
			L2 += *(P2 + 407)*(Lv2); // Lhv
			L2 += *(P2 + 428)*(Lb2); // Lhb
			*(Ls3 + 8) = L2;
			
			// L(I)
			L2 = *(P2 + 9)*(La2); // Lia
			L2 += *(P2 + 30)*(Lr2); // Lir
			L2 += *(P2 + 51)*(Ln2); // Lin
			L2 += *(P2 + 72)*(Ld2); // Lid
			L2 += *(P2 + 93)*(Lc2); // Lic
			L2 += *(P2 + 114)*(Lq2); // Liq
			L2 += *(P2 + 135)*(Le2); // Lie
			L2 += *(P2 + 156)*(Lg2); // Lig
			L2 += *(P2 + 177)*(Lh2); // Lih
			L2 += *(P2 + 198)*(Li2); // Lii
			L2 += *(P2 + 219)*(Ll2); // Lil
			L2 += *(P2 + 240)*(Lk2); // Lik
			L2 += *(P2 + 261)*(Lm2); // Lim
			L2 += *(P2 + 282)*(Lf2); // Lif
			L2 += *(P2 + 303)*(Lp2); // Lip
			L2 += *(P2 + 324)*(Lz2); // Lis
			L2 += *(P2 + 345)*(Lt2); // Lit
			L2 += *(P2 + 366)*(Lw2); // Liw
			L2 += *(P2 + 387)*(Ly2); // Liy
			L2 += *(P2 + 408)*(Lv2); // Liv
			L2 += *(P2 + 429)*(Lb2); // Lib
			*(Ls3 + 9) = L2;
			
			// L(L)
			L2 = *(P2 + 10)*(La2); // Lla
			L2 += *(P2 + 31)*(Lr2); // Llr
			L2 += *(P2 + 52)*(Ln2); // Lln
			L2 += *(P2 + 73)*(Ld2); // Lld
			L2 += *(P2 + 94)*(Lc2); // Llc
			L2 += *(P2 + 115)*(Lq2); // Llq
			L2 += *(P2 + 136)*(Le2); // Lle
			L2 += *(P2 + 157)*(Lg2); // Llg
			L2 += *(P2 + 178)*(Lh2); // Llh
			L2 += *(P2 + 199)*(Li2); // Lli
			L2 += *(P2 + 220)*(Ll2); // Lll
			L2 += *(P2 + 241)*(Lk2); // Llk
			L2 += *(P2 + 262)*(Lm2); // Llm
			L2 += *(P2 + 283)*(Lf2); // Llf
			L2 += *(P2 + 304)*(Lp2); // Llp
			L2 += *(P2 + 325)*(Lz2); // Lls
			L2 += *(P2 + 346)*(Lt2); // Llt
			L2 += *(P2 + 367)*(Lw2); // Llw
			L2 += *(P2 + 388)*(Ly2); // Lly
			L2 += *(P2 + 409)*(Lv2); // Llv
			L2 += *(P2 + 430)*(Lb2); // Llb
			*(Ls3 + 10) = L2;
			
			// L(K)
			L2 = *(P2 + 11)*(La2); // Lka
			L2 += *(P2 + 32)*(Lr2); // Lkr
			L2 += *(P2 + 53)*(Ln2); // Lkn
			L2 += *(P2 + 74)*(Ld2); // Lkd
			L2 += *(P2 + 95)*(Lc2); // Lkc
			L2 += *(P2 + 116)*(Lq2); // Lkq
			L2 += *(P2 + 137)*(Le2); // Lke
			L2 += *(P2 + 158)*(Lg2); // Lkg
			L2 += *(P2 + 179)*(Lh2); // Lkh
			L2 += *(P2 + 200)*(Li2); // Lki
			L2 += *(P2 + 221)*(Ll2); // Lkl
			L2 += *(P2 + 242)*(Lk2); // Lkk
			L2 += *(P2 + 263)*(Lm2); // Lkm
			L2 += *(P2 + 284)*(Lf2); // Lkf
			L2 += *(P2 + 305)*(Lp2); // Lkp
			L2 += *(P2 + 326)*(Lz2); // Lks
			L2 += *(P2 + 347)*(Lt2); // Lkt
			L2 += *(P2 + 368)*(Lw2); // Lkw
			L2 += *(P2 + 389)*(Ly2); // Lky
			L2 += *(P2 + 410)*(Lv2); // Lkv
			L2 += *(P2 + 431)*(Lb2); // Lkb
			*(Ls3 + 11) = L2;
			
			// L(M)
			L2 = *(P2 + 12)*(La2); // Lma
			L2 += *(P2 + 33)*(Lr2); // Lmr
			L2 += *(P2 + 54)*(Ln2); // Lmn
			L2 += *(P2 + 75)*(Ld2); // Lmd
			L2 += *(P2 + 96)*(Lc2); // Lmc
			L2 += *(P2 + 117)*(Lq2); // Lmq
			L2 += *(P2 + 138)*(Le2); // Lme
			L2 += *(P2 + 159)*(Lg2); // Lmg
			L2 += *(P2 + 180)*(Lh2); // Lmh
			L2 += *(P2 + 201)*(Li2); // Lmi
			L2 += *(P2 + 222)*(Ll2); // Lml
			L2 += *(P2 + 243)*(Lk2); // Lmk
			L2 += *(P2 + 264)*(Lm2); // Lmm
			L2 += *(P2 + 285)*(Lf2); // Lmf
			L2 += *(P2 + 306)*(Lp2); // Lmp
			L2 += *(P2 + 327)*(Lz2); // Lms
			L2 += *(P2 + 348)*(Lt2); // Lmt
			L2 += *(P2 + 369)*(Lw2); // Lmw
			L2 += *(P2 + 390)*(Ly2); // Lmy
			L2 += *(P2 + 411)*(Lv2); // Lmv
			L2 += *(P2 + 432)*(Lb2); // Lmb
			*(Ls3 + 12) = L2;
			
			// L(F)
			L2 = *(P2 + 13)*(La2); // Lfa
			L2 += *(P2 + 34)*(Lr2); // Lfr
			L2 += *(P2 + 55)*(Ln2); // Lfn
			L2 += *(P2 + 76)*(Ld2); // Lfd
			L2 += *(P2 + 97)*(Lc2); // Lfc
			L2 += *(P2 + 118)*(Lq2); // Lfq
			L2 += *(P2 + 139)*(Le2); // Lfe
			L2 += *(P2 + 160)*(Lg2); // Lfg
			L2 += *(P2 + 181)*(Lh2); // Lfh
			L2 += *(P2 + 202)*(Li2); // Lfi
			L2 += *(P2 + 223)*(Ll2); // Lfl
			L2 += *(P2 + 244)*(Lk2); // Lfk
			L2 += *(P2 + 265)*(Lm2); // Lfm
			L2 += *(P2 + 286)*(Lf2); // Lff
			L2 += *(P2 + 307)*(Lp2); // Lfp
			L2 += *(P2 + 328)*(Lz2); // Lfs
			L2 += *(P2 + 349)*(Lt2); // Lft
			L2 += *(P2 + 370)*(Lw2); // Lfw
			L2 += *(P2 + 391)*(Ly2); // Lfy
			L2 += *(P2 + 412)*(Lv2); // Lfv
			L2 += *(P2 + 433)*(Lb2); // Lfb
			*(Ls3 + 13) = L2;
			
			// L(P)
			L2 = *(P2 + 14)*(La2); // Lpa
			L2 += *(P2 + 35)*(Lr2); // Lpr
			L2 += *(P2 + 56)*(Ln2); // Lpn
			L2 += *(P2 + 77)*(Ld2); // Lpd
			L2 += *(P2 + 98)*(Lc2); // Lpc
			L2 += *(P2 + 119)*(Lq2); // Lpq
			L2 += *(P2 + 140)*(Le2); // Lpe
			L2 += *(P2 + 161)*(Lg2); // Lpg
			L2 += *(P2 + 182)*(Lh2); // Lph
			L2 += *(P2 + 203)*(Li2); // Lpi
			L2 += *(P2 + 224)*(Ll2); // Lpl
			L2 += *(P2 + 245)*(Lk2); // Lpk
			L2 += *(P2 + 266)*(Lm2); // Lpm
			L2 += *(P2 + 287)*(Lf2); // Lpf
			L2 += *(P2 + 308)*(Lp2); // Lpp
			L2 += *(P2 + 329)*(Lz2); // Lps
			L2 += *(P2 + 350)*(Lt2); // Lpt
			L2 += *(P2 + 371)*(Lw2); // Lpw
			L2 += *(P2 + 392)*(Ly2); // Lpy
			L2 += *(P2 + 413)*(Lv2); // Lpv
			L2 += *(P2 + 434)*(Lb2); // Lpb
			*(Ls3 + 14) = L2;
			
			// L(S)
			L2 = *(P2 + 15)*(La2); // Lsa
			L2 += *(P2 + 36)*(Lr2); // Lsr
			L2 += *(P2 + 57)*(Ln2); // Lsn
			L2 += *(P2 + 78)*(Ld2); // Lsd
			L2 += *(P2 + 99)*(Lc2); // Lsc
			L2 += *(P2 + 120)*(Lq2); // Lsq
			L2 += *(P2 + 141)*(Le2); // Lse
			L2 += *(P2 + 162)*(Lg2); // Lsg
			L2 += *(P2 + 183)*(Lh2); // Lsh
			L2 += *(P2 + 204)*(Li2); // Lsi
			L2 += *(P2 + 225)*(Ll2); // Lsl
			L2 += *(P2 + 246)*(Lk2); // Lsk
			L2 += *(P2 + 267)*(Lm2); // Lsm
			L2 += *(P2 + 288)*(Lf2); // Lsf
			L2 += *(P2 + 309)*(Lp2); // Lsp
			L2 += *(P2 + 330)*(Lz2); // Lss
			L2 += *(P2 + 351)*(Lt2); // Lst
			L2 += *(P2 + 372)*(Lw2); // Lsw
			L2 += *(P2 + 393)*(Ly2); // Lsy
			L2 += *(P2 + 414)*(Lv2); // Lsv
			L2 += *(P2 + 435)*(Lb2); // Lsb
			*(Ls3 + 15) = L2;
			
			// L(T)
			L2 = *(P2 + 16)*(La2); // Lta
			L2 += *(P2 + 37)*(Lr2); // Ltr
			L2 += *(P2 + 58)*(Ln2); // Ltn
			L2 += *(P2 + 79)*(Ld2); // Ltd
			L2 += *(P2 + 100)*(Lc2); // Ltc
			L2 += *(P2 + 121)*(Lq2); // Ltq
			L2 += *(P2 + 142)*(Le2); // Lte
			L2 += *(P2 + 163)*(Lg2); // Ltg
			L2 += *(P2 + 184)*(Lh2); // Lth
			L2 += *(P2 + 205)*(Li2); // Lti
			L2 += *(P2 + 226)*(Ll2); // Ltl
			L2 += *(P2 + 247)*(Lk2); // Ltk
			L2 += *(P2 + 268)*(Lm2); // Ltm
			L2 += *(P2 + 289)*(Lf2); // Ltf
			L2 += *(P2 + 310)*(Lp2); // Ltp
			L2 += *(P2 + 331)*(Lz2); // Lts
			L2 += *(P2 + 352)*(Lt2); // Ltt
			L2 += *(P2 + 373)*(Lw2); // Ltw
			L2 += *(P2 + 394)*(Ly2); // Lty
			L2 += *(P2 + 415)*(Lv2); // Ltv
			L2 += *(P2 + 436)*(Lb2); // Ltb
			*(Ls3 + 16) = L2;
			
			// L(W)
			L2 = *(P2 + 17)*(La2); // Lwa
			L2 += *(P2 + 38)*(Lr2); // Lwr
			L2 += *(P2 + 59)*(Ln2); // Lwn
			L2 += *(P2 + 80)*(Ld2); // Lwd
			L2 += *(P2 + 101)*(Lc2); // Lwc
			L2 += *(P2 + 122)*(Lq2); // Lwq
			L2 += *(P2 + 143)*(Le2); // Lwe
			L2 += *(P2 + 164)*(Lg2); // Lwg
			L2 += *(P2 + 185)*(Lh2); // Lwh
			L2 += *(P2 + 206)*(Li2); // Lwi
			L2 += *(P2 + 227)*(Ll2); // Lwl
			L2 += *(P2 + 248)*(Lk2); // Lwk
			L2 += *(P2 + 269)*(Lm2); // Lwm
			L2 += *(P2 + 290)*(Lf2); // Lwf
			L2 += *(P2 + 311)*(Lp2); // Lwp
			L2 += *(P2 + 332)*(Lz2); // Lws
			L2 += *(P2 + 353)*(Lt2); // Lwt
			L2 += *(P2 + 374)*(Lw2); // Lww
			L2 += *(P2 + 395)*(Ly2); // Lwy
			L2 += *(P2 + 416)*(Lv2); // Lwv
			L2 += *(P2 + 437)*(Lb2); // Lwb
			*(Ls3 + 17) = L2;
			
			// L(Y)
			L2 = *(P2 + 18)*(La2); // Lya
			L2 += *(P2 + 39)*(Lr2); // Lyr
			L2 += *(P2 + 60)*(Ln2); // Lyn
			L2 += *(P2 + 81)*(Ld2); // Lyd
			L2 += *(P2 + 102)*(Lc2); // Lyc
			L2 += *(P2 + 123)*(Lq2); // Lyq
			L2 += *(P2 + 144)*(Le2); // Lye
			L2 += *(P2 + 165)*(Lg2); // Lyg
			L2 += *(P2 + 186)*(Lh2); // Lyh
			L2 += *(P2 + 207)*(Li2); // Lyi
			L2 += *(P2 + 228)*(Ll2); // Lyl
			L2 += *(P2 + 249)*(Lk2); // Lyk
			L2 += *(P2 + 270)*(Lm2); // Lym
			L2 += *(P2 + 291)*(Lf2); // Lyf
			L2 += *(P2 + 312)*(Lp2); // Lyp
			L2 += *(P2 + 333)*(Lz2); // Lys
			L2 += *(P2 + 354)*(Lt2); // Lyt
			L2 += *(P2 + 375)*(Lw2); // Lyw
			L2 += *(P2 + 396)*(Ly2); // Lyy
			L2 += *(P2 + 417)*(Lv2); // Lyv
			L2 += *(P2 + 438)*(Lb2); // Lyb
			*(Ls3 + 18) = L2;
			
			// L(V)
			L2 = *(P2 + 19)*(La2); // Lva
			L2 += *(P2 + 40)*(Lr2); // Lvr
			L2 += *(P2 + 61)*(Ln2); // Lvn
			L2 += *(P2 + 82)*(Ld2); // Lvd
			L2 += *(P2 + 103)*(Lc2); // Lvc
			L2 += *(P2 + 124)*(Lq2); // Lvq
			L2 += *(P2 + 145)*(Le2); // Lve
			L2 += *(P2 + 166)*(Lg2); // Lvg
			L2 += *(P2 + 187)*(Lh2); // Lvh
			L2 += *(P2 + 208)*(Li2); // Lvi
			L2 += *(P2 + 229)*(Ll2); // Lvl
			L2 += *(P2 + 250)*(Lk2); // Lvk
			L2 += *(P2 + 271)*(Lm2); // Lvm
			L2 += *(P2 + 292)*(Lf2); // Lvf
			L2 += *(P2 + 313)*(Lp2); // Lvp
			L2 += *(P2 + 334)*(Lz2); // Lvs
			L2 += *(P2 + 355)*(Lt2); // Lvt
			L2 += *(P2 + 376)*(Lw2); // Lvw
			L2 += *(P2 + 397)*(Ly2); // Lvy
			L2 += *(P2 + 418)*(Lv2); // Lvv
			L2 += *(P2 + 439)*(Lb2); // Lvb
			*(Ls3 + 19) = L2;
			
			// L(Indels)
			L2 = *(P2 + 20)*(La2); // Lba
			L2 += *(P2 + 41)*(Lr2); // Lbr
			L2 += *(P2 + 62)*(Ln2); // Lbn
			L2 += *(P2 + 83)*(Ld2); // Lbd
			L2 += *(P2 + 104)*(Lc2); // Lbc
			L2 += *(P2 + 125)*(Lq2); // Lbq
			L2 += *(P2 + 146)*(Le2); // Lbe
			L2 += *(P2 + 167)*(Lg2); // Lbg
			L2 += *(P2 + 188)*(Lh2); // Lbh
			L2 += *(P2 + 209)*(Li2); // Lbi
			L2 += *(P2 + 230)*(Ll2); // Lbl
			L2 += *(P2 + 251)*(Lk2); // Lbk
			L2 += *(P2 + 272)*(Lm2); // Lbm
			L2 += *(P2 + 293)*(Lf2); // Lbf
			L2 += *(P2 + 314)*(Lp2); // Lbp
			L2 += *(P2 + 335)*(Lz2); // Lbs
			L2 += *(P2 + 356)*(Lt2); // Lbt
			L2 += *(P2 + 377)*(Lw2); // Lbw
			L2 += *(P2 + 398)*(Ly2); // Lby
			L2 += *(P2 + 419)*(Lv2); // Lbv
			L2 += *(P2 + 440)*(Lb2); // Lbb
			*(Ls3 + 20) = L2;
			
			if (root &&
				(La1 != 0 ||
				Lr1 != 0 ||
				Ln1 != 0 ||
				Ld1 != 0 ||
				Lc1 != 0 ||
				Lq1 != 0 ||
				Le1 != 0 ||
				Lg1 != 0 ||
				Lh1 != 0 ||
				Li1 != 0 ||
				Ll1 != 0 ||
				Lk1 != 0 ||
				Lm1 != 0 ||
				Lf1 != 0 ||
				Lp1 != 0 ||
				Lz1 != 0 ||
				Lt1 != 0 ||
				Lw1 != 0 ||
				Ly1 != 0 ||
				Lv1 != 0 ||
				Lb1 != 0)) {
				*(Ls3) *= La1;
				*(Ls3 + 1) *= Lr1;
				*(Ls3 + 2) *= Ln1;
				*(Ls3 + 3) *= Ld1;
				*(Ls3 + 4) *= Lc1;
				*(Ls3 + 5) *= Lq1;
				*(Ls3 + 6) *= Le1;
				*(Ls3 + 7) *= Lg1;
				*(Ls3 + 8) *= Lh1;
				*(Ls3 + 9) *= Li1;
				*(Ls3 + 10) *= Ll1;
				*(Ls3 + 11) *= Lk1;
				*(Ls3 + 12) *= Lm1;
				*(Ls3 + 13) *= Lf1;
				*(Ls3 + 14) *= Lp1;
				*(Ls3 + 15) *= Lz1;
				*(Ls3 + 16) *= Lt1;
				*(Ls3 + 17) *= Lw1;
				*(Ls3 + 18) *= Ly1;
				*(Ls3 + 19) *= Lv1;
				*(Ls3 + 20) *= Lb1;
				*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			} else {
				*(Ls3 + 21) = *(Ls2 + 21);
			}
			
			if ((*(Ls3) > 0 && *(Ls3) < inv_epsilon) ||
				(*(Ls3 + 1) > 0 && *(Ls3 + 1) < inv_epsilon) ||
				(*(Ls3 + 2) > 0 && *(Ls3 + 2) < inv_epsilon) ||
				(*(Ls3 + 3) > 0 && *(Ls3 + 3) < inv_epsilon) ||
				(*(Ls3 + 4) > 0 && *(Ls3 + 4) < inv_epsilon) ||
				(*(Ls3 + 5) > 0 && *(Ls3 + 5) < inv_epsilon) ||
				(*(Ls3 + 6) > 0 && *(Ls3 + 6) < inv_epsilon) ||
				(*(Ls3 + 7) > 0 && *(Ls3 + 7) < inv_epsilon) ||
				(*(Ls3 + 8) > 0 && *(Ls3 + 8) < inv_epsilon) ||
				(*(Ls3 + 9) > 0 && *(Ls3 + 9) < inv_epsilon) ||
				(*(Ls3 + 10) > 0 && *(Ls3 + 10) < inv_epsilon) ||
				(*(Ls3 + 11) > 0 && *(Ls3 + 11) < inv_epsilon) ||
				(*(Ls3 + 12) > 0 && *(Ls3 + 12) < inv_epsilon) ||
				(*(Ls3 + 13) > 0 && *(Ls3 + 13) < inv_epsilon) ||
				(*(Ls3 + 14) > 0 && *(Ls3 + 14) < inv_epsilon) ||
				(*(Ls3 + 15) > 0 && *(Ls3 + 15) < inv_epsilon) ||
				(*(Ls3 + 16) > 0 && *(Ls3 + 16) < inv_epsilon) ||
				(*(Ls3 + 17) > 0 && *(Ls3 + 17) < inv_epsilon) ||
				(*(Ls3 + 18) > 0 && *(Ls3 + 18) < inv_epsilon) ||
				(*(Ls3 + 19) > 0 && *(Ls3 + 19) < inv_epsilon) ||
				(*(Ls3 + 20) > 0 && *(Ls3 + 20) < inv_epsilon)) {
				*(Ls3) *= epsilon;
				*(Ls3 + 1) *= epsilon;
				*(Ls3 + 2) *= epsilon;
				*(Ls3 + 3) *= epsilon;
				*(Ls3 + 4) *= epsilon;
				*(Ls3 + 5) *= epsilon;
				*(Ls3 + 6) *= epsilon;
				*(Ls3 + 7) *= epsilon;
				*(Ls3 + 8) *= epsilon;
				*(Ls3 + 9) *= epsilon;
				*(Ls3 + 10) *= epsilon;
				*(Ls3 + 11) *= epsilon;
				*(Ls3 + 12) *= epsilon;
				*(Ls3 + 13) *= epsilon;
				*(Ls3 + 14) *= epsilon;
				*(Ls3 + 15) *= epsilon;
				*(Ls3 + 16) *= epsilon;
				*(Ls3 + 17) *= epsilon;
				*(Ls3 + 18) *= epsilon;
				*(Ls3 + 19) *= epsilon;
				*(Ls3 + 20) *= epsilon;
				*(Ls3 + 21) += 1;
			}
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lr1;
			*(Ls3 + 2) = Ln1;
			*(Ls3 + 3) = Ld1;
			*(Ls3 + 4) = Lc1;
			*(Ls3 + 5) = Lq1;
			*(Ls3 + 6) = Le1;
			*(Ls3 + 7) = Lg1;
			*(Ls3 + 8) = Lh1;
			*(Ls3 + 9) = Li1;
			*(Ls3 + 10) = Ll1;
			*(Ls3 + 11) = Lk1;
			*(Ls3 + 12) = Lm1;
			*(Ls3 + 13) = Lf1;
			*(Ls3 + 14) = Lp1;
			*(Ls3 + 15) = Lz1;
			*(Ls3 + 16) = Lt1;
			*(Ls3 + 17) = Lw1;
			*(Ls3 + 18) = Ly1;
			*(Ls3 + 19) = Lv1;
			*(Ls3 + 20) = Lb1;
			*(Ls3 + 21) = *(Ls1 + 21);
		}
	}
}

static void L_unknown_SIMD(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	typedef double __attribute__((vector_size(32))) v4d;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lc1 = *(Ls1 + 1);
	const double Lg1 = *(Ls1 + 2);
	const double Lt1 = *(Ls1 + 3);
	const double La2 = *(Ls2);
	const double Lc2 = *(Ls2 + 1);
	const double Lg2 = *(Ls2 + 2);
	const double Lt2 = *(Ls2 + 3);
	
	unsigned int flag = 0;
	
	if (root == 0 &&
		(La1 != 0 ||
		Lc1 != 0 ||
		Lg1 != 0 ||
		Lt1 != 0)) {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0) {
			// neither branch can be disregarded
			
			v4d *L1 = (v4d*)(Ls3);
			v4d *vP11 = (v4d*)(P1);
			v4d *vP12 = (v4d*)(P1 + 5);
			v4d *vP13 = (v4d*)(P1 + 10);
			v4d *vP14 = (v4d*)(P1 + 15);
			// broadcast likelihoods
			const v4d LA1 = {La1, La1, La1, La1};
			const v4d LC1 = {Lc1, Lc1, Lc1, Lc1};
			const v4d LG1 = {Lg1, Lg1, Lg1, Lg1};
			const v4d LT1 = {Lt1, Lt1, Lt1, Lt1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LC1;
			*L1 += *vP13*LG1;
			*L1 += *vP14*LT1;
			
			v4d L2;
			v4d *vP21 = (v4d*)(P2);
			v4d *vP22 = (v4d*)(P2 + 5);
			v4d *vP23 = (v4d*)(P2 + 10);
			v4d *vP24 = (v4d*)(P2 + 15);
			// broadcast likelihoods
			const v4d LA2 = {La2, La2, La2, La2};
			const v4d LC2 = {Lc2, Lc2, Lc2, Lc2};
			const v4d LG2 = {Lg2, Lg2, Lg2, Lg2};
			const v4d LT2 = {Lt2, Lt2, Lt2, Lt2};
			L2 = *vP21*LA2;
			L2 += *vP22*LC2;
			L2 += *vP23*LG2;
			L2 += *vP24*LT2;
			
			*L1 *= L2;
			for (int i = 0; i < 4; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
		} else {
			// second branch can be disregarded
			
			v4d *L1 = (v4d*)(Ls3);
			v4d *vP11 = (v4d*)(P1);
			v4d *vP12 = (v4d*)(P1 + 5);
			v4d *vP13 = (v4d*)(P1 + 10);
			v4d *vP14 = (v4d*)(P1 + 15);
			// broadcast likelihoods
			const v4d LA1 = {La1, La1, La1, La1};
			const v4d LC1 = {Lc1, Lc1, Lc1, Lc1};
			const v4d LG1 = {Lg1, Lg1, Lg1, Lg1};
			const v4d LT1 = {Lt1, Lt1, Lt1, Lt1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LC1;
			*L1 += *vP13*LG1;
			*L1 += *vP14*LT1;
			
			for (int i = 0; i < 4; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 5) = *(Ls1 + 5);
		}
	} else {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0) {
			// first branch can be disregarded
			
			v4d *L2 = (v4d*)(Ls3);
			v4d *vP21 = (v4d*)(P2);
			v4d *vP22 = (v4d*)(P2 + 5);
			v4d *vP23 = (v4d*)(P2 + 10);
			v4d *vP24 = (v4d*)(P2 + 15);
			// broadcast likelihoods
			const v4d LA2 = {La2, La2, La2, La2};
			const v4d LC2 = {Lc2, Lc2, Lc2, Lc2};
			const v4d LG2 = {Lg2, Lg2, Lg2, Lg2};
			const v4d LT2 = {Lt2, Lt2, Lt2, Lt2};
			*L2 = *vP21*LA2;
			*L2 += *vP22*LC2;
			*L2 += *vP23*LG2;
			*L2 += *vP24*LT2;
			
			if (root &&
				(La1 != 0 ||
				Lc1 != 0 ||
				Lg1 != 0 ||
				Lt1 != 0)) {
				const v4d L1 = {La1, Lc1, Lg1, Lt1};
				*L2 *= L1;
				*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
			} else {
				*(Ls3 + 5) = *(Ls2 + 5);
			}
			
			for (int i = 0; i < 4; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lc1;
			*(Ls3 + 2) = Lg1;
			*(Ls3 + 3) = Lt1;
			*(Ls3 + 5) = *(Ls1 + 5);
		}
	}
	
	*(Ls3 + 4) = 0;
	if (flag) {
		*(Ls3) *= epsilon;
		*(Ls3 + 1) *= epsilon;
		*(Ls3 + 2) *= epsilon;
		*(Ls3 + 3) *= epsilon;
		*(Ls3 + 5) += 1;
	}
}
/*
static void L_unknown_AA_SIMD(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	typedef double __attribute__((vector_size(160))) v20d;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lr1 = *(Ls1 + 1);
	const double Ln1 = *(Ls1 + 2);
	const double Ld1 = *(Ls1 + 3);
	const double Lc1 = *(Ls1 + 4);
	const double Lq1 = *(Ls1 + 5);
	const double Le1 = *(Ls1 + 6);
	const double Lg1 = *(Ls1 + 7);
	const double Lh1 = *(Ls1 + 8);
	const double Li1 = *(Ls1 + 9);
	const double Ll1 = *(Ls1 + 10);
	const double Lk1 = *(Ls1 + 11);
	const double Lm1 = *(Ls1 + 12);
	const double Lf1 = *(Ls1 + 13);
	const double Lp1 = *(Ls1 + 14);
	const double Lz1 = *(Ls1 + 15);
	const double Lt1 = *(Ls1 + 16);
	const double Lw1 = *(Ls1 + 17);
	const double Ly1 = *(Ls1 + 18);
	const double Lv1 = *(Ls1 + 19);
	const double La2 = *(Ls2);
	const double Lr2 = *(Ls2 + 1);
	const double Ln2 = *(Ls2 + 2);
	const double Ld2 = *(Ls2 + 3);
	const double Lc2 = *(Ls2 + 4);
	const double Lq2 = *(Ls2 + 5);
	const double Le2 = *(Ls2 + 6);
	const double Lg2 = *(Ls2 + 7);
	const double Lh2 = *(Ls2 + 8);
	const double Li2 = *(Ls2 + 9);
	const double Ll2 = *(Ls2 + 10);
	const double Lk2 = *(Ls2 + 11);
	const double Lm2 = *(Ls2 + 12);
	const double Lf2 = *(Ls2 + 13);
	const double Lp2 = *(Ls2 + 14);
	const double Lz2 = *(Ls2 + 15);
	const double Lt2 = *(Ls2 + 16);
	const double Lw2 = *(Ls2 + 17);
	const double Ly2 = *(Ls2 + 18);
	const double Lv2 = *(Ls2 + 19);
	
	unsigned int flag = 0;
	
	if (root == 0 &&
		(La1 != 0 ||
		Lr1 != 0 ||
		Ln1 != 0 ||
		Ld1 != 0 ||
		Lc1 != 0 ||
		Lq1 != 0 ||
		Le1 != 0 ||
		Lg1 != 0 ||
		Lh1 != 0 ||
		Li1 != 0 ||
		Ll1 != 0 ||
		Lk1 != 0 ||
		Lm1 != 0 ||
		Lf1 != 0 ||
		Lp1 != 0 ||
		Lz1 != 0 ||
		Lt1 != 0 ||
		Lw1 != 0 ||
		Ly1 != 0 ||
		Lv1 != 0)) {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0) {
			// neither branch can be disregarded
			
			v20d *L1 = (v20d*)(Ls3);
			v20d *vP11 = (v20d*)(P1);
			v20d *vP12 = (v20d*)(P1 + 21);
			v20d *vP13 = (v20d*)(P1 + 42);
			v20d *vP14 = (v20d*)(P1 + 63);
			v20d *vP15 = (v20d*)(P1 + 84);
			v20d *vP16 = (v20d*)(P1 + 105);
			v20d *vP17 = (v20d*)(P1 + 126);
			v20d *vP18 = (v20d*)(P1 + 147);
			v20d *vP19 = (v20d*)(P1 + 168);
			v20d *vP110 = (v20d*)(P1 + 189);
			v20d *vP111 = (v20d*)(P1 + 210);
			v20d *vP112 = (v20d*)(P1 + 231);
			v20d *vP113 = (v20d*)(P1 + 252);
			v20d *vP114 = (v20d*)(P1 + 273);
			v20d *vP115 = (v20d*)(P1 + 294);
			v20d *vP116 = (v20d*)(P1 + 315);
			v20d *vP117 = (v20d*)(P1 + 336);
			v20d *vP118 = (v20d*)(P1 + 357);
			v20d *vP119 = (v20d*)(P1 + 378);
			v20d *vP120 = (v20d*)(P1 + 399);
			// broadcast likelihoods
			const v20d LA1 = {La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1};
			const v20d LR1 = {Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1};
			const v20d LN1 = {Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1};
			const v20d LD1 = {Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1};
			const v20d LC1 = {Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1};
			const v20d LQ1 = {Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1};
			const v20d LE1 = {Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1};
			const v20d LG1 = {Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1};
			const v20d LH1 = {Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1};
			const v20d LI1 = {Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1};
			const v20d LL1 = {Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1};
			const v20d LK1 = {Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1};
			const v20d LM1 = {Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1};
			const v20d LF1 = {Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1};
			const v20d LP1 = {Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1};
			const v20d LZ1 = {Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1};
			const v20d LT1 = {Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1};
			const v20d LW1 = {Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1};
			const v20d LY1 = {Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1};
			const v20d LV1 = {Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LR1;
			*L1 += *vP13*LN1;
			*L1 += *vP14*LD1;
			*L1 += *vP15*LC1;
			*L1 += *vP16*LQ1;
			*L1 += *vP17*LE1;
			*L1 += *vP18*LG1;
			*L1 += *vP19*LH1;
			*L1 += *vP110*LI1;
			*L1 += *vP111*LL1;
			*L1 += *vP112*LK1;
			*L1 += *vP113*LM1;
			*L1 += *vP114*LF1;
			*L1 += *vP115*LP1;
			*L1 += *vP116*LZ1;
			*L1 += *vP117*LT1;
			*L1 += *vP118*LW1;
			*L1 += *vP119*LY1;
			*L1 += *vP120*LV1;
			
			v20d L2;
			v20d *vP21 = (v20d*)(P2);
			v20d *vP22 = (v20d*)(P2 + 21);
			v20d *vP23 = (v20d*)(P2 + 42);
			v20d *vP24 = (v20d*)(P2 + 63);
			v20d *vP25 = (v20d*)(P2 + 84);
			v20d *vP26 = (v20d*)(P2 + 105);
			v20d *vP27 = (v20d*)(P2 + 126);
			v20d *vP28 = (v20d*)(P2 + 147);
			v20d *vP29 = (v20d*)(P2 + 168);
			v20d *vP210 = (v20d*)(P2 + 189);
			v20d *vP211 = (v20d*)(P2 + 210);
			v20d *vP212 = (v20d*)(P2 + 231);
			v20d *vP213 = (v20d*)(P2 + 252);
			v20d *vP214 = (v20d*)(P2 + 273);
			v20d *vP215 = (v20d*)(P2 + 294);
			v20d *vP216 = (v20d*)(P2 + 315);
			v20d *vP217 = (v20d*)(P2 + 336);
			v20d *vP218 = (v20d*)(P2 + 357);
			v20d *vP219 = (v20d*)(P2 + 378);
			v20d *vP220 = (v20d*)(P2 + 399);
			// broadcast likelihoods
			const v20d LA2 = {La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2};
			const v20d LR2 = {Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2};
			const v20d LN2 = {Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2};
			const v20d LD2 = {Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2};
			const v20d LC2 = {Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2};
			const v20d LQ2 = {Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2};
			const v20d LE2 = {Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2};
			const v20d LG2 = {Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2};
			const v20d LH2 = {Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2};
			const v20d LI2 = {Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2};
			const v20d LL2 = {Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2};
			const v20d LK2 = {Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2};
			const v20d LM2 = {Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2};
			const v20d LF2 = {Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2};
			const v20d LP2 = {Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2};
			const v20d LZ2 = {Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2};
			const v20d LT2 = {Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2};
			const v20d LW2 = {Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2};
			const v20d LY2 = {Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2};
			const v20d LV2 = {Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2};
			L2 = *vP21*LA2;
			L2 += *vP22*LR2;
			L2 += *vP23*LN2;
			L2 += *vP24*LD2;
			L2 += *vP25*LC2;
			L2 += *vP26*LQ2;
			L2 += *vP27*LE2;
			L2 += *vP28*LG2;
			L2 += *vP29*LH2;
			L2 += *vP210*LI2;
			L2 += *vP211*LL2;
			L2 += *vP212*LK2;
			L2 += *vP213*LM2;
			L2 += *vP214*LF2;
			L2 += *vP215*LP2;
			L2 += *vP216*LZ2;
			L2 += *vP217*LT2;
			L2 += *vP218*LW2;
			L2 += *vP219*LY2;
			L2 += *vP220*LV2;
			
			*L1 *= L2;
			for (int i = 0; i < 20; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
		} else {
			// second branch can be disregarded
			
			v20d *L1 = (v20d*)(Ls3);
			v20d *vP11 = (v20d*)(P1);
			v20d *vP12 = (v20d*)(P1 + 21);
			v20d *vP13 = (v20d*)(P1 + 42);
			v20d *vP14 = (v20d*)(P1 + 63);
			v20d *vP15 = (v20d*)(P1 + 84);
			v20d *vP16 = (v20d*)(P1 + 105);
			v20d *vP17 = (v20d*)(P1 + 126);
			v20d *vP18 = (v20d*)(P1 + 147);
			v20d *vP19 = (v20d*)(P1 + 168);
			v20d *vP110 = (v20d*)(P1 + 189);
			v20d *vP111 = (v20d*)(P1 + 210);
			v20d *vP112 = (v20d*)(P1 + 231);
			v20d *vP113 = (v20d*)(P1 + 252);
			v20d *vP114 = (v20d*)(P1 + 273);
			v20d *vP115 = (v20d*)(P1 + 294);
			v20d *vP116 = (v20d*)(P1 + 315);
			v20d *vP117 = (v20d*)(P1 + 336);
			v20d *vP118 = (v20d*)(P1 + 357);
			v20d *vP119 = (v20d*)(P1 + 378);
			v20d *vP120 = (v20d*)(P1 + 399);
			// broadcast likelihoods
			const v20d LA1 = {La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1};
			const v20d LR1 = {Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1};
			const v20d LN1 = {Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1};
			const v20d LD1 = {Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1};
			const v20d LC1 = {Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1};
			const v20d LQ1 = {Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1};
			const v20d LE1 = {Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1};
			const v20d LG1 = {Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1};
			const v20d LH1 = {Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1};
			const v20d LI1 = {Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1};
			const v20d LL1 = {Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1};
			const v20d LK1 = {Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1};
			const v20d LM1 = {Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1};
			const v20d LF1 = {Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1};
			const v20d LP1 = {Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1};
			const v20d LZ1 = {Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1};
			const v20d LT1 = {Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1};
			const v20d LW1 = {Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1};
			const v20d LY1 = {Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1};
			const v20d LV1 = {Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LR1;
			*L1 += *vP13*LN1;
			*L1 += *vP14*LD1;
			*L1 += *vP15*LC1;
			*L1 += *vP16*LQ1;
			*L1 += *vP17*LE1;
			*L1 += *vP18*LG1;
			*L1 += *vP19*LH1;
			*L1 += *vP110*LI1;
			*L1 += *vP111*LL1;
			*L1 += *vP112*LK1;
			*L1 += *vP113*LM1;
			*L1 += *vP114*LF1;
			*L1 += *vP115*LP1;
			*L1 += *vP116*LZ1;
			*L1 += *vP117*LT1;
			*L1 += *vP118*LW1;
			*L1 += *vP119*LY1;
			*L1 += *vP120*LV1;
			
			for (int i = 0; i < 20; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 21) = *(Ls1 + 21);
		}
	} else {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0) {
			// first branch can be disregarded
			
			v20d *L2 = (v20d*)(Ls3);
			v20d *vP21 = (v20d*)(P2);
			v20d *vP22 = (v20d*)(P2 + 21);
			v20d *vP23 = (v20d*)(P2 + 42);
			v20d *vP24 = (v20d*)(P2 + 63);
			v20d *vP25 = (v20d*)(P2 + 84);
			v20d *vP26 = (v20d*)(P2 + 105);
			v20d *vP27 = (v20d*)(P2 + 126);
			v20d *vP28 = (v20d*)(P2 + 147);
			v20d *vP29 = (v20d*)(P2 + 168);
			v20d *vP210 = (v20d*)(P2 + 189);
			v20d *vP211 = (v20d*)(P2 + 210);
			v20d *vP212 = (v20d*)(P2 + 231);
			v20d *vP213 = (v20d*)(P2 + 252);
			v20d *vP214 = (v20d*)(P2 + 273);
			v20d *vP215 = (v20d*)(P2 + 294);
			v20d *vP216 = (v20d*)(P2 + 315);
			v20d *vP217 = (v20d*)(P2 + 336);
			v20d *vP218 = (v20d*)(P2 + 357);
			v20d *vP219 = (v20d*)(P2 + 378);
			v20d *vP220 = (v20d*)(P2 + 399);
			// broadcast likelihoods
			const v20d LA2 = {La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2};
			const v20d LR2 = {Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2};
			const v20d LN2 = {Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2};
			const v20d LD2 = {Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2};
			const v20d LC2 = {Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2};
			const v20d LQ2 = {Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2};
			const v20d LE2 = {Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2};
			const v20d LG2 = {Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2};
			const v20d LH2 = {Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2};
			const v20d LI2 = {Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2};
			const v20d LL2 = {Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2};
			const v20d LK2 = {Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2};
			const v20d LM2 = {Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2};
			const v20d LF2 = {Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2};
			const v20d LP2 = {Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2};
			const v20d LZ2 = {Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2};
			const v20d LT2 = {Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2};
			const v20d LW2 = {Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2};
			const v20d LY2 = {Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2};
			const v20d LV2 = {Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2};
			*L2 = *vP21*LA2;
			*L2 += *vP22*LR2;
			*L2 += *vP23*LN2;
			*L2 += *vP24*LD2;
			*L2 += *vP25*LC2;
			*L2 += *vP26*LQ2;
			*L2 += *vP27*LE2;
			*L2 += *vP28*LG2;
			*L2 += *vP29*LH2;
			*L2 += *vP210*LI2;
			*L2 += *vP211*LL2;
			*L2 += *vP212*LK2;
			*L2 += *vP213*LM2;
			*L2 += *vP214*LF2;
			*L2 += *vP215*LP2;
			*L2 += *vP216*LZ2;
			*L2 += *vP217*LT2;
			*L2 += *vP218*LW2;
			*L2 += *vP219*LY2;
			*L2 += *vP220*LV2;
			
			if (root &&
				(La1 != 0 ||
				Lr1 != 0 ||
				Ln1 != 0 ||
				Ld1 != 0 ||
				Lc1 != 0 ||
				Lq1 != 0 ||
				Le1 != 0 ||
				Lg1 != 0 ||
				Lh1 != 0 ||
				Li1 != 0 ||
				Ll1 != 0 ||
				Lk1 != 0 ||
				Lm1 != 0 ||
				Lf1 != 0 ||
				Lp1 != 0 ||
				Lz1 != 0 ||
				Lt1 != 0 ||
				Lw1 != 0 ||
				Ly1 != 0 ||
				Lv1 != 0)) {
				const v20d L1 = {La1, Lr1, Ln1, Ld1, Lc1, Lq1, Le1, Lg1, Lh1, Li1, Ll1, Lk1, Lm1, Lf1, Lp1, Lz1, Lt1, Lw1, Ly1, Lv1};
				*L2 *= L1;
				*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			} else {
				*(Ls3 + 21) = *(Ls2 + 21);
			}
			
			for (int i = 0; i < 20; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lr1;
			*(Ls3 + 2) = Ln1;
			*(Ls3 + 3) = Ld1;
			*(Ls3 + 4) = Lc1;
			*(Ls3 + 5) = Lq1;
			*(Ls3 + 6) = Le1;
			*(Ls3 + 7) = Lg1;
			*(Ls3 + 8) = Lh1;
			*(Ls3 + 9) = Li1;
			*(Ls3 + 10) = Ll1;
			*(Ls3 + 11) = Lk1;
			*(Ls3 + 12) = Lm1;
			*(Ls3 + 13) = Lf1;
			*(Ls3 + 14) = Lp1;
			*(Ls3 + 15) = Lz1;
			*(Ls3 + 16) = Lt1;
			*(Ls3 + 17) = Lw1;
			*(Ls3 + 18) = Ly1;
			*(Ls3 + 19) = Lv1;
			*(Ls3 + 21) = *(Ls1 + 21);
		}
	}
	
	*(Ls3 + 20) = 0;
	if (flag) {
		*(Ls3) *= epsilon;
		*(Ls3 + 1) *= epsilon;
		*(Ls3 + 2) *= epsilon;
		*(Ls3 + 3) *= epsilon;
		*(Ls3 + 4) *= epsilon;
		*(Ls3 + 5) *= epsilon;
		*(Ls3 + 6) *= epsilon;
		*(Ls3 + 7) *= epsilon;
		*(Ls3 + 8) *= epsilon;
		*(Ls3 + 9) *= epsilon;
		*(Ls3 + 10) *= epsilon;
		*(Ls3 + 11) *= epsilon;
		*(Ls3 + 12) *= epsilon;
		*(Ls3 + 13) *= epsilon;
		*(Ls3 + 14) *= epsilon;
		*(Ls3 + 15) *= epsilon;
		*(Ls3 + 16) *= epsilon;
		*(Ls3 + 17) *= epsilon;
		*(Ls3 + 18) *= epsilon;
		*(Ls3 + 19) *= epsilon;
		*(Ls3 + 21) += 1;
	}
}
*/
static void L_unknown_5_SIMD(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	typedef double __attribute__((vector_size(32))) v4d;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lc1 = *(Ls1 + 1);
	const double Lg1 = *(Ls1 + 2);
	const double Lt1 = *(Ls1 + 3);
	const double Li1 = *(Ls1 + 4);
	const double La2 = *(Ls2);
	const double Lc2 = *(Ls2 + 1);
	const double Lg2 = *(Ls2 + 2);
	const double Lt2 = *(Ls2 + 3);
	const double Li2 = *(Ls2 + 4);
	
	unsigned int flag = 0;
	
	if (root == 0 &&
		(La1 != 0 ||
		Lc1 != 0 ||
		Lg1 != 0 ||
		Lt1 != 0 ||
		Li1 != 0)) {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0 ||
			Li2 != 0) {
			// neither branch can be disregarded
			
			v4d *L1 = (v4d*)(Ls3);
			v4d *vP11 = (v4d*)(P1);
			v4d *vP12 = (v4d*)(P1 + 5);
			v4d *vP13 = (v4d*)(P1 + 10);
			v4d *vP14 = (v4d*)(P1 + 15);
			v4d *vP15 = (v4d*)(P1 + 20);
			// broadcast likelihoods
			const v4d LA1 = {La1, La1, La1, La1};
			const v4d LC1 = {Lc1, Lc1, Lc1, Lc1};
			const v4d LG1 = {Lg1, Lg1, Lg1, Lg1};
			const v4d LT1 = {Lt1, Lt1, Lt1, Lt1};
			const v4d LI1 = {Li1, Li1, Li1, Li1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LC1;
			*L1 += *vP13*LG1;
			*L1 += *vP14*LT1;
			*L1 += *vP15*LI1;
			
			v4d L2;
			v4d *vP21 = (v4d*)(P2);
			v4d *vP22 = (v4d*)(P2 + 5);
			v4d *vP23 = (v4d*)(P2 + 10);
			v4d *vP24 = (v4d*)(P2 + 15);
			v4d *vP25 = (v4d*)(P2 + 20);
			// broadcast likelihoods
			const v4d LA2 = {La2, La2, La2, La2};
			const v4d LC2 = {Lc2, Lc2, Lc2, Lc2};
			const v4d LG2 = {Lg2, Lg2, Lg2, Lg2};
			const v4d LT2 = {Lt2, Lt2, Lt2, Lt2};
			const v4d LI2 = {Li2, Li2, Li2, Li2};
			L2 = *vP21*LA2;
			L2 += *vP22*LC2;
			L2 += *vP23*LG2;
			L2 += *vP24*LT2;
			L2 += *vP25*LI2;
			
			*L1 *= L2;
			
			// L(-)
			double L3, L4;
			L3 = *(P1 + 4)*(La1); // L-a
			L3 += *(P1 + 9)*(Lc1); // L-c
			L3 += *(P1 + 14)*(Lg1); // L-g
			L3 += *(P1 + 19)*(Lt1); // L-t
			L3 += *(P1 + 24)*(Li1); // L--
			L4 = *(P2 + 4)*(La2); // L-a
			L4 += *(P2 + 9)*(Lc2); // L-c
			L4 += *(P2 + 14)*(Lg2); // L-g
			L4 += *(P2 + 19)*(Lt2); // L-t
			L4 += *(P2 + 24)*(Li2); // L--
			*(Ls3 + 4) = L3*L4;
			
			for (int i = 0; i < 5; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
		} else {
			// second branch can be disregarded
			
			v4d *L1 = (v4d*)(Ls3);
			v4d *vP11 = (v4d*)(P1);
			v4d *vP12 = (v4d*)(P1 + 5);
			v4d *vP13 = (v4d*)(P1 + 10);
			v4d *vP14 = (v4d*)(P1 + 15);
			v4d *vP15 = (v4d*)(P1 + 20);
			// broadcast likelihoods
			const v4d LA1 = {La1, La1, La1, La1};
			const v4d LC1 = {Lc1, Lc1, Lc1, Lc1};
			const v4d LG1 = {Lg1, Lg1, Lg1, Lg1};
			const v4d LT1 = {Lt1, Lt1, Lt1, Lt1};
			const v4d LI1 = {Li1, Li1, Li1, Li1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LC1;
			*L1 += *vP13*LG1;
			*L1 += *vP14*LT1;
			*L1 += *vP15*LI1;
			
			// L(-)
			double L3;
			L3 = *(P1 + 4)*(La1); // L-a
			L3 += *(P1 + 9)*(Lc1); // L-c
			L3 += *(P1 + 14)*(Lg1); // L-g
			L3 += *(P1 + 19)*(Lt1); // L-t
			L3 += *(P1 + 24)*(Li1); // L--
			*(Ls3 + 4) = L3;
			
			for (int i = 0; i < 5; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 5) = *(Ls1 + 5);
		}
	} else {
		if (La2 != 0 ||
			Lc2 != 0 ||
			Lg2 != 0 ||
			Lt2 != 0 ||
			Li2 != 0) {
			// first branch can be disregarded
			
			v4d *L2 = (v4d*)(Ls3);
			v4d *vP21 = (v4d*)(P2);
			v4d *vP22 = (v4d*)(P2 + 5);
			v4d *vP23 = (v4d*)(P2 + 10);
			v4d *vP24 = (v4d*)(P2 + 15);
			v4d *vP25 = (v4d*)(P2 + 20);
			// broadcast likelihoods
			const v4d LA2 = {La2, La2, La2, La2};
			const v4d LC2 = {Lc2, Lc2, Lc2, Lc2};
			const v4d LG2 = {Lg2, Lg2, Lg2, Lg2};
			const v4d LT2 = {Lt2, Lt2, Lt2, Lt2};
			const v4d LI2 = {Li2, Li2, Li2, Li2};
			*L2 = *vP21*LA2;
			*L2 += *vP22*LC2;
			*L2 += *vP23*LG2;
			*L2 += *vP24*LT2;
			*L2 += *vP25*LI2;
			
			// L(-)
			double L4;
			L4 = *(P2 + 4)*(La2); // L-a
			L4 += *(P2 + 9)*(Lc2); // L-c
			L4 += *(P2 + 14)*(Lg2); // L-g
			L4 += *(P2 + 19)*(Lt2); // L-t
			L4 += *(P2 + 24)*(Li2); // L--
			*(Ls3 + 4) = L4;
			
			if (root &&
				(La1 != 0 ||
				Lc1 != 0 ||
				Lg1 != 0 ||
				Lt1 != 0 ||
				Li1 != 0)) {
				const v4d L1 = {La1, Lc1, Lg1, Lt1};
				*L2 *= L1;
				*(Ls3 + 4) *= Li1;
				*(Ls3 + 5) = *(Ls1 + 5) + *(Ls2 + 5);
			} else {
				*(Ls3 + 5) = *(Ls2 + 5);
			}
			
			for (int i = 0; i < 5; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lc1;
			*(Ls3 + 2) = Lg1;
			*(Ls3 + 3) = Lt1;
			*(Ls3 + 4) = Li1;
			*(Ls3 + 5) = *(Ls1 + 5);
		}
	}
	
	if (flag) {
		*(Ls3) *= epsilon;
		*(Ls3 + 1) *= epsilon;
		*(Ls3 + 2) *= epsilon;
		*(Ls3 + 3) *= epsilon;
		*(Ls3 + 4) *= epsilon;
		*(Ls3 + 5) += 1;
	}
}
/*
static void L_unknown_AA_5_SIMD(double *__restrict Ls, const int i3, const int i1, const int i2, const double *P1, const double *P2, const double epsilon, const double inv_epsilon, const int root)
{
	typedef double __attribute__((vector_size(160))) v20d;
	
	double *Ls1 = Ls + i1;
	double *Ls2 = Ls + i2;
	double *Ls3 = Ls + i3;
	
	const double La1 = *(Ls1);
	const double Lr1 = *(Ls1 + 1);
	const double Ln1 = *(Ls1 + 2);
	const double Ld1 = *(Ls1 + 3);
	const double Lc1 = *(Ls1 + 4);
	const double Lq1 = *(Ls1 + 5);
	const double Le1 = *(Ls1 + 6);
	const double Lg1 = *(Ls1 + 7);
	const double Lh1 = *(Ls1 + 8);
	const double Li1 = *(Ls1 + 9);
	const double Ll1 = *(Ls1 + 10);
	const double Lk1 = *(Ls1 + 11);
	const double Lm1 = *(Ls1 + 12);
	const double Lf1 = *(Ls1 + 13);
	const double Lp1 = *(Ls1 + 14);
	const double Lz1 = *(Ls1 + 15);
	const double Lt1 = *(Ls1 + 16);
	const double Lw1 = *(Ls1 + 17);
	const double Ly1 = *(Ls1 + 18);
	const double Lv1 = *(Ls1 + 19);
	const double Lb1 = *(Ls1 + 20);
	const double La2 = *(Ls2);
	const double Lr2 = *(Ls2 + 1);
	const double Ln2 = *(Ls2 + 2);
	const double Ld2 = *(Ls2 + 3);
	const double Lc2 = *(Ls2 + 4);
	const double Lq2 = *(Ls2 + 5);
	const double Le2 = *(Ls2 + 6);
	const double Lg2 = *(Ls2 + 7);
	const double Lh2 = *(Ls2 + 8);
	const double Li2 = *(Ls2 + 9);
	const double Ll2 = *(Ls2 + 10);
	const double Lk2 = *(Ls2 + 11);
	const double Lm2 = *(Ls2 + 12);
	const double Lf2 = *(Ls2 + 13);
	const double Lp2 = *(Ls2 + 14);
	const double Lz2 = *(Ls2 + 15);
	const double Lt2 = *(Ls2 + 16);
	const double Lw2 = *(Ls2 + 17);
	const double Ly2 = *(Ls2 + 18);
	const double Lv2 = *(Ls2 + 19);
	const double Lb2 = *(Ls2 + 20);
	
	unsigned int flag = 0;
	
	if (root == 0 &&
		(La1 != 0 ||
		Lr1 != 0 ||
		Ln1 != 0 ||
		Ld1 != 0 ||
		Lc1 != 0 ||
		Lq1 != 0 ||
		Le1 != 0 ||
		Lg1 != 0 ||
		Lh1 != 0 ||
		Li1 != 0 ||
		Ll1 != 0 ||
		Lk1 != 0 ||
		Lm1 != 0 ||
		Lf1 != 0 ||
		Lp1 != 0 ||
		Lz1 != 0 ||
		Lt1 != 0 ||
		Lw1 != 0 ||
		Ly1 != 0 ||
		Lv1 != 0 ||
		Lb1 != 0)) {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0 ||
			Lb2 != 0) {
			// neither branch can be disregarded
			
			v20d *L1 = (v20d*)(Ls3);
			v20d *vP11 = (v20d*)(P1);
			v20d *vP12 = (v20d*)(P1 + 21);
			v20d *vP13 = (v20d*)(P1 + 42);
			v20d *vP14 = (v20d*)(P1 + 63);
			v20d *vP15 = (v20d*)(P1 + 84);
			v20d *vP16 = (v20d*)(P1 + 105);
			v20d *vP17 = (v20d*)(P1 + 126);
			v20d *vP18 = (v20d*)(P1 + 147);
			v20d *vP19 = (v20d*)(P1 + 168);
			v20d *vP110 = (v20d*)(P1 + 189);
			v20d *vP111 = (v20d*)(P1 + 210);
			v20d *vP112 = (v20d*)(P1 + 231);
			v20d *vP113 = (v20d*)(P1 + 252);
			v20d *vP114 = (v20d*)(P1 + 273);
			v20d *vP115 = (v20d*)(P1 + 294);
			v20d *vP116 = (v20d*)(P1 + 315);
			v20d *vP117 = (v20d*)(P1 + 336);
			v20d *vP118 = (v20d*)(P1 + 357);
			v20d *vP119 = (v20d*)(P1 + 378);
			v20d *vP120 = (v20d*)(P1 + 399);
			v20d *vP121 = (v20d*)(P1 + 420);
			// broadcast likelihoods
			const v20d LA1 = {La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1};
			const v20d LR1 = {Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1};
			const v20d LN1 = {Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1};
			const v20d LD1 = {Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1};
			const v20d LC1 = {Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1};
			const v20d LQ1 = {Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1};
			const v20d LE1 = {Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1};
			const v20d LG1 = {Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1};
			const v20d LH1 = {Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1};
			const v20d LI1 = {Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1};
			const v20d LL1 = {Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1};
			const v20d LK1 = {Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1};
			const v20d LM1 = {Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1};
			const v20d LF1 = {Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1};
			const v20d LP1 = {Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1};
			const v20d LZ1 = {Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1};
			const v20d LT1 = {Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1};
			const v20d LW1 = {Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1};
			const v20d LY1 = {Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1};
			const v20d LV1 = {Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1};
			const v20d LB1 = {Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LR1;
			*L1 += *vP13*LN1;
			*L1 += *vP14*LD1;
			*L1 += *vP15*LC1;
			*L1 += *vP16*LQ1;
			*L1 += *vP17*LE1;
			*L1 += *vP18*LG1;
			*L1 += *vP19*LH1;
			*L1 += *vP110*LI1;
			*L1 += *vP111*LL1;
			*L1 += *vP112*LK1;
			*L1 += *vP113*LM1;
			*L1 += *vP114*LF1;
			*L1 += *vP115*LP1;
			*L1 += *vP116*LZ1;
			*L1 += *vP117*LT1;
			*L1 += *vP118*LW1;
			*L1 += *vP119*LY1;
			*L1 += *vP120*LV1;
			*L1 += *vP121*LB1;
			
			v20d L2;
			v20d *vP21 = (v20d*)(P2);
			v20d *vP22 = (v20d*)(P2 + 21);
			v20d *vP23 = (v20d*)(P2 + 42);
			v20d *vP24 = (v20d*)(P2 + 63);
			v20d *vP25 = (v20d*)(P2 + 84);
			v20d *vP26 = (v20d*)(P2 + 105);
			v20d *vP27 = (v20d*)(P2 + 126);
			v20d *vP28 = (v20d*)(P2 + 147);
			v20d *vP29 = (v20d*)(P2 + 168);
			v20d *vP210 = (v20d*)(P2 + 189);
			v20d *vP211 = (v20d*)(P2 + 210);
			v20d *vP212 = (v20d*)(P2 + 231);
			v20d *vP213 = (v20d*)(P2 + 252);
			v20d *vP214 = (v20d*)(P2 + 273);
			v20d *vP215 = (v20d*)(P2 + 294);
			v20d *vP216 = (v20d*)(P2 + 315);
			v20d *vP217 = (v20d*)(P2 + 336);
			v20d *vP218 = (v20d*)(P2 + 357);
			v20d *vP219 = (v20d*)(P2 + 378);
			v20d *vP220 = (v20d*)(P2 + 399);
			v20d *vP221 = (v20d*)(P2 + 420);
			// broadcast likelihoods
			const v20d LA2 = {La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2};
			const v20d LR2 = {Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2};
			const v20d LN2 = {Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2};
			const v20d LD2 = {Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2};
			const v20d LC2 = {Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2};
			const v20d LQ2 = {Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2};
			const v20d LE2 = {Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2};
			const v20d LG2 = {Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2};
			const v20d LH2 = {Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2};
			const v20d LI2 = {Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2};
			const v20d LL2 = {Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2};
			const v20d LK2 = {Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2};
			const v20d LM2 = {Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2};
			const v20d LF2 = {Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2};
			const v20d LP2 = {Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2};
			const v20d LZ2 = {Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2};
			const v20d LT2 = {Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2};
			const v20d LW2 = {Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2};
			const v20d LY2 = {Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2};
			const v20d LV2 = {Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2};
			const v20d LB2 = {Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2};
			L2 = *vP21*LA2;
			L2 += *vP22*LR2;
			L2 += *vP23*LN2;
			L2 += *vP24*LD2;
			L2 += *vP25*LC2;
			L2 += *vP26*LQ2;
			L2 += *vP27*LE2;
			L2 += *vP28*LG2;
			L2 += *vP29*LH2;
			L2 += *vP210*LI2;
			L2 += *vP211*LL2;
			L2 += *vP212*LK2;
			L2 += *vP213*LM2;
			L2 += *vP214*LF2;
			L2 += *vP215*LP2;
			L2 += *vP216*LZ2;
			L2 += *vP217*LT2;
			L2 += *vP218*LW2;
			L2 += *vP219*LY2;
			L2 += *vP220*LV2;
			L2 += *vP221*LB2;
			
			*L1 *= L2;
			
			// L(Indels)
			double L3, L4;
			L3 = *(P1 + 20)*(La1); // Lba
			L3 += *(P1 + 41)*(Lr1); // Lbr
			L3 += *(P1 + 62)*(Ln1); // Lbn
			L3 += *(P1 + 83)*(Ld1); // Lbd
			L3 += *(P1 + 104)*(Lc1); // Lbc
			L3 += *(P1 + 125)*(Lq1); // Lbq
			L3 += *(P1 + 146)*(Le1); // Lbe
			L3 += *(P1 + 167)*(Lg1); // Lbg
			L3 += *(P1 + 188)*(Lh1); // Lbh
			L3 += *(P1 + 209)*(Li1); // Lbi
			L3 += *(P1 + 230)*(Ll1); // Lbl
			L3 += *(P1 + 251)*(Lk1); // Lbk
			L3 += *(P1 + 272)*(Lm1); // Lbm
			L3 += *(P1 + 293)*(Lf1); // Lbf
			L3 += *(P1 + 314)*(Lp1); // Lbp
			L3 += *(P1 + 335)*(Lz1); // Lbs
			L3 += *(P1 + 356)*(Lt1); // Lbt
			L3 += *(P1 + 377)*(Lw1); // Lbw
			L3 += *(P1 + 398)*(Ly1); // Lby
			L3 += *(P1 + 419)*(Lv1); // Lbv
			L3 += *(P1 + 440)*(Lb1); // Lbb
			L4 = *(P2 + 20)*(La2); // Lba
			L4 += *(P2 + 41)*(Lr2); // Lbr
			L4 += *(P2 + 62)*(Ln2); // Lbn
			L4 += *(P2 + 83)*(Ld2); // Lbd
			L4 += *(P2 + 104)*(Lc2); // Lbc
			L4 += *(P2 + 125)*(Lq2); // Lbq
			L4 += *(P2 + 146)*(Le2); // Lbe
			L4 += *(P2 + 167)*(Lg2); // Lbg
			L4 += *(P2 + 188)*(Lh2); // Lbh
			L4 += *(P2 + 209)*(Li2); // Lbi
			L4 += *(P2 + 230)*(Ll2); // Lbl
			L4 += *(P2 + 251)*(Lk2); // Lbk
			L4 += *(P2 + 272)*(Lm2); // Lbm
			L4 += *(P2 + 293)*(Lf2); // Lbf
			L4 += *(P2 + 314)*(Lp2); // Lbp
			L4 += *(P2 + 335)*(Lz2); // Lbs
			L4 += *(P2 + 356)*(Lt2); // Lbt
			L4 += *(P2 + 377)*(Lw2); // Lbw
			L4 += *(P2 + 398)*(Ly2); // Lby
			L4 += *(P2 + 419)*(Lv2); // Lbv
			L4 += *(P2 + 440)*(Lb2); // Lbb
			*(Ls3 + 20) = L3*L4;
			
			for (int i = 0; i < 21; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
		} else {
			// second branch can be disregarded
			
			v20d *L1 = (v20d*)(Ls3);
			v20d *vP11 = (v20d*)(P1);
			v20d *vP12 = (v20d*)(P1 + 21);
			v20d *vP13 = (v20d*)(P1 + 42);
			v20d *vP14 = (v20d*)(P1 + 63);
			v20d *vP15 = (v20d*)(P1 + 84);
			v20d *vP16 = (v20d*)(P1 + 105);
			v20d *vP17 = (v20d*)(P1 + 126);
			v20d *vP18 = (v20d*)(P1 + 147);
			v20d *vP19 = (v20d*)(P1 + 168);
			v20d *vP110 = (v20d*)(P1 + 189);
			v20d *vP111 = (v20d*)(P1 + 210);
			v20d *vP112 = (v20d*)(P1 + 231);
			v20d *vP113 = (v20d*)(P1 + 252);
			v20d *vP114 = (v20d*)(P1 + 273);
			v20d *vP115 = (v20d*)(P1 + 294);
			v20d *vP116 = (v20d*)(P1 + 315);
			v20d *vP117 = (v20d*)(P1 + 336);
			v20d *vP118 = (v20d*)(P1 + 357);
			v20d *vP119 = (v20d*)(P1 + 378);
			v20d *vP120 = (v20d*)(P1 + 399);
			v20d *vP121 = (v20d*)(P1 + 420);
			// broadcast likelihoods
			const v20d LA1 = {La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1, La1};
			const v20d LR1 = {Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1, Lr1};
			const v20d LN1 = {Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1, Ln1};
			const v20d LD1 = {Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1, Ld1};
			const v20d LC1 = {Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1, Lc1};
			const v20d LQ1 = {Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1, Lq1};
			const v20d LE1 = {Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1, Le1};
			const v20d LG1 = {Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1, Lg1};
			const v20d LH1 = {Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1, Lh1};
			const v20d LI1 = {Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1, Li1};
			const v20d LL1 = {Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1, Ll1};
			const v20d LK1 = {Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1, Lk1};
			const v20d LM1 = {Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1, Lm1};
			const v20d LF1 = {Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1, Lf1};
			const v20d LP1 = {Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1, Lp1};
			const v20d LZ1 = {Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1, Lz1};
			const v20d LT1 = {Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1, Lt1};
			const v20d LW1 = {Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1, Lw1};
			const v20d LY1 = {Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1, Ly1};
			const v20d LV1 = {Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1, Lv1};
			const v20d LB1 = {Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1, Lb1};
			*L1 = *vP11*LA1;
			*L1 += *vP12*LR1;
			*L1 += *vP13*LN1;
			*L1 += *vP14*LD1;
			*L1 += *vP15*LC1;
			*L1 += *vP16*LQ1;
			*L1 += *vP17*LE1;
			*L1 += *vP18*LG1;
			*L1 += *vP19*LH1;
			*L1 += *vP110*LI1;
			*L1 += *vP111*LL1;
			*L1 += *vP112*LK1;
			*L1 += *vP113*LM1;
			*L1 += *vP114*LF1;
			*L1 += *vP115*LP1;
			*L1 += *vP116*LZ1;
			*L1 += *vP117*LT1;
			*L1 += *vP118*LW1;
			*L1 += *vP119*LY1;
			*L1 += *vP120*LV1;
			*L1 += *vP121*LB1;
			
			// L(Indels)
			double L3;
			L3 = *(P1 + 20)*(La1); // Lba
			L3 += *(P1 + 41)*(Lr1); // Lbr
			L3 += *(P1 + 62)*(Ln1); // Lbn
			L3 += *(P1 + 83)*(Ld1); // Lbd
			L3 += *(P1 + 104)*(Lc1); // Lbc
			L3 += *(P1 + 125)*(Lq1); // Lbq
			L3 += *(P1 + 146)*(Le1); // Lbe
			L3 += *(P1 + 167)*(Lg1); // Lbg
			L3 += *(P1 + 188)*(Lh1); // Lbh
			L3 += *(P1 + 209)*(Li1); // Lbi
			L3 += *(P1 + 230)*(Ll1); // Lbl
			L3 += *(P1 + 251)*(Lk1); // Lbk
			L3 += *(P1 + 272)*(Lm1); // Lbm
			L3 += *(P1 + 293)*(Lf1); // Lbf
			L3 += *(P1 + 314)*(Lp1); // Lbp
			L3 += *(P1 + 335)*(Lz1); // Lbs
			L3 += *(P1 + 356)*(Lt1); // Lbt
			L3 += *(P1 + 377)*(Lw1); // Lbw
			L3 += *(P1 + 398)*(Ly1); // Lby
			L3 += *(P1 + 419)*(Lv1); // Lbv
			L3 += *(P1 + 440)*(Lb1); // Lbb
			*(Ls3 + 20) = L3;
			
			for (int i = 0; i < 21; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
			*(Ls3 + 21) = *(Ls1 + 21);
		}
	} else {
		if (La2 != 0 ||
			Lr2 != 0 ||
			Ln2 != 0 ||
			Ld2 != 0 ||
			Lc2 != 0 ||
			Lq2 != 0 ||
			Le2 != 0 ||
			Lg2 != 0 ||
			Lh2 != 0 ||
			Li2 != 0 ||
			Ll2 != 0 ||
			Lk2 != 0 ||
			Lm2 != 0 ||
			Lf2 != 0 ||
			Lp2 != 0 ||
			Lz2 != 0 ||
			Lt2 != 0 ||
			Lw2 != 0 ||
			Ly2 != 0 ||
			Lv2 != 0 ||
			Lb2 != 0) {
			// first branch can be disregarded
			
			v20d *L2 = (v20d*)(Ls3);
			v20d *vP21 = (v20d*)(P2);
			v20d *vP22 = (v20d*)(P2 + 21);
			v20d *vP23 = (v20d*)(P2 + 42);
			v20d *vP24 = (v20d*)(P2 + 63);
			v20d *vP25 = (v20d*)(P2 + 84);
			v20d *vP26 = (v20d*)(P2 + 105);
			v20d *vP27 = (v20d*)(P2 + 126);
			v20d *vP28 = (v20d*)(P2 + 147);
			v20d *vP29 = (v20d*)(P2 + 168);
			v20d *vP210 = (v20d*)(P2 + 189);
			v20d *vP211 = (v20d*)(P2 + 210);
			v20d *vP212 = (v20d*)(P2 + 231);
			v20d *vP213 = (v20d*)(P2 + 252);
			v20d *vP214 = (v20d*)(P2 + 273);
			v20d *vP215 = (v20d*)(P2 + 294);
			v20d *vP216 = (v20d*)(P2 + 315);
			v20d *vP217 = (v20d*)(P2 + 336);
			v20d *vP218 = (v20d*)(P2 + 357);
			v20d *vP219 = (v20d*)(P2 + 378);
			v20d *vP220 = (v20d*)(P2 + 399);
			v20d *vP221 = (v20d*)(P2 + 420);
			// broadcast likelihoods
			const v20d LA2 = {La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2, La2};
			const v20d LR2 = {Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2, Lr2};
			const v20d LN2 = {Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2, Ln2};
			const v20d LD2 = {Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2, Ld2};
			const v20d LC2 = {Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2, Lc2};
			const v20d LQ2 = {Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2, Lq2};
			const v20d LE2 = {Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2, Le2};
			const v20d LG2 = {Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2, Lg2};
			const v20d LH2 = {Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2, Lh2};
			const v20d LI2 = {Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2, Li2};
			const v20d LL2 = {Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2, Ll2};
			const v20d LK2 = {Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2, Lk2};
			const v20d LM2 = {Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2, Lm2};
			const v20d LF2 = {Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2, Lf2};
			const v20d LP2 = {Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2, Lp2};
			const v20d LZ2 = {Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2, Lz2};
			const v20d LT2 = {Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2, Lt2};
			const v20d LW2 = {Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2, Lw2};
			const v20d LY2 = {Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2, Ly2};
			const v20d LV2 = {Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2, Lv2};
			const v20d LB2 = {Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2, Lb2};
			*L2 = *vP21*LA2;
			*L2 += *vP22*LR2;
			*L2 += *vP23*LN2;
			*L2 += *vP24*LD2;
			*L2 += *vP25*LC2;
			*L2 += *vP26*LQ2;
			*L2 += *vP27*LE2;
			*L2 += *vP28*LG2;
			*L2 += *vP29*LH2;
			*L2 += *vP210*LI2;
			*L2 += *vP211*LL2;
			*L2 += *vP212*LK2;
			*L2 += *vP213*LM2;
			*L2 += *vP214*LF2;
			*L2 += *vP215*LP2;
			*L2 += *vP216*LZ2;
			*L2 += *vP217*LT2;
			*L2 += *vP218*LW2;
			*L2 += *vP219*LY2;
			*L2 += *vP220*LV2;
			*L2 += *vP221*LB2;
			
			// L(Indels)
			double L4;
			L4 = *(P2 + 20)*(La2); // Lba
			L4 += *(P2 + 41)*(Lr2); // Lbr
			L4 += *(P2 + 62)*(Ln2); // Lbn
			L4 += *(P2 + 83)*(Ld2); // Lbd
			L4 += *(P2 + 104)*(Lc2); // Lbc
			L4 += *(P2 + 125)*(Lq2); // Lbq
			L4 += *(P2 + 146)*(Le2); // Lbe
			L4 += *(P2 + 167)*(Lg2); // Lbg
			L4 += *(P2 + 188)*(Lh2); // Lbh
			L4 += *(P2 + 209)*(Li2); // Lbi
			L4 += *(P2 + 230)*(Ll2); // Lbl
			L4 += *(P2 + 251)*(Lk2); // Lbk
			L4 += *(P2 + 272)*(Lm2); // Lbm
			L4 += *(P2 + 293)*(Lf2); // Lbf
			L4 += *(P2 + 314)*(Lp2); // Lbp
			L4 += *(P2 + 335)*(Lz2); // Lbs
			L4 += *(P2 + 356)*(Lt2); // Lbt
			L4 += *(P2 + 377)*(Lw2); // Lbw
			L4 += *(P2 + 398)*(Ly2); // Lby
			L4 += *(P2 + 419)*(Lv2); // Lbv
			L4 += *(P2 + 440)*(Lb2); // Lbb
			*(Ls3 + 20) = L4;
			
			if (root &&
				(La1 != 0 ||
				Lr1 != 0 ||
				Ln1 != 0 ||
				Ld1 != 0 ||
				Lc1 != 0 ||
				Lq1 != 0 ||
				Le1 != 0 ||
				Lg1 != 0 ||
				Lh1 != 0 ||
				Li1 != 0 ||
				Ll1 != 0 ||
				Lk1 != 0 ||
				Lm1 != 0 ||
				Lf1 != 0 ||
				Lp1 != 0 ||
				Lz1 != 0 ||
				Lt1 != 0 ||
				Lw1 != 0 ||
				Ly1 != 0 ||
				Lv1 != 0 ||
				Lb1 != 0)) {
				const v20d L1 = {La1, Lr1, Ln1, Ld1, Lc1, Lq1, Le1, Lg1, Lh1, Li1, Ll1, Lk1, Lm1, Lf1, Lp1, Lz1, Lt1, Lw1, Ly1, Lv1};
				*L2 *= L1;
				*(Ls3 + 20) *= Lb1;
				*(Ls3 + 21) = *(Ls1 + 21) + *(Ls2 + 21);
			} else {
				*(Ls3 + 21) = *(Ls2 + 21);
			}
			
			for (int i = 0; i < 21; i++)
				flag |= (Ls3[i] < inv_epsilon && Ls3[i] > 0);
		} else {
			*(Ls3) = La1;
			*(Ls3 + 1) = Lr1;
			*(Ls3 + 2) = Ln1;
			*(Ls3 + 3) = Ld1;
			*(Ls3 + 4) = Lc1;
			*(Ls3 + 5) = Lq1;
			*(Ls3 + 6) = Le1;
			*(Ls3 + 7) = Lg1;
			*(Ls3 + 8) = Lh1;
			*(Ls3 + 9) = Li1;
			*(Ls3 + 10) = Ll1;
			*(Ls3 + 11) = Lk1;
			*(Ls3 + 12) = Lm1;
			*(Ls3 + 13) = Lf1;
			*(Ls3 + 14) = Lp1;
			*(Ls3 + 15) = Lz1;
			*(Ls3 + 16) = Lt1;
			*(Ls3 + 17) = Lw1;
			*(Ls3 + 18) = Ly1;
			*(Ls3 + 19) = Lv1;
			*(Ls3 + 20) = Lb1;
			*(Ls3 + 21) = *(Ls1 + 21);
		}
	}
	
	if (flag) {
		*(Ls3) *= epsilon;
		*(Ls3 + 1) *= epsilon;
		*(Ls3 + 2) *= epsilon;
		*(Ls3 + 3) *= epsilon;
		*(Ls3 + 4) *= epsilon;
		*(Ls3 + 5) *= epsilon;
		*(Ls3 + 6) *= epsilon;
		*(Ls3 + 7) *= epsilon;
		*(Ls3 + 8) *= epsilon;
		*(Ls3 + 9) *= epsilon;
		*(Ls3 + 10) *= epsilon;
		*(Ls3 + 11) *= epsilon;
		*(Ls3 + 12) *= epsilon;
		*(Ls3 + 13) *= epsilon;
		*(Ls3 + 14) *= epsilon;
		*(Ls3 + 15) *= epsilon;
		*(Ls3 + 16) *= epsilon;
		*(Ls3 + 17) *= epsilon;
		*(Ls3 + 18) *= epsilon;
		*(Ls3 + 19) *= epsilon;
		*(Ls3 + 20) *= epsilon;
		*(Ls3 + 21) += 1;
	}
}
*/
static void ProbChangeExp(double *m, double *E, double v)
{
	v = (v < 1e-6) ? 1e-6 : v;
	double A = *(m), C = *(m + 1), G = *(m + 2), T = *(m + 3), I = *(m + 4);
	double k1 = *(m + 5), k2 = *(m + 6), k3 = *(m + 7), k4 = *(m + 8), k5 = *(m + 9), k6 = *(m + 10);
	v = v/(2*(A*I*k6 + C*I*k6 + G*I*k6 + I*T*k6 + G*T + A*C*k3 + A*G*k1 + C*G*k5 + A*T*k4 + C*T*k2));
	
	double *Q = (double *) calloc(25, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *F = (double *) calloc(25, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *H = (double *) calloc(25, sizeof(double)); // initialized to zero (thread-safe on Windows)
	
	*(Q + 1) = C*k3*v;
	*(Q + 2) = G*k1*v;
	*(Q + 3) = T*k4*v;
	*(Q + 4) = I*k6*v;
	*(Q + 0) = -1*(*(Q + 1) + *(Q + 2) + *(Q + 3) + *(Q + 4));
	
	*(Q + 5) = A*k3*v;
	*(Q + 7) = G*k5*v;
	*(Q + 8) = T*k2*v;
	*(Q + 9) = I*k6*v;
	*(Q + 6) = -1*(*(Q + 5) + *(Q + 7) + *(Q + 8) + *(Q + 9));
	
	*(Q + 10) = A*k1*v;
	*(Q + 11) = C*k5*v;
	*(Q + 13) = T*v;
	*(Q + 14) = I*k6*v;
	*(Q + 12) = -1*(*(Q + 10) + *(Q + 11) + *(Q + 13) + *(Q + 14));
	
	*(Q + 15) = A*k4*v;
	*(Q + 16) = C*k2*v;
	*(Q + 17) = G*v;
	*(Q + 19) = I*k6*v;
	*(Q + 18) = -1*(*(Q + 15) + *(Q + 16) + *(Q + 17) + *(Q + 19));
	
	*(Q + 20) = A*k6*v;
	*(Q + 21) = C*k6*v;
	*(Q + 22) = G*k6*v;
	*(Q + 23) = T*k6*v;
	*(Q + 24) = -1*(*(Q + 20) + *(Q + 21) + *(Q + 22) + *(Q + 23));
	
	*(F + 0) = 1;
	*(F + 6) = 1;
	*(F + 12) = 1;
	*(F + 18) = 1;
	*(F + 24) = 1;
	
	// calculate infinity norm of the matrix Q
	int i, j, k, l;
	double r, x = 0, y;
	for (i = 0; i < 5; i++) {
		r = 0;
		for (j = i; j < 25; j += 5) {
			if (*(Q + j) >= 0) {
				r += *(Q + j);
			} else {
				r -= *(Q + j);
			}
		}
		if (r > x)
			x = r; // maximum absolute row sum
	}
	
	// E = exp(Q) = exp(Q/x)^x
	x = ceil(log2(x));
	if (x > 0) {
		double m = exp2(x);
		for (i = 0; i < 25; i++)
			*(Q + i) /= m;
	}
	
	k = 0;
	do {
		k++;
		
		// E = E + F
		for (i = 0; i < 25; i++)
			*(E + i) += *(F + i);
		
		// F = (A*F)/k
		for (i = 0; i < 5; i++) {
			for (j = 0; j < 5; j++) {
				*(H + 5*j + i) = 0;
				for (l = 0; l < 5; l++)
					*(H + 5*j + i) += *(Q + 5*l + i) * *(F + 5*j + l);
			}
		}
		for (i = 0; i < 25; i++)
			*(F + i) = *(H + i)/k;
		
		// G = E + F - E gives machine precision
		for (i = 0; i < 25; i++) {
			*(H + i) = *(E + i);
			*(H + i) += *(F + i);
			*(H + i) -= *(E + i);
		}
		
		// calculate the one norm of the matrix G
		y = 0;
		for (i = 0; i < 25;) {
			r = 0;
			for (j = i; j < i + 5; j++) {
				if (*(H + j) >= 0) {
					r += *(H + j);
				} else {
					r -= *(H + j);
				}
			}
			if (r > y)
				y = r; // maximum absolute column sum
			i = j;
		}
	} while (y > 0);
	
	if (x > 0) {
		for (k = 1; k <= x; k++) {
			for (i = 0; i < 5; i++) {
				for (j = 0; j < 5; j++) {
					*(H + 5*j + i) = 0;
					for (l = 0; l < 5; l++)
						*(H + 5*j + i) += *(E + 5*l + i) * *(E + 5*j + l);
				}
			}
			for (i = 0; i < 25; i++)
				*(E + i) = *(H + i);
		}
	}
	
	free(Q);
	free(F);
	free(H);
}

static void ProbChangeExpAA(double *m, double *E, double v)
{
	int i, j, k, l;
	double r, x = 0, y;
	
	double *Q = (double *) calloc(441, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *F = (double *) calloc(441, sizeof(double)); // initialized to zero (thread-safe on Windows)
	double *H = (double *) calloc(441, sizeof(double)); // initialized to zero (thread-safe on Windows)
	
	// fill Q matrix
	k = 0; // starting index in m
	for (j = 1; j < 20; j++) {
		for (i = 0; i < j; i++) {
			*(Q + 21*j + i) = *(m + k) * *(m + 190 + i);
			*(Q + 21*i + j) = *(m + k) * *(m + 190 + j);
			k++;
		}
	}
	for (i = 0; i < 20; i++) { // j = 20 (indels)
		*(Q + 21*j + i) = *(m + 211) * *(m + 190 + i);
		*(Q + 21*i + j) = *(m + 211) * *(m + 210);
	}
	
	r = 0; // sum of Q*freqs
	for (j = 0; j < 21; j++)
		for (i = 0; i < 21; i++)
			r += *(Q + 21*j + i) * *(m + 190 + j);
	
	v = (v < 1e-6) ? 1e-6 : v;
	v /= r;
	for (j = 0; j < 21; j++) {
		for (i = 0; i < 21; i++) {
			if (i != j) {
				*(Q + 21*j + i) *= v;
				*(Q + 21*j + j) -= *(Q + 21*j + i);
			}
		}
	}
	
	*(F + 0) = 1;
	*(F + 22) = 1;
	*(F + 44) = 1;
	*(F + 66) = 1;
	*(F + 88) = 1;
	*(F + 110) = 1;
	*(F + 132) = 1;
	*(F + 154) = 1;
	*(F + 176) = 1;
	*(F + 198) = 1;
	*(F + 220) = 1;
	*(F + 242) = 1;
	*(F + 264) = 1;
	*(F + 286) = 1;
	*(F + 308) = 1;
	*(F + 330) = 1;
	*(F + 352) = 1;
	*(F + 374) = 1;
	*(F + 396) = 1;
	*(F + 418) = 1;
	*(F + 440) = 1;
	
	// calculate infinity norm of the matrix Q
	for (i = 0; i < 21; i++) {
		r = 0;
		for (j = i; j < 441; j += 21) {
			if (*(Q + j) >= 0) {
				r += *(Q + j);
			} else {
				r -= *(Q + j);
			}
		}
		if (r > x)
			x = r; // maximum absolute row sum
	}
	
	// E = exp(Q) = exp(Q/x)^x
	x = ceil(log2(x));
	if (x > 0) {
		double m = exp2(x);
		for (i = 0; i < 441; i++)
			*(Q + i) /= m;
	}
	
	k = 0;
	do {
		k++;
		
		// E = E + F
		for (i = 0; i < 441; i++)
			*(E + i) += *(F + i);
		
		// F = (A*F)/k
		for (i = 0; i < 21; i++) {
			for (j = 0; j < 21; j++) {
				*(H + 21*j + i) = 0;
				for (l = 0; l < 21; l++)
					*(H + 21*j + i) += *(Q + 21*l + i) * *(F + 21*j + l);
			}
		}
		for (i = 0; i < 441; i++)
			*(F + i) = *(H + i)/k;
		
		// G = E + F - E gives machine precision
		for (i = 0; i < 441; i++) {
			*(H + i) = *(E + i);
			*(H + i) += *(F + i);
			*(H + i) -= *(E + i);
		}
		
		// calculate the one norm of the matrix G
		y = 0;
		for (i = 0; i < 441;) {
			r = 0;
			for (j = i; j < i + 21; j++) {
				if (*(H + j) >= 0) {
					r += *(H + j);
				} else {
					r -= *(H + j);
				}
			}
			if (r > y)
				y = r; // maximum absolute column sum
			i = j;
		}
	} while (y > 0);
	
	if (x > 0) {
		for (k = 1; k <= x; k++) {
			for (i = 0; i < 21; i++) {
				for (j = 0; j < 21; j++) {
					*(H + 21*j + i) = 0;
					for (l = 0; l < 21; l++)
						*(H + 21*j + i) += *(E + 21*l + i) * *(E + 21*j + l);
				}
			}
			for (i = 0; i < 441; i++)
				*(E + i) = *(H + i);
		}
	}
	
	free(Q);
	free(F);
	free(H);
}

static void Transpose(double *E, int n)
{
	double temp;
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++) {
			temp = *(E + n*i + j);
			*(E + n*i + j) = *(E + n*j + i);
			*(E + n*j + i) = temp;
		}
	}
}

SEXP clusterML(SEXP x, SEXP y, SEXP model, SEXP branches, SEXP lengths, SEXP states, SEXP type, SEXP weights, SEXP nThreads)
{
	// initialize variables
	XStringSet_holder y_set;
	Chars_holder y_i;
	y_set = hold_XStringSet(y);
	int length = get_length_from_XStringSet_holder(&y_set);
	int i, j, k, o, p, numRates, indels, row, params;
	int *Up;
	double *I;
	double *T = REAL(x); // Tree Topology
	double *m = REAL(model); // Substitution Model
	double s = asReal(states); // > 0 for reconstruct
	int t = asInteger(type);
	const int s1 = (t == 3) ? 21 : 5;
	const int width = s1 + 1;
	const int s2 = s1*s1;
	const int s22 = s2 + s2;
	const int s23 = s22 + s2;
	const int s24 = s23 + s2;
	const int s25 = s24 + s2;
	int *W = INTEGER(weights);
	int nthreads = asInteger(nThreads);
	
	SEXP dims;
	PROTECT(dims = GET_DIM(x));
	int l1 = INTEGER(dims)[0]; // number of rows in x (at most the length of x - 1)
	UNPROTECT(1);
	
	// alternative branch lengths
	int altL = length(lengths); // number of altered lengths
	double *ls = REAL(lengths); // new branch lengths
	int altB = length(branches); // number of altered branches
	int *bs = INTEGER(branches); // altered branches
	// If altL != altB then NNI mode, in which case,
	// bs is the index of the central branch of the NNI,
	// and ls are the five associated branch lengths:
	// center, opposite, down-left, down-right, up
	// If bs > 0 then swap down-left with opposite
	// If bs < 0 then swap down-right with opposite
	if (s != 0 && altL != 0)
		error("branches or lengths specified when states is not zero.");
	
	// calculate a vector of sequence lengths
	y_i = get_elt_from_XStringSet_holder(&y_set, 0);
	int maxWidth = y_i.length; // assume all sequences are maxWidth (aligned)
	
	double *node;
	// matrix of base probabilities at each internal node (including the root)
	// [node base site]
	if (s > 0)
		node = (double *) calloc(l1*s1*maxWidth, sizeof(double)); // initialized to zero (thread-safe on Windows)
	
	if (s > 0 || altB > 0) { // initialize pointers up the tree
		Up = Calloc(l1 - 1, int); // root node is untouched
		for (j = 1; j < l1; j++) { // for each node
			if ((int)T[6*l1 + j] > 0) { // first branch is a node
				row = (int)T[6*l1 + j] - 1;
				if (s > 0) {
					Up[row] = j;
				} else {
					Up[row] = -j;
				}
			}
			if ((int)T[7*l1 + j] > 0) { // second branch is a node
				row = (int)T[7*l1 + j] - 1;
				if (s > 0) {
					Up[row] = j;
				} else {
					Up[row] = -j;
				}
			}
		}
		for (j = 0; j < altB; j++) {
			p = *(bs + j);
			if (p < 0) // only applicable in NNI mode
				p *= -1;
			if (p > l1)
				p -= l1;
			p--;
			while (p < l1 - 1 && Up[p] < 0) {
				Up[p] *= -1;
				p = Up[p];
			}
		}
		if (t == 3) {
			I = Calloc(441, double); // initialized to zero
			I[0] = 1;
			I[22] = 1;
			I[44] = 1;
			I[66] = 1;
			I[88] = 1;
			I[110] = 1;
			I[132] = 1;
			I[154] = 1;
			I[176] = 1;
			I[198] = 1;
			I[220] = 1;
			I[242] = 1;
			I[264] = 1;
			I[286] = 1;
			I[308] = 1;
			I[330] = 1;
			I[352] = 1;
			I[374] = 1;
			I[396] = 1;
			I[418] = 1;
			I[440] = 1;
		} else {
			I = Calloc(25, double); // initialized to zero
			I[0] = 1;
			I[6] = 1;
			I[12] = 1;
			I[18] = 1;
			I[24] = 1;
		}
	}
	
	static void (*L_unknown_)(double *, const int, const int, const int, const double *, const double *, const double, const double, const int);
	static void (*L_known_)(const char *, double *, int);
	if (t == 3) {
		params = 212;
		numRates = (length(model) - params)/2; // number of bins for the gamma distribution
		indels = (*(m + 210) == 0) ? 0 : 1;
		if (indels) {
			L_unknown_ = &_L_unknown_AA_5;
		} else {
			L_unknown_ = &_L_unknown_AA;
		}
		L_known_ = &L_known_AA;
	} else {
		params = 11;
		numRates = (length(model) - params)/2; // number of bins for the gamma distribution
		indels = (*(m + 4) == 0) ? 0 : 1;
		if (indels) {
			L_unknown_ = &_L_unknown_5;
		} else {
			L_unknown_ = &_L_unknown;
		}
		L_known_ = &L_known;
	}
	
	const double epsilon = pow(2, DBL_MAX_EXP/4);
	const double inv_epsilon = 1/epsilon;
	const double log_epsilon = log(epsilon);
	int size = l1*s22 + altL*s2;
	double *P = Calloc(size*numRates, double); // initialized to zero
	for (k = 0; k < numRates; k++) { // for each bin of the gamma distribution determined by alpha
		// P = expM(Q*v)
		// transpose for cache efficiency
		#pragma omp parallel for num_threads(nthreads)
		for (i = 0; i < l1; i++) {
			if (t == 3) {
				ProbChangeExpAA(m, (P + i*s22 + k*size), T[3*l1 + i] * *(m + k + params));
			} else {
				ProbChangeExp(m, (P + i*s22 + k*size), T[3*l1 + i] * *(m + k + params));
			}
			Transpose(P + i*s22 + k*size, s1);
			if (t == 3) {
				ProbChangeExpAA(m, (P + i*s22 + s2 + k*size), T[4*l1 + i] * *(m + k + params));
			} else {
				ProbChangeExp(m, (P + i*s22 + s2 + k*size), T[4*l1 + i] * *(m + k + params));
			}
			Transpose(P + i*s22 + s2 + k*size, s1);
		}
		#pragma omp parallel for num_threads(nthreads)
		for (i = 0; i < altL; i++) {
			if (t == 3) {
				ProbChangeExpAA(m, (P + l1*s22 + i*s2 + k*size), *(ls + i) * *(m + k + params));
			} else {
				ProbChangeExp(m, (P + l1*s22 + i*s2 + k*size), *(ls + i) * *(m + k + params));
			}
			Transpose(P + l1*s22 + i*s2 + k*size, s1);
		}
	}
	
	double *sumL = Calloc(maxWidth*(altB + 1), double);
	#pragma omp parallel for private(j,k,o,p,y_i,row) num_threads(nthreads)
	for (i = 0; i < maxWidth; i++) { // for each position
		int weight;
		if (s > 0) { // reconstruct ancestral states
			weight = 1;
		} else {
			weight = *(W + i);
		}
		if (weight > 0) {
			double *Ls;
			if (altL != altB) {
				Ls = (double *) calloc((length + 1)*3*width, sizeof(double)); // initialized to zero (thread-safe on Windows)
			} else {
				Ls = (double *) calloc(l1*3*width, sizeof(double)); // initialized to zero (thread-safe on Windows)
			}
			// Ls[row, 0] = likelihood left of node
			// Ls[row, width] = likelihood right of node
			// Ls[row, 2*width] = likelihood top of node
			double *mins = (double *) calloc(altB + 1, sizeof(double)); // initialized to zero (thread-safe on Windows)
			
			for (k = 0; k < numRates; k++) { // for each bin of the gamma distribution determined by alpha
				for (j = 0; j < l1; j++) { // for each node
					// if first branch is leaf then its base L is 1
					if ((int)T[6*l1 + j] < 0) { // first branch is a leaf
						if (k == 0) {
							y_i = get_elt_from_XStringSet_holder(&y_set, (-1*(int)T[6*l1 + j] - 1));
							(*L_known_)(&y_i.ptr[i], (Ls + 0 + j*3*width), indels);
						}
					} else  { // first branch is a node
						// L is probability(branch lengths) * L(previous nodes)
						row = (int)T[6*l1 + j] - 1;
						(*L_unknown_)(Ls, j*3*width + 0, row*3*width + 0, row*3*width + width, (P + row*s22 + k*size), (P + row*s22 + s2 + k*size), epsilon, inv_epsilon, 0);
					}
					// if second branch is leaf then its base L is 1
					if ((int)T[7*l1 + j] < 0) { // second branch is a leaf
						if (k == 0) {
							y_i = get_elt_from_XStringSet_holder(&y_set, (-1*(int)T[7*l1 + j] - 1));
							(*L_known_)(&y_i.ptr[i], (Ls + width + j*3*width), indels);
						}
					} else { // second branch is a node
						// L is probability(branch lengths) * L(previous nodes)
						row = (int)T[7*l1 + j] - 1;
						(*L_unknown_)(Ls, j*3*width + width, row*3*width + 0, row*3*width + width, (P + row*s22 + k*size), (P + row*s22 + s2 + k*size), epsilon, inv_epsilon, 0);
					}
				}
				
				if (k > 0) // clear Ls[root, 2*width]
					for (j = 2*width; j < 3*width; j++)
						*(Ls + (l1 - 1)*3*width + j) = 0;
				
				if (s > 0 || altL > 0) { // descend tree
					for (j = l1 - 2; j >= 0; j--) { // for each node below the root
						row = Up[j]; // calculate likelihood at node above
						if (row > 0) {
							int side = ((int)T[6*l1 + row] == j + 1) ? width : 0; // edge opposite j
							int off1 = (side == 0) ? 0 : s2;
							if (row == (l1 - 1)) { // root is above
								// no node above so ignore by giving Ls[root, 2*width] = {0} currently
								(*L_unknown_)(Ls, j*3*width + 2*width, row*3*width + 2*width, row*3*width + side, P, (P + row*s22 + off1 + k*size), epsilon, inv_epsilon, 0);
							} else {
								int off2 = ((int)T[6*l1 + Up[row]] == row + 1) ? 0 : s2; // edge of row
								(*L_unknown_)(Ls, j*3*width + 2*width, row*3*width + 2*width, row*3*width + side, (P + Up[row]*s22 + off2 + k*size), (P + row*s22 + off1 + k*size), epsilon, inv_epsilon, 0);
							}
						}
					}
				}
				
				// calculate Likelihood of each base at final node
				row = l1 - 1;
				(*L_unknown_)(Ls, row*3*width + 2*width, row*3*width + 0, row*3*width + width, (P + row*s22 + k*size), (P + row*s22 + s2 + k*size), epsilon, inv_epsilon, 0);
				
				if (s > 0) { // final node
					if (t == 3) {
						if (!((*(Ls + 0 + row*3*width)==0 &&
							*(Ls + 1 + row*3*width)==0 &&
							*(Ls + 2 + row*3*width)==0 &&
							*(Ls + 3 + row*3*width)==0 &&
							*(Ls + 4 + row*3*width)==0 &&
							*(Ls + 5 + row*3*width)==0 &&
							*(Ls + 6 + row*3*width)==0 &&
							*(Ls + 7 + row*3*width)==0 &&
							*(Ls + 8 + row*3*width)==0 &&
							*(Ls + 9 + row*3*width)==0 &&
							*(Ls + 10 + row*3*width)==0 &&
							*(Ls + 11 + row*3*width)==0 &&
							*(Ls + 12 + row*3*width)==0 &&
							*(Ls + 13 + row*3*width)==0 &&
							*(Ls + 14 + row*3*width)==0 &&
							*(Ls + 15 + row*3*width)==0 &&
							*(Ls + 16 + row*3*width)==0 &&
							*(Ls + 17 + row*3*width)==0 &&
							*(Ls + 18 + row*3*width)==0 &&
							*(Ls + 19 + row*3*width)==0) ||
							(*(Ls + width + 0 + row*3*width)==0 &&
							*(Ls + width + 1 + row*3*width)==0 &&
							*(Ls + width + 2 + row*3*width)==0 &&
							*(Ls + width + 3 + row*3*width)==0 &&
							*(Ls + width + 4 + row*3*width)==0 &&
							*(Ls + width + 5 + row*3*width)==0 &&
							*(Ls + width + 6 + row*3*width)==0 &&
							*(Ls + width + 7 + row*3*width)==0 &&
							*(Ls + width + 8 + row*3*width)==0 &&
							*(Ls + width + 9 + row*3*width)==0 &&
							*(Ls + width + 10 + row*3*width)==0 &&
							*(Ls + width + 11 + row*3*width)==0 &&
							*(Ls + width + 12 + row*3*width)==0 &&
							*(Ls + width + 13 + row*3*width)==0 &&
							*(Ls + width + 14 + row*3*width)==0 &&
							*(Ls + width + 15 + row*3*width)==0 &&
							*(Ls + width + 16 + row*3*width)==0 &&
							*(Ls + width + 17 + row*3*width)==0 &&
							*(Ls + width + 18 + row*3*width)==0 &&
							*(Ls + width + 19 + row*3*width)==0))) { // neither branch is a gap
							for (j = 0; j < s1; j++)
								node[row*maxWidth*s1 + i*s1 + j] += *(Ls + 2*width + j + row*3*width) * *(m + numRates + k + params);
						}
					} else {
						if (!((*(Ls + 0 + row*3*width)==0 &&
							*(Ls + 1 + row*3*width)==0 &&
							*(Ls + 2 + row*3*width)==0 &&
							*(Ls + 3 + row*3*width)==0) ||
							(*(Ls + width + 0 + row*3*width)==0 &&
							*(Ls + width + 1 + row*3*width)==0 &&
							*(Ls + width + 2 + row*3*width)==0 &&
							*(Ls + width + 3 + row*3*width)==0))) { // neither branch is a gap
							for (j = 0; j < s1; j++)
								node[row*maxWidth*s1 + i*s1 + j] += *(Ls + 2*width + j + row*3*width) * *(m + numRates + k + params);
						}
					}
				}
				
				// calculate overall Likelihood
				int count = 0;
				for (o = 0; o <= altB; o++) {
					if (o > 0) {
						if (altL != altB) { // NNI mode
							int down, flip, side;
							p = *(bs + o - 1);
							if (p < 0) {
								flip = -1; // swap down-right with opposite
								p *= -1;
							} else {
								flip = 1; // swap down-left with opposite
							}
							if (p > l1) {
								side = 0;
								p -= l1;
								p--;
								down = (int)T[7*l1 + p] - 1; // center is right branch
							} else {
								side = width;
								p--;
								down = (int)T[6*l1 + p] - 1; // center is left branch
							}
							
							if (p == row) { // root node
								int opposite;
								if (side == 0) {
									opposite = (int)T[6*l1 + p] - 1;
								} else {
									opposite = (int)T[7*l1 + p] - 1;
								}
								if (flip > 0) {
									// opposite-left merging with down-right
									(*L_unknown_)(Ls, l1*3*width + 0, opposite*3*width + 0, down*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s22), (P + l1*s22 + (o - 1)*s25 + k*size + s23), epsilon, inv_epsilon, 0);
									// opposite-right merging with down-left
									(*L_unknown_)(Ls, l1*3*width + width, opposite*3*width + width, down*3*width + 0, (P + l1*s22 + (o - 1)*s25 + k*size + s24), (P + l1*s22 + (o - 1)*s25 + k*size + s2), epsilon, inv_epsilon, 0);
								} else {
									// opposite-left merging with down-left
									(*L_unknown_)(Ls, l1*3*width + 0, opposite*3*width + 0, down*3*width + 0, (P + l1*s22 + (o - 1)*s25 + k*size + s23), (P + l1*s22 + (o - 1)*s25 + k*size + s22), epsilon, inv_epsilon, 0);
									// opposite-right merging with down-right
									(*L_unknown_)(Ls, l1*3*width + width, opposite*3*width + width, down*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s24), (P + l1*s22 + (o - 1)*s25 + k*size + s2), epsilon, inv_epsilon, 0);
								}
								// merge both across center
								if (side == 0) { // center is right branch
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, (P + row*s22 + k*size), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
								} else { // center is left branch
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, (P + row*s22 + s2 + k*size), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
								}
							} else {
								if (count == 11) {
									count = 1;
								} else {
									count++;
								}
								if (count == 1) {
									if (flip > 0) {
										// opposite merging with down-right
										(*L_unknown_)(Ls, l1*3*width + 0, p*3*width + side, down*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s22), (P + l1*s22 + (o - 1)*s25 + k*size + s23), epsilon, inv_epsilon, 0);
										// up merging with down-left
										(*L_unknown_)(Ls, l1*3*width + width, p*3*width + 2*width, down*3*width + 0, (P + l1*s22 + (o - 1)*s25 + k*size + s24), (P + l1*s22 + (o - 1)*s25 + k*size + s2), epsilon, inv_epsilon, 0);
									} else {
										// opposite merging with down-left
										(*L_unknown_)(Ls, l1*3*width + 0, p*3*width + side, down*3*width + 0, (P + l1*s22 + (o - 1)*s25 + k*size + s23), (P + l1*s22 + (o - 1)*s25 + k*size + s22), epsilon, inv_epsilon, 0);
										// up merging with down-right
										(*L_unknown_)(Ls, l1*3*width + width, p*3*width + 2*width, down*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s24), (P + l1*s22 + (o - 1)*s25 + k*size + s2), epsilon, inv_epsilon, 0);
									}
									// merge both across center
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, I, (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 1);
									
									// up merging with below
									(*L_unknown_)(Ls, length*3*width + 0, p*3*width + 2*width, l1*3*width + 0, (P + l1*s22 + (o - 1)*s25 + k*size + s24), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
									if (flip > 0) {
										// down-left merging with below
										(*L_unknown_)(Ls, l1*3*width + 2*width, down*3*width + 0, l1*3*width + 0, (P + l1*s22 + (o - 1)*s25 + k*size + s2), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
										// down-right merging with above
										(*L_unknown_)(Ls, length*3*width + width, down*3*width + width, l1*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s23), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
										// opposite merging with above
										(*L_unknown_)(Ls, length*3*width + 2*width, p*3*width + side, l1*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s22), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
									} else {
										// down-right merging with below
										(*L_unknown_)(Ls, l1*3*width + 2*width, down*3*width + width, l1*3*width + 0, (P + l1*s22 + (o - 1)*s25 + k*size + s2), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
										// down-left merging with above
										(*L_unknown_)(Ls, length*3*width + width, down*3*width + 0, l1*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s22), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
										// opposite merging with above
										(*L_unknown_)(Ls, length*3*width + 2*width, p*3*width + side, l1*3*width + width, (P + l1*s22 + (o - 1)*s25 + k*size + s23), (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 0);
									}
								} else if (count == 2 || count == 7) {
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 0, l1*3*width + width, I, (P + l1*s22 + (o - 1)*s25 + k*size + 0), epsilon, inv_epsilon, 1);
								} else if (count == 3 || count == 8) {
									if (flip > 0) {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 0, down*3*width + 0, I, (P + l1*s22 + (o - 1)*s25 + k*size + s2), epsilon, inv_epsilon, 1);
									} else {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 0, down*3*width + width, I, (P + l1*s22 + (o - 1)*s25 + k*size + s2), epsilon, inv_epsilon, 1);
									}
								} else if (count == 4 || count == 9) {
									if (flip > 0) {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + width, p*3*width + side, I, (P + l1*s22 + (o - 1)*s25 + k*size + s22), epsilon, inv_epsilon, 1);
									} else {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 2*width, down*3*width + 0, I, (P + l1*s22 + (o - 1)*s25 + k*size + s22), epsilon, inv_epsilon, 1);
									}
								} else if (count == 5 || count == 10) {
									if (flip > 0) {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + 2*width, down*3*width + width, I, (P + l1*s22 + (o - 1)*s25 + k*size + s23), epsilon, inv_epsilon, 1);
									} else {
										(*L_unknown_)(Ls, row*3*width + 2*width, length*3*width + width, p*3*width + side, I, (P + l1*s22 + (o - 1)*s25 + k*size + s23), epsilon, inv_epsilon, 1);
									}
								} else { // count == 6 || count == 11
									(*L_unknown_)(Ls, row*3*width + 2*width, l1*3*width + 2*width, p*3*width + 2*width, I, (P + l1*s22 + (o - 1)*s25 + k*size + s24), epsilon, inv_epsilon, 1);
								}
							}
						} else {
							p = *(bs + o - 1);
							int side;
							if (p > l1) {
								side = width;
								p -= l1;
							} else {
								side = 0;
							}
							
							if (p == l1) { // root node
								if (side == width) {
									(*L_unknown_)(Ls, row*3*width + 2*width, row*3*width + 0, row*3*width + width, (P + row*s22 + k*size), (P + l1*s22 + (o - 1)*s2 + k*size), epsilon, inv_epsilon, 0);
								} else {
									(*L_unknown_)(Ls, row*3*width + 2*width, row*3*width + 0, row*3*width + width, (P + l1*s22 + (o - 1)*s2 + k*size), (P + row*s22 + s2 + k*size), epsilon, inv_epsilon, 0);
								}
							} else { // root at p - 1
								// calculate likelihood at node
								int off1 = ((int)T[6*l1 + Up[p - 1]] == p) ? 0 : s2; // branch of row
								if (side == 0) {
									(*L_unknown_)(Ls, row*3*width + 2*width, (p - 1)*3*width + width, (p - 1)*3*width + 2*width, (P + (p - 1)*s22 + s2 + k*size), (P + Up[p - 1]*s22 + off1 + k*size), epsilon, inv_epsilon, 0);
								} else {
									(*L_unknown_)(Ls, row*3*width + 2*width, (p - 1)*3*width + 0, (p - 1)*3*width + 2*width, (P + (p - 1)*s22 + k*size), (P + Up[p - 1]*s22 + off1 + k*size), epsilon, inv_epsilon, 0);
								}
								(*L_unknown_)(Ls, row*3*width + 2*width, (p - 1)*3*width + side, row*3*width + 2*width, I, (P + l1*s22 + (o - 1)*s2 + k*size), epsilon, inv_epsilon, 1);
							}
						}
					}
					
					if (k > 0) {
						double delta = *(Ls + 2*width + s1 + row*3*width) - mins[o];
						if (delta < 0) {
							for (j = 0; j > delta; j--)
								*(sumL + i + o*maxWidth) *= inv_epsilon;
							mins[o] = *(Ls + 2*width + s1 + row*3*width);
						} else if (delta > 0) {
							while (delta > 0) {
								for (j = 0; j < s1; j++)
									*(Ls + 2*width + j + row*3*width) *= inv_epsilon;
								delta--;
							}
						}
					} else {
						mins[o] = *(Ls + 2*width + s1 + row*3*width);
					}
					
					if (t == 3) {
						for (j = 0; j < s1; j++)
							*(sumL + i + o*maxWidth) += (*(m + j + 190) * *(Ls + 2*width + j + row*3*width)) * *(m + numRates + k + params);
					} else {
						for (j = 0; j < s1; j++)
							*(sumL + i + o*maxWidth) += (*(m + j) * *(Ls + 2*width + j + row*3*width)) * *(m + numRates + k + params);
					}
				}
				
				if (s > 0) {
					for (j = l1 - 2; j >= 0; j--) { // for each node below the root
						// apply three-way parsimony to resolve gaps
						int c1 = 0, c2 = 0;
						for (o = 0; o < s1 - 1; o++) {
							c1 |= *(Ls + o + j*3*width) != 0;
							c2 |= *(Ls + width + o + j*3*width) != 0;
						}
						int c = c1 + c2;
						if (c >= 1) { // at least one non-gap
							int side = ((int)T[6*l1 + Up[j]] == j + 1) ? 0 : width; // edge of row
							if (c == 2 || c2) {
								int off1 = (side == 0) ? 0 : s2;
								(*L_unknown_)(Ls, row*3*width + 2*width, Up[j]*3*width + side, j*3*width + 2*width, I, (P + Up[j]*s22 + off1 + k*size), epsilon, inv_epsilon, 1);
								for (o = 0; o < s1; o++)
									node[j*maxWidth*s1 + i*s1 + o] += *(Ls + 2*width + o + row*3*width) * *(m + numRates + k + params);
							}
						}
					}
				}
			}
			
			for (o = 0; o <= altB; o++)
				if (*(sumL + i + o*maxWidth) > 0)
					*(sumL + i + o*maxWidth) = log(*(sumL + i + o*maxWidth)) - mins[o]*log_epsilon;
			
			free(Ls);
			free(mins);
		}
	}
	Free(P);
	
	if (altL > 0) {
		Free(Up);
		Free(I);
	}
	
	SEXP ans;
	if (s > 0) { // return character states
		SEXP ans1, ans2;
		PROTECT(ans = allocVector(VECSXP, 2));
		PROTECT(ans1 = allocVector(STRSXP, l1));
		s = (s >= 1) ? 0.9999999 : s; // machine precision
		if (t == 3) {
			for (j = 0; j < l1; j++) {
				char *rans1 = Calloc(maxWidth + 1, char);
				for (i = 0; i < maxWidth; i++) {
					double La = *(m) * node[j*maxWidth*s1 + i*s1 + 0];
					double Lr = *(m + 1) * node[j*maxWidth*s1 + i*s1 + 1];
					double Ln = *(m + 2) * node[j*maxWidth*s1 + i*s1 + 2];
					double Ld = *(m + 3) * node[j*maxWidth*s1 + i*s1 + 3];
					double Lc = *(m + 4) * node[j*maxWidth*s1 + i*s1 + 4];
					double Lq = *(m + 5) * node[j*maxWidth*s1 + i*s1 + 5];
					double Le = *(m + 6) * node[j*maxWidth*s1 + i*s1 + 6];
					double Lg = *(m + 7) * node[j*maxWidth*s1 + i*s1 + 7];
					double Lh = *(m + 8) * node[j*maxWidth*s1 + i*s1 + 8];
					double Li = *(m + 9) * node[j*maxWidth*s1 + i*s1 + 9];
					double Ll = *(m + 10) * node[j*maxWidth*s1 + i*s1 + 10];
					double Lk = *(m + 11) * node[j*maxWidth*s1 + i*s1 + 11];
					double Lm = *(m + 12) * node[j*maxWidth*s1 + i*s1 + 12];
					double Lf = *(m + 13) * node[j*maxWidth*s1 + i*s1 + 13];
					double Lp = *(m + 14) * node[j*maxWidth*s1 + i*s1 + 14];
					double Lz = *(m + 15) * node[j*maxWidth*s1 + i*s1 + 15];
					double Lt = *(m + 16) * node[j*maxWidth*s1 + i*s1 + 16];
					double Lw = *(m + 17) * node[j*maxWidth*s1 + i*s1 + 17];
					double Ly = *(m + 18) * node[j*maxWidth*s1 + i*s1 + 18];
					double Lv = *(m + 19) * node[j*maxWidth*s1 + i*s1 + 19];
					double L = *(m + 20) * node[j*maxWidth*s1 + i*s1 + 20];
					if (s*L >= La && s*L >= Lr && s*L >= Ln && s*L >= Ld && s*L >= Lc && s*L >= Lq && s*L >= Le && s*L >= Lg && s*L >= Lh && s*L >= Li && s*L >= Ll && s*L >= Lk && s*L >= Lm && s*L >= Lf && s*L >= Lp && s*L >= Lz && s*L >= Lt && s*L >= Lw && s*L >= Ly && s*L >= Lv) {
						rans1[i] = '-';
					} else if (s*La >= Lr && s*La >= Ln && s*La >= Ld && s*La >= Lc && s*La >= Lq && s*La >= Le && s*La >= Lg && s*La >= Lh && s*La >= Li && s*La >= Ll && s*La >= Lk && s*La >= Lm && s*La >= Lf && s*La >= Lp && s*La >= Lz && s*La >= Lt && s*La >= Lw && s*La >= Ly && s*La >= Lv) {
						rans1[i] = 'A';
					} else if (s*Lr >= La && s*Lr >= Ln && s*Lr >= Ld && s*Lr >= Lc && s*Lr >= Lq && s*Lr >= Le && s*Lr >= Lg && s*Lr >= Lh && s*Lr >= Li && s*Lr >= Ll && s*Lr >= Lk && s*Lr >= Lm && s*Lr >= Lf && s*Lr >= Lp && s*Lr >= Lz && s*Lr >= Lt && s*Lr >= Lw && s*Lr >= Ly && s*Lr >= Lv) {
						rans1[i] = 'R';
					} else if (s*Ln >= La && s*Ln >= Lr && s*Ln >= Ld && s*Ln >= Lc && s*Ln >= Lq && s*Ln >= Le && s*Ln >= Lg && s*Ln >= Lh && s*Ln >= Li && s*Ln >= Ll && s*Ln >= Lk && s*Ln >= Lm && s*Ln >= Lf && s*Ln >= Lp && s*Ln >= Lz && s*Ln >= Lt && s*Ln >= Lw && s*Ln >= Ly && s*Ln >= Lv) {
						rans1[i] = 'N';
					} else if (s*Ld >= La && s*Ld >= Lr && s*Ld >= Ln && s*Ld >= Lc && s*Ld >= Lq && s*Ld >= Le && s*Ld >= Lg && s*Ld >= Lh && s*Ld >= Li && s*Ld >= Ll && s*Ld >= Lk && s*Ld >= Lm && s*Ld >= Lf && s*Ld >= Lp && s*Ld >= Lz && s*Ld >= Lt && s*Ld >= Lw && s*Ld >= Ly && s*Ld >= Lv) {
						rans1[i] = 'D';
					} else if (s*Lc >= La && s*Lc >= Lr && s*Lc >= Ln && s*Lc >= Ld && s*Lc >= Lq && s*Lc >= Le && s*Lc >= Lg && s*Lc >= Lh && s*Lc >= Li && s*Lc >= Ll && s*Lc >= Lk && s*Lc >= Lm && s*Lc >= Lf && s*Lc >= Lp && s*Lc >= Lz && s*Lc >= Lt && s*Lc >= Lw && s*Lc >= Ly && s*Lc >= Lv) {
						rans1[i] = 'C';
					} else if (s*Lq >= La && s*Lq >= Lr && s*Lq >= Ln && s*Lq >= Ld && s*Lq >= Lc && s*Lq >= Le && s*Lq >= Lg && s*Lq >= Lh && s*Lq >= Li && s*Lq >= Ll && s*Lq >= Lk && s*Lq >= Lm && s*Lq >= Lf && s*Lq >= Lp && s*Lq >= Lz && s*Lq >= Lt && s*Lq >= Lw && s*Lq >= Ly && s*Lq >= Lv) {
						rans1[i] = 'Q';
					} else if (s*Le >= La && s*Le >= Lr && s*Le >= Ln && s*Le >= Ld && s*Le >= Lc && s*Le >= Lq && s*Le >= Lg && s*Le >= Lh && s*Le >= Li && s*Le >= Ll && s*Le >= Lk && s*Le >= Lm && s*Le >= Lf && s*Le >= Lp && s*Le >= Lz && s*Le >= Lt && s*Le >= Lw && s*Le >= Ly && s*Le >= Lv) {
						rans1[i] = 'E';
					} else if (s*Lg >= La && s*Lg >= Lr && s*Lg >= Ln && s*Lg >= Ld && s*Lg >= Lc && s*Lg >= Lq && s*Lg >= Le && s*Lg >= Lh && s*Lg >= Li && s*Lg >= Ll && s*Lg >= Lk && s*Lg >= Lm && s*Lg >= Lf && s*Lg >= Lp && s*Lg >= Lz && s*Lg >= Lt && s*Lg >= Lw && s*Lg >= Ly && s*Lg >= Lv) {
						rans1[i] = 'G';
					} else if (s*Lh >= La && s*Lh >= Lr && s*Lh >= Ln && s*Lh >= Ld && s*Lh >= Lc && s*Lh >= Lq && s*Lh >= Le && s*Lh >= Lg && s*Lh >= Li && s*Lh >= Ll && s*Lh >= Lk && s*Lh >= Lm && s*Lh >= Lf && s*Lh >= Lp && s*Lh >= Lz && s*Lh >= Lt && s*Lh >= Lw && s*Lh >= Ly && s*Lh >= Lv) {
						rans1[i] = 'H';
					} else if (s*Li >= La && s*Li >= Lr && s*Li >= Ln && s*Li >= Ld && s*Li >= Lc && s*Li >= Lq && s*Li >= Le && s*Li >= Lg && s*Li >= Lh && s*Li >= Ll && s*Li >= Lk && s*Li >= Lm && s*Li >= Lf && s*Li >= Lp && s*Li >= Lz && s*Li >= Lt && s*Li >= Lw && s*Li >= Ly && s*Li >= Lv) {
						rans1[i] = 'I';
					} else if (s*Ll >= La && s*Ll >= Lr && s*Ll >= Ln && s*Ll >= Ld && s*Ll >= Lc && s*Ll >= Lq && s*Ll >= Le && s*Ll >= Lg && s*Ll >= Lh && s*Ll >= Li && s*Ll >= Lk && s*Ll >= Lm && s*Ll >= Lf && s*Ll >= Lp && s*Ll >= Lz && s*Ll >= Lt && s*Ll >= Lw && s*Ll >= Ly && s*Ll >= Lv) {
						rans1[i] = 'L';
					} else if (s*Lk >= La && s*Lk >= Lr && s*Lk >= Ln && s*Lk >= Ld && s*Lk >= Lc && s*Lk >= Lq && s*Lk >= Le && s*Lk >= Lg && s*Lk >= Lh && s*Lk >= Li && s*Lk >= Ll && s*Lk >= Lm && s*Lk >= Lf && s*Lk >= Lp && s*Lk >= Lz && s*Lk >= Lt && s*Lk >= Lw && s*Lk >= Ly && s*Lk >= Lv) {
						rans1[i] = 'K';
					} else if (s*Lm >= La && s*Lm >= Lr && s*Lm >= Ln && s*Lm >= Ld && s*Lm >= Lc && s*Lm >= Lq && s*Lm >= Le && s*Lm >= Lg && s*Lm >= Lh && s*Lm >= Li && s*Lm >= Ll && s*Lm >= Lk && s*Lm >= Lf && s*Lm >= Lp && s*Lm >= Lz && s*Lm >= Lt && s*Lm >= Lw && s*Lm >= Ly && s*Lm >= Lv) {
						rans1[i] = 'M';
					} else if (s*Lf >= La && s*Lf >= Lr && s*Lf >= Ln && s*Lf >= Ld && s*Lf >= Lc && s*Lf >= Lq && s*Lf >= Le && s*Lf >= Lg && s*Lf >= Lh && s*Lf >= Li && s*Lf >= Ll && s*Lf >= Lk && s*Lf >= Lm && s*Lf >= Lp && s*Lf >= Lz && s*Lf >= Lt && s*Lf >= Lw && s*Lf >= Ly && s*Lf >= Lv) {
						rans1[i] = 'F';
					} else if (s*Lp >= La && s*Lp >= Lr && s*Lp >= Ln && s*Lp >= Ld && s*Lp >= Lc && s*Lp >= Lq && s*Lp >= Le && s*Lp >= Lg && s*Lp >= Lh && s*Lp >= Li && s*Lp >= Ll && s*Lp >= Lk && s*Lp >= Lm && s*Lp >= Lf && s*Lp >= Lz && s*Lp >= Lt && s*Lp >= Lw && s*Lp >= Ly && s*Lp >= Lv) {
						rans1[i] = 'P';
					} else if (s*Lz >= La && s*Lz >= Lr && s*Lz >= Ln && s*Lz >= Ld && s*Lz >= Lc && s*Lz >= Lq && s*Lz >= Le && s*Lz >= Lg && s*Lz >= Lh && s*Lz >= Li && s*Lz >= Ll && s*Lz >= Lk && s*Lz >= Lm && s*Lz >= Lf && s*Lz >= Lp && s*Lz >= Lt && s*Lz >= Lw && s*Lz >= Ly && s*Lz >= Lv) {
						rans1[i] = 'S';
					} else if (s*Lt >= La && s*Lt >= Lr && s*Lt >= Ln && s*Lt >= Ld && s*Lt >= Lc && s*Lt >= Lq && s*Lt >= Le && s*Lt >= Lg && s*Lt >= Lh && s*Lt >= Li && s*Lt >= Ll && s*Lt >= Lk && s*Lt >= Lm && s*Lt >= Lf && s*Lt >= Lp && s*Lt >= Lz && s*Lt >= Lw && s*Lt >= Ly && s*Lt >= Lv) {
						rans1[i] = 'T';
					} else if (s*Lw >= La && s*Lw >= Lr && s*Lw >= Ln && s*Lw >= Ld && s*Lw >= Lc && s*Lw >= Lq && s*Lw >= Le && s*Lw >= Lg && s*Lw >= Lh && s*Lw >= Li && s*Lw >= Ll && s*Lw >= Lk && s*Lw >= Lm && s*Lw >= Lf && s*Lw >= Lp && s*Lw >= Lz && s*Lw >= Lt && s*Lw >= Ly && s*Lw >= Lv) {
						rans1[i] = 'W';
					} else if (s*Ly >= La && s*Ly >= Lr && s*Ly >= Ln && s*Ly >= Ld && s*Ly >= Lc && s*Ly >= Lq && s*Ly >= Le && s*Ly >= Lg && s*Ly >= Lh && s*Ly >= Li && s*Ly >= Ll && s*Ly >= Lk && s*Ly >= Lm && s*Ly >= Lf && s*Ly >= Lp && s*Ly >= Lz && s*Ly >= Lt && s*Ly >= Lw && s*Ly >= Lv) {
						rans1[i] = 'Y';
					} else if (s*Lv >= La && s*Lv >= Lr && s*Lv >= Ln && s*Lv >= Ld && s*Lv >= Lc && s*Lv >= Lq && s*Lv >= Le && s*Lv >= Lg && s*Lv >= Lh && s*Lv >= Li && s*Lv >= Ll && s*Lv >= Lk && s*Lv >= Lm && s*Lv >= Lf && s*Lv >= Lp && s*Lv >= Lz && s*Lv >= Lt && s*Lv >= Lw && s*Lv >= Ly) {
						rans1[i] = 'V';
					} else if (s*Ln >= La && s*Ln >= Lr && s*Ln >= Lc && s*Ln >= Lq && s*Ln >= Le && s*Ln >= Lg && s*Ln >= Lh && s*Ln >= Li && s*Ln >= Ll && s*Ln >= Lk && s*Ln >= Lm && s*Ln >= Lf && s*Ln >= Lp && s*Ln >= Lz && s*Ln >= Lt && s*Ln >= Lw && s*Ln >= Ly && s*Ln >= Lv &&
						s*Ld >= La && s*Ld >= Lr && s*Ld >= Lc && s*Ld >= Lq && s*Ld >= Le && s*Ld >= Lg && s*Ld >= Lh && s*Ld >= Li && s*Ld >= Ll && s*Ld >= Lk && s*Ld >= Lm && s*Ld >= Lf && s*Ld >= Lp && s*Ld >= Lz && s*Ld >= Lt && s*Ld >= Lw && s*Ld >= Ly && s*Ld >= Lv) {
						rans1[i] = 'B';
					} else if (s*Lq >= La && s*Lq >= Lr && s*Lq >= Ln && s*Lq >= Ld && s*Lq >= Lc && s*Lq >= Lg && s*Lq >= Lh && s*Lq >= Li && s*Lq >= Ll && s*Lq >= Lk && s*Lq >= Lm && s*Lq >= Lf && s*Lq >= Lp && s*Lq >= Lz && s*Lq >= Lt && s*Lq >= Lw && s*Lq >= Ly && s*Lq >= Lv &&
						s*Le >= La && s*Le >= Lr && s*Le >= Ln && s*Le >= Ld && s*Le >= Lc && s*Le >= Lg && s*Le >= Lh && s*Le >= Li && s*Le >= Ll && s*Le >= Lk && s*Le >= Lm && s*Le >= Lf && s*Le >= Lp && s*Le >= Lz && s*Le >= Lt && s*Le >= Lw && s*Le >= Ly && s*Le >= Lv) {
						rans1[i] = 'Z';
					} else if (s*Li >= La && s*Li >= Lr && s*Li >= Ln && s*Li >= Ld && s*Li >= Lc && s*Li >= Lq && s*Li >= Le && s*Li >= Lg && s*Li >= Lh && s*Li >= Lk && s*Li >= Lm && s*Li >= Lf && s*Li >= Lp && s*Li >= Lz && s*Li >= Lt && s*Li >= Lw && s*Li >= Ly && s*Li >= Lv &&
						s*Ll >= La && s*Ll >= Lr && s*Ll >= Ln && s*Ll >= Ld && s*Ll >= Lc && s*Ll >= Lq && s*Ll >= Le && s*Ll >= Lg && s*Ll >= Lh && s*Ll >= Lk && s*Ll >= Lm && s*Ll >= Lf && s*Ll >= Lp && s*Ll >= Lz && s*Ll >= Lt && s*Ll >= Lw && s*Ll >= Ly && s*Ll >= Lv) {
						rans1[i] = 'J';
					} else {
						rans1[i] = 'X';
					}
				}
				rans1[i] = '\0'; // end (null terminate) the string
				SET_STRING_ELT(ans1, j, mkChar(rans1));
				Free(rans1);
			}
		} else {
			char TU = t == 1 ? 'T' : 'U';
			for (j = 0; j < l1; j++) {
				char *rans1 = Calloc(maxWidth + 1, char);
				for (i = 0; i < maxWidth; i++) {
					double La = *(m) * node[j*maxWidth*s1 + i*s1 + 0];
					double Lc = *(m + 1) * node[j*maxWidth*s1 + i*s1 + 1];
					double Lg = *(m + 2) * node[j*maxWidth*s1 + i*s1 + 2];
					double Lt = *(m + 3) * node[j*maxWidth*s1 + i*s1 + 3];
					double Li = *(m + 4) * node[j*maxWidth*s1 + i*s1 + 4];
					if (s*Li >= La && s*Li >= Lc && s*Li >= Lg && s*Li >= Lt) {
						rans1[i] = '-';
					} else if (s*La > Lc && s*La > Lg && s*La > Lt) {
						rans1[i] = 'A';
					} else if (s*Lc > La && s*Lc > Lg && s*Lc > Lt) {
						rans1[i] = 'C';
					} else if (s*Lg > La && s*Lg > Lc && s*Lg > Lt) {
						rans1[i] = 'G';
					} else if (s*Lt > La && s*Lt > Lc && s*Lt > Lg) {
						rans1[i] = TU;
					} else if (s*La > Lg && s*La > Lt && s*Lc > Lg && s*Lc > Lt) { // M = A or C
						rans1[i] = 'M';
					} else if (s*La > Lc && s*La > Lt && s*Lg > Lc && s*Lg > Lt) { // R = A or G
						rans1[i] = 'R';
					} else if (s*La > Lc && s*La > Lg && s*Lt > Lc && s*Lt > Lg) { // W = A or T
						rans1[i] = 'W';
					} else if (s*Lc > La && s*Lc > Lt && s*Lg > La && s*Lg > Lt) { // S = C or G
						rans1[i] = 'S';
					} else if (s*Lc > La && s*Lc > Lg && s*Lt > La && s*Lt > Lg) { // Y = C or T
						rans1[i] = 'Y';
					} else if (s*Lg > La && s*Lg > Lc && s*Lt > La && s*Lt > Lc) { // K = G or T
						rans1[i] = 'K';
					} else if (s*La > Lt && s*Lc > Lt && s*Lg > Lt) { // V = A or C or G
						rans1[i] = 'V';
					} else if (s*La > Lg && s*Lc > Lg && s*Lt > Lg) { // H = A or C or T
						rans1[i] = 'H';
					} else if (s*La > Lc && s*Lg > Lc && s*Lt > Lc) { // D = A or G or T
						rans1[i] = 'D';
					} else if (s*Lc > La && s*Lg > La && s*Lt > La) { // B = C or G or T
						rans1[i] = 'B';
					} else { // N = A or C or G or T
						rans1[i] = 'N';
					}
				}
				rans1[i] = '\0'; // end (null terminate) the string
				SET_STRING_ELT(ans1, j, mkChar(rans1));
				Free(rans1);
			}
		}
		free(node);
		PROTECT(ans2 = allocVector(REALSXP, maxWidth));
		double *rans2 = REAL(ans2);
		for (i = 0; i < maxWidth; i++) {
			if (*(sumL + i) < 0) {
				rans2[i] = *(sumL + i);
			} else {
				rans2[i] = 0;
			}
		}
		
		SET_VECTOR_ELT(ans, 0, ans1);
		SET_VECTOR_ELT(ans, 1, ans2);
		UNPROTECT(2);
	} else { // return -LnL
		PROTECT(ans = allocVector(REALSXP, altB + 1));
		double *rans = REAL(ans);
		for (o = 0; o <= altB; o++) {
			rans[o] = 0;
			for (i = 0; i < maxWidth; i++) {
				int weight = *(W + i);
				if (weight > 0)
					rans[o] -= ((double)weight) * *(sumL + i + o*maxWidth);
			}
		}
	}
	
	Free(sumL);
	UNPROTECT(1);
	
	return ans;
}

SEXP expM(SEXP x, SEXP model, SEXP type)
{
	// initialize variables
	double l = asReal(x);
	double *m = REAL(model); // Substitution Model
	int t = asInteger(type);
	int size = (t == 3) ? 21 : 5;
	
	SEXP ans;
	PROTECT(ans = allocMatrix(REALSXP, size, size));
	double *rans = REAL(ans);
	size *= size;
	for (int i = 0; i < size; i++)
		rans[i] = 0;
	
	if (t == 3) {
		ProbChangeExpAA(m, rans, l);
	} else {
		ProbChangeExp(m, rans, l);
	}
	
	UNPROTECT(1);
	
	return ans;
}
