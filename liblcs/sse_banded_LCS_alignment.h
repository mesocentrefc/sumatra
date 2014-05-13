/*
 * sse_banded_LCS_alignment.h
 *
 *  Created on: november 29, 2012
 *      Author: mercier
 */

#ifndef SSE_BANDED_LCS_ALIGNMENT_H_
#define SSE_BANDED_LCS_ALIGNMENT_H_
#include <stdint.h>

double sse_banded_lcs_align(int16_t* seq1, int16_t* seq2, int l1, int l2, BOOL normalize, int reference, BOOL lcsmode, int16_t* address, int LCSmin);
int calculateSizeToAllocate(int maxLen, int minLen, int LCSmin);
void calculateThresholdFromErrorNumber(int error, int length, double* threshold);
void iniSeq(int16_t* seq, int size, int16_t iniValue);
void putSeqInSeq(int16_t* seq, char* s, int l, BOOL reverse);
double generic_sse_banded_lcs_align(char* seq1, char* seq2, double threshold, BOOL normalize, int reference, BOOL lcsmode, int16_t** address, int* buffer_size, int16_t** iseq1,
		int16_t** iseq2, int* buffer_sizeS);
int prepareTablesForSumathings(int lmax, int lmin, double threshold, BOOL normalize, int reference, BOOL lcsmode,
		int16_t** address, int16_t** iseq1, int16_t** iseq2);
double alignForSumathings(char* seq1, int16_t* iseq1, char* seq2, int16_t* iseq2, int l1, int l2, BOOL normalize,
		int reference, BOOL lcsmode, int16_t* address, int sizeForSeqs, int LCSmin);
int calculateLCSmin(int l1, int l2, double threshold, BOOL normalize, int reference, BOOL lcsmode);
#endif
