/*
 * sse_banded_LCS_alignment.c
 *
 *  Created on: 7 nov. 2012
 *      Author: celine mercier
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "../libutils/utilities.h"
#include "../libsse/_sse.h"



/*static void  printreg(__m128i r)
{
	int16_t a0,a1,a2,a3,a4,a5,a6,a7;

	a0= _MM_EXTRACT_EPI16(r,0);
	a1= _MM_EXTRACT_EPI16(r,1);
	a2= _MM_EXTRACT_EPI16(r,2);
	a3= _MM_EXTRACT_EPI16(r,3);
	a4= _MM_EXTRACT_EPI16(r,4);
	a5= _MM_EXTRACT_EPI16(r,5);
	a6= _MM_EXTRACT_EPI16(r,6);
	a7= _MM_EXTRACT_EPI16(r,7);

fprintf(stderr, "a00 :-> %7d  %7d  %7d  %7d "
		" %7d  %7d  %7d  %7d "
		"\n"
		, a0,a1,a2,a3,a4,a5,a6,a7
		);
}
*/

static inline int extract_reg(__m128i r, int p)
{
	switch (p) {
	case 0: return(_MM_EXTRACT_EPI16(r,0));
	case 1: return(_MM_EXTRACT_EPI16(r,1));
	case 2: return(_MM_EXTRACT_EPI16(r,2));
	case 3: return(_MM_EXTRACT_EPI16(r,3));
	case 4: return(_MM_EXTRACT_EPI16(r,4));
	case 5: return(_MM_EXTRACT_EPI16(r,5));
	case 6: return(_MM_EXTRACT_EPI16(r,6));
	case 7: return(_MM_EXTRACT_EPI16(r,7));
	}
	return(0);
}


static inline __m128i insert_reg(__m128i r, int16_t v, int p)
{
	switch (p) {
	case 0: return(_MM_INSERT_EPI16(r,v,0));
	case 1: return(_MM_INSERT_EPI16(r,v,1));
	case 2: return(_MM_INSERT_EPI16(r,v,2));
	case 3: return(_MM_INSERT_EPI16(r,v,3));
	case 4: return(_MM_INSERT_EPI16(r,v,4));
	case 5: return(_MM_INSERT_EPI16(r,v,5));
	case 6: return(_MM_INSERT_EPI16(r,v,6));
	case 7: return(_MM_INSERT_EPI16(r,v,7));
	}
	return(_MM_SETZERO_SI128());
}


void sse_banded_align_lcs_and_ali_len(int16_t* seq1, int16_t* seq2, int l1, int l2, int bandLengthLeft, int bandLengthTotal, int16_t* address, double* lcs_length, int* ali_length)
{
	register int j;
	int k1, k2;
	int max, diff;
	int l_reg, l_loc;
	int line;
	int numberOfRegistersPerLine;
	int numberOfRegistersFor3Lines;

	BOOL even_line;
	BOOL odd_line;
	BOOL even_BLL;
	BOOL odd_BLL;

	um128*  SSEregisters;
	um128*  p_diag;
	um128*  p_gap1;
	um128*  p_gap2;
	um128*  p_diag_j;
	um128*  p_gap1_j;
	um128*  p_gap2_j;
	um128   current;

	um128*  l_ali_SSEregisters;
	um128*  p_l_ali_diag;
	um128*  p_l_ali_gap1;
	um128*  p_l_ali_gap2;
	um128*  p_l_ali_diag_j;
	um128*  p_l_ali_gap1_j;
	um128*  p_l_ali_gap2_j;
	um128   l_ali_current;

	um128  nucs1;
	um128  nucs2;
	um128  scores;

	um128 boolean_reg;

	// Initialisations

	odd_BLL = bandLengthLeft & 1;
	even_BLL  = !odd_BLL;

	max = INT16_MAX - l1;

	numberOfRegistersPerLine = bandLengthTotal / 8;
	numberOfRegistersFor3Lines   = 3 * numberOfRegistersPerLine;

	SSEregisters = (um128*) malloc(numberOfRegistersFor3Lines * sizeof(um128));
	l_ali_SSEregisters = (um128*) malloc(numberOfRegistersFor3Lines * sizeof(um128));

	// preparer registres SSE

	for (j=0; j<numberOfRegistersFor3Lines; j++)
	{
		(*(SSEregisters+j)).i       = _MM_SETZERO_SI128();
		(*(l_ali_SSEregisters+j)).i = _MM_LOAD_SI128(address+j*8);
	}

	p_diag    = SSEregisters;
	p_gap1    = SSEregisters+numberOfRegistersPerLine;
	p_gap2    = SSEregisters+2*numberOfRegistersPerLine;

	p_l_ali_diag    = l_ali_SSEregisters;
	p_l_ali_gap1    = l_ali_SSEregisters+numberOfRegistersPerLine;
	p_l_ali_gap2    = l_ali_SSEregisters+2*numberOfRegistersPerLine;

	// Loop on diagonals = 'lines' :
	for (line=2; line <= l1+l2; line++)
	{
		odd_line = line & 1;
		even_line  = !odd_line;

		// loop on the registers of a line :
		for (j=0; j < numberOfRegistersPerLine; j++)
		{
			p_diag_j = p_diag+j;
			p_gap1_j = p_gap1+j;
			p_gap2_j = p_gap2+j;
			p_l_ali_diag_j = p_l_ali_diag+j;
			p_l_ali_gap1_j = p_l_ali_gap1+j;
			p_l_ali_gap2_j = p_l_ali_gap2+j;

		// comparing nucleotides for diagonal scores :

			// k1 = position of the 1st nucleotide to align for seq1 and k2 = position of the 1st nucleotide to align for seq2
			if (odd_line && odd_BLL)
				k1 = (line / 2) + ((bandLengthLeft+1) / 2) - j*8;
			else
				k1 = (line / 2) + (bandLengthLeft/2) - j*8;

			k2 = line - k1 - 1;

			nucs1.i = _MM_LOADU_SI128(seq1+l1-k1);
			nucs2.i = _MM_LOADU_SI128(seq2+k2);

/*			fprintf(stderr, "\nnucs, r %d, k1 = %d, k2 = %d\n", j, k1, k2);
			printreg(nucs1.i);
			printreg(nucs2.i);
*/

		// computing diagonal score :
			scores.i = _MM_AND_SI128(_MM_CMPEQ_EPI16(nucs1.i, nucs2.i), _MM_SET1_EPI16(1));
			current.i = _MM_ADDS_EPU16((*(p_diag_j)).i, scores.i);

	// Computing alignment length

			l_ali_current.i = (*(p_l_ali_diag_j)).i;
			boolean_reg.i = _MM_CMPGT_EPI16((*(p_gap1_j)).i, current.i);
			l_ali_current.i = _MM_OR_SI128(
								_MM_AND_SI128((*(p_l_ali_gap1_j)).i, boolean_reg.i),
								_MM_ANDNOT_SI128(boolean_reg.i, l_ali_current.i));
			current.i = _MM_OR_SI128(
							_MM_AND_SI128((*(p_gap1_j)).i, boolean_reg.i),
							_MM_ANDNOT_SI128(boolean_reg.i, current.i));
			boolean_reg.i = _MM_AND_SI128(
								_MM_CMPEQ_EPI16((*(p_gap1_j)).i, current.i),
								_MM_CMPLT_EPI16((*(p_l_ali_gap1_j)).i, l_ali_current.i));
			l_ali_current.i = _MM_OR_SI128(
								_MM_AND_SI128((*(p_l_ali_gap1_j)).i, boolean_reg.i),
								_MM_ANDNOT_SI128(boolean_reg.i, l_ali_current.i));
			current.i = _MM_OR_SI128(
							_MM_AND_SI128((*(p_gap1_j)).i, boolean_reg.i),
							_MM_ANDNOT_SI128(boolean_reg.i, current.i));
			boolean_reg.i = _MM_CMPGT_EPI16((*(p_gap2_j)).i, current.i);
			l_ali_current.i = _MM_OR_SI128(
								_MM_AND_SI128((*(p_l_ali_gap2_j)).i, boolean_reg.i),
								_MM_ANDNOT_SI128(boolean_reg.i, l_ali_current.i));
			current.i = _MM_OR_SI128(
							_MM_AND_SI128((*(p_gap2_j)).i, boolean_reg.i),
							_MM_ANDNOT_SI128(boolean_reg.i, current.i));
			boolean_reg.i = _MM_AND_SI128(
								_MM_CMPEQ_EPI16((*(p_gap2_j)).i, current.i),
								_MM_CMPLT_EPI16((*(p_l_ali_gap2_j)).i, l_ali_current.i));
			l_ali_current.i = _MM_OR_SI128(
								_MM_AND_SI128((*(p_l_ali_gap2_j)).i, boolean_reg.i),
								_MM_ANDNOT_SI128(boolean_reg.i, l_ali_current.i));
			current.i = _MM_OR_SI128(
							_MM_AND_SI128((*(p_gap2_j)).i, boolean_reg.i),
							_MM_ANDNOT_SI128(boolean_reg.i, current.i));


/*			fprintf(stderr, "\nline = %d", line);
			fprintf(stderr, "\nDiag, r %d : ", j);
			printreg((*(p_diag_j)).i);
			fprintf(stderr, "Gap1      : ");
			printreg((*(p_gap1_j)).i);
			fprintf(stderr, "Gap2      : ");
			printreg((*(p_gap2_j)).i);
			fprintf(stderr, "current   : ");
			printreg(current.i);
			fprintf(stderr, "L ALI\nDiag  r %d : ", j);
			printreg((*(p_l_ali_diag_j)).i);
			fprintf(stderr, "Gap1      : ");
			printreg((*(p_l_ali_gap1_j)).i);
			fprintf(stderr, "Gap2      : ");
			printreg((*(p_l_ali_gap2_j)).i);
			fprintf(stderr, "current   : ");
			printreg(l_ali_current.i);
*/

		// diag = gap1 and gap1 = current
			(*(p_diag_j)).i = (*(p_gap1_j)).i;
			(*(p_gap1_j)).i = current.i;

		// l_ali_diag = l_ali_gap1 and l_ali_gap1 = l_ali_current+1
			(*(p_l_ali_diag_j)).i = (*(p_l_ali_gap1_j)).i;
			(*(p_l_ali_gap1_j)).i = _MM_ADD_EPI16(l_ali_current.i, _MM_SET1_EPI16(1));
		}

		// shifts for gap2, to do only once all the registers of a line have been computed     Copier gap2 puis le charger depuis la copie?

			for (j=0; j < numberOfRegistersPerLine; j++)
			{
				if ((odd_line && even_BLL) || (even_line && odd_BLL))
				{
					(*(p_gap2+j)).i = _MM_LOADU_SI128(((*(p_gap1+j)).s16)-1);
					(*(p_l_ali_gap2+j)).i = _MM_LOADU_SI128(((*(p_l_ali_gap1+j)).s16)-1);
					if (j == 0)
					{
						(*(p_gap2+j)).i = insert_reg((*(p_gap2+j)).i, 0, 0);
						(*(p_l_ali_gap2+j)).i = insert_reg((*(p_l_ali_gap2+j)).i, max, 0);
					}
				}
				else
				{
					(*(p_gap2+j)).i = _MM_LOADU_SI128(((*(p_gap1+j)).s16)+1);
					(*(p_l_ali_gap2+j)).i = _MM_LOADU_SI128(((*(p_l_ali_gap1+j)).s16)+1);
					if (j == numberOfRegistersPerLine - 1)
					{
						(*(p_gap2+j)).i = insert_reg((*(p_gap2+j)).i, 0, 7);
						(*(p_l_ali_gap2+j)).i = insert_reg((*(p_l_ali_gap2+j)).i, max, 7);
					}
				}
			}
		// end shifts for gap2

	}

/*  /// Recovering LCS and alignment lengths  \\\  */

	// finding the location of the results in the registers :
	diff = l1-l2;
	if ((diff & 1) && odd_BLL)
		l_loc = (int) floor((double)(bandLengthLeft) / (double)2) - floor((double)(diff) / (double)2);
	else
		l_loc = (int) floor((double)(bandLengthLeft) / (double)2) - ceil((double)(diff) / (double)2);

	l_reg = (int)floor((double)l_loc/(double)8.0);
	//fprintf(stderr, "\nl_reg = %d, l_loc = %d\n", l_reg, l_loc);
	l_loc = l_loc - l_reg*8;

	// extracting the results from the registers :
	*(lcs_length) = extract_reg((*(p_gap1+l_reg)).i, l_loc);
	*(ali_length) = extract_reg((*(p_l_ali_gap1+l_reg)).i, l_loc) - 1;

	// freeing the registers
	free(SSEregisters);
	free(l_ali_SSEregisters);
}


double sse_banded_align_just_lcs(int16_t* seq1, int16_t* seq2, int l1, int l2, int bandLengthLeft, int bandLengthTotal)
{
	register int j;
	int k1, k2;
	int diff;
	int l_reg, l_loc;
	int16_t l_lcs;
	int line;
	int numberOfRegistersPerLine;
	int numberOfRegistersFor3Lines;

	BOOL even_line;
	BOOL odd_line;
	BOOL even_BLL;
	BOOL odd_BLL;

	um128*  SSEregisters;
	um128*  p_diag;
	um128*  p_gap1;
	um128*  p_gap2;
	um128*  p_diag_j;
	um128*  p_gap1_j;
	um128*  p_gap2_j;
	um128   current;

	um128  nucs1;
	um128  nucs2;
	um128  scores;

	// Initialisations

	odd_BLL = bandLengthLeft & 1;
	even_BLL  = !odd_BLL;

	numberOfRegistersPerLine = bandLengthTotal / 8;
	numberOfRegistersFor3Lines   = 3 * numberOfRegistersPerLine;

	SSEregisters = malloc(numberOfRegistersFor3Lines * sizeof(um128));

	// preparer registres SSE

	for (j=0; j<numberOfRegistersFor3Lines; j++)
		(*(SSEregisters+j)).i       = _MM_SETZERO_SI128();

	p_diag    = SSEregisters;
	p_gap1    = SSEregisters+numberOfRegistersPerLine;
	p_gap2    = SSEregisters+2*numberOfRegistersPerLine;

	// Loop on diagonals = 'lines' :
	for (line=2; line <= l1+l2; line++)
	{
		odd_line = line & 1;
		even_line  = !odd_line;

		// loop on the registers of a line :
		for (j=0; j < numberOfRegistersPerLine; j++)
		{
			p_diag_j = p_diag+j;
			p_gap1_j = p_gap1+j;
			p_gap2_j = p_gap2+j;

		// comparing nucleotides for diagonal scores :

			// k1 = position of the 1st nucleotide to align for seq1 and k2 = position of the 1st nucleotide to align for seq2
			if (odd_line && odd_BLL)
				k1 = (line / 2) + ((bandLengthLeft+1) / 2) - j*8;
			else
				k1 = (line / 2) + (bandLengthLeft/2) - j*8;

			k2 = line - k1 - 1;

			nucs1.i = _MM_LOADU_SI128(seq1+l1-k1);
			nucs2.i = _MM_LOADU_SI128(seq2+k2);

		// computing diagonal score :
			scores.i = _MM_AND_SI128(_MM_CMPEQ_EPI16(nucs1.i, nucs2.i), _MM_SET1_EPI16(1));
			current.i = _MM_ADDS_EPU16((*(p_diag_j)).i, scores.i);

		// current = max(gap1, current)
			current.i = _MM_MAX_EPI16((*(p_gap1_j)).i, current.i);

		// current  = max(gap2, current)
			current.i = _MM_MAX_EPI16((*(p_gap2_j)).i, current.i);

		// diag = gap1 and gap1 = current
			(*(p_diag_j)).i = (*(p_gap1_j)).i;
			(*(p_gap1_j)).i = current.i;
		}

		// shifts for gap2, to do only once all the registers of a line have been computed

			for (j=0; j < numberOfRegistersPerLine; j++)
			{
				if ((odd_line && even_BLL) || (even_line && odd_BLL))
				{
					(*(p_gap2+j)).i = _MM_LOADU_SI128(((*(p_gap1+j)).s16)-1);
					if (j == 0)
					{
						(*(p_gap2+j)).i = insert_reg((*(p_gap2+j)).i, 0, 0);
					}
				}
				else
				{
					(*(p_gap2+j)).i = _MM_LOADU_SI128(((*(p_gap1+j)).s16)+1);
					if (j == numberOfRegistersPerLine - 1)
					{
						(*(p_gap2+j)).i = insert_reg((*(p_gap2+j)).i, 0, 7);
					}
				}
			}
		// end shifts for gap2

	}

/*  /// Recovering LCS and alignment lengths  \\\  */

	// finding the location of the results in the registers :
	diff = l1-l2;
	if ((diff & 1) && odd_BLL)
		l_loc = (int) floor((double)(bandLengthLeft) / (double)2) - floor((double)(diff) / (double)2);
	else
		l_loc = (int) floor((double)(bandLengthLeft) / (double)2) - ceil((double)(diff) / (double)2);

	l_reg = (int)floor((double)l_loc/(double)8.0);
	//fprintf(stderr, "\nl_reg = %d, l_loc = %d\n", l_reg, l_loc);
	l_loc = l_loc - l_reg*8;

	// extracting LCS from the registers :
	l_lcs = extract_reg((*(p_gap1+l_reg)).i, l_loc);

	// freeing the registers
	free(SSEregisters);

	return((double) l_lcs);
}


void calculateBandLengths(int l1, int l2, int* bandLengthRight, int* bandLengthLeft, int LCSmin)
{
	(*bandLengthLeft)  = l1 - LCSmin;
	(*bandLengthRight) = l2 - LCSmin;
}


int calculateLCSmin(int l1, int l2, double threshold, BOOL normalize, int reference, BOOL lcsmode)
{
	int LCSmin;

	if (threshold > 0)
	{
		if (normalize)
		{
			if (reference == MINLEN)
				LCSmin = threshold*l2;
			else 		// ref = maxlen or alilen
				LCSmin = threshold*l1;
		}
		else if (lcsmode)
			LCSmin = threshold;
		else if ((reference == MINLEN)) // not lcsmode
			LCSmin = l2 - threshold;
		else	// not lcsmode and ref = maxlen or alilen
			LCSmin = l1 - threshold;
	}
	else
		LCSmin = 0;

	return(LCSmin);
}


void calculateSSEBandLength(int bandLengthRight, int bandLengthLeft, int* bandLengthTotal)
{
	(*bandLengthTotal)= (double) floor(bandLengthRight + bandLengthLeft) / 2.0 + 1;

    while (((*bandLengthTotal)%8) != 0)
    	(*bandLengthTotal)+=1;
}


int calculateSizeToAllocate(int maxLen, int minLen, int LCSmin)
{
	int size;
	int notUsed;

	calculateBandLengths(maxLen, minLen, &notUsed, &size, LCSmin);		// max size = max left band length * 2

	//fprintf(stderr, "\nsize for address before %8 = %d", size);

	size = size*2;

	while ((size%8) != 0)
		size+=1;
	size*=3;
	size+=16;

	//fprintf(stderr, "\nsize for address = %d", size);

	return(size*sizeof(int16_t));
}


void iniSeq(int16_t* seq, int size, int16_t iniValue)
{
	int i;
	for (i=0; i<size; i++)
	{
		seq[i] = iniValue;
	}
}


void putSeqInSeq(int16_t* seq, char* s, int l, BOOL reverse)
{
	int i,j;

	if (reverse == FALSE)
	{
		for (i=0; i<l; i++)
			seq[i] = s[i];
	}
	else
	{
		for (i=0,j=l-1; i<l; i++,j--)
			seq[i] = s[j];
	}
}


void initializeAddressWithGaps(int16_t* address, int bandLengthTotal, int bandLengthLeft, int l1)
{
	int i;
	int address_00, x_address_10, address_01, address_01_shifted;
	int numberOfRegistersPerLine;
	int bm;

	numberOfRegistersPerLine   = bandLengthTotal / 8;
	bm = bandLengthLeft%2;

	for (i=0; i < (3*numberOfRegistersPerLine*8); i++)
		address[i] = INT16_MAX-l1;

 // 0,0 set to 1 and 0,1 and 1,0 set to 2

	address_00   = bandLengthLeft / 2;

	x_address_10 = address_00 + bm - 1;
	address_01   = numberOfRegistersPerLine*8 + x_address_10;

	address_01_shifted = numberOfRegistersPerLine*16 + address_00 - bm;

	// fill address_00, address_01,+1, address_01_shifted,+1

	address[address_00] = 1;
	address[address_01] = 2;
	address[address_01+1] = 2;
	address[address_01_shifted] = 2;
	address[address_01_shifted+1] = 2;
}


double sse_banded_lcs_align(int16_t* seq1, int16_t* seq2, int l1, int l2, BOOL normalize, int reference, BOOL lcsmode, int16_t* address, int LCSmin)
{
	double id;
	int bandLengthRight, bandLengthLeft, bandLengthTotal;
	int ali_length;

	//fprintf(stderr, "\nl1 = %d, l2 = %d\n", l1, l2);

	calculateBandLengths(l1, l2, &bandLengthRight, &bandLengthLeft, LCSmin);

	//fprintf(stderr, "\nBLL = %d, BLR = %d, LCSmin = %d\n", bandLengthLeft, bandLengthRight, LCSmin);

	calculateSSEBandLength(bandLengthRight, bandLengthLeft, &bandLengthTotal);

	//fprintf(stderr, "\nBLT = %d\n", bandLengthTotal);

	if ((reference == ALILEN) && (normalize || !lcsmode))
	{
		initializeAddressWithGaps(address, bandLengthTotal, bandLengthLeft, l1);
		sse_banded_align_lcs_and_ali_len(seq1, seq2, l1, l2, bandLengthLeft, bandLengthTotal, address, &id, &ali_length);
	}

	else
	{
		id = sse_banded_align_just_lcs(seq1, seq2, l1, l2, bandLengthLeft, bandLengthTotal);
	}
	//fprintf(stderr, "\nid before normalizations = %f", id);

	//fprintf(stderr, "\nlcs = %f, ali = %d\n", id, ali_length);

	if (!lcsmode && !normalize)
	{
		if (reference == ALILEN)
			id = ali_length - id;
		else if (reference == MAXLEN)
			id = l1 - id;
		else if (reference == MINLEN)
			id = l2 - id;
	}
	//fprintf(stderr, "\n2>>> %f, %d\n", id, ali_length);
	if (normalize)
	{
		if (reference == ALILEN)
		{
			//fprintf(stderr, "\n3>>> %f, %d\n", id, ali_length);
			id = (double) (id / (double) ali_length);
		}
		else if (reference == MAXLEN)
			id = (double) (id / (double) l1);
		else if (reference == MINLEN)
			id = (double) (id / (double) l2);
	}

	//fprintf(stderr, "\nid = %f\n", id);
	return(id);
}


double generic_sse_banded_lcs_align(char* seq1, char* seq2, double threshold, BOOL normalize, int reference, BOOL lcsmode, int16_t** address, int* buffer_size, int16_t** iseq1,
		int16_t** iseq2, int* buffer_sizeS)
{
	double id;
	int l1;
	int l2;
	int lmax, lmin;
	int sizeToAllocateForBand;
	int	maxBLL, notUsed;
	int sizeToAllocateForSeqs;
	int LCSmin;

	l1 = strlen(seq1);
	l2 = strlen(seq2);

	if (l2 > l1)
	{
		lmax = l1;
		lmin = l2;
	}
	else
	{
		lmax = l2;
		lmin = l1;
	}

	if (!lcsmode && (normalize==TRUE))
	{
		threshold = 1.0 - threshold;
	}

	LCSmin = calculateLCSmin(lmax, lmin, threshold, normalize, reference, lcsmode);

// Allocating space for matrix band if the alignment must be computed

	if ((reference == ALILEN) && ((lcsmode && normalize) || (!lcsmode))) // checking if alignment must be computed
	{
		sizeToAllocateForBand = calculateSizeToAllocate(lmax, lmin, LCSmin);

		if (sizeToAllocateForBand > (*buffer_size))
		{
			// reallocating if needed
			*address = reallocA16Address(*address, sizeToAllocateForBand);
		}
	}

// Allocating space for the int16_t arrays representing the sequences

	calculateBandLengths(lmax, lmin, &notUsed, &maxBLL, LCSmin);

	sizeToAllocateForSeqs = 2*maxBLL+lmax;

	if (sizeToAllocateForSeqs > *buffer_sizeS)
	{
		(*(iseq1)) = realloc((*(iseq1)), sizeToAllocateForSeqs*sizeof(int16_t));
		(*(iseq2)) = realloc((*(iseq2)), sizeToAllocateForSeqs*sizeof(int16_t));
	}

	iniSeq(*(iseq1), maxBLL, 0);
	iniSeq(*(iseq2), maxBLL, 255);
	*(iseq1) = *(iseq1)+maxBLL;
	*(iseq2) = *(iseq2)+maxBLL;

	// longest seq must be first argument of sse_align function
	if (l2 > l1)
	{
		putSeqInSeq((*(iseq1)), seq2, l2, TRUE);
		putSeqInSeq((*(iseq2)), seq1, l1, FALSE);
		id = sse_banded_lcs_align(*(iseq1), *(iseq2), l2, l1, normalize, reference, lcsmode, *address, LCSmin);
	}
	else
	{
		putSeqInSeq((*(iseq1)), seq1, l1, TRUE);
		putSeqInSeq((*(iseq2)), seq2, l2, FALSE);
		id = sse_banded_lcs_align(*(iseq1), *(iseq2), l1, l2, normalize, reference, lcsmode, *address, LCSmin);
	}

	return(id);
}


int prepareTablesForSumathings(int lmax, int lmin, double threshold, BOOL normalize, int reference, BOOL lcsmode,
		int16_t** address, int16_t** iseq1, int16_t** iseq2)
{
	int sizeToAllocateForBand;
	int	maxBLL;
	int notUsed;
	int sizeToAllocateForSeqs;
	int LCSmin;

	LCSmin = calculateLCSmin(lmax, lmin, threshold, normalize, reference, lcsmode);

	// Allocating space for matrix band if the alignment must be computed

	if ((reference == ALILEN) && (normalize || !lcsmode)) // checking if alignment must be computed
	{
		sizeToAllocateForBand = calculateSizeToAllocate(lmax, lmin, LCSmin);
		(*(address)) = getA16Address(sizeToAllocateForBand);
	}

	// Allocating space for the int16_t arrays representing the sequences

	calculateBandLengths(lmax, lmin, &notUsed, &maxBLL, LCSmin);

	sizeToAllocateForSeqs = 2*maxBLL+lmax;
	(*(iseq1)) = malloc(sizeToAllocateForSeqs*sizeof(int16_t));
	(*(iseq2)) = malloc(sizeToAllocateForSeqs*sizeof(int16_t));

	iniSeq(*(iseq1), maxBLL, 0);
	iniSeq(*(iseq2), maxBLL, 255);
	*(iseq1) = *(iseq1)+maxBLL;
	*(iseq2) = *(iseq2)+maxBLL;

	return(maxBLL+lmax);
}


double alignForSumathings(char* seq1, int16_t* iseq1, char* seq2, int16_t* iseq2, int l1, int l2,
		BOOL normalize, int reference, BOOL lcsmode, int16_t* address, int sizeForSeqs, int LCSmin)
{
	double id;

	iniSeq(iseq1, sizeForSeqs, 0);
	iniSeq(iseq2, sizeForSeqs, 255);

	if (l2 > l1)
	{
		putSeqInSeq(iseq1, seq2, l2, TRUE);
		putSeqInSeq(iseq2, seq1, l1, FALSE);
		id = sse_banded_lcs_align(iseq1, iseq2, l2, l1, normalize, reference, lcsmode, address, LCSmin);
	}
	else
	{
		putSeqInSeq(iseq1, seq1, l1, TRUE);
		putSeqInSeq(iseq2, seq2, l2, FALSE);
		id = sse_banded_lcs_align(iseq1, iseq2, l1, l2, normalize, reference, lcsmode, address, LCSmin);
	}

	return(id);
}

