/*
 * banded_LCS_alignment.c
 *
 *  Created on: 7 nov. 2012
 *      Author: merciece
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../libutils/utilities.h"


typedef struct {
	int		score;
	int		l_path;
}infos;


int calculateScore(char nuc1, char nuc2)
{
	return(nuc1 == nuc2);
}

infos** banded_align(char *seq1, char *seq2, int l1, int l2, int bandLengthRight, int bandLengthLeft)
{
	int i, j;
	//int c;
	//double id;
	int start, end;
	int diag_score, delete, insert, mismatch;
	int l_path, l_path_i, l_path_d;
	int bestScore;
	int mismatch_margin;
	int stop;
	int diag_index;
	infos **matrix;

	l1++;
	l2++;
	mismatch_margin = bandLengthLeft;		// the biggest one
	diag_index = l1-l2;						// diagonal index
	stop=0;

	//fprintf(stderr,"\nseq1 = %s, seq2=%s, bandLengthR = %d, bandLengthL = %d", seq1, seq2, bandLengthRight, bandLengthLeft);

	// Matrix initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	matrix = (infos**) malloc(l1 * sizeof(infos*));
	for (i = 0; i < l1; i++)
		matrix[i] = (infos*) malloc(l2 * sizeof(infos));

	for (i = 0; i < l1; i++)
		for (j = 0; j < l2; j++)
		{
			matrix[i][j].score  = 0;
			matrix[i][j].l_path = 0;
		}

	for (i = 0; i < l1; i++)
		matrix[i][0].l_path = i;

	for (j = 0; j < l2; j++)
		matrix[0][j].l_path = j;

	// Matrix initialized~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	for (i = 1; i < l1; i++)
	{
		start = i - bandLengthLeft;
		if (start < 1)
			start = 1;
		end = i+bandLengthRight+1;
		if (end > l2)
			end = l2;

		for (j = start; j < end; j++)
		{
			delete = matrix[i-1][j].score;
			l_path_d = matrix[i-1][j].l_path + 1;
			insert = matrix[i][j-1].score;
			l_path_i = matrix[i][j-1].l_path + 1;
			mismatch = 0;

			diag_score = calculateScore(seq1[i-1], seq2[j-1]);
			bestScore = matrix[i-1][j-1].score + diag_score;
			l_path = matrix[i-1][j-1].l_path + 1;
			if (diag_score == 0)   // mismatch
				mismatch = 1;

			if ((insert > bestScore) || ((insert == bestScore) && (l_path_i < l_path)))
			{
				bestScore = matrix[i][j-1].score;
				l_path = l_path_i;
				mismatch = 0;
			}

			if ((delete > bestScore) || ((delete == bestScore) && (l_path_d < l_path)))
			{
				bestScore = delete;
				l_path = l_path_d;
				mismatch = 0;
			}

			/*if (((i-j) - diag_index == 0) && (mismatch == 1))
			{
				//fprintf(stderr, "\nR = %d, L = %d\n", bandLengthRight, bandLengthLeft);
				if (bandLengthRight+bandLengthLeft == 0)
				{
					stop = 1;
					//fprintf(stderr, "\nBREAKING LOOPS\n");
					break;
				}
				if (bandLengthRight != 0)
					bandLengthRight = bandLengthRight - 1;
				if (bandLengthLeft != 0)
					bandLengthLeft = bandLengthLeft - 1;
			}*/

			(matrix[i][j]).score = bestScore;
			(matrix[i][j]).l_path = l_path;
		}

		//if ((bandLengthRight + bandLengthLeft == 0) && ((matrix[i][j].l_path - matrix[i][j].score) > mismatch_margin))
		if (stop==1)
			break;
	}
	return(matrix);
}


void calculateBandLength(int l1, int l2, double threshold, int* bandLengthRight, int* bandLengthLeft)
{
	(*bandLengthLeft)  = round(-l1 * threshold + l1);
	(*bandLengthRight) = round(-l1 * threshold + l2);

//	fprintf(stderr,"\nR=%d, L=%d", (*bandLengthRight), (*bandLengthLeft));
}


double calculateId(infos** matrix, int len1, int len2)
{
	double id;
	int l_ali;
	int l_lcs;

	l_lcs = matrix[len1][len2].score;
	l_ali = matrix[len1][len2].l_path;

	if (l_lcs == 0)
		id = 0.0;
	else
		id = (double) l_lcs / (double) l_ali;

	//fprintf(stderr, "\n%d, %d\n", l_lcs, l_ali);
	return(id);
}


double banded_lcs_align(int16_t* seq1, int16_t* seq2, int l1, int l2, double threshold, BOOL n, int ref, BOOL lcsmode, int16_t* address)
{
	double id;
	int bandLengthRight, bandLengthLeft;
	int i,j;

	char* s1;
	char* s2;

	s1 = (char*) malloc(l1*sizeof(char)+1);
	s2 = (char*) malloc(l2*sizeof(char)+1);

	for (i=l1-1, j=0; i>=0, j<l1; i--, j++)
		*(s1+i) = (char) *(seq1+j);

	for (i=0; i<l2; i++)
		*(s2+i) = (char) *(seq2+i);

	*(s1+l1) = 0;
	*(s2+l2) = 0;

	//fprintf(stderr, "\nl1=%d, %s\nl2=%d, %s\n", l1, s1, l2, s2);

	infos** matrix;

	calculateBandLength(l1, l2, threshold, &bandLengthRight, &bandLengthLeft);

	matrix = banded_align(s1, s2, l1, l2, bandLengthRight, bandLengthLeft);

	/*fprintf(stderr, "\n");
	for (i = 0; i <= l1; i++)
	{
		fprintf(stderr, "\n");
		for (j = 0; j <= l2; j++)
			fprintf(stderr, "%d/%d\t", matrix[i][j].score, matrix[i][j].l_path); //matrix[i][j].stop);
	}
	fprintf(stderr, "\n");*/

	id = calculateId(matrix, l1, l2);

	for (i = 0; i <= l1; i++)
		free(matrix[i]);

	free(matrix);

	free(s1);
	free(s2);

	//fprintf(stderr, "\nscore = %lf\n", id);

	return(id);
}
