#include "_lcs.h"

#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <stdio.h>



// Allocate a band allowing to align sequences of length : 'length'

column_t* allocateColumn(int length,column_t *column, bool mode8bits)
{
	int size;
	bool newc = false;

			// The band length should be equal to the length
			// of the sequence + 7 for taking into account its
			// shape

	size = (length+1) * ((mode8bits) ? sizeof(int8_t):sizeof(int16_t));


			// If the pointer to the old column is NULL we allocate
			// a new column

	if (column==NULL)
	{

		column = malloc(sizeof(column_t));
		if (!column)
			return NULL;

		column->size = 0;
		column->data.shrt=NULL;
		column->score.shrt=NULL;
		newc = true;
	}

			// Otherwise we check if its size is sufficient
			// or if it should be extended

	if (size > column->size)
	{
		int16_t *old = column->data.shrt;
		int16_t *olds= column->score.shrt;

		column->data.shrt = malloc(size);
		column->score.shrt= malloc(size);

		if (column->data.shrt==NULL || column->score.shrt==NULL)
		{
			fprintf(stderr,"Allocation Error on column for a size of %d\n" , size);
			column->data.shrt = old;
			column->score.shrt= olds;

			if (newc)
			{
				free(column);
				column=NULL;
				return NULL;
			}
			return NULL;
		}
		else
			column->size = size;
	}

	return column;
}

void freeColumn(column_p column)
{
	if (column)
	{
		if (column->data.shrt)
			free(column->data.shrt);

		if (column->score.shrt)
			free(column->score.shrt);

		free(column);
	}
}

int fastLCSScore(const char* seq1, const char* seq2,column_pp column,int32_t* lpath)
{
	return fastLCSScore16(seq1,seq2,column,lpath);
}

int simpleLCS(const char* seq1, const char* seq2,column_pp ppcolumn,int32_t* lpath)
{
	int lseq1,lseq2;		// length of the both sequences
	int lcs;
	int itmp;				// tmp variables for swap
	const char* stmp;		//
	int32_t *score;
	int32_t *path;
	column_t *column;
	int32_t i,j;
	int32_t sl,su,sd;
	int32_t pl,pu,pd;

		// Made seq1 the longest sequences
	lseq1=strlen(seq1);
	lseq2=strlen(seq2);

	if (lseq1 < lseq2)
	{
		itmp=lseq1;
		lseq1=lseq2;
		lseq2=itmp;

		stmp=seq1;
		seq1=seq2;
		seq2=stmp;
	}

	lseq1++;
	lseq2++;

							// a band sized to the smallest sequence is allocated

	if (ppcolumn)
		column = *ppcolumn;
	else
		column=NULL;

	column = allocateColumn(lseq1*2,column,0);
	score = (int32_t*) column->score.shrt;
	path = (int32_t*) column->data.shrt;

	memset(score,0,lseq1 * sizeof(int32_t));

	for (j=0; j < lseq1; j++)
		path[j]=j;

	for (i=1; i< lseq2; i++)
	{
		sl=0;
		pl=i;
		for (j=1; j < lseq1; j++)
		{
			sd=score[j-1] + (seq2[i-1]==seq1[j-1] ? 1:0);
			pd=path[j-1]  + 1;

			su=score[j];
			pu=path[j] + 1;

			score[j-1]=sl;

			if (su > sl) sl=su, pl=pu;
			if (sd > sl) sl=sd, pl=pd;
		}
	}

	lcs = sl;
	if(lpath) *lpath=pl;

	if (ppcolumn)
		*ppcolumn=column;
	else
		freeColumn(column);

	return lcs;
}

