/*
 * mtcompare.c
 *
 *  Created on: 17 aožt 2010
 *      Authors: Eric Coissac, Celine Mercier
 */

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sumatra.h"
#include "../libfasta/sequence.h"
#include "../libutils/utilities.h"
#include "../liblcs/upperband.h"
#include "../liblcs/sse_banded_LCS_alignment.h"


typedef struct {
	int32_t         tocompute;    // next line to compute
    int64_t         computed;     // count of alignment computed
    int32_t         thread_count; // count of active thread
    int32_t         tokill;
    pthread_t        *threads;
	pthread_mutex_t m_tocompute;
	pthread_mutex_t m_thread_count;
	pthread_mutex_t m_incharge;
	pthread_mutex_t m_print;
	pthread_mutex_t m_tokill;
	pthread_mutex_t m_finished;

    int32_t         lineset;
    int64_t 		pairs;
	fastaSeqCount   *db1;
	fastaSeqCount   *db2;
	BOOL            normalize;
	BOOL            lcsmode;
	int             reference;
	BOOL            extradata;
	BOOL            always;
	double			threshold;

	int16_t**		addresses;
	int16_t**		iseqs1;
	int16_t**		iseqs2;
	int				sizeForSeqs;
} thread_control_t;


void printLine(fastaSeqPtr seq1, fastaSeqPtr seq2, int32_t seqcount, double *score, BOOL extradata, int64_t pairs)
{
	int j;

	for (j=0; j < seqcount; j++,seq2++,score++)
		printResults(seq1,seq2,*score,extradata,pairs,(*score >=0));
}


void *computeline(void* c)
{
    thread_control_t *control=(thread_control_t*)c;
    int32_t line;
    int32_t start;
    int32_t end;
    int32_t threadid;
    int32_t i,j;
    double  *score;
	double  *scores;
	fastaSeqPtr   db2;
	fastaSeqPtr   db1;
	fastaSeqPtr   sq2;
	fastaSeqPtr   sq1;
	char*	s1;
	int	l1;
	int LCSmin;

    threadid = control->thread_count++;
    pthread_mutex_unlock(&(control->m_thread_count));

    db1 = control->db1->fastaSeqs;
    if (control->db2)
    {
    	end = control->db2->count;
    	db2 = control->db2->fastaSeqs;
    }
    else
	{
    	end = control->db1->count;
    	db2 = control->db1->fastaSeqs;
    }

    scores = (double*)malloc(end * sizeof(double));

    while (control->tocompute < control->db1->count)
    {
    	pthread_mutex_lock(&(control->m_tocompute));
    	line = control->tocompute;
    	pthread_mutex_unlock(&(control->m_incharge));

    	if (line < control->db1->count)
    	{
    		for (i=0,sq1=db1+line; i < control->lineset && (i+line) < control->db1->count; i++,sq1++)
    		{
    			s1 = sq1->sequence;
    			l1 = sq1->length;
				start = (control->db2) ? 0:line+i+1;

				for (score=scores,j=start; j < end ; j++,score++)
				{
					sq2=db2+j;
					filtersSumatra(sq1, sq2, control->threshold, control->normalize, control->reference, control->lcsmode, score, &LCSmin);
					if (control->always || (*score == -1.0))
					{
						*score = alignForSumathings(s1, *((control->iseqs1)+threadid), sq2->sequence, *((control->iseqs2)+threadid), l1, sq2->length,
								control->normalize, control->reference, control->lcsmode, *((control->addresses)+threadid), control->sizeForSeqs, LCSmin);

						if (!control->always && (((*score < control->threshold) && (control->lcsmode || control->normalize)) || ((!control->lcsmode && !control->normalize) && (*score > control->threshold))))
							*score = -1.0;
						else if (!control->lcsmode && control->normalize)
							*score = 1.0 - *score;
					}
					else if (*score == -2.0)
						*score = -1.0;
				}

				pthread_mutex_lock(&(control->m_print));
				printLine(sq1,db2+start,end-start,scores,control->extradata,control->pairs);
				pthread_mutex_unlock(&(control->m_print));
    		}
    	}
    }

	pthread_mutex_unlock(&(control->m_incharge));
    free(scores);
    pthread_mutex_lock(&(control->m_finished));
    control->tokill=threadid;
	pthread_mutex_unlock(&(control->m_tokill));

    return (void*)threadid;
}


int mt_compare_sumatra(fastaSeqCount *db1, fastaSeqCount *db2, double threshold, BOOL normalize, int reference, BOOL lcsmode, BOOL extradata, int n)
{
	int64_t  pairs;
    thread_control_t control;
	int32_t  i;
	int	lmax, lmax1;
	int lmin, lmin1;

	if (db2==NULL)
	{
		fprintf(stderr,"Pairwise alignment of one database against itself\n");
		pairs = (int64_t)(db1->count - 1) * (int64_t)db1->count /2;
	}
	else
	{
		fprintf(stderr,"Pairwise alignment of two databases\n");
		pairs = (int64_t)db1->count * (int64_t)db2->count;
	}

	fprintf(stderr,"Count of alignment to do : %lld\n",pairs);

    control.addresses  = (int16_t**) malloc(n*sizeof(int16_t*));
    control.iseqs1  = (int16_t**) malloc(n*sizeof(int16_t*));
    control.iseqs2  = (int16_t**) malloc(n*sizeof(int16_t*));

 	calculateMaxAndMinLenDB(*db1, &lmax, &lmin);

	if (!(db2==NULL))
	{
		calculateMaxAndMinLenDB(*db2, &lmax1, &lmin1);
		if (lmax1 > lmax)
			lmax = lmax1;
		if (lmin1 < lmin)
			lmin = lmin1;
	}

    for (i=0; i < n; i++)
    	control.sizeForSeqs = prepareTablesForSumathings(lmax, lmin, threshold, normalize, reference, lcsmode, (control.addresses)+i, (control.iseqs1)+i, (control.iseqs2)+i);

	control.db1 = db1;
	control.db2 = db2;
	control.tocompute = 0;
	control.normalize = normalize;
	control.reference = reference;
	control.extradata = extradata;
	control.threshold = threshold;
	control.lcsmode = lcsmode;
	control.computed = 0;
	control.thread_count = 0;
	control.pairs = pairs;

	if (n > control.db1->count/2)
		n = control.db1->count/2;

	control.lineset = control.db1->count / n / 2;

	if (threshold > 0)
	{
		fprintf(stderr,"Compute exact LCS only for score > %lf\n", threshold);
		control.always = FALSE;
	}

	else
		control.always=TRUE;

	if (pthread_mutex_init(&(control.m_thread_count),NULL))
	{
		fprintf(stderr,"m_thread_count mutex init error\n");
		exit(1);
	}

	if (pthread_mutex_init(&(control.m_tocompute),NULL))
	{
		fprintf(stderr,"m_tocompute mutex init error\n");
		exit(1);
	}

	if (pthread_mutex_init(&(control.m_incharge),NULL))
	{
		fprintf(stderr,"m_incharge mutex init error\n");
		exit(1);
	}

	if (pthread_mutex_init(&(control.m_print),NULL))
	{
		fprintf(stderr,"m_print mutex init error\n");
		exit(1);
	}

	if (pthread_mutex_init(&(control.m_tokill),NULL))
	{
		fprintf(stderr,"m_tokill mutex init error\n");
		exit(1);
	}

	if (pthread_mutex_init(&(control.m_finished),NULL))
	{
		fprintf(stderr,"m_finished mutex init error\n");
		exit(1);
	}

	control.threads = (pthread_t*)malloc(n * sizeof(pthread_t));

	if (!control.threads)
	{
		fprintf(stderr,"Cannot allocate memory for threads\n");
		exit(2);
	}

	pthread_mutex_lock(&(control.m_thread_count));
	pthread_mutex_lock(&(control.m_tocompute));
	pthread_mutex_lock(&(control.m_incharge));
	pthread_mutex_unlock(&(control.m_print));
	pthread_mutex_lock(&(control.m_tokill));
	pthread_mutex_unlock(&(control.m_finished));

	fprintf(stderr,"\n");

	for (i=0; i < n; i++)
	{
		fprintf(stderr,"Initializing thread...");
		if (pthread_create(control.threads+i,NULL,computeline,&control))
		{
			fprintf(stderr," : thread %d Error\n",i);
			exit(3);
		}

		pthread_mutex_lock(&(control.m_thread_count));
		fprintf(stderr," : thread %d Ok\n",control.thread_count);
	}

	pthread_mutex_unlock(&(control.m_thread_count));

	for (control.tocompute=0;
		 control.tocompute < db1->count+1;
		 control.tocompute+=control.lineset)
	{
		pthread_mutex_unlock(&(control.m_tocompute));
		pthread_mutex_lock(&(control.m_incharge));
	}

	pthread_mutex_unlock(&(control.m_tocompute));

	fprintf(stderr,"\n");

	for (i=0; i < n; i++)
	{
		pthread_mutex_lock(&(control.m_tokill));
		fprintf(stderr,"Joining thread %d...",control.tokill);
		pthread_mutex_unlock(&(control.m_tocompute));
		pthread_join(control.threads[control.tokill],NULL);
		fprintf(stderr," : Ok\n");
		pthread_mutex_unlock(&(control.m_finished));
	}

	// Freeing

	for (i=0; i < n; i++)
	{
		free((*((control.iseqs1)+i))-(control.sizeForSeqs)+lmax);
		free((*((control.iseqs2)+i))-(control.sizeForSeqs)+lmax);
	}

	free(control.iseqs1);
	free(control.iseqs2);

	if ((reference == ALILEN) && ((lcsmode && normalize) || (!lcsmode)))
	{
		for (i=0; i < n; i++)
			free(*((control.addresses)+i));
		free(control.addresses);
	}

	return 0;
}
