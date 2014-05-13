/*
 * mtcompare_cumaclust.c
 *
 *  Author: Celine Mercier
 *
 */


#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include "../libfasta/sequence.h"
#include "../libutils/utilities.h"
#include "../liblcs/upperband.h"
#include "../liblcs/sse_banded_LCS_alignment.h"
#include <math.h>


typedef struct {
    int32_t			next;
    int32_t			threads_number;
    int*			potential_nexts_list;
	fastaSeqPtr*    db;
	int				n;
	int             normalize;
	int 			reference;
	BOOL			lcsmode;
	BOOL			fast;
	double			threshold;
    BOOL 			stop;
    int				sizeForSeqs;
    int16_t**		addresses;
    int16_t**		iseqs1;
    int16_t**		iseqs2;
    int				seeds_counter;
    double			worstscore;
} thread_control_t;


void computeScore(void* c, double* score, int32_t seed_number, int32_t i, int thread_number)
{
    thread_control_t *control=(thread_control_t*)c;

	*score = control->worstscore;
	int LCSmin;

    filters((*((control->db)+i)), (*((control->db)+seed_number)), control->threshold, control->normalize, control->reference, control->lcsmode, score, &LCSmin);

    if (*score == -1.0)
    {
		*score = alignForSumathings((*((control->db)+seed_number))->sequence, *((control->iseqs1)+thread_number),
				(*((control->db)+i))->sequence, *((control->iseqs2)+thread_number), (*((control->db)+seed_number))->length, (*((control->db)+i))->length,
				control->normalize, control->reference, control->lcsmode, *((control->addresses)+thread_number), control->sizeForSeqs, LCSmin);
    }
}


void putSeqInClusterMT(void *c, int32_t center_idx, int32_t seq, double score)
{
	// saves a sequence as belonging to a cluster and its score with the seed

    thread_control_t *control=(thread_control_t*)c;

    (*((control->db)+seq))->center         = (control->db)+center_idx;
    (*((control->db)+seq))->center_index   = center_idx;	// saves cluster
    (*((control->db)+seq))->score          = score;			    // saves score with the seed
	(*((control->db)+seq))->cluster_center = FALSE;
}


void getNext(void* c)
{
    thread_control_t *control=(thread_control_t*)c;
    int32_t 		i;

    control->next = control->n;
    for (i=0;i<control->threads_number;i++)
    {
    	if (((*(control->potential_nexts_list+i)) < control->next) && ((*(control->potential_nexts_list+i)) != 0))
    		control->next = (*(control->potential_nexts_list+i));
    }
    if (control->next < (control->n)-1)
    	(control->seeds_counter)++;
    else if (control->next == (control->n)-1)
	{
		control->stop = TRUE;
		(control->seeds_counter)++;
	}
	else if (control->next == control->n)
		control->stop = TRUE;
}


void computeOneSeed(void* c)
{
    thread_control_t *control=(thread_control_t*)c;
    int32_t 		i;
    BOOL 		  	found;
    int32_t current_seed;
    int32_t seed_number;
    double  score;
    int		thread_id;

    seed_number = control->next;
    found = FALSE;

    //printf("\n seed = %d, n = %d", seed_number, control->n);

    omp_set_num_threads(control->threads_number);

	#pragma omp parallel for firstprivate(found) private(score) private(current_seed) private(thread_id) schedule(static)

    	for (i=seed_number+1; i < control->n; i++)
		{
    		thread_id = omp_get_thread_num();

			current_seed = (*((control->db)+i))->center_index;
    		if ((control->fast == FALSE) || (current_seed == i))
    		{
				computeScore(c, &score, seed_number, i, thread_id);		// computes LCS score or 0 if k-mer filter not passed
    			if (((score < control->threshold) && (control->lcsmode || control->normalize)) || ((!control->lcsmode && !control->normalize) && (score > control->threshold))) // similarity under threshold
				{
					if ((found == FALSE) && (current_seed == i))
					{
						(*((control->potential_nexts_list)+thread_id)) = i;	// saves potential next seed
						found = TRUE;		// saves the fact that a next seed has been found for this thread
					}
				}
				else if ((current_seed == i) || ((control->fast == FALSE) && ((*((control->db)+i))->score < score)))
				{ // if seq matching with current seed : clustering with seed if seq doesn't belong to any cluster yet OR not in fast mode and the score is better with this seed
					if (!control->lcsmode && control->normalize)
						score = 1.0 - score;
					putSeqInClusterMT(c, seed_number, i, score);		// saves new seed for this seq
				}
    		}
		}

    getNext(c);
}


void initializeCentersAndScores(void *c)
{
	// Initializing the scores table for each seq :

    thread_control_t *control = (thread_control_t*) c;
    int32_t  i;
    int scoremax;

	if (control->normalize && control->lcsmode)
		scoremax = 1.0;
	else if (!control->lcsmode)
		scoremax = 0.0;
	else
		scoremax = (*(control->db))->length;

	for (i=0; i <= control->n-1; i++)
	{
		(*((control->db)+i))->center         = (control->db)+i;
		(*((control->db)+i))->center_index   = i;
		(*((control->db)+i))->score          = scoremax;
		(*((control->db)+i))->cluster_center = TRUE;
	}
}


void freeEverything(void *c)
{
    thread_control_t *control=(thread_control_t*)c;
    int i;

	free(control->potential_nexts_list);
	if ((control->reference == ALILEN) && (control->normalize || !control->lcsmode))
	{
		for (i=0; i < control->threads_number; i++)
			free(*((control->addresses)+i));
		free(control->addresses);
	}
	free(control->iseqs1);
	free(control->iseqs2);
}


int mt_compare_sumaclust(fastaSeqPtr* db, int n, BOOL fast, double threshold, BOOL normalize, int reference, BOOL lcsmode, int threads_number)
{
    thread_control_t control;
	int32_t  i;
	int lmax, lmin;

	if (lcsmode || normalize)
		fprintf(stderr,"Clustering sequences when similarity >= %lf\n", threshold);
	else
		fprintf(stderr,"Clustering sequences when distance <= %lf\n", threshold);

	fprintf(stderr,"Aligning and clustering... \n");

    control.threads_number = omp_get_max_threads();
    if (threads_number < control.threads_number)
    	control.threads_number = threads_number;

	calculateMaxAndMinLen(db, n, &lmax, &lmin);

    control.addresses  = (int16_t**) malloc(control.threads_number*sizeof(int16_t*));
    control.iseqs1  = (int16_t**) malloc(control.threads_number*sizeof(int16_t*));
    control.iseqs2  = (int16_t**) malloc(control.threads_number*sizeof(int16_t*));

    for (i=0; i < control.threads_number; i++)
    	control.sizeForSeqs = prepareTablesForSumathings(lmax, lmin, threshold, normalize, reference, lcsmode, (control.addresses)+i, (control.iseqs1)+i, (control.iseqs2)+i);

	control.db = db;
	control.next = 0;
    control.normalize = normalize;
    control.reference = reference;
    control.threshold = threshold;
    control.lcsmode = lcsmode;
    control.stop = FALSE;
    control.fast = fast;
    control.seeds_counter = 1;
    control.potential_nexts_list = (int*) calloc(control.threads_number, sizeof(int));
    control.n = n;

    if (lcsmode || normalize)
		control.worstscore = 0.0;
	else
		control.worstscore = lmax;

    fprintf(stderr, "%d threads running\n", control.threads_number);

    // initialize scores table :
    initializeCentersAndScores(&control);

	while (control.stop == FALSE)
	{
		if ((control.next)%10 == 0)
		{
			float p = ((float)(control.next)/(float)n)*100;
			fprintf(stderr,"\rDone : %f %%      ",p);
		}
		computeOneSeed(&control);

		for (i=0;i<control.threads_number;i++)
			(*((control.potential_nexts_list)+i)) = 0;
	}

	for (i=0; i < control.threads_number; i++)
	{
		free((*((control.iseqs1)+i))-(control.sizeForSeqs)+lmax);
		free((*((control.iseqs2)+i))-(control.sizeForSeqs)+lmax);
	}

	freeEverything(&control);
    fprintf(stderr,"\rDone : 100 %%       %d clusters created.                        \n\n", control.seeds_counter);
    return(control.seeds_counter);
}
