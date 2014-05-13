/**
 * FileName:    sumatra.c
 * Authors:      Eric Coissac, Celine Mercier
 * Description: computation of pairwise similarities of DNA sequences
 * **/

#include "sumatra.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#include "../libfasta/sequence.h"
#include "../liblcs/upperband.h"
#include "../liblcs/sse_banded_LCS_alignment.h"
#include "../libutils/utilities.h"
#include "mtcompare_sumatra.h"

#define VERSION "1.0"


/* ----------------------------------------------- */
/* printout help                                   */
/* ----------------------------------------------- */
#define PP fprintf(stdout,

static void PrintHelp()
{
        PP      "-----------------------------------------------------------------------------------------------------------------------------\n");
        PP      " SUMATRA Version %s\n", VERSION);
        PP      "-----------------------------------------------------------------------------------------------------------------------------\n");
        PP      " Synopsis : sumatra computes all the pairwise LCS (Longest Common Subsequence) scores\n");
        PP		" of one nucleotide dataset or between two nucleotide datasets.\n");
        PP      " Usage: sumatra [options] <dataset1> [dataset2]\n");
        PP      "-----------------------------------------------------------------------------------------------------------------------------\n");
        PP      " Options:\n\n");
        PP      " -h       : [H]elp - print <this> help\n\n");
        PP      " -l       : Reference sequence length is the shortest. \n\n");
        PP      " -L       : Reference sequence length is the largest. \n\n");
        PP      " -a       : Reference sequence length is the alignment length (default). \n\n");
        PP      " -n       : Score is normalized by reference sequence length (default).\n\n");
        PP      " -r       : Raw score, not normalized. \n\n");
        PP      " -d       : Score is expressed in distance (default: score is expressed in similarity). \n\n");
        PP      " -t ##.## : Score threshold. If the score is normalized and expressed in similarity (default),\n");
        PP		"            it is an identity, e.g. 0.95 for an identity of 95%%. If the score is normalized\n");
        PP		"            and expressed in distance, it is (1.0 - identity), e.g. 0.05 for an identity of 95%%.\n");
        PP		"            If the score is not normalized and expressed in similarity, it is the length of the\n");
        PP		"            Longest Common Subsequence. If the score is not normalized and expressed in distance,\n");
        PP		"            it is (reference length - LCS length).\n");
        PP		"            Only sequence pairs with a similarity above ##.## are printed. Default: 0.00 \n");
        PP		"            (no threshold).\n\n");
        PP      " -p ##    : Number of threads used for computation (default=1).\n\n");
        PP      " -g       : n's are replaced with a's (default: sequences with n's are discarded).\n");
        PP      " -x       : Adds four extra columns with the count and length of both sequences.\n");
        PP      "-----------------------------------------------------------------------------------------------------------------------------\n");
        PP      " First argument  : the nucleotide dataset to analyze\n\n");
        PP      " Second argument : optionally the second nucleotide dataset\n");
        PP      "-----------------------------------------------------------------------------------------------------------------------------\n");
        PP      " Results table description : \n");
        PP      " column 1 : Identifier sequence 1\n");
        PP      " column 2 : Identifier sequence 2\n");
        PP      " column 3 : Score\n");
        PP      " column 4 : Count of sequence 1  (only with option -x)\n");
        PP      " column 5 : Count of sequence 2  (only with option -x)\n");
        PP      " column 6 : Length of sequence 1 (only with option -x)\n");
        PP      " column 7 : Length of sequence 2 (only with option -x)\n");
        PP      "-----------------------------------------------------------------------------------------------------------------------------\n");
        PP		" http://metabarcoding.org/sumatra\n");
        PP      "-----------------------------------------------------------------------------------------------------------------------------\n\n");
}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit                         */
/* ----------------------------------------------- */

#define PP fprintf(stderr,

static void ExitUsage(stat)
        int stat;
{
        PP      "usage: sumatra [-l|L|a|n|r|d|g|x] [-t threshold_value] [-p number of threads] dataset1 [dataset2]\n");
        PP      "type \"sumatra -h\" for help\n");

        if (stat)
            exit(stat);
}

#undef  PP


void printResults(fastaSeqPtr seq1,fastaSeqPtr seq2,
			      double score,
			      BOOL extradata,
			      int64_t pairs,
			      BOOL print)
{
	static struct timeval start;
	static struct timeval lastprint;
	static BOOL first=TRUE;
	static uint64_t aligned=0;

	struct timeval current;
	double  fraction;
	time_t fulltime;
	time_t remaintime;
	double elapsedtime;
	int32_t day;
	int32_t hour;
	int32_t minute;
	int32_t seconde;


	aligned++;

	if (first)
	{
		first=FALSE;
		gettimeofday(&start,NULL);
		lastprint=start;
	}

	gettimeofday(&current,NULL);

	if (current.tv_sec!=lastprint.tv_sec)
	{
		lastprint=current;
		fraction = (double)aligned/(double)pairs;
		elapsedtime = difftime(current.tv_sec,start.tv_sec);
		fulltime = elapsedtime / fraction;
		remaintime = (time_t)difftime(fulltime,(time_t)elapsedtime);


		fprintf(stderr,
				"Computed %lld / %lld -> %5.2lf%%",
				aligned, pairs, fraction*100.
				);
		seconde = fulltime % 60;
		minute  = fulltime / 60;
		hour    = minute / 60;
		minute  = minute % 60;
		day     = hour / 24;
		hour    = hour % 24;
		if (day)
		fprintf(stderr,
				", estimated computation time = %3d days %02d:%02d:%02d",
				day,
				hour,
				minute,
				seconde
				);
		else
			fprintf(stderr,
					", estimated computation time = %02d:%02d:%02d",
					hour,
					minute,
					seconde
					);

		seconde = remaintime % 60;
		minute  = remaintime / 60;
		hour    = minute / 60;
		minute  = minute % 60;
		day     = hour / 24;
		hour    = hour % 24;
		if (day)
		fprintf(stderr,
				", about %3d days %02d:%02d:%02d remaining                  \r",
				day,
				hour,
				minute,
				seconde
				);
		else
			fprintf(stderr,
					", about %02d:%02d:%02d remaining                  \r",

					hour,
					minute,
					seconde
					);

	}

	if (print)
	{

	if (extradata)
		printf("%s\t%s\t%lf\t%d\t%d\t%d\t%d\n", seq1->accession_id,
										seq2->accession_id,
										score,
										seq1->count,
										seq2->count,
										seq1->length,
										seq2->length
		      );
	else
		printf("%s\t%s\t%lf\n", seq1->accession_id,
				seq2->accession_id,
				score);
	}
}


int compare1(fastaSeqCount db1, double threshold, BOOL normalize, int reference, BOOL lcsmode, BOOL extradata)
{
	BOOL     always = TRUE;
	int64_t  pairs = (int64_t)(db1.count - 1) * (int64_t)db1.count /2;
	BOOL     print;
	double   score;
	int32_t  i,j;
	char* s1;
	char* s2;
	int l1;
	int l2;
	int16_t* iseq1;
	int16_t* iseq2;
	int16_t* address;
	int sizeForSeqs;
	int lmax, lmin;
	int LCSmin;

	fprintf(stderr,"Pairwise alignments of one dataset against itself\n");
	fprintf(stderr,"Count of alignments to do : %lld\n",pairs);

	if (threshold > 0)
	{
		fprintf(stderr,"Computing for scores > %lf\n",threshold);
		always = FALSE;
	}

	calculateMaxAndMinLenDB(db1, &lmax, &lmin);
	sizeForSeqs = prepareTablesForSumathings(lmax, lmin, threshold, normalize, reference, lcsmode, &address, &iseq1, &iseq2);

	for (i=0; i < db1.count; i++) // ...??
		for (j=i+1; j < db1.count; j++)
		{
			print = FALSE;
			filtersSumatra(db1.fastaSeqs+i, db1.fastaSeqs+j, threshold, normalize, reference, lcsmode, &score, &LCSmin);
			if (score >= 0) // identical sequences
				print = TRUE;
			else if (always || (score == -1.0))	// filter passed or no threshold
			{
				s1 = (db1.fastaSeqs+i)->sequence;
				l1 = (db1.fastaSeqs+i)->length;
				s2 = (db1.fastaSeqs+j)->sequence;
				l2 = (db1.fastaSeqs+j)->length;
				score = alignForSumathings(s1, iseq1, s2, iseq2, l1, l2, normalize, reference, lcsmode, address, sizeForSeqs, LCSmin);
				print = always || (((normalize || lcsmode) && (score >= threshold)) || ((!lcsmode && !normalize) && (score <= threshold)));
				if (print && !lcsmode && normalize)
					score = 1.0 - score;
			}
			printResults(db1.fastaSeqs+i, db1.fastaSeqs+j, score, extradata, pairs, print);
		}

	free(iseq1-sizeForSeqs+lmax);
	free(iseq2-sizeForSeqs+lmax);
	return 0;
}


int compare2(fastaSeqCount db1, fastaSeqCount db2, double threshold, BOOL normalize, int reference, BOOL lcsmode, BOOL extradata)
{
	BOOL     always = TRUE;
	int64_t  pairs = (int64_t)(db1.count) * (int64_t)(db2.count);
	BOOL     print;
	double   score;
	int32_t  i,j;
	char* s1;
	char* s2;
	int l1;
	int l2;
	int16_t* iseq1;
	int16_t* iseq2;
	int16_t* address;
	int sizeForSeqs;
	int lmax;
	int lmax1;
	int lmin;
	int lmin1;
	int LCSmin;

	fprintf(stderr,"Pairwise alignments of two datasets\n");
	fprintf(stderr,"Count of alignments to do : %lld\n",pairs);

	if (threshold > 0)
	{
		fprintf(stderr,"Computing for scores > %lf\n",threshold);
		always = FALSE;
	}

	calculateMaxAndMinLenDB(db1, &lmax, &lmin);
	calculateMaxAndMinLenDB(db2, &lmax1, &lmin1);

	if (lmax1 > lmax)
		lmax = lmax1;
	if (lmin1 < lmin)
		lmin = lmin1;

	sizeForSeqs = prepareTablesForSumathings(lmax, lmin, threshold, normalize, reference, lcsmode, &address, &iseq1, &iseq2);

	for (i=0; i < db1.count; i++)
		for (j=0; j < db2.count; j++)
		{
			print = FALSE;
			filtersSumatra(db1.fastaSeqs+i, db2.fastaSeqs+j, threshold, normalize, reference, lcsmode, &score, &LCSmin);
			if (score >= 0) // identical sequences
				print = TRUE;
			else if (always || (score == -1.0))	// filter passed or no threshold
			{
				s1 = (db1.fastaSeqs+i)->sequence;
				l1 = (db1.fastaSeqs+i)->length;
				s2 = (db2.fastaSeqs+j)->sequence;
				l2 = (db2.fastaSeqs+j)->length;
				score = alignForSumathings(s1, iseq1, s2, iseq2, l1, l2, normalize, reference, lcsmode, address, sizeForSeqs, LCSmin);
				print = always || (((normalize || lcsmode) && (score >= threshold)) || ((!lcsmode && !normalize) && (score <= threshold)));
				if (print && !lcsmode && normalize)
					score = 1.0 - score;
			}
			printResults(db1.fastaSeqs+i, db2.fastaSeqs+j, score, extradata, pairs, print);
		}

//	for (i=0; i < db1.count; i++)
//	{
//		print = FALSE;
//		filtersSumatra(db1.fastaSeqs+i, db2.fastaSeqs+i, threshold, normalize, reference, lcsmode, &score, &LCSmin);
//		if (score >= 0) // identical sequences
//			print = TRUE;
//		else if (always || (score == -1.0))	// filter passed or no threshold
//		{
//			s1 = (db1.fastaSeqs+i)->sequence;
//			l1 = (db1.fastaSeqs+i)->length;
//			s2 = (db2.fastaSeqs+i)->sequence;
//			l2 = (db2.fastaSeqs+i)->length;
//			score = alignForSumathings(s1, iseq1, s2, iseq2, l1, l2, normalize, reference, lcsmode, address, sizeForSeqs, LCSmin);
//			print = always || (((normalize || lcsmode) && (score >= threshold)) || ((!lcsmode && !normalize) && (score <= threshold)));
//			if (print && !lcsmode && normalize)
//				score = 1.0 - score;
//		}
//		printResults(db1.fastaSeqs+i, db2.fastaSeqs+i, score, extradata, pairs, print);
//	}

	free(iseq1-sizeForSeqs+lmax);
	free(iseq2-sizeForSeqs+lmax);
	return 0;
}


int main(int argc, char **argv)
{

	int32_t     	carg		= 0;
	int32_t         errflag     = 0;
	BOOL            normalize   = TRUE;
	BOOL            lcsmode     = TRUE;
	int             reference   = ALILEN;
	BOOL            extradata   = FALSE;
	BOOL			onlyATGC	   = TRUE;
	double          threshold   = 0.0;
	int             ndb         = 0;
	int				nproc       = 1;
	fastaSeqCount   db1;
	fastaSeqCount   db2;

	while ((carg = getopt(argc, argv, "hlLanrdt:p:gx")) != -1) {
	    switch (carg) {
								   /* -------------------- */
	    case 'h':                  /* help                 */
								   /* -------------------- */
	    	PrintHelp();
	    	exit(0);
	    	break;

									  /* -------------------------------------------------- */
	    case 'l':               	  /* Normalize LCS/Error by the shortest sequence length*/
									  /* -------------------------------------------------- */
	    	reference=MINLEN;
	        break;

									  /* -------------------------------------------------- */
	    case 'L':               	  /* Normalize LCS/Error by the largest sequence length */
									  /* -------------------------------------------------- */
	    	reference=MAXLEN;
	        break;

	          	  	  	  	  	  	  /* -------------------------------------------------- */
	    case 'a':               	  /* Normalize LCS/Error by the alignment length        */
	    	   	   	   	   	   	   	  /* -------------------------------------------------- */
	    	reference=ALILEN;
	    	break;

									  /* -------------------------------------------------- */
	    case 'n':               	  /* Normalize LCS by the reference length              */
									  /* -------------------------------------------------- */
	    	normalize=TRUE;
	    	break;

	    								/* -------------------------------------------------- */
		case 'r':   					/* No normalization     				              */
					  	  	  	  	  	/* -------------------------------------------------- */
			normalize=FALSE;
			break;

									  /* -------------------------------------------------- */
	    case 'd':               	  /* Score is expressed in distance                  */
								      /* -------------------------------------------------- */
	    	lcsmode=FALSE;
	    	break;

						/* ------------------------------------------------------------------- */
		case 't':   	/* Prints only pairs with similarity higher than (threshold)         */
						/* ------------------------------------------------------------------- */
	    	sscanf(optarg,"%lf",&threshold);
	    	break;

							 /* -------------------------------------------------- */
	    case 'p':            /* number of processors to use                        */
							 /* -------------------------------------------------- */
	    	sscanf(optarg,"%d",&nproc);
			break;

								  /* -------------------------------------------------- */
	    case 'x':              	  /* Print extra data (node weight, lseq1, lseq2)       */
								  /* -------------------------------------------------- */
	    	extradata=TRUE;
			break;


									/* -------------------------------------------------- */
		case 'g':   				/* replace n's with a's in sequences                  */
									/* -------------------------------------------------- */
			onlyATGC=FALSE;
			break;


	    case '?':               	/* bad option   	        */
	    	errflag++;
	    	break;
	     }
	}

	ndb = argc - optind;
	if (ndb < 1)
        errflag++;

	if (errflag)
		ExitUsage(errflag);

    fprintf(stderr,"===============================================================\n");
	fprintf(stderr," SUMATRA version %s\n",VERSION);
#ifdef __SSE2__
	fprintf(stderr,"Alignment using SSE2 code\n");
#else
	fprintf(stderr,"Alignment using standard code, SSE2 unavailable\n");
#endif
	fprintf(stderr,"===============================================================\n");


	if (normalize && (threshold > 1.0))
	{
		fprintf(stderr, "\nERROR : Please specify a threshold included between 0.0 and 1.0 when normalizing scores.\n\n");
		exit(1);
	}

	fprintf(stderr,"Reading dataset 1...");
	db1 = seq_readAllSeq2(argv[optind], TRUE, onlyATGC);
	fprintf(stderr,"\n%d sequences\n",db1.count);

	if (!onlyATGC)
		(void)cleanDB(db1);

	if (threshold>0)
		(void)hashDB(db1);

	optind++;

	if (ndb == 2)
	{
	    fprintf(stderr,"Reading dataset 2...");
		db2 = seq_readAllSeq2(argv[optind], TRUE, onlyATGC);
		fprintf(stderr,"\n%d sequences\n",db2.count);

		if (!onlyATGC)
			(void)cleanDB(db2);

		if (threshold>0)
			(void)hashDB(db2);
	}

	if (!lcsmode && normalize && (threshold > 0))
		threshold = 1.0 - threshold;

	if (ndb==1)
	{
		if (nproc==1)
			compare1(db1, threshold, normalize, reference, lcsmode, extradata);
		else
		    mt_compare_sumatra(&db1, NULL, threshold, normalize, reference, lcsmode, extradata, nproc);
	}
	else
	{
		if (nproc==1)
			compare2(db1, db2, threshold, normalize, reference, lcsmode, extradata);
		else
		    mt_compare_sumatra(&db1, &db2, threshold, normalize, reference, lcsmode, extradata, nproc);
	}

	fprintf(stderr,"\nDone\n\n");
	return 0;

}
