/**
 * FileName:    sumaclust.c
 * Author:      Celine Mercier
 * Description: star clustering of DNA sequences
 * **/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#include "../libutils/utilities.h"
#include "../libfasta/sequence.h"
#include "../libfasta/fasta_header_parser.h"
#include "../libfasta/fasta_header_handler.h"
#include "../libfasta/fasta_seq_writer.h"
#include "mtcompare_sumaclust.h"
#include "../liblcs/upperband.h"
#include "../liblcs/sse_banded_LCS_alignment.h"
#include "sumaclust.h"

#define VERSION "1.0"


/* ----------------------------------------------- */
/* printout help                                   */
/* ----------------------------------------------- */

#define PP fprintf(stdout,


static void PrintHelp()
{
        PP      "------------------------------------------------------------\n");
		PP      " SUMACLUST Version %s\n", VERSION);
        PP      "------------------------------------------------------------\n");
        PP      " Synopsis : star clustering of sequences.\n");
        PP      " Usage: sumaclust [options] <dataset>\n");
        PP      "------------------------------------------------------------\n");
        PP      " Options:\n");
        PP      " -h       : [H]elp - print <this> help\n\n");
        PP      " -l       : Reference sequence length is the shortest. \n\n");
        PP      " -L       : Reference sequence length is the largest. \n\n");
        PP      " -a       : Reference sequence length is the alignment length (default). \n\n");
        PP      " -n       : Score is normalized by reference sequence length (default).\n\n");
        PP      " -r       : Raw score, not normalized. \n\n");
        PP      " -d       : Score is expressed in distance (default : score is expressed in similarity). \n\n");
        PP      " -t ##.## : Score threshold for clustering. If the score is normalized and expressed in similarity (default),\n");
        PP		"            it is an identity, e.g. 0.95 for an identity of 95%%. If the score is normalized\n");
        PP		"            and expressed in distance, it is (1.0 - identity), e.g. 0.05 for an identity of 95%%.\n");
        PP		"            If the score is not normalized and expressed in similarity, it is the length of the\n");
        PP		"            Longest Common Subsequence. If the score is not normalized and expressed in distance,\n");
        PP		"            it is (reference length - LCS length).\n");
        PP		"            Only sequences with a similarity above ##.## with the center sequence of a cluster\n");
        PP		"            are assigned to that cluster. Default: 0.97.\n\n");
        PP      " -e       : Exact option : A sequence is assigned to the cluster with the center sequence presenting the\n");
        PP		"            highest similarity score > threshold, as opposed to the default 'fast' option where a sequence is\n");
        PP		"            assigned to the first cluster found with a center sequence presenting a score > threshold.\n\n");
		PP      " -p ##    : Multithreading with ## threads using openMP.\n\n");
        PP      " -s ####  : Sorting by ####. Must be 'None' for no sorting, or a key in the fasta header of each sequence,\n");
        PP		"            except for the count that can be computed (default : sorting by count).\n\n");
        PP      " -o       : Sorting is in ascending order (default : descending).\n\n");
        PP      " -g       : n's are replaced with a's (default: sequences with n's are discarded).\n\n");
        PP      " -B ###   : Output of the OTU table in BIOM format is activated, and written to file ###.\n\n");
        PP      " -O ###   : Output of the OTU map (observation map) is activated, and written to file ###.\n\n");
        PP      " -F       : Output in FASTA format is deactivated.\n");
        PP      "\n");
        PP      "------------------------------------------------------------\n");
        PP      " Argument : the nucleotide dataset to cluster\n");
        PP      "------------------------------------------------------------\n");
        PP		" http://metabarcoding.org/sumatra\n");
        PP      "------------------------------------------------------------\n\n");
}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit                         */
/* ----------------------------------------------- */

#define PP fprintf(stderr,


static void ExitUsage(stat)
        int stat;
{
        PP      "usage: sumaclust [-l|L|a|n|r|d|e|o|g|F] [-t threshold_value] [-s sorting_key] [-p number_of_threads]\n");
        PP		"[-B file_name_for_BIOM-formatted_output] [-O file_name_for_OTU_table-formatted_output] dataset\n");
        PP      "type \"sumaclust -h\" for help\n");

        if (stat)
            exit(stat);
}

#undef  PP


static char* sortingKey="count";

static int sortSeqsP(const void **s1, const void **s2)
{
	int res;
	double r1;
	double r2;

	r1 = atof(getItemFromHeader(sortingKey, ((fastaSeqPtr) *s2)->header));
	r2 = atof(getItemFromHeader(sortingKey, ((fastaSeqPtr) *s2)->header));
	if (r2 > r1)
		res = 1;
	else if (r2 < r1)
		res = -1;
	else
		res = 0;

	return(res);
}


static int reverseSortSeqsP(const void **s1, const void **s2)
{
	int res;
	double r1;
	double r2;

	r1 = atof(getItemFromHeader(sortingKey, ((fastaSeqPtr) *s2)->header));
	r2 = atof(getItemFromHeader(sortingKey, ((fastaSeqPtr) *s2)->header));

	if (r1 > r2)
		res = 1;
	else if (r1 < r2)
		res = -1;
	else
		res = 0;

	return(res);
}


int uniqSeqsDoubleSortFunction(const void *s1, const void *s2)
{
	int c;
	char* str_r1;
	double r1;
	double r2;

	c = strcmp(((fastaSeqPtr) s1)->sequence, ((fastaSeqPtr) s2)->sequence);
	if (c == 0)
	{
		str_r1 = getItemFromHeader(sortingKey, ((fastaSeqPtr) s1)->header);
		if (str_r1 == NULL)
		{
			fprintf(stderr, "\nERROR: '%s' not in sequence header(s).\n\n", sortingKey);
			exit(1);
		}
		r1 = atof(str_r1);
		r2 = atof(getItemFromHeader(sortingKey, ((fastaSeqPtr) s2)->header));

		if (r2 > r1)
			c = 1;
		else if (r2 < r1)
			c = -1;
		else
			c = 0;
	}
	return(c);
}


int uniqSeqsDoubleReverseSortFunction(const void *s1, const void *s2)
{
	int c;
	char* str_r1;
	double r1;
	double r2;

	c = strcmp(((fastaSeqPtr) s1)->sequence, ((fastaSeqPtr) s2)->sequence);
	if (c == 0)
	{
		str_r1 = getItemFromHeader(sortingKey, ((fastaSeqPtr) s1)->header);
		if (str_r1 == NULL)
		{
			fprintf(stderr, "\nERROR: '%s' not in sequence header(s).\n\n", sortingKey);
			exit(1);
		}
		r1 = atof(str_r1);
		r2 = atof(getItemFromHeader(sortingKey, ((fastaSeqPtr) s2)->header));

		if (r1 > r2)
			c = 1;
		else if (r1 < r2)
			c = -1;
		else
			c = 0;
	}
	return(c);
}


void printInBIOMformat(fastaSeqPtr* uniqSeqs, int count, int numberOfCenters, char* biomFile_name)
{
	int i, j, n;
	FILE* biomFile;
	struct tm* tm_info;
	time_t timer;
	char buffer_date[20];
	fastaSeqPtr* c;
	fastaSeqPtr* seq;
	int id_len;
	int row_number;
	BOOL first_center = TRUE;

	int buffer_col_rows;
	int buffer_col_rows_1;
	int buffer_col_rows_2;

	buffer_col_rows = 29;
	buffer_col_rows_1 = 9;
	buffer_col_rows_2 = 20;

	n = 0;

	biomFile = fopen(biomFile_name, "w");
	if (biomFile == NULL)
	  fprintf(stderr, "\nCan't open BIOM output file.\n"); //, %s outputFilename);

    for (i=0; i<count; i++)    // Loop to store columns
    {
    	seq = uniqSeqs+i;
    	id_len = strlen((*seq)->accession_id);
    	j=0;

    	if ((*seq)->cluster_center) 	// center sequence
    	{
    		n++;
    		(*seq)->cluster_weight_unique_ids = 1;

    		if (first_center)
    		{
    			(*seq)->columns_BIOM_size = id_len + buffer_col_rows;
    			(*seq)->columns_BIOM = (char*) malloc(((*seq)->columns_BIOM_size)*sizeof(char));
    			strcpy((*seq)->columns_BIOM, "{\"id\": \"");
    			first_center = FALSE;
    		}
    		else
    		{
    			(*seq)->columns_BIOM_size = id_len + buffer_col_rows + 1;
    			(*seq)->columns_BIOM = (char*) malloc(((*seq)->columns_BIOM_size)*sizeof(char));
    			strcpy((*seq)->columns_BIOM, ",{\"id\": \"");
    		}

    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - id_len - buffer_col_rows_2 - 1, (*seq)->accession_id, id_len);
    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - buffer_col_rows_2 - 1, "\", \"metadata\": null}", buffer_col_rows_2+1);

    		if ((*seq)->next != NULL)	// not last sequence
    		{
				for (j=1; ((((*seq)+j)->next != NULL) && (((*seq)+j)->uniqHead == FALSE)); j++)		// identical sequences
				{
					id_len = strlen((*(seq)+j)->accession_id);
					n++;

					(*seq)->cluster_weight_unique_ids++;
					(*seq)->columns_BIOM_size = (*seq)->columns_BIOM_size + id_len + buffer_col_rows;
					(*seq)->columns_BIOM = realloc((*seq)->columns_BIOM, ((*seq)->columns_BIOM_size) * sizeof(char));
		    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - buffer_col_rows - id_len - 1, ",{\"id\": \"", buffer_col_rows_1);
		    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - id_len - buffer_col_rows_2 - 1, (*(seq)+j)->accession_id, id_len);
		    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - buffer_col_rows_2 - 1, "\", \"metadata\": null}", buffer_col_rows_2+1);
				}
				if ((((*seq)+j)->next == NULL) && (((*seq)+j)->uniqHead == FALSE))	// last sequence
				{
					id_len = strlen((*(seq)+j)->accession_id);
					n++;

					(*seq)->cluster_weight_unique_ids++;
					(*seq)->columns_BIOM_size = (*seq)->columns_BIOM_size + id_len + buffer_col_rows;
					(*seq)->columns_BIOM = realloc((*seq)->columns_BIOM, ((*seq)->columns_BIOM_size) * sizeof(char));
		    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - buffer_col_rows - id_len - 1, ",{\"id\": \"", buffer_col_rows_1);
		    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - id_len - buffer_col_rows_2 - 1, (*(seq)+j)->accession_id, id_len);
		    		memcpy((*seq)->columns_BIOM + (*seq)->columns_BIOM_size - buffer_col_rows_2 - 1, "\", \"metadata\": null}", buffer_col_rows_2+1);
				}
    		}
    	}
    	else	// not a center sequence
    	{
    		n++;

			c = (*seq)->center;

			id_len = strlen((*seq)->accession_id);
			n++;

			(*c)->cluster_weight_unique_ids++;
			(*c)->columns_BIOM_size = (*c)->columns_BIOM_size + id_len + buffer_col_rows;
			(*c)->columns_BIOM = realloc((*c)->columns_BIOM, ((*c)->columns_BIOM_size) * sizeof(char));
			memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - buffer_col_rows - id_len - 1, ",{\"id\": \"", buffer_col_rows_1);
			memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - id_len - buffer_col_rows_2 - 1, (*seq)->accession_id, id_len);
			memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - buffer_col_rows_2 - 1, "\", \"metadata\": null}", buffer_col_rows_2+1);

			if ((*seq)->next != NULL)	// not last sequence
			{
				for (j=1; ((((*seq)+j)->next != NULL) && (((*seq)+j)->uniqHead == FALSE)); j++)		// identical sequences
				{
					id_len = strlen((*(seq)+j)->accession_id);
					n++;

					(*c)->cluster_weight_unique_ids++;
					(*c)->columns_BIOM_size = (*c)->columns_BIOM_size + id_len + buffer_col_rows;
					(*c)->columns_BIOM = realloc((*c)->columns_BIOM, ((*c)->columns_BIOM_size) * sizeof(char));
					memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - buffer_col_rows - id_len - 1, ",{\"id\": \"", buffer_col_rows_1);
					memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - id_len - buffer_col_rows_2 - 1, (*(seq)+j)->accession_id, id_len);
					memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - buffer_col_rows_2 - 1, "\", \"metadata\": null}", buffer_col_rows_2+1);
				}

				if ((((*seq)+j)->next == NULL) && (((*seq)+j)->uniqHead == FALSE))	// last sequence
				{
					id_len = strlen((*(seq)+j)->accession_id);
					n++;

					(*c)->cluster_weight_unique_ids++;
					(*c)->columns_BIOM_size = (*c)->columns_BIOM_size + id_len + buffer_col_rows;
					(*c)->columns_BIOM = realloc((*c)->columns_BIOM, ((*c)->columns_BIOM_size) * sizeof(char));
					memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - buffer_col_rows - id_len - 1, ",{\"id\": \"", buffer_col_rows_1);
					memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - id_len - buffer_col_rows_2 - 1, (*(seq)+j)->accession_id, id_len);
					memcpy((*c)->columns_BIOM + (*c)->columns_BIOM_size - buffer_col_rows_2 - 1, "\", \"metadata\": null}", buffer_col_rows_2+1);
				}
			}
		}
    }

	time(&timer);
	tm_info = localtime(&timer);
	strftime(buffer_date, 20, "%Y-%m-%dT%H:%M:%S", tm_info);

	fprintf(biomFile, "{\"id\": \"None\",\"format\": \"Biological Observation Matrix 1.0.0\","
			"\"format_url\": \"http://biom-format.org\",\"type\": \"OTU table\","
			"\"generated_by\": \"SUMACLUST %s\",\"date\": \"%s\",\"matrix_type\": \"sparse\","
			"\"matrix_element_type\": \"int\",\"shape\": [%d, %d],",
			VERSION, buffer_date, numberOfCenters, n);

	// print data

	row_number = 0;
	n = 0;

	fprintf(biomFile, "\"data\": [");

	for (i=0; i<count; i++)
	{
		seq = uniqSeqs+i;
	    if ((*seq)->cluster_center) 	// center sequence
	    {
	    	for (j=0; j<(*seq)->cluster_weight_unique_ids; j++)
	    	{
	    		if ((row_number == (numberOfCenters - 1)) && (j == ((*seq)->cluster_weight_unique_ids - 1)))	// last seq to print
	    			fprintf(biomFile, "[%d,%d,1]],", row_number, n);
	    		else
	    			fprintf(biomFile, "[%d,%d,1],", row_number, n);
	    		n++;
	    	}
	    	row_number++;
	    }
	}
	// end data

	// Print rows

	first_center = TRUE;

    for (i=0; i<count; i++)
    {
    	seq = uniqSeqs+i;
    	if ((*seq)->cluster_center) 	// center sequence
    	{
    		if (first_center)
    		{
    			fprintf(biomFile, "\"rows\": [{\"id\": \"%s\", \"metadata\": null}", (*seq)->accession_id);
    			first_center = FALSE;
    		}
    		else
    			fprintf(biomFile, ",{\"id\": \"%s\", \"metadata\": null}", (*seq)->accession_id);
    	}
     }

    // Print columns

	fprintf(biomFile, "],\"columns\": [");
	for (i=0; i<count; i++)
	{
		seq = uniqSeqs+i;
	    if ((*seq)->cluster_center) 	// center sequence
	    	fprintf(biomFile, (*seq)->columns_BIOM);
	}
	fprintf(biomFile, "]}");

	fclose(biomFile);
}


void printInOTUtableFormat(fastaSeqPtr* uniqSeqs, int count, char* OTUtableFile_name)
{
	int i, j;
	FILE* OTUtableFile;
	fastaSeqPtr* c;
	fastaSeqPtr* seq;
	int id_len;

	OTUtableFile = fopen(OTUtableFile_name, "w");
	if (OTUtableFile == NULL)
	  fprintf(stderr, "\nCan't open OTU table output file.\n"); //, %s outputFilename);

    for (i=0; i<count; i++)
    {
    	seq = uniqSeqs+i;
    	id_len = strlen((*seq)->accession_id);
    	j=0;


    	if ((*seq)->cluster_center) 	// center sequence
    	{
   			(*seq)->line_OTU_table_size = id_len*2 + 2;
   			(*seq)->line_OTU_table = (char*) malloc(((*seq)->line_OTU_table_size)*sizeof(char));
   			strcpy((*seq)->line_OTU_table, (*seq)->accession_id);
    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - id_len - 2, "\t", 1);
    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - id_len - 1, (*seq)->accession_id, id_len);
    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - 1, "\0", 1);

    		if ((*seq)->next != NULL)	// not last sequence
    		{
				for (j=1; ((((*seq)+j)->next != NULL) && (((*seq)+j)->uniqHead == FALSE)); j++)		// identical sequences
				{
					id_len = strlen((*(seq)+j)->accession_id);

					(*seq)->line_OTU_table_size = (*seq)->line_OTU_table_size + id_len + 1;
					(*seq)->line_OTU_table = realloc((*seq)->line_OTU_table, ((*seq)->line_OTU_table_size) * sizeof(char));
		    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - id_len - 2, "\t", 1);
		    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - id_len - 1, (*(seq)+j)->accession_id, id_len);
		    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - 1, "\0", 1);
				}

				if ((((*seq)+j)->next == NULL) && (((*seq)+j)->uniqHead == FALSE))	// last sequence
				{
					id_len = strlen((*(seq)+j)->accession_id);

					(*seq)->line_OTU_table_size = (*seq)->line_OTU_table_size + id_len + 1;
					(*seq)->line_OTU_table = realloc((*seq)->line_OTU_table, ((*seq)->line_OTU_table_size) * sizeof(char));
		    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - id_len - 2, "\t", 1);
		    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - id_len - 1, (*(seq)+j)->accession_id, id_len);
		    		memcpy((*seq)->line_OTU_table + (*seq)->line_OTU_table_size - 1, "\0", 1);
				}
    		}
    	}
    	else	// not a center sequence
    	{
			c = (*seq)->center;

			(*c)->line_OTU_table_size = (*c)->line_OTU_table_size + id_len + 1;
			(*c)->line_OTU_table = realloc((*c)->line_OTU_table, ((*c)->line_OTU_table_size) * sizeof(char));
    		memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - id_len - 2, "\t", 1);
    		memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - id_len - 1, (*seq)->accession_id, id_len);
    		memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - 1, "\0", 1);

    		if ((*seq)->next != NULL)	// not last sequence
    		{
				for (j=1; ((((*seq)+j)->next != NULL) && (((*seq)+j)->uniqHead == FALSE)); j++)		// identical sequences
				{
					id_len = strlen((*(seq)+j)->accession_id);

					(*c)->line_OTU_table_size = (*c)->line_OTU_table_size + id_len + 1;
					(*c)->line_OTU_table = realloc((*c)->line_OTU_table, ((*c)->line_OTU_table_size) * sizeof(char));
					memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - id_len - 2, "\t", 1);
					memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - id_len - 1, (*(seq)+j)->accession_id, id_len);
					memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - 1, "\0", 1);
				}

				if ((((*seq)+j)->next == NULL) && (((*seq)+j)->uniqHead == FALSE))	// last sequence
				{
					id_len = strlen((*(seq)+j)->accession_id);

					(*c)->line_OTU_table_size = (*c)->line_OTU_table_size + id_len + 1;
					(*c)->line_OTU_table = realloc((*c)->line_OTU_table, ((*c)->line_OTU_table_size) * sizeof(char));
					memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - id_len - 2, "\t", 1);
					memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - id_len - 1, (*(seq)+j)->accession_id, id_len);
					memcpy((*c)->line_OTU_table + (*c)->line_OTU_table_size - 1, "\0", 1);
				}
    		}
		}
    }

	// Print rows

    for (i=0; i<count; i++)
    {
    	seq = uniqSeqs+i;
    	if ((*seq)->cluster_center) 	// center sequence
    	{
   			fprintf(OTUtableFile, (*seq)->line_OTU_table);
   			fprintf(OTUtableFile, "\n");
    	}
     }

	fclose(OTUtableFile);
}


void printSeq(fastaSeqPtr* seq, fastaSeqPtr* center, double score)
{
	int i;

	char* score_n;
	char* score_v;
	char* cluster_n;
	char* cluster_v;
	char* center_n;
	char* center_true;
	char* center_false;
	int id_size;

	score_n = (char*) malloc(14*sizeof(char));
	score_v = (char*) malloc(20*sizeof(char));

	strcpy(score_n, "cluster_score");
	sprintf(score_v,"%f", score);

	id_size = strlen((*center)->accession_id);

	cluster_n = (char*) malloc(8*sizeof(char));
	cluster_v = (char*) malloc((id_size+1)*sizeof(char));

	strcpy(cluster_n, "cluster");
	strcpy(cluster_v, (*center)->accession_id);

	center_n = (char*) malloc(15*sizeof(char));
	strcpy(center_n, "cluster_center");
	center_true = (char*) malloc(5*sizeof(char));
	strcpy(center_true, "True");
	center_false = (char*) malloc(6*sizeof(char));
	strcpy(center_false, "False");

	(*seq)->header = table_header_add_field((*seq)->header, cluster_n, cluster_v);
	(*seq)->header = table_header_add_field((*seq)->header, score_n, score_v);
	if ((*seq)->cluster_center)
		(*seq)->header = table_header_add_field((*seq)->header, center_n, center_true);
	else
		(*seq)->header = table_header_add_field((*seq)->header, center_n, center_false);

	printOnlyHeaderFromTable((*seq)->header);
	printOnlySeqFromFastaSeqPtr((*seq));

	if ((*seq)->next != NULL)
	{
		for (i=1; ((((*seq)+i)->next != NULL) && (((*seq)+i)->uniqHead == FALSE)); i++)
		{
			((*seq)+i)->header = table_header_add_field(((*seq)+i)->header, cluster_n, cluster_v);
			((*seq)+i)->header = table_header_add_field(((*seq)+i)->header, score_n, score_v);
			((*seq)+i)->header = table_header_add_field(((*seq)+i)->header, center_n, center_false);

			printOnlyHeaderFromTable(((*seq)+i)->header);
			printOnlySeqFromFastaSeqPtr(((*seq)+i));
		}

		if ((((*seq)+i)->next == NULL) && (((*seq)+i)->uniqHead == FALSE))	// last sequence
		{
			((*seq)+i)->header = table_header_add_field(((*seq)+i)->header, cluster_n, cluster_v);
			((*seq)+i)->header = table_header_add_field(((*seq)+i)->header, score_n, score_v);
			((*seq)+i)->header = table_header_add_field(((*seq)+i)->header, center_n, center_false);

			printOnlyHeaderFromTable(((*seq)+i)->header);
			printOnlySeqFromFastaSeqPtr(((*seq)+i));
		}
	}
}


void putSeqInCluster(fastaSeqPtr* seq, fastaSeqPtr* center, double score)
{
	(*seq)->center = center;
	(*seq)->score  = score;
}


int compare(fastaSeqPtr* db,  int n, BOOL fastOption, double threshold, BOOL normalize, int reference, BOOL lcsmode)
{
	double       score;
	double	     scoremax;
	double	     worstscore;
	BOOL	     toCluster;
	static BOOL  first=TRUE;
	int32_t      i,j,k;
	int          center;
	float        p;
	BOOL         found;
	int			 lmax, lmin;
	int16_t*	 address;
	int16_t*	 iseq1;
	int16_t*	 iseq2;
	int			 l1;
	int			 l2;
	char* 		 s1;
	char*		 s2;
	int			 sizeForSeqs;
	int			 LCSmin;

	if (lcsmode || normalize)
		fprintf(stderr,"Clustering sequences when similarity >= %lf\n", threshold);
	else
		fprintf(stderr,"Clustering sequences when distance <= %lf\n", threshold);

	fprintf(stderr,"Aligning and clustering... \n");

	int* centers = (int*) malloc(n * sizeof(int));

	for (i=0; i < n; i++)
		centers[i] = -1;

	k=0;
	found = FALSE;

	calculateMaxAndMinLen(db, n, &lmax, &lmin);

	sizeForSeqs = prepareTablesForSumathings(lmax, lmin, threshold, normalize, reference, lcsmode, &address, &iseq1, &iseq2);

	if (lcsmode || normalize)
		worstscore = 0.0;
	else
		worstscore = lmax;

	for (i=0; i < n; i++)
	{
		if (i%100 == 0)
		{
			p = (i/(float)n)*100;
			fprintf(stderr,"\rDone : %f %%       %d clusters created",p,k);
		}

		if (first)
		{
			first = FALSE;
			if (normalize && lcsmode)
				score = 1.0;
			else if (!lcsmode)
				score = 0.0;
			else
				score = (*(db+i))->length;
			(*(db+i))->cluster_center = TRUE;
			putSeqInCluster(db+i, db+i, score);
			centers[k] = i;
			k++;
		}

		else
		{
			scoremax = worstscore;
			center = 0;
			found = FALSE;
			toCluster = FALSE;
			j=0;

			s1 = (*(db+i))->sequence;
			l1 = (*(db+i))->length;

			while (((found == FALSE) && (centers[j] != -1) && (fastOption == TRUE)) || ((fastOption == FALSE) && (centers[j] != -1)))
			{
				score = worstscore;
				filters((*(db+i)), (*(db+centers[j])), threshold, normalize, reference, lcsmode, &score, &LCSmin);
				if (score == -1.0)
				{
					s2 = (*(db+centers[j]))->sequence;
					l2 = (*(db+centers[j]))->length;

					score = alignForSumathings(s1, iseq1, s2, iseq2, l1, l2, normalize, reference, lcsmode, address, sizeForSeqs, LCSmin);
				}

				if (((score >= threshold) && (lcsmode || normalize) && (score > scoremax)) || ((!lcsmode && !normalize) && (score <= threshold) && (score < scoremax)))
				{
					toCluster = TRUE;
					scoremax = score;
					center = centers[j];
					if (fastOption == TRUE)
						found = TRUE;
				}
				j++;
			}

			if (toCluster)
			{
				if (!lcsmode && normalize)
					scoremax = 1.0 - scoremax;
				(*(db+i))->cluster_center = FALSE;
				putSeqInCluster(db+i, db+center, scoremax);
			}
			else
			{
				if (normalize && lcsmode)
					score = 1.0;
				else if (!lcsmode)
					score = 0.0;
				else
					score = (*(db+i))->length;
				(*(db+i))->cluster_center = TRUE;
				putSeqInCluster(db+i, db+i, score);
				centers[k] = i;
				k++;
			}
		}
	}
	fprintf(stderr,"\rDone : 100 %%       %d clusters created.                        \n",k);

	free(centers);

	free(iseq1-sizeForSeqs+lmax);
	free(iseq2-sizeForSeqs+lmax);

	if (normalize && reference == ALILEN)
		free(address);

	return(k);

}


void computeClusterWeights(fastaSeqPtr* uniqSeqs, int n)
{
	int i,j;
	fastaSeqPtr* seq;
	fastaSeqPtr* center;
	char* cluster_weight_n;
	char* cluster_weight_v;
	int cluster_weight;

    for (i=0; i<n; i++)
    {
    	seq = uniqSeqs+i;

    	if ((*seq)->cluster_center)
    		(*seq)->cluster_weight = (*seq)->count;
    	else
    	{
    		center = (*seq)->center;
    		((*center)->cluster_weight)+=(*seq)->count;
    	}
    }

    for (i=0; i<n; i++)
    {
    	seq = uniqSeqs+i;

    	if ((*seq)->cluster_center)
    		cluster_weight = (*seq)->cluster_weight;
    	else
    	{
    		center = (*seq)->center;
    		cluster_weight = (*center)->cluster_weight;
    	}
		cluster_weight_n = (char*) malloc(15*sizeof(char));
		cluster_weight_v = (char*) malloc(20*sizeof(char));
		strcpy(cluster_weight_n, "cluster_weight");
		sprintf(cluster_weight_v,"%d", cluster_weight);
		(*seq)->header = table_header_add_field((*seq)->header, cluster_weight_n, cluster_weight_v);

		if ((*seq)->next != NULL)		// not the last sequence
		{
			for (j=1; ((((*seq)+j)->next != NULL) && (((*seq)+j)->uniqHead == FALSE)); j++)
				(*(seq)+j)->header = table_header_add_field((*(seq)+j)->header, cluster_weight_n, cluster_weight_v);

			if ((((*seq)+j)->next == NULL) && (((*seq)+j)->uniqHead == FALSE))	// last sequence
				(*(seq)+j)->header = table_header_add_field((*(seq)+j)->header, cluster_weight_n, cluster_weight_v);
		}
    }
}


int main(int argc, char** argv)
{

	int32_t     	carg		   = 0;
	int32_t         errflag        = 0;
	char*			sort;
	double          threshold      = 0.97;
	BOOL			lcsmode		   = TRUE;
	BOOL 			fastOption     = TRUE;
	BOOL     		normalize      = TRUE;
	BOOL			reverse		   = FALSE;
	BOOL			onlyATGC	   = TRUE;
	int             reference      = ALILEN;
	int             ndb            = 0;
	int				nproc          = 1;
	BOOL			printBIOM      = FALSE;
	BOOL			printOTUtable      = FALSE;
	BOOL			printFASTA	   = TRUE;
	fastaSeqCount   db;
	int				i,n;
	fastaSeqPtr*	uniqSeqs;
	char* 			biomFile_name;
	char* 			OTUtableFile_name;
	int 			numberOfCenters;

	sort = malloc(1024*sizeof(char));
	biomFile_name = malloc(1024*sizeof(char));
	OTUtableFile_name = malloc(1024*sizeof(char));
	strcpy(sort, "count");

	while ((carg = getopt(argc, argv, "hlLanrdet:p:s:ogB:O:F")) != -1) {
		switch (carg) {
									/* -------------------- */
		case 'h':                   /* help                 */
									/* -------------------- */
			PrintHelp();
			exit(0);
			break;

					  /* -------------------------------------------------- */
		case 'l':     /* Normalize LCS/Error by the shortest sequence length*/
					  /* -------------------------------------------------- */
			reference=MINLEN;
			break;

					  /* -------------------------------------------------- */
		case 'L':     /* Normalize LCS/Error by the largest sequence length */
					  /* -------------------------------------------------- */
			reference=MAXLEN;
			break;

						  /* -------------------------------------------------- */
		case 'a':         /* Normalize LCS/Error by the alignment length        */
						  /* -------------------------------------------------- */
			reference=ALILEN;
			break;

					  /* -------------------------------------------------- */
		case 'n':     /* Normalize LCS by the reference length              */
					  /* -------------------------------------------------- */
			normalize=TRUE;
			break;

					  /* -------------------------------------------------- */
		case 'r':     /* No normalization					                */
					  /* -------------------------------------------------- */
			normalize=FALSE;
			break;

					  /* -------------------------------------------------- */
		case 'd':     /* Score is expressed in distance                  */
					  /* -------------------------------------------------- */
			lcsmode=FALSE;
			break;

								/* ---------------------------------------------------------------------------------------------------------- */
		case 'e':               /* center with the best score > threshold is chosen, otherwise first center with a score > threshold     */
								/* ---------------------------------------------------------------------------------------------------------- */
			fastOption=FALSE;
			break;

						/* ------------------------------------------------------------------- */
		case 't':   	/* Clusters only pairs with similarity higher than (threshold)         */
						/* ------------------------------------------------------------------- */
			sscanf(optarg,"%lf",&threshold);
			break;

					 	 	 /* -------------------------------------------------- */
		case 'p':            /* number of processors to use                        */
							 /* -------------------------------------------------- */
			sscanf(optarg,"%d",&nproc);
			break;

									/* -------------------------------------------------- */
		case 's':   				/* Sorting option								      */
									/* -------------------------------------------------- */
			sscanf(optarg, "%s", sort);
			sortingKey = sort;
			break;

									/* -------------------------------------------------- */
		case 'o':   				/* reverse sorting                                   */
									/* -------------------------------------------------- */
			reverse=TRUE;
			break;

									/* -------------------------------------------------- */
		case 'g':   				/* replace n's with a's in sequences                  */
									/* -------------------------------------------------- */
			onlyATGC=FALSE;
			break;

									/* -------------------------------------------------- */
		case 'B':   				/* file name to print results in BIOM format          */
									/* -------------------------------------------------- */
			sscanf(optarg, "%s", biomFile_name);
			printBIOM=TRUE;
			break;

									/* -------------------------------------------------- */
		case 'O':   				/* file name to print results in OTU table format     */
									/* -------------------------------------------------- */
			sscanf(optarg, "%s", OTUtableFile_name);
			printOTUtable=TRUE;
			break;

									/* -------------------------------------------------- */
		case 'F':   				/* don't print results in FASTA format    	          */
									/* -------------------------------------------------- */
			printFASTA=FALSE;
			break;


		case '?':               	/* invalid option   	        */
			errflag++;
			break;
		}
	}

	ndb = argc - optind;
	if (ndb != 1)
        errflag++;

	if (errflag)
		ExitUsage(errflag);

    fprintf(stderr,"===========================================================\n");
	fprintf(stderr," SUMACLUST version %s\n",VERSION);
#ifdef __SSE2__
	fprintf(stderr," Alignment using SSE2 instructions.\n");
#else
	fprintf(stderr," Alignment using standard code, SSE2 unavailable.\n");
#endif
	fprintf(stderr,"===========================================================\n");

	if ((threshold == 0.0) || (normalize && (threshold > 1.0)))
	{
		fprintf(stderr, "\nERROR: Please specify a threshold > 0, and < 1 when scores are normalized.\n\n");
		exit(1);
	}

	fprintf(stderr,"Reading dataset...");
	db = seq_readAllSeq2(argv[optind], TRUE, onlyATGC);
	fprintf(stderr,"\n%d sequences\n",db.count);

	if (!onlyATGC)
		(void)cleanDB(db);

	if (!lcsmode && normalize)
		threshold = 1.0 - threshold;

	if (threshold > 0)
		(void)hashDB(db);

	addCounts(&db);

	// first sorting of sequences to have good unique heads

	if ((strcmp(sortingKey, "None") != 0) && (strcmp(sortingKey, "none") != 0))
	{
		if (reverse == FALSE)
			qsort((void*) db.fastaSeqs, db.count, sizeof(fastaSeq), uniqSeqsDoubleSortFunction);
		else
			qsort((void*) db.fastaSeqs, db.count, sizeof(fastaSeq), uniqSeqsDoubleReverseSortFunction);
	}

	// getting the vector of unique seqs
	uniqSeqs = (fastaSeqPtr*) malloc((db.count)*sizeof(fastaSeqPtr));
	n = uniqSeqsVector(&db, &uniqSeqs);
	uniqSeqs = realloc(uniqSeqs, n*sizeof(fastaSeqPtr));

	// putting a flag on the last sequence
	for (i=0; i<(db.count-1); i++)
		((db.fastaSeqs)+i)->next = (db.fastaSeqs)+i-1;
    ((db.fastaSeqs)+(db.count)-1)->next = NULL;

	// sorting unique sequences
	if (strcmp(sortingKey, "count") == 0)
	{
		fprintf(stderr,"Sorting sequences by count...\n", n);
		if (reverse == FALSE)
			qsort((void*) uniqSeqs, n, sizeof(fastaSeqPtr), sortSeqsWithCounts);
		else
			qsort((void*) uniqSeqs, n, sizeof(fastaSeqPtr), reverseSortSeqsWithCounts);
	}
	else if ((strcmp(sortingKey, "None") != 0) && (strcmp(sortingKey, "none") != 0))
	{
		fprintf(stderr,"Sorting sequences by %s...\n", sortingKey);
		if (reverse == FALSE)
			qsort((void*) uniqSeqs, n, sizeof(fastaSeqPtr), sortSeqsP);
		else
			qsort((void*) uniqSeqs, n, sizeof(fastaSeqPtr), reverseSortSeqsP);
	}

	// Computing
	if (nproc==1)
		numberOfCenters = compare(uniqSeqs, n, fastOption, threshold, normalize, reference, lcsmode);

	else
		numberOfCenters = mt_compare_sumaclust(uniqSeqs, n, fastOption, threshold, normalize, reference, lcsmode, nproc);

	// Computing cluster weights
	computeClusterWeights(uniqSeqs, n);

	// Printing results

	// FASTA file
	if (printFASTA)
	{
		fprintf(stderr,"Printing results in FASTA format...\n");
		for (i=0; i<n; i++)
		{
			printSeq(uniqSeqs+i, (*(uniqSeqs+i))->center, (*(uniqSeqs+i))->score);
		}
	}
	fprintf(stderr,"Done.\n");

    // BIOM file
	if (printBIOM)
	{
		fprintf(stderr,"Printing results in BIOM format...\n");
		printInBIOMformat(uniqSeqs, n, numberOfCenters, biomFile_name);
	}
	fprintf(stderr,"Done.\n");

    // OTU table file
	if (printOTUtable)
	{
		fprintf(stderr,"Printing results in OTU table format...\n");
		printInOTUtableFormat(uniqSeqs, n, OTUtableFile_name);
	}
	fprintf(stderr,"Done.\n");

	// Freeing
	for (i=0; i < db.count; i++)
	{
		free(((db.fastaSeqs)[i]).table);
		free_header_table(((db.fastaSeqs)[i]).header);
	}
	free(db.fastaSeqs);
	free(sort);
	free(uniqSeqs);

	return(0);

}
