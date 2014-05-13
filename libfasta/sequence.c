/**
 * FileName:    sequence.c
 * Authors:      Tiayyba Riaz, Celine Mercier
 * Description: C file for sequence reading and parsing
 * **/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "../libutils/utilities.h"
#include "sequence.h"
#include "../libfile/fileHandling.h"
#include "fasta_header_handler.h"
#include "fasta_header_parser.h"
#include "dic_parser.h"


/*
 * Function Name: seq_getNext(FILE *fp, char *fieldDelim)
 * Description:   Gets the next sequence from file by calling another function, passes the sequence
 * to other function  to get the header elements and nucleotide suquence into a strcuture of 
 * type fastaSeq and returns a pointer to this newly populated structure. 
 */

fastaSeqPtr seq_getNext(FILE *fp, char *fieldDelim,  BOOL isStandardSeq, BOOL onlyATGC)
{
	char *seq;
	char *header;
	char *strTemp;
	fastaSeqPtr seqElem;
	int seqLen;

	seq = seq_readNextFromFilebyLine(fp);
	if (seq == NULL) return NULL;
	
	/* Find header separator \n, if not found return NULL */
	strTemp = strchr(seq, '\n');
	if(strTemp == NULL)
		return NULL;
	
	seqLen = strlen(strTemp);
	header = (char*) util_malloc(1+(strlen(seq) - seqLen)*sizeof(char), __FILE__, __LINE__);

	/* Separate header in header variable */
	strncpy(header, seq, strTemp - seq);
	header[strTemp - seq] = '\0';
	/* Get memory for new sequence structure element */
	seqElem = (fastaSeqPtr) util_malloc(sizeof(fastaSeq), __FILE__, __LINE__);
	/* Parse header and assign values to structure fields */
	seq_fillHeader(header, fieldDelim, seqElem);
	/* Get clean sequence and assign to structure field */
	if (isStandardSeq)
		if (onlyATGC)
			seq_fillSeqOnlyATGC(strTemp, seqElem, seqLen);
		else
			seq_fillSeq(strTemp, seqElem, seqLen);
	else
		seq_fillDigitSeq(strTemp, seqElem, seqLen);
	/* Type cast the char * seq to void pointer and deallocate the memory pointed by this */
	util_free((void *)seq);
	/* Return new sequence structure element */
	return seqElem;
} 


char *seq_readNextFromFilebyLine(FILE* fp)
{
	char newc = '\0';
	BOOL seqCompleted = FALSE;
	int length = 500;
	int32_t len;
	char tempstr[length];
	char* buffer;
	
	if (feof(fp)) return NULL;
	newc = file_nextChar(fp);
	if (newc != '>') ungetc(newc, fp);
		
	buffer = util_malloc(1*sizeof(char), __FILE__, __LINE__);
	buffer[0] = '\0';
	
	while(!seqCompleted)
	{
		newc = file_nextChar(fp);
		if(newc == '>' || newc == '\0')
		{
			seqCompleted = TRUE;
			if (newc == '>')
				ungetc(newc, fp); // Make sure next time we start from sequence delimiter >
		}
		else
		{
			ungetc(newc, fp);
			if(file_nextLine( fp, tempstr, length) != NULL)
			{
				len = strlen(tempstr) + strlen(buffer) + 1;
				buffer = util_realloc(buffer, len, __FILE__, __LINE__);
				strcat(buffer, tempstr);
			}
			else
			{
				seqCompleted = TRUE;
			}
		}
	}
	return buffer;
}


/*
 * Function Name: seq_fillHeader(char* header, char *fieldDelim, fastaSeqPtr seqElem)
 */
void seq_fillHeader(char* header, char *fieldDelim, fastaSeqPtr seqElem)
{
	char* IdEnd;
	int IdSize;
	
	seqElem->rawheader = strdup(header);

	IdEnd = strchr(header, ' ');
	if (IdEnd == NULL)
		IdSize = strlen(header);
	else
		IdSize = strlen(header) - strlen(IdEnd);

	seqElem->accession_id = (char*) util_malloc(1+IdSize*sizeof(char), __FILE__, __LINE__);

	strncpy(seqElem->accession_id, header, IdSize);

	(seqElem->accession_id)[IdSize] = '\0';
}


/*
 * Function Name: seq_fillSeq(char *seq, fastaSeqPtr seqElem)
 * Description:   Parses the whole sequences for actual nucleotide sequences and stores that
 * sequence in the field of structure 'seqElem' . 
 */
void seq_fillSeq(char *seq, fastaSeqPtr seqElem, int seqLen)
{
	char* seqTemp;
	char c;
	int32_t index = 0, seqIndex = 0, len = strlen(seq);
	char* seqAlphabets = "acgtACGT-nN";
	
	seqTemp = (char*) util_malloc(seqLen*sizeof(char), __FILE__, __LINE__);

	while (index < len)
	{
		c = seq[index++];
		if (strchr(seqAlphabets, c) != NULL)
			seqTemp[seqIndex++] = tolower(c);
	}
	seqTemp[seqIndex] = '\0';
	seqElem->length=seqIndex;
	seqElem->sequence = strdup(seqTemp);
}


void seq_fillSeqOnlyATGC(char *seq, fastaSeqPtr seqElem, int seqLen)
{
	char* seqTemp;
	char c;
	int32_t index = 0, seqIndex = 0, len = strlen(seq);
	char* seqAlphabets = "acgtACGT";
	int notAllATGC = 0;

	seqTemp = (char*) util_malloc(seqLen*sizeof(char), __FILE__, __LINE__);

	while (index < len)
	{
		c = seq[index++];
		if (strchr(seqAlphabets, c) != NULL)
			seqTemp[seqIndex++] = tolower(c);
		else if (c != '\n')
			notAllATGC = 1;
	}

	if (notAllATGC)
		seqTemp[0] = '\0';
	else
	{
		seqTemp[seqIndex] = '\0';
		seqElem->length=seqIndex;
	}
	seqElem->sequence = strdup(seqTemp);
}


void seq_fillDigitSeq(char *seq, fastaSeqPtr seqElem, int seqLen)
{
	char* seqTemp;
	char c;
	int32_t index = 0, seqIndex = 0, len = strlen(seq);
	
	seqTemp = (char*) util_malloc(seqLen*sizeof(char), __FILE__, __LINE__);

	while (index < len)
	{
		c = seq[index++];
		if ((c >= '0' && c <= '9') || c == ' ')
			seqTemp[seqIndex++] = c;
		/*else
		{
			printf("Error in input file");
			exit(0);			
		}*/
	}
	seqTemp[seqIndex] = '\0';  
	seqElem->sequence = strdup(seqTemp);
}


fastaSeqCount seq_readAllSeq2(char *fileName, BOOL isStandardSeq, BOOL onlyATGC)
{
    FILE* fp;
    fastaSeqPtr seqPtr;
    fastaSeqPtr seqPtrAr;
    
    int32_t counter = 0;
    int32_t slots = 1000;
    fastaSeqCount allseqs;
    int32_t discarded=0;

    fp = file_open(fileName, TRUE);

    if (fp == NULL)
    {
            fprintf(stderr, "\nCould not open file.\n");
            exit(1);
    }

    exitIfEmptyFile(fp);

    seqPtrAr = (fastaSeqPtr) util_malloc(slots*sizeof(fastaSeq), __FILE__, __LINE__);

    seqPtr = seq_getNext(fp, " ", isStandardSeq, onlyATGC);

    while (seqPtr != NULL)
    {
    	if (counter == slots)
    	{
    		slots += 1000;
    		seqPtrAr = (fastaSeqPtr)util_realloc(seqPtrAr, slots*sizeof(fastaSeq), __FILE__, __LINE__);
    	}
    	
    	if ((seqPtr->sequence)[0] != '\0')
    		seqPtrAr[counter++] = *seqPtr;
    	else
    		discarded++;

    	util_free((void *)seqPtr);
    	seqPtr = seq_getNext(fp, " ", isStandardSeq, onlyATGC);
    }
    fclose(fp);

    if (counter != slots)
    	seqPtrAr = (fastaSeqPtr)util_realloc(seqPtrAr, counter*sizeof(fastaSeq), __FILE__, __LINE__);

    allseqs.count = counter;
    allseqs.fastaSeqs = seqPtrAr;
    
    if (discarded)
    	fprintf(stderr, "\nDiscarded %d sequences that did not contain only 'AaTtGgCc' characters.", discarded);

    return allseqs;
}


int32_t seq_findSeqByAccId (char *accid, fastaSeqCountPtr allseqs)
{
	int32_t i;
	
	for (i = 0; i < allseqs->count; i++)
	{
		if (strcmp (accid, allseqs->fastaSeqs[i].accession_id) == 0)
			return i;
	}
	return -1;
}


void seq_printSeqs (fastaSeqCountPtr allseq)
{
	int32_t i;
	
	for (i = 0; i < allseq->count; i++)
	//for (i = 0; i < 4; i++)
	{
		if (allseq->fastaSeqs[i].sequence == NULL) continue;
		if (allseq->fastaSeqs[i].rawheader)
			printf (">%s\n", allseq->fastaSeqs[i].rawheader);
		else
			printf (">%s\n", allseq->fastaSeqs[i].accession_id);
		printf ("%s\n", allseq->fastaSeqs[i].sequence);
	}
}


int cleanDB(fastaSeqCount db) // replace not a/t/g/c with a's
{
	int32_t i;
	char    *seq;
	BOOL    changed;
	int32_t seqchanged=0;
	int32_t nucchanged=0;

	fprintf(stderr,"Cleaning dataset...");

	for (i=0; i < db.count;i++)
	{

		changed=FALSE;
		for (seq = db.fastaSeqs[i].sequence; *seq!=0; seq++)
		{
			if (*seq!='a' && *seq!='c' && *seq!='g' && *seq!='t')
			{
				changed=TRUE;
				nucchanged++;
				*seq='a';
			}
		}
		if (changed)
			seqchanged++;
	}

	if (seqchanged)
		fprintf(stderr," : %d nucleotides substituted in %d sequences\n",nucchanged,seqchanged);
	else
		fprintf(stderr," : Done\n");

	return(db.count);
}


void addCounts(fastaSeqCount* db)
{
	int s;
	char* count;
	element_from_header* header;
	char* count_n;
	char* count_v;

	count_n = (char*) malloc(6*sizeof(char));
	count_v = (char*) malloc(2*sizeof(char));

	strcpy(count_n, "count");
	strcpy(count_v, "1");

	for (s=0; s < db->count; s++)
	{
		header = header_parser_main(db->fastaSeqs[s].rawheader);
		count = getItemFromHeader("count", header);
		if (count == 0)	 // no count field
		{
			header = table_header_add_field(header, count_n, count_v);
			db->fastaSeqs[s].count = 1;
		}
		else
			db->fastaSeqs[s].count = atoi(count);
		db->fastaSeqs[s].header = header;
	}
}


void readSampleCounts(fastaSeqCount* db, char* key_name)
{
	int s;
	element_from_header* header;
	char* to_parse;

	for (s=0; s < db->count; s++)
	{
		header = header_parser_main(db->fastaSeqs[s].rawheader);
		db->fastaSeqs[s].header = header;
		to_parse = getItemFromHeader(key_name, header);
		if (to_parse == 0)	 // no field named key_name
		{
			fprintf(stderr, "\nERROR: '%s' not in sequence header(s).\n\n", key_name);
			exit(1);
		}
		else
		{
			db->fastaSeqs[s].sample_counts = parseDictionary(to_parse);
		}
	}
}


int uniqSeqsVector(fastaSeqCount* db, fastaSeqPtr** uniqSeqs)
{
	int i, j, k;
	*(*(uniqSeqs)) = db->fastaSeqs;
	db->fastaSeqs[0].uniqHead = TRUE;

	i = 0;
	k = 1;

	for (j=1; j < db->count; j++)
	{
		if (strcmp(db->fastaSeqs[i].sequence, db->fastaSeqs[j].sequence) == 0)
		{
			db->fastaSeqs[i].count += db->fastaSeqs[j].count;
			db->fastaSeqs[j].uniqHead = FALSE;
		}
		else
		{
			db->fastaSeqs[j].uniqHead = TRUE;
			*(*(uniqSeqs)+k) = (db->fastaSeqs)+j;
			k++;
			i = j;
		}
	}
	return(k);
}


void calculateMaxAndMinLen(fastaSeqPtr* db, int n, int* lmax, int* lmin)
{
	int i;
	int l;

	*lmax = 0;
	for (i=0; i < n; i++)
	{
		l = (*(db+i))->length;
		if (l > *lmax)
			*lmax = l;
	}

	*lmin = *lmax;
	for (i=0; i < n; i++)
	{
		l = (*(db+i))->length;
		if (l < *lmin)
			*lmin = l;
	}
}


void calculateMaxAndMinLenDB(fastaSeqCount db, int* lmax, int* lmin)
{
	int i;
	int l;

	*lmax = 0;
	for (i=0; i < db.count; i++)
	{
		l = ((db.fastaSeqs)+i)->length;
		if (l > *lmax)
			*lmax = l;
	}

	*lmin = *lmax;
	for (i=0; i < db.count; i++)
	{
		l = ((db.fastaSeqs)+i)->length;;
		if (l < *lmin)
			*lmin = l;
	}
}


int sortSeqsWithCounts(const void **s1, const void **s2)
{
	return(((fastaSeqPtr) *s2)->count - ((fastaSeqPtr) *s1)->count);
}


int reverseSortSeqsWithCounts(const void **s1, const void **s2)
{
	return(((fastaSeqPtr) *s1)->count - ((fastaSeqPtr) *s2)->count);
}
