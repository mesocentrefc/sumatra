/**
 * FileName:    sequence.h
 * Authors:      Tiayyba Riaz, Celine Mercier
 * Description: Prototypes and other declarations for sequences
 * **/
#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <stdint.h>
#include <stdio.h>
#include "../libutils/utilities.h"
#include "fasta_header_parser.h"


typedef struct {
	char* accession_id;					// identifier
	char *rawheader;					// not parsed header
	element_from_header* header;		// parsed header
	char *sequence;						// DNA sequence itself
	int32_t length;						// DNA sequence's length
	int32_t count;						// abundance of the sequence
	unsigned char *table;      			// 4mer occurrence table build using function buildTable
	int32_t over;              			// count of 4mer with occurrences greater than 255 (overflow)
	struct fastaSeqPtr* next;			// next unique sequence for example
	BOOL cluster_center;				// whether the sequence is a cluster center or not
	int32_t cluster_weight;				// cluster weight when sequence is cluster center
	int32_t cluster_weight_unique_ids;	// cluster weight when sequence is cluster center, counting the number sequence records
	double score;						// score with cluster center for example
	struct fastaSeqPtr* center;			// pointer to the sequence's cluster center
	int32_t center_index;				// index of the sequence's cluster center
	BOOL uniqHead;						// whether the sequence is a unique head or not
	char* columns_BIOM;					// to print in BIOM format
	int   columns_BIOM_size;			// size allocated for columns_BIOM
	char* line_OTU_table;				// to print in OTU table format
	int	  line_OTU_table_size;			// size allocated for line_OTU_table
	struct hashtable *sample_counts;	// sample counts for sumaclean
}fastaSeq,*fastaSeqPtr;


typedef struct {
	int32_t count;
	fastaSeqPtr fastaSeqs;
}fastaSeqCount, *fastaSeqCountPtr;


fastaSeqPtr seq_getNext(FILE *fp, char *fieldDelim, BOOL isStandardSeq, BOOL onlyATGC);
char *seq_readNextFromFilebyLine(FILE* fp);
void seq_fillSeq(char *seq, fastaSeqPtr seqElem, int seqLen);
void seq_fillSeqOnlyATGC(char *seq, fastaSeqPtr seqElem, int seqLen);
void seq_fillDigitSeq(char *seq, fastaSeqPtr seqElem, int seqLen);
void seq_fillHeader(char* header, char *fieldDelim, fastaSeqPtr seqElem);
fastaSeqCount seq_readAllSeq2(char *fileName, BOOL isStandardSeq, BOOL onlyATGC);
int32_t seq_findSeqByAccId (char *accid, fastaSeqCountPtr allseqs);
void seq_printSeqs (fastaSeqCountPtr allseq);
int cleanDB(fastaSeqCount);
void addCounts(fastaSeqCount* db);
int uniqSeqsVector(fastaSeqCount* db, fastaSeqPtr** uniqSeqs);
void calculateMaxAndMinLen(fastaSeqPtr* db, int n, int* lmax, int* lmin);
void calculateMaxAndMinLenDB(fastaSeqCount db, int* lmax, int* lmin);
int sortSeqsWithCounts(const void **s1, const void **s2);
int reverseSortSeqsWithCounts(const void **s1, const void **s2);
void readSampleCounts(fastaSeqCount* db, char* key_name);

#endif /*SEQUENCE_H_*/
