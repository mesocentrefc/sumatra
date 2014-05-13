
#ifndef UPPERBAND_H_
#define UPPERBAND_H_


int buildTable(const char *sequence, unsigned char *table, int *count);
int compareTable(unsigned char *t1, int over1, unsigned char* t2,  int over2);
int threshold4(int wordcount,double identity);
int thresholdLCS4(int32_t reflen,int32_t lcs);


int hashDB(fastaSeqCount);
BOOL isPossible(fastaSeqPtr, fastaSeqPtr, BOOL, int, double, BOOL);
BOOL isPossibleSumathings(fastaSeqPtr seq1, fastaSeqPtr seq2, int l1, int l2, double threshold, BOOL normalize, int reference, BOOL lcsmode);
void filters(fastaSeqPtr seq1, fastaSeqPtr seq2, double threshold, BOOL normalize, int reference, BOOL lcsmode, double* score, int* LCSmin);
void filtersSumatra(fastaSeqPtr seq1, fastaSeqPtr seq2, double threshold, BOOL normalize, int reference, BOOL lcsmode, double* score, int* LCSmin);
#endif

