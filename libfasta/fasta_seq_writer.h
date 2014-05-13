
#ifndef FASTA_SEQ_WRITER_H_
#define FASTA_SEQ_WRITER_H_


void printOnlySeqFromFastaSeqPtr(fastaSeqPtr);

void printOnlySeqFromChar(char*);

void printOnlyHeaderFromFastaSeqPtr(fastaSeqPtr);

void printOnlyHeaderFromTable(element_from_header*);

void printHeaderAndSeqFromFastaSeqPtr(fastaSeqPtr);


#endif
