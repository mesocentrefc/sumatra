#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sequence.h"
#include "fasta_header_parser.h"


void printOnlySeqFromFastaSeqPtr(fastaSeqPtr seq)
{
	char nuc;
	int n=60;
	int l = strlen(seq->sequence);
	for (n=60; n<l; n+=60)
	{
		nuc = seq->sequence[n];
		seq->sequence[n]=0;
		printf("%s\n",seq->sequence+n-60);
		seq->sequence[n]=nuc;
	}
	printf("%s\n",seq->sequence+n-60);
}


void printOnlySeqFromChar(char* seq)
{
	char nuc;
	int n=60;
	int l = strlen(seq);
	for (n=60; n<l; n+=60)
	{
		nuc = seq[n];
		seq[n]=0;
		printf("%s\n",seq+n-60);
		seq[n]=nuc;
	}
	printf("%s\n",seq+n-60);
}


void printOnlyHeaderFromFastaSeqPtr(fastaSeqPtr seq)
{
	printf(">%s\n",seq->rawheader);
}


void printOnlyHeaderFromTable(element_from_header* header)
{
	int i;
	int nbf;

	nbf = atoi(header[0].value);

	printf(">%s ",header[1].value);

	for (i = 2; i <= nbf; i++)
	{
		if (strcmp(header[i].name, "definition") != 0)
		{
			printf("%s",header[i].name);
			printf("=");
		    printf("%s; ",header[i].value);
		}
	}

	if (strcmp(header[nbf].name, "definition") == 0)
		printf("%s; ",header[nbf].value);

	printf("\n");
}


void printHeaderAndSeqFromFastaSeqPtr(fastaSeqPtr seq)
{
	printOnlyHeaderFromFastaSeqPtr(seq);
	printOnlySeqFromFastaSeqPtr(seq);
}
