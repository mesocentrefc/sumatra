#ifndef HEADER_MEM_HANDLER_H_
#define HEADER_MEM_HANDLER_H_

#include "fasta_header_parser.h"

char* malloc_field(int*);

int check_mem_field(int);

char* realloc_field(int, char*);

char* check_and_realloc_field(char*, int, int*);

char* store_in_field(char*, char*, int*, int*);

char* store_in_header_table(char*, char**, int*, int*);

element_from_header** check_and_realloc_mem_in_header_table(element_from_header**, int*, int*);

void end_header_table(element_from_header** p_header, int nbf);

#endif
