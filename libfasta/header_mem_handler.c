#include <stdio.h>
#include <stdlib.h>
#include "header_mem_handler.h"
#include <string.h>

#define FIELD_BUFFER 1024


char* malloc_field(int *free_size)
{
	char* field = (char*) malloc(sizeof(char) * FIELD_BUFFER);
	field[0] = 0;
	(*free_size) = FIELD_BUFFER;
	return field;
}

int check_mem_field(int size_needed)
{
	int number_of_chunks_to_alloc;
	number_of_chunks_to_alloc = size_needed / FIELD_BUFFER + 1;
	return number_of_chunks_to_alloc;
}

char* realloc_field(int number_of_chunks_to_alloc, char* field)
{
	int size_needed;
	size_needed = number_of_chunks_to_alloc * FIELD_BUFFER;
	field = realloc(field, (size_needed)*sizeof(char));
	return field;
}

char* check_and_realloc_field(char* field, int size_needed, int* free_size)
{
	size_needed = size_needed + strlen(field);
	int number_of_chunks_to_alloc = check_mem_field(size_needed);
	if (strlen(field)>0)
		field = realloc_field(number_of_chunks_to_alloc, field);
	else
	{
		free(field);
		field = malloc(number_of_chunks_to_alloc * FIELD_BUFFER);
	}
	(*free_size) = number_of_chunks_to_alloc*FIELD_BUFFER - size_needed + 1;
	return field;
}


char* store_in_field(char* field, char* yytext, int* free_size, int* i)
{
	int size_needed;
	size_needed = strlen(yytext)+1;
	if (size_needed > (*free_size))
		field = check_and_realloc_field(field, size_needed, free_size);
	else
		(*free_size) = (*free_size) - size_needed + 1;
	strcpy(&(field[(*i)]),yytext);
	(*i) = (*i)+size_needed-1;
	return field;
}


char* store_in_header_table(char* field, char** storing_place, int* free_size, int* i)
{
	int size_needed;
	size_needed = strlen(field)+1;
	*storing_place = (char*) malloc(size_needed*sizeof(char));
	strcpy(*storing_place,field);
	(*i)=0;
	free(field);
	field = malloc_field(free_size);
	return field;
}


element_from_header** check_and_realloc_mem_in_header_table(element_from_header** p_header, int* nbf, int* memory_allocated)
{
	(*nbf)++;

	if (*nbf == *memory_allocated)
	{
		(*memory_allocated)++;
		*p_header = (element_from_header*) realloc(*p_header, (*memory_allocated) * sizeof(element_from_header));
	}

	return p_header;
}

void end_header_table(element_from_header** p_header, int nbf)
{
	nbf = nbf - 1;
	//fprintf(stderr, "nbf = %d", nbf);
	sprintf((*p_header)->value, "%d", nbf);
}
