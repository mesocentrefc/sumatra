/**
 * FileName:    utilities.c
 * Author:      Tiayyba Riaz
 * Description: C file for miscellenious functions and macros
 * **/

#include "utilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Function Name: errorAbort(int errorCode, char* errorMsg, char* fileName, int lineNumber)
 * Description:   Reports an error on standard error and aborts
 */
void errorAbort(int32_t errorCode, char* errorMsg, char* fileName, int32_t lineNumber)
{
	fprintf(stderr,"Error %d in file %s line %d : %s\n",
			errorCode,
			fileName,
	        lineNumber,
	        errorMsg);
	
	abort();
}

void   *util_malloc(size_t chunksize, const char *filename, int32_t line)
{
	void * chunk;
	
	chunk = calloc(1,chunksize);
	
	if (!chunk)
		errorAbort(MEM_ALLOC_ERROR,"Could not allocate memory.",filename,line);
	
	return chunk;
}

/*
 * Function Name: util_realloc(void *chunk, int32_t newsize, const char *filename, int32_t line)
 * Description:   Overloading realloc funstion, changes the size of the memory object pointed to by chunk 
 * to the size specified by newsize. If memory cannot be allocated, gives the error on stderr and aborts.
 */
void   *util_realloc(void *chunk, size_t newsize, const char *filename, int32_t    line)
{
	void *newchunk;
	
	newchunk = realloc(chunk,newsize);
	
	if (!newchunk)
	{
		errorAbort(MEM_ALLOC_ERROR,"Could not allocate memory.",filename,line);
	}

	return newchunk;	
}

/*
 * Function Name: util_free(void *chunk)
 * Description:   Returns the memory specified by chunk back to operating syste.
 */
void    util_free(void *chunk)
{
	free(chunk);
}

BOOL util_findInArr(int32_t tempArr[], int seqNo, int32_t noOfSeqs)
{
	int index;
	
	for(index = 0; index < noOfSeqs; index++)
	{
		if(tempArr[index] == seqNo) return TRUE;
	}
	
	return FALSE;
}


/**
 * 
 * String handling utilities
 * 
 **/

/*
 * Function Name: str_chopAtDelim(char *dest, char *src, char *delim, BOOL includeDelim)
 * Description:   chops the string startig from source to the delimeter specified.
 */
char *str_chopAtDelim(char *dest, char *src, char *delim, BOOL includeDelim)
{
	char *temp;
	int32_t len;
	
	/* returns a pointer to the first occurance of delim in src*/
	temp = strstr(src, delim); 
	
	if (temp == NULL)
		return NULL;
	
	if (includeDelim)
	{
		/* temp - src + strlen(delim) -> a string between src and delimeter including delimeter*/
		len = temp - src + strlen(delim);
		strncpy(dest, src, len);
	}
	else
	{
		len = temp - src;
		strncpy(dest, src, temp - src);
	}
	dest[len] = '\0';
	
	return dest;
}

/*
 * Function Name: str_sepNameValue(char *name, char *value, char *src, char *delim)
 * Description:   .
 */
void str_sepNameValue(char *name, char *value, char *src, char *delim)
{
	char *temp;
	
	temp = strstr(src, delim);
	
	if(temp != NULL)
	{
		strncpy(name, src, temp - src);
		strcpy(value, temp + strlen(delim));
	}
	else
	{
		strcpy(name, src);
		strcpy(value, "");
	}
}

/*
 * Function Name: str_removeSpaces(char *src)
 * Description:   Removes the spaces from the start and end of the string.
 */
int str_isSpace (char ch)
{
	switch (ch)
	{
	case ' ':
	case '\t':
	case '\n':
		return 1;
	}
	return 0;
}

void str_removeSpaces(char *src)
{
	int32_t start = 0, end = strlen(src) - 1;
	int32_t index = 0;
	
	if (src == NULL || end < 0) return;
	
	while(str_isSpace(src[start]) && start < end) start++;
	while(str_isSpace(src[end]) && end > start) end--;
	
	if ( start == end && src[start] == ' ')
	{
		src[0] = '\0';
		return;
	}
	if (start > 0)
	{
		while(start <= end)
		{
			src[index] = src[start];
			index++;
			start++;
		}
		src[index] = '\0';
		return;
	}
	src[end+1] = '\0';
}

/*
 * Function Name: str_strrstr(char *src, char *delim)
 * Description:   Searches the position of last occurence of string delim in src.
 */
char *str_strrstr(char *src, char *delim)
{
	char *last, *next;
	next = strstr(src, delim);
	last = next;
	while(next != NULL)
	{
		last = next;
		next = strstr(last + 1, delim);
	}
	return last;
}


void* getA16Address(int size)
{
	void* address;
	address = (void*) malloc(size);
	while ((((long long unsigned int) (address))%16) != 0)
		address++;
	return(address);
}


void** reallocA16Address(void** address, int size)
{
	if (*(address) == NULL)
		*(address) = malloc(size);
	*(address) = realloc(address, size);
	while ((((long long unsigned int) (*(address)))%16) != 0)
		(*(address))++;
	return(address);
}










