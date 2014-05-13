/**
 * FileName:    utilities.h
 * Author:      Tiayyba Riaz
 * Description: Header file for miscellenious functions and macros
 * **/

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <stdint.h>
#include <stdio.h>
#include <time.h>



//static char *basecodes = "00100020000000000003000000";

//#define BASEIDX(ch) basecodes[ch - 'a'] - 48

#ifndef MAX
#define MAX(x,y) (((x)>(y)) ? (x):(y))
#define MIN(x,y) (((x)<(y)) ? (x):(y))
#endif

typedef char BOOL;
#define TRUE (3==3)
#define FALSE (!TRUE)
#define ALILEN (0)
#define MAXLEN (1)
#define MINLEN (2)


/* Error Codes */
#define FILE_OPENING_ERROR (1)
#define MEM_ALLOC_ERROR (2)

/* Prototypes */
void errorAbort(int32_t code, char* errorMsg, char* fileName, int32_t lineNumber);
char *str_strrstr(char *src, char *delim);
void str_removeSpaces(char *src);
void str_sepNameValue(char *name, char *value, char *src, char *delim);
char *str_chopAtDelim(char *dest, char *src, char *delim, BOOL includeDelim);
void util_free(void *chunk);
void *util_realloc(void *chunk, size_t newsize, const char *filename, int32_t    line);
void *util_malloc(size_t chunksize, const char *filename, int32_t line);
BOOL util_findInArr(int32_t tempArr[], int seqNo, int32_t noOfSeqs);
void* getA16Address(int size);
void** reallocA16Address(void** address, int size);

/* Macros */
#define ERRORABORT(code, msg)  errorAbort((code), (msg), __FILE__, __LINE__)

#endif /*UTILITIES_H_*/



