/**
 * FileName:    fileHandling.h
 * Authors:      Tiayyba Riaz, Celine Mercier
 * Description: Header file for file handling functions
 * **/


#ifndef FILEHANDLING_H_
#define FILEHANDLING_H_

#include "../libutils/utilities.h"
/* Prototypes */

FILE *file_open(char* fileName, BOOL abortOnError);
char file_nextChar(FILE* fp);
char *file_nextLine(FILE *fp, char *buffer, int32_t bufferSize);
FILE *file_openrw(char* fileName, BOOL abortOnError);
void exitIfEmptyFile(FILE *file);

#endif /*FILEHANDLING_H_*/
