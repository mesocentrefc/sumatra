/**
 * FileName:    fileHandling.c
 * Authors:      Tiayyba Riaz, Celine Mercier
 * Description: C file for file handling functions
 * **/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../libutils/utilities.h"

/*
 * Function Name: fileOpen(char* fileName, BOOL abortOnError)
 * Description:   Opens the file and returns the pointer to file object
 */
FILE *file_open(char* fileName, BOOL abortOnError)
{
	FILE* fp;
	
	if (fileName == NULL && abortOnError)
		ERRORABORT(FILE_OPENING_ERROR, "File name not given.");
	
	if (fileName == NULL)
		return NULL;
	
	fp = fopen(fileName, "r");
	return fp;
}

FILE *file_openrw(char* fileName, BOOL abortOnError)
{
	FILE* fp;
	
	if (fileName == NULL && abortOnError)
		ERRORABORT(FILE_OPENING_ERROR, "File name not given.");
	
	if (fileName == NULL)
		return NULL;
	
	fp = fopen(fileName, "w+");
	 return fp;
}

/*
 * Function Name: fileNextChar(FILE* fp)
 * Description:   Reads the file and returns next character, if file is null or its end of file, returns \¯.
 */
char file_nextChar(FILE* fp)
{
	if (fp == NULL)
		return '\0';
	
	if(feof(fp))
		return '\0';
	
	return (char) fgetc(fp);
}

/*
 * Function Name: *fileNextLine(FILE *fp, char *buffer, int32_t bufferSize)
 * Description:   Reads the file and returns next line, if file is null or its end of file, returns \¯.
 */
char *file_nextLine(FILE *fp, char *buffer, int32_t bufferSize)
{
	if(fp == NULL)
		return NULL;

	if(feof(fp))
		return NULL;

	return fgets(buffer, bufferSize, fp);
}


void exitIfEmptyFile(FILE *file)
{
    long savedOffset = ftell(file);
    fseek(file, 0, SEEK_END);

    if (ftell(file) == 0)
    {
        fprintf(stderr, "\nInput file is empty.\n");
        exit(1);
    }
    fseek(file, savedOffset, SEEK_SET);
}

