/*
 * Add -ll in Makefile if you modify this file to convert to .c
 */

%x REGKEY
%x REGVAL

%{

#include <stdlib.h>
#include <string.h>
#include "header_mem_handler.h"
#include "fasta_header_handler.h"

	
#define YY_DECL int dic_parser(struct hashtable* dic)

	
%}

WORD				[[:alnum:]\-._()\[\]\|\&\/\%]+
SEP                 ,
COLON				:
OPENING_BRACE		"{"
CLOSING_BRACE		"}"
NOT_KEEPING			[[:blank:]'"]

%%

	int size_needed;
	char* key;
	char* value;
	struct hashtable* temp;

	
<INITIAL>{OPENING_BRACE}	{
								//fprintf(stderr,"\n<INITIAL>{OPENING_BRACE} **%s**",yytext);
								BEGIN(REGKEY);
							}

<REGKEY>{WORD}				{
								//fprintf(stderr,"\n<REGKEY>{WORD} **%s**",yytext);
								size_needed = strlen(yytext)+1;
								key = (char*) malloc(size_needed*sizeof(char));
								strcpy(key, yytext);
							}

<REGKEY>{COLON}				{
								//fprintf(stderr,"\n<REGKEY>{COLON} **%s**",yytext);
								BEGIN(REGVAL);
							}

<REGVAL>{WORD}				{
								//fprintf(stderr,"\n<REGVAL>{WORD} **%s**",yytext);
								size_needed = strlen(yytext)+1;
								value = (char*) malloc(size_needed*sizeof(char));
								strcpy(value, yytext);
							}

<REGVAL>{SEP}       		{
								//fprintf(stderr,"\n<REGVAL>{SEP} **%s**",yytext);
								install(dic, key, value);
								temp = lookup(dic, key);
								BEGIN(REGKEY);
				    		}

<INITIAL,REGKEY,REGVAL>{NOT_KEEPING}	{;}	

<REGVAL>{CLOSING_BRACE}		{
								//fprintf(stderr,"\n<REGVAL>{CLOSING_BRACE} **%s**",yytext);
								install(dic, key, value);
								temp = lookup(dic, key);
								BEGIN(INITIAL);
								return 0;
							}
%%


int hashtable_yywrap()
{
	return 1;
}


struct hashtable* parseDictionary(char* to_parse)
{
	int i;
	int l;
	int nb_elements = 0;
	struct hashtable* dic;
	
	l = strlen(to_parse);

	// count number of ':' or ',' (useless for now but eventually replae HASHSIZE with it)
	for (i=0; i<l; i++)
	{
		if (to_parse[i] == ':')
			nb_elements++;
	}
	
	dic = (struct hashtable*) malloc(HASHSIZE*sizeof(struct hashtable));
	
	for (i=0; i<HASHSIZE; i++)
	{
		(dic+i)->next = NULL;
		(dic+i)->name = "";
		(dic+i)->defn = NULL;
	}
	
	YY_BUFFER_STATE state;
	
	state=yy_scan_string(to_parse);

    dic_parser(dic);
          
    yy_delete_buffer(state);

    return dic;
}

