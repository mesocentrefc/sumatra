/*
 * Add -ll in Makefile if you modify this file to convert to .c
 */

%x REGID
%x REGNAME
%x REGVAL

%{

#include <stdlib.h>
#include <string.h>
#include "header_mem_handler.h"
#include "fasta_header_handler.h"

#define MEMALLOCATED 10
#define BUFFER 5
	
#define YY_DECL int header_parser(int *nbf, int *memory_allocated, element_from_header **p_header)

	
%}

WORD				[[:alnum:]:\-.{},'_()\#\[\]\|\&\"\'\/\%\+]+
SUP					>
EOL                	\n
SEP                 ;
SPACE				[[:blank:]]+
EQUAL				=

%%

	int i;
	int size_needed;
	int free_size;
	char* field;
	
	
<INITIAL>{SUP}	{	
					/*printf("\n<INITIAL>{SUP},%s",yytext);*/
					BEGIN(REGID);
				}
				
<INITIAL,REGID>{WORD}     {
								i=0;						
								
								field = malloc_field(&free_size);
								(*p_header)[*nbf].name = (char*) malloc(3*sizeof(char));
								strcpy(((*p_header)[*nbf]).name,"id");
										
								size_needed = strlen(yytext)+1;
								(*p_header)[*nbf].value = (char*) malloc(sizeof(char)*size_needed);
								strcpy(((*p_header)[*nbf]).value,yytext);
										
								(*nbf)++;
						  }

<INITIAL,REGID>{EQUAL}     {
								fprintf(stderr, "\nERROR: '=' symbol not allowed in sequence identifiers.\n\n");
								exit(1);
						  }

<INITIAL,REGID>{SPACE}     {
						BEGIN(REGNAME);
					}

<REGNAME>{WORD}      {
						/*fprintf(stderr,"\n<REGNAME>{WORD} **%s**",yytext);*/
						field = store_in_field(field,yytext,&free_size,&i);
				    }

<REGNAME>{SPACE}     {
						/*fprintf(stderr,"\n<REGNAME>{SPACE} **%s**",yytext);*/
						if (i != 0)
							field = store_in_field(field,yytext,&free_size,&i);
				    }

<REGNAME>{EQUAL}     {
						/*fprintf(stderr,"\n<REGNAME>{EQUAL},%s",yytext);*/
						field = store_in_header_table(field, &((*p_header)[*nbf].name), &free_size, &i);
						BEGIN(REGVAL);
				    }

<REGNAME>{SEP}       {
						/*fprintf(stderr,"\n<REGNAME>{SEP},%s",yytext);*/
						(*p_header)[*nbf].name = (char*) malloc(19*sizeof(char));
						strcpy((*p_header)[*nbf].name,"definition");
						field = store_in_header_table(field, &((*p_header)[*nbf].value), &free_size, &i);
						p_header = check_and_realloc_mem_in_header_table(p_header, nbf, memory_allocated);
						BEGIN(REGNAME);
				    }

<REGVAL>{WORD}      {
					/*fprintf(stderr,"\n<REGVAL>{WORD} **%s**\n",yytext);*/
					field = store_in_field(field,yytext,&free_size,&i);
				    }

<REGVAL>{SPACE}     {
					/*fprintf(stderr,"\n<REGVAL>{SPACE} **%s**\n",yytext);*/
					field = store_in_field(field,yytext,&free_size,&i);
				    }

<REGVAL>{SEP}      {
					/*fprintf(stderr,"\n<REGVAL>{SEP},%s\n",yytext);*/
					
					field = store_in_header_table(field, &((*p_header)[*nbf].value), &free_size, &i);
					p_header = check_and_realloc_mem_in_header_table(p_header, nbf, memory_allocated);
					BEGIN(REGNAME);
				    }


<REGVAL>{EQUAL}     {
						/*fprintf(stderr, "\nWarning : separator ';' probably missing in header after %s",(*p_header)[*nbf].name);*/
				    }

<REGVAL><<EOF>>  {
					field = store_in_header_table(field, &((*p_header)[*nbf].value), &free_size, &i);
					p_header = check_and_realloc_mem_in_header_table(p_header, nbf, memory_allocated);
					end_header_table(p_header, *nbf);
					
					free(field);
					BEGIN(INITIAL);
                	return 0;
                 }

<REGNAME><<EOF>>      {
						/*(*p_header)[*nbf].name = (char*) malloc(sizeof(char)*19);
						strcpy((*p_header)[*nbf].name,"other_informations");
						field = store_in_header_table(field, &((*p_header)[*nbf].value), &free_size, &i);
						p_header = check_and_realloc_mem_in_header_table(p_header, nbf, memory_allocated);						
						*/
						end_header_table(p_header, *nbf);
						
						free(field);
						BEGIN(INITIAL);
						return 0;
                	}

%%

int header_yywrap()
{
	return 1;
}

element_from_header* header_parser_main(char *h)
{
	int nbfields,memory_allocated;
	element_from_header* header;
	char* nbfields_n;
	char* nbfields_v;
	
	nbfields_n = (char*) malloc(9*sizeof(char));
	nbfields_v = (char*) malloc(5*sizeof(char));
	
	memory_allocated=MEMALLOCATED;
	
	nbfields=1;

	strcpy(nbfields_n, "nbfields");
	strcpy(nbfields_v, "1");

    header = (element_from_header*) malloc(memory_allocated * sizeof(element_from_header));
    
    header[0].name = nbfields_n;
    header[0].value = nbfields_v;
    
	YY_BUFFER_STATE state;

	state=yy_scan_string(h);

    header_parser(&nbfields, &memory_allocated, &header);
      
    yy_delete_buffer(state);
        
    return header;
}




