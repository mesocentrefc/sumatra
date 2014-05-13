#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sequence.h"
#include "fasta_header_parser.h"
#include "dic_parser.h"
#include "fasta_header_handler.h"


char* char_header_add_field(char* header, char* name, char* value)
{
	int lheader = strlen(header);
	header = (char*) realloc(header, (lheader+strlen(name)+strlen(value)+4)*sizeof(char));
	if (header[lheader-1] == '.')
	{
		strcpy(header+lheader-1,";");
		strcpy(header+lheader," ");
		strcpy(header+lheader+1,name);
		strcpy(header+lheader+1+strlen(name),"=");
		strcpy(header+lheader+1+strlen(name)+1,value);
	}
	else
	{
		strcpy(header+lheader,";");
		strcpy(header+lheader+1," ");
		strcpy(header+lheader+2,name);
		strcpy(header+lheader+2+strlen(name),"=");
		strcpy(header+lheader+2+strlen(name)+1,value);
	}
	return header;
}


char* fastaSeqPtr_header_add_field(fastaSeqPtr seq, char* name, char* value)
{
	int lheader = strlen(seq->rawheader);
	int i;
	char* buffer;
	char* rawheader;

	rawheader = (char*) malloc((lheader+strlen(name)+strlen(value)+5)*sizeof(char));
	strcpy(rawheader, seq->rawheader);

	buffer = calloc(lheader, sizeof(char));

	i=0;

	while ((rawheader[i] != ' ') && (rawheader[i] != 0))
		i++;

	if (rawheader[i] == ' ')
		strcpy(buffer, rawheader+i);
	else
		strcpy(rawheader+i, " ");

	i++;

	strcpy(rawheader+i,name);
	strcpy(rawheader+i+strlen(name),"=");
	strcpy(rawheader+i+strlen(name)+1,value);
	strcpy(rawheader+i+strlen(name)+1+strlen(value),";");
	strcpy(rawheader+i+strlen(name)+1+strlen(value)+1, buffer);

	free(buffer);

	return(rawheader);
}


element_from_header* table_header_add_field(element_from_header* header, char* name, char* value)
{
	int nbf;
	nbf = atoi(header[0].value);
	nbf++;
	header = (element_from_header*) realloc(header, (nbf+1)*sizeof(element_from_header));
	header[nbf].name = (char*) malloc((1+strlen(name))*sizeof(char));
	strcpy(header[nbf].name, name);
	header[nbf].value = (char*) malloc((1+strlen(value))*sizeof(char));
	strcpy(header[nbf].value, value);
	sprintf(header[0].value, "%d", nbf);
	return(header);
}


void free_header_table(element_from_header* header)
{
	int i;
	int nbf = atoi(header[0].value);

	for (i = 0; i <= nbf; i++)
	{
		free((header[i]).name);
	    free((header[i]).value);
	}
	free(header);
}


char* getItemFromHeader(char* name, element_from_header* header)
{
	char* value = 0;
	int nbf;
	int i;
	nbf = atoi(header[0].value);
	for (i = 1; i <= nbf; i++)
	{
		if (strcmp(header[i].name,name)==0)
			value = header[i].value;
	}
	return value;
}


void changeValue(element_from_header* header, char* name, char* newValue)
{
	int i;
	int nbf = atoi(header[0].value);

	for (i = 1; i <= nbf; i++)
	{
		if (strcmp(header[i].name, name)==0)
		{
			header[i].value = realloc(header[i].value, (1+strlen(newValue))*sizeof(char));
			strcpy(header[i].value, newValue);
		}
	}
}


/// For hash table ///

/* hash: form hash value for string s */
unsigned hash(char *s)
{
    unsigned hashval;
    for (hashval = 0; *s != '\0'; s++)
      hashval = *s + 31 * hashval;
    return hashval % HASHSIZE;
}

/* lookup: look for s in hashtab */
struct hashtable* lookup(struct hashtable *hashtab, char *s)
{
    struct hashtable* np;
    for (np = hashtab+hash(s); np != NULL; np = np->next)
    {
        if (strcmp(s, np->name) == NULL)
        {
        	return np; /* found */
        }
    }
    return NULL; /* not found */
}

/* install: put (name, defn) in hashtab */
struct hashtable* install(struct hashtable *hashtab, char *name, char *defn)
{
    struct hashtable* np;
    unsigned hashval;

    if ((np = lookup(hashtab, name)) == NULL)
    { /* not found */
        np = (struct hashtable*) malloc(sizeof(*np));
        if (np == NULL || (np->name = strdup(name)) == NULL)
          return NULL;
        hashval = hash(name);
        np->next = &hashtab[hashval];
        np->defn = defn;
        hashtab[hashval] = *np;
    }
    else /* already there */
        free((void *) np->defn); /*free previous defn */

    if ((np->defn = strdup(defn)) == NULL)
       return NULL;

    return np;
}
/// End for hash table ///


//char* getValueFromHashTable(element_from_header* header, char* dic_name, char* key_name)
//{
//	char* value = 0;
//	int nbf;
//	int i;
//	struct hashtable* dic;
//
//	nbf = atoi(header[0].value);
//	for (i = 1; i <= nbf; i++)
//	{
//		if (strcmp(header[i].name, dic_name)==0)
//		{
//			dic = (struct nlist*) header[i].value;
//			value = (lookup(dic, key_name))->defn;
//		}
//	}
//	return value;
//}


//element_from_header* table_header_add_dic(element_from_header* header, char* name, struct hashtable *hashtab)
//{
//	int nbf;
//	nbf = atoi(header[0].value);
//	nbf++;
//	header = (element_from_header*) realloc(header, (nbf+1)*sizeof(element_from_header));
//	header[nbf].name = (char*) malloc((1+strlen(name))*sizeof(char));
//	strcpy(header[nbf].name, name);
//	header[nbf].value = hashtab;
//	sprintf(header[0].value, "%d", nbf);
//	return(header);
//}


