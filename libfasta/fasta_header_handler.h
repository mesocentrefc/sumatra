
#ifndef FASTA_HEADER_HANDLER_H_
#define FASTA_HEADER_HANDLER_H_

#define HASHSIZE 101

#include "sequence.h"

struct hashtable { /* table entry: */
    struct hashtable *next; /* next entry in chain */
    char *name; /* defined name */
    char *defn; /* replacement text */
};

struct hashtable* lookup(struct hashtable *hashtab, char *s);

struct hashtable* install(struct hashtable *hashtab, char *name, char *defn);

char* char_header_add_field(char*,char*,char*);

char* fastaSeqPtr_header_add_field(fastaSeqPtr seq, char* name, char* value);

element_from_header* table_header_add_dic(element_from_header* header, char* name, struct hashtable *hashtab);

element_from_header* table_header_add_field(element_from_header* header, char* name, char* value);

void free_header_table(element_from_header*);

char* getItemFromHeader(char*, element_from_header*);

void changeValue(element_from_header* header, char* name, char* newValue);

#endif
