
#ifndef FASTA_HEADER_PARSER_H_
#define FASTA_HEADER_PARSER_H_

typedef struct {
	char *name;
	void *value;
}element_from_header;

element_from_header* header_parser_main(char*);


#endif
