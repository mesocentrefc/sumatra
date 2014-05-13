
SOURCES = fasta_header_parser.c \
		  fasta_seq_writer.c \
		  fasta_header_handler.c \
		  header_mem_handler.c \
          sequence.c \
          dic_parser.c
         
SRCS=$(SOURCES)

         
OBJECTS= $(patsubst %.c,%.o,$(SOURCES))

LIBFILE = libfasta.a
RANLIB = ranlib


include ../global.mk

all: $(LIBFILE)

fasta_header_parser.c: fasta_header_parser.l
	lex -Pheader_yy -t $< > $@

dic_parser.c: dic_parser.l
	lex -Phashtable_yy -t $< > $@

clean:
	rm -rf $(OBJECTS) $(LIBFILE)
	rm -f *.a

$(LIBFILE): $(OBJECTS)
	ar -cr $@ $?
	$(RANLIB) $@
