EXEC=sumatra

SUMATRA_SRC= sumatra.c \
			 mtcompare_sumatra.c
			              
SUMATRA_OBJ= $(patsubst %.c,%.o,$(SUMATRA_SRC))


SRCS= $(SUMATRA_SRC)

LIB= -lz -lm -lfasta -llcs -lfile -lutils

LIBFASTA = ../libfasta/libfasta.a
LIBLCS   = ../liblcs/liblcs.a
LIBFILE  = ../libfile/libfile.a
LIBUTILS = ../libutils/libutils.a

include ../global.mk

all: $(EXEC)

########
#
# libraries compilation
#
########

../libfasta/libfasta.a:
	$(MAKE) -C ../libfasta

../liblcs/liblcs.a:
	$(MAKE) -C ../liblcs

../libfile/libfile.a:
	$(MAKE) -C ../libfile

../libutils/libutils.a:
	$(MAKE) -C ../libutils
	
	
########
#
# sumatra compilation
#
########

# executable compilation and link

sumatra: $(SUMATRA_OBJ) $(LIBFASTA) $(LIBLCS) $(LIBFILE) $(LIBUTILS)
	$(CC) $(LDFLAGS) -o $@ -pthread $(SUMATRA_OBJ) $(LIBFASTAPATH) $(LIBLCSPATH) $(LIBFILEPATH) $(LIBUTILSPATH) $(LIB)
	
########
#
# project management
#
########

clean:
	rm -f *.o	
	rm -f *.P
	rm -f $(EXEC)
	$(MAKE) -C ../libfasta clean
	$(MAKE) -C ../liblcs clean
	$(MAKE) -C ../libfile clean
	$(MAKE) -C ../libutils clean

		
