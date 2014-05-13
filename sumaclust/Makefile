EXEC=sumaclust

SUMACLUST_SRC= sumaclust.c \
               mtcompare_sumaclust.c

             
SUMACLUST_OBJ= $(patsubst %.c,%.o,$(SUMACLUST_SRC))


SRCS= $(SUMACLUST_SRC)

LIB= -lfasta -llcs -lfile -lutils -lm #-ll

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
# sumaclust compilation
#
########

# executable compilation and link

sumaclust: $(SUMACLUST_OBJ) $(LIBFASTA) $(LIBLCS) $(LIBFILE) $(LIBUTILS)
	$(CC) $(LDFLAGS) -o $@ -fopenmp $(SUMACLUST_OBJ) $(LIB) $(LIBFASTAPATH) $(LIBLCSPATH) $(LIBFILEPATH) $(LIBUTILSPATH)
	
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

		
