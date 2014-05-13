
LIBFASTAPATH = -L../libfasta
LIBLCSPATH   = -L../liblcs
LIBFILEPATH  = -L../libfile
LIBUTILSPATH = -L../libutils

CC=gcc
CFLAGS= -O3 -s -fopenmp -w
LDFLAGS=

default: all

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<
