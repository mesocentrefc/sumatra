
SOURCES = sse_banded_LCS_alignment.c  \
		  upperband.c
         
SRCS=$(SOURCES)

         
OBJECTS= $(patsubst %.c,%.o,$(SOURCES))

LIBFILE= liblcs.a
RANLIB=ranlib


include ../global.mk

all: $(LIBFILE)

clean:
	rm -rf $(OBJECTS) $(LIBFILE)
	rm -f *.P
	rm -f *.a

$(LIBFILE): $(OBJECTS)
	ar -cr $@ $?
	$(RANLIB) $@
