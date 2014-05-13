
SOURCES = utilities.c  \
          debug.c      
         
SRCS=$(SOURCES)

         
OBJECTS= $(patsubst %.c,%.o,$(SOURCES))

LIBFILE= libutils.a
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
