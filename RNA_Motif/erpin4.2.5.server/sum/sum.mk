SHELL = /bin/sh

#
# Fichier 'sum.mk', cree le 03/05/04                    A.Lambert
# Usage:  make -f sum.mk ;           recompile les programmes du rep 'sum'.
#
#==============================================================================

INCLUDE = ../include
LIBDIR = ../lib

CC = cc
CFLAGS = -O3 -Wall
INCL = -I${INCLUDE}
LIBS = -L${LIBDIR} -lrnaIV -lm

#==============================================================================

exes: tview2 tstrip2 mksum

tview2: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a tview2.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

tstrip2: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a tstrip2.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

mksum: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a mksum.c stsum.c hlxsum.c
	cc -O2 -Wall ${INCL} -c stsum.c hlxsum.c ; \
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c stsum.o hlxsum.o ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ; \
	rm *.o ;
