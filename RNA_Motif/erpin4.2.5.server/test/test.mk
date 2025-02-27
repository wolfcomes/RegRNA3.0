SHELL = /bin/sh

#
# Fichier 'test.mk', cree le 03/05/04                    A.Lambert
# Usage:  make -f test.mk ;           recompile les programmes du rep 'test'.
#
#==============================================================================

INCLUDE = ../include
LIBDIR = ../lib

CC = cc
CFLAGS = -O3 -Wall
INCL = -I${INCLUDE}
LIBS = -L${LIBDIR} -lrnaIV -lm

#==============================================================================

exes: convsstat convhstat

convsstat: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a convsstat.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

convhstat: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a convhstat.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;
