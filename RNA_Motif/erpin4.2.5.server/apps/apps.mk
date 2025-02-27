SHELL = /bin/sh

# Fichier 'apps.mk'      cree le 02/01/02      A.Lambert  CPT
# Usage:    make -f apps.mk ;                  recompile les programmes, et
#                                              place les binaires dans ../bin
#
#------------------------------------------------------------------------------

INCLUDE = ../include
LIBDIR = ../lib
BINDIR = ../bin

CC = cc
CFLAGS = -O3 -Wall
INCL = -I${INCLUDE}
LIBS = -L${LIBDIR} -lrnaIV -lm

BINS = erpin tview tstat mstat cfgs sview tstrip mhistview pview ev frandseq \
       epnstat

#------------------------------------------------------------------------------

exes:  erpin tview tstat mstat cfgs sview tstrip mhistview pview ev frandseq \
       epnstat install

erpin: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  erpin.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

tview: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  tview.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

tstat: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  tstat.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

mstat: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  mstat.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

cfgs:  ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  cfgs.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

sview:  ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  sview.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

tstrip: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  tstrip.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

mhistview: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  mhistview.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

pview: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  pview.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

ev: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  ev.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

epnstat: ${INCLUDE}/rnaIV.h ${LIBDIR}/librnaIV.a  epnstat.c
	${CC} ${CFLAGS} ${INCL} -o $@ $@.c ${LIBS} ; \
	strip $@ ; \
	chmod 755 $@ ;

frandseq: frandseq.c
	${CC} ${CFLAGS} -o $@ $@.c ; \
	strip $@ ; \
	chmod 755 $@ ;

install:
	mv ${BINS} ${BINDIR} ;
