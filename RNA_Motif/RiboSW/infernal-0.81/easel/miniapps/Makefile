################################################################
# Makefile.in for Easel's miniapps
# SVN $Id$
SHELL = /bin/sh

## your compiler and compiler flags.
#
CC       = gcc
CFLAGS   = -g -O2
DEFS     = -DHAVE_CONFIG_H
LDFLAGS  = 
LIBS     = -leasel  -lm

## where is Easel?
#
ESLLIBDIR = -L..
ESLINCDIR = -I..


## list of the miniapps to compile.
#
PROGS = esl-compstruct\
	esl-histplot\
	esl-reformat

## Targets defining how to make miniapp executables.
#
all: 	$(PROGS) 

${PROGS}: %: %.o  ../libeasel.a
	${CC} ${CFLAGS} ${DEFS} ${LDFLAGS} ${ESLLIBDIR} -o $@ $@.o ${OBJS} ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${DEFS} ${ESLINCDIR} -c $<

clean:
	-rm *.o *~

distclean:
	make clean
	-rm -f ${PROGS}
	-rm Makefile

################################################################
#    This copyrighted source code is freely distributed 
#    under the terms of the GNU General Public License. See
#    the files COPYRIGHT and LICENSE for details.
################################################################