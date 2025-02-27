#################################################################
# Infernal
# testsuite Makefile
# SVN $Id: Makefile.in 1717 2006-11-29 11:36:24Z nawrockie $
#################################################################
#    This copyrighted source code is freely distributed 
#    under the terms of the GNU General Public License. See
#    the files COPYRIGHT and LICENSE for details.
#################################################################

SHELL  = /bin/sh

## your compiler and compiler flags;
#  and other configuration of the compilation and linking
#  environment.
#  DEFS contains stuff that autoconf set up.
#
CC       = gcc
CFLAGS   = -g -O2
CPPFLAGS = 
LDFLAGS  = 
DEFS     = -DHAVE_CONFIG_H
LIBS     =  -lm 
MYLIBS   = -linfernal -lsquid -leasel

# Where the build arenas are
#
MYLIBDIR = -L../squid -L../easel/ -L../src
MYINCDIR = -I../squid -I../easel/ -I../src 

TESTPROGS = \
	bandcyk-montecarlo-test\
	bandcyk-truncation-test\
	cp9-test\
	sub_cm-test

#################################################################
## Building the test programs.
##

.c.o: 
	${CC} ${CFLAGS} ${DEFS} ${MYINCDIR} -c $<

all:	${TESTPROGS}

${TESTPROGS}: %: %.o 
	${CC} ${CFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIR} -o $@ $@.o $(MYLIBS) $(LIBS)



#################################################################
## Miscellaneous targets.
##
clean:
	-rm -f *.o *~ Makefile.bak core ${TESTPROGS} TAGS gmon.out 

distclean:
	-rm -f *.o *~ Makefile.bak core ${TESTPROGS} TAGS gmon.out 
	-rm -f Makefile 

TAGS:
	etags -t *.c *.h Makefile.in

