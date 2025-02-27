#################################################################
# Infernal
# src Makefile
# SVN $Id: Makefile.in 1970 2007-04-27 17:57:09Z nawrockie $
#################################################################
#    This copyrighted source code is freely distributed 
#    under the terms of the GNU General Public License. See
#    the files COPYRIGHT and LICENSE for details.
#################################################################

SHELL  = /bin/sh

## your compiler and compiler flags
#
CC     = gcc
CFLAGS = -g -O2

## other defined flags. 
#  DEFS contains stuff that autoconf decides on. 
#  LIBS contains system libraries that the configure script decides we need.
#
DEFS  = -DHAVE_CONFIG_H
LIBS  =  -lm 

## archiving command, and ranlib command if you need one.
#  In general, you shouldn't need to change these, and they're
#  only used for building the testsuite anyway... e.g. we
#  make a "libinfernal.a" library for building the testsuite.
#
AR     = /usr/bin/ar rcv 
RANLIB = ranlib

# configuration for optional MPI functionality
MPICC     = 
MPILIBS   = 
MPIPROGS  = 

# Where the build arena is (module .h and .o files)
#
MYLIBS   = -lsquid -leasel 
MYLIBDIR = -L../squid -L../easel 
MYINCDIR = -I../squid -I../easel 

PROGS = cmalign cmbuild cmscore cmsearch cmemit

OBJS  = alphabet.o\
	bandcyk.o\
	cm.o\
	cmio.o\
	display.o\
	emit.o\
	emitmap.o\
	globals.o\
	modelconfig.o\
	modelmaker.o\
	parsetree.o\
	prior.o\
	scancyk.o\
	smallcyk.o\
	truncyk.o\
	cm_eweight.o\
	hmmband.o\
	hbandcyk.o\
	cplan9.o\
	CP9_cm2wrhmm.o\
	CP9_hmmio.o\
	CP9_scan.o\
	cm_postprob.o\
	scaninside.o\
	cm_masks.o\
	sub_cm.o\
	stats.o\
	mpifuncs.o\
	cm_dispatch.o\
	rnamat.o\
	cm_cluster.o

HDRS  = funcs.h\
	prior.h\
	structs.h\
	hmmband.h\
	cplan9.h\
	cm_postprob.h\
	stats.h\
	mpifuncs.h

.c.o: 
	$(CC) $(CFLAGS) $(DEFS) $(MYINCDIR) -c $<

#################################################################
## Targets defining how to make INFERNAL executables.
##
all: 	$(PROGS) $(MPIPROGS)

$(PROGS): %: %.o $(OBJS) ${HDRS}
	$(CC) $(CFLAGS) $(DEFS) $(MYLIBDIR) -o $@ $@.o $(OBJS) $(MYLIBS) $(LIBS)

mpi-cmalign: cmalign.c $(OBJS) ${HDRS}
	$(CC) -DMPI_EXECUTABLE $(CFLAGS) $(DEFS) $(MYINCDIR) -c cmalign.c
	$(CC) -DMPI_EXECUTABLE $(CFLAGS) $(DEFS) $(MYLIBDIR) -o mpi-cmalign cmalign.o $(OBJS) $(MYLIBS) $(LIBS)

mpi-cmsearch: cmsearch.c $(OBJS) ${HDRS}
	$(CC) -DMPI_EXECUTABLE $(CFLAGS) $(DEFS) $(MYINCDIR) -c cmsearch.c
	$(CC) -DMPI_EXECUTABLE $(CFLAGS) $(DEFS) $(MYLIBDIR) -o mpi-cmsearch cmsearch.o $(OBJS) $(MYLIBS) $(LIBS)

#################################################################
## Targets used in making INFERNAL module for testsuite compilation.
##
module: libinfernal.a

libinfernal.a: $(OBJS) ${HDRS}
	$(AR) libinfernal.a $(OBJS) 
	$(RANLIB) libinfernal.a
	chmod 644 libinfernal.a


#################################################################
## Miscellaneous targets.
##
clean:
	-rm -f *.o *~ Makefile.bak core $(PROGS) $(MPIPROGS) TAGS gmon.out libinfernal.a

distclean:
	-rm -f *.o *~ Makefile.bak core $(PROGS) $(MPIPROGS) TAGS gmon.out libinfernal.a
	-rm -f Makefile config.h

TAGS:
	etags -t *.c *.h Makefile.in

