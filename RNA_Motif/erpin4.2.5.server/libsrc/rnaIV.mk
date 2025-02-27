SHELL = /bin/sh

# Fichier 'rnaIV.mk'   cree le 02/01/02  A.Lambert CPT  revu le 12/02/04
# Usage:   make -f rnaIV.mk ;            recompile la librairie 'librnaIV.a'
#          make -f rnaIV.mk install ;    mv 'librnaIV.a' dans ./lib
#          make -f rnaIV.mk clean ;      efface les ".o"

#------------------------------------------------------------------------------

SRCDIR = ./
INCLUDE = ../include
LIBDIR = ../lib

CC = cc -O2 -Wall
INCL = -I${INCLUDE}

SRCS = \
env.c \
tab1.c \
tab2.c \
io.c \
trset.c \
atom.c \
pattern.c \
helix.c \
strand.c \
align.c \
profs.c \
profsSM.c \
sum.c \
scores.c \
sstat.c \
fstat.c \
cfg.c \
tscores.c \
cfgstr.c \
mask.c \
maskcfg.c \
masks.c \
mscores.c \
outputs.c \
Seqs.c \
thresholds.c \
msearch.c \
args.c \
list.c \
dmp.c \
ctrlcfgs.c \
maps.c \
dhisto.c \
histools.c \
hshisto.c \
mhisto.c \
conv.c \
cdf.c \
lfit.c \
bkgstat.c \
Eval.c \
ntcode.c

OBJS = \
env.o \
tab1.o \
tab2.o \
io.o \
trset.o \
atom.o \
pattern.o \
helix.o \
strand.o \
align.o \
profs.o \
profsSM.o \
sum.o \
scores.o \
sstat.o \
fstat.o \
cfg.o \
tscores.o \
cfgstr.o \
mask.o \
maskcfg.o \
masks.o \
mscores.o \
outputs.o \
Seqs.o \
thresholds.o \
msearch.o \
args.o \
list.o \
dmp.o \
ctrlcfgs.o \
maps.o \
dhisto.o \
histools.o \
hshisto.o \
mhisto.o \
conv.o \
cdf.o \
lfit.o \
bkgstat.o \
Eval.o \
ntcode.o

#-------------- regles de dependances et commandes -----------

librnaIV.a: ${SRCS} ${INCLUDE}/rnaIV.h
	${CC} -c ${SRCS} ${INCL} ; \
	ar -rs ./$@ ${OBJS} ; \
	chmod 755 ./$@ ;

install:
	mv librnaIV.a ${LIBDIR} ;

clean:
	rm ${OBJS} ;




