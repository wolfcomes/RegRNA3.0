SHELL = /bin/sh

#
# Fichier 'erpin.mk'           cree le 03/05/04, A.Lambert
# Usage:  make -f erpin.mk ;   recompile la distribution 'erpin'.
#

erpin:
	cd libsrc ; \
	make -f rnaIV.mk ; \
	make -f rnaIV.mk install ; \
	make -f rnaIV.mk clean ;  \
	cd ../apps ; \
	make -f apps.mk ; \
	cd ../sum ; \
	make -f sum.mk ; \
	cd ../test ; \
	make -f test.mk ;
