######################################################################
#	Makefile for PSIpred4                                        #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= cc

CFLAGS		= -O
LIBS		= -lm

all:		chkparse solvpred alnstats psicov metapsicov metapsicovp2

install:
	/bin/cp chkparse ../bin
	/bin/cp solvpred ../bin
	/bin/cp alnstats ../bin
	/bin/cp psicov ../bin
	/bin/cp metapsicov ../bin
	/bin/cp metapsicovp2 ../bin

clean:
	/bin/rm -f chkparse
	/bin/rm -f solvpred
	/bin/rm -f alnstats
	/bin/rm -f psicov
	/bin/rm -f metapsicov
	/bin/rm -f metapsicovp2

chkparse:
	$(CC) $(CFLAGS) chkparse.c $(LIBS) -o chkparse

solvpred:
	$(CC) $(CFLAGS) solvpred_avpred.c $(LIBS) -o solvpred

alnstats:
	$(CC) $(CFLAGS) alnstats.c $(LIBS) -o alnstats

psicov:
	$(CC) $(CFLAGS) psicov23.c $(LIBS) -o psicov

metapsicov:
	$(CC) $(CFLAGS) metapsicov2_pred.c $(LIBS) -o metapsicov

metapsicovp2:
	$(CC) $(CFLAGS) metapsicov_pass2_pred.c $(LIBS) -o metapsicovp2
