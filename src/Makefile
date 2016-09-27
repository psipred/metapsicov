######################################################################
#	Makefile for PSIpred4                                        #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= cc

CFLAGS		= -O
LIBS		= -lm

all:		chkparse solvpred alnstats

install:
	/bin/cp chkparse ../bin
	/bin/cp solvpred ../bin
	/bin/cp alnstats ../bin

clean:
	/bin/rm -f chkparse
	/bin/rm -f solvpred
	/bin/rm -f alnstats

chkparse:
	$(CC) $(CFLAGS) chkparse.c $(LIBS) -o chkparse

solvpred:
	$(CC) $(CFLAGS) solvpred_avpred.c $(LIBS) -o solvpred

alnstats:
	$(CC) $(CFLAGS) alnstats.c $(LIBS) -o alnstats
