######################################################################
#	Makefile for MetaPSICOV                                      #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= cc

CFLAGS		= -O
LIBS		= -lm

all:		metapsicov metapsicovp2 metapsicovhb alnstats solvpred

install:
		/bin/cp -f metapsicov metapsicovp2 metapsicovhb alnstats solvpred ../bin

clean:
		/bin/rm -f metapsicov metapsicovp2 metapsicovhb alnstats solvpred *.o

solvpred:	solvpred.c solvation_net.h
		$(CC) $(CFLAGS) solvpred.c $(LIBS) -o solvpred

metapsicov:	metapsicov2_pred.c metapsicov2_net.h
		$(CC) $(CFLAGS) metapsicov2_pred.c $(LIBS) -o metapsicov

metapsicovp2:	metapsicov_pass2_pred.c metapsicov_pass2_net.h
		$(CC) $(CFLAGS) metapsicov_pass2_pred.c $(LIBS) -o metapsicovp2

metapsicovhb:	metapsicovhb.c metapsicov_pass2_net.h
		$(CC) $(CFLAGS) metapsicovhb.c $(LIBS) -o metapsicovhb

alnstats:	alnstats.c
		$(CC) $(CFLAGS) alnstats.c $(LIBS) -o alnstats
