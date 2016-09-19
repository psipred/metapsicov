######################################################################
#	Makefile for PSIpred4                                        #
######################################################################

#### CC		Command to execute C compiler
#### CFLAGS	Flags to pass to C compiler.

INCDIR		= .
CC		= cc

CFLAGS		= -O
LIBS		= -lm

all:		chkparse

install:
		/bin/cp chkparse ../bin

clean:
		/bin/rm -f chkparse

chkparse:	chkparse.c
		$(CC) $(CFLAGS) chkparse.c $(LIBS) -o chkparse
