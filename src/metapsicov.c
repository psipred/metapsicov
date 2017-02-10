/* METAPSICOV - Neural Network Meta-prediction of Residue Contacts */

/* Written by David T. Jones */

/* Copyright (C) 2014 University College London - Created : January 2014 */
/* Original Neural Network code Copyright (C) 1990 David T. Jones */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#define MAXSEQLEN 5000

/* logistic 'squashing' function (+/- 1.0) */
#define logistic(x) (1.0F / (1.0F + (float)exp(-(x))))

#define MAX(x,y) ((x)>(y)?(x):(y))

#include "metapsicov_net.h"

char           *wtfnm;

int             nwtsum, fwt_to[TOTAL], lwt_to[TOTAL];
float           activation[TOTAL], bias[TOTAL], *weight[TOTAL];

float **cmutmat;

int             profile[MAXSEQLEN][20];

char seq[MAXSEQLEN];

struct entry
{
    char           *id, *seq, **contmap;
    float         **mimat, **mipmat, **potmat, **psicov, **evfold, **ccmpred, *cprob, *hprob, *eprob, *entropy, *psolv, **profile, effnseq;
    float         *aacomp, cmean, hmean, emean, smean, entropmean;
    int           length, nseq;
} target;

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};

void
err(char *s)
{
    fprintf(stderr, "%s\n", s);
}

void
fail(char *s)
{
    err(s);
    exit(1);
}

/* Allocate matrix */
void           *allocmat(int rows, int columns, int size)
{
    int             i;
    void          **p, *rp;

    rp = malloc(rows * sizeof(void *) + sizeof(int));

    if (rp == NULL)
	fail("allocmat: malloc [] failed!");

    *((int *)rp) = rows;

    p = rp + sizeof(int);

    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}

/* Allocate vector */
void           *allocvec(int columns, int size)
{
    void          *p;

    p = calloc(columns, size);

    if (p == NULL)
	fail("allocvec: calloc failed!");

    return p;
}

void
compute_output(void)
{
    int             i, j;
    float           netinp;

    for (i = NUM_IN; i < TOTAL; i++)
    {
	netinp = bias[i];

	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	    netinp += activation[j] * weight[i][j];

	/* Trigger neuron */
	activation[i] = logistic(netinp);
    }
}

/*
 * load weights - load all link weights from a disk file
 */
void
load_wts(char *fname)
{
    int             i, j;
    double          t;
    FILE           *ifp;

    if (!(ifp = fopen(fname, "r")))
	fail("Cannot open weight file!\n");

    /* Load input units to hidden layer weights */
    for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	{
	    fscanf(ifp, "%lf", &t);
	    weight[i][j] = t;
	}

    /* Load hidden layer to output units weights */
    for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
	for (j = fwt_to[i]; j < lwt_to[i]; j++)
	{
	    fscanf(ifp, "%lf", &t);
	    weight[i][j] = t;
	}

    /* Load bias weights */
    for (j = NUM_IN; j < TOTAL; j++)
    {
	fscanf(ifp, "%lf", &t);
	bias[j] = t;
    }

    fclose(ifp);
}

/* Initialize network */
void
init(void)
{
    int             i, j;

    for (i = NUM_IN; i < TOTAL; i++)
	if (!(weight[i] = calloc(TOTAL - NUM_OUT, sizeof(float))))
	  fail("init: Out of Memory!");

    /* Connect input units to hidden layer */
    for (i = NUM_IN; i < NUM_IN + NUM_HID; i++)
    {
	fwt_to[i] = 0;
	lwt_to[i] = NUM_IN;
    }

    /* Connect hidden units to output layer */
    for (i = NUM_IN + NUM_HID; i < TOTAL; i++)
    {
	fwt_to[i] = NUM_IN;
	lwt_to[i] = NUM_IN + NUM_HID;
    }
}

/* Convert AA letter to numeric code (0-20) */
int
aanum(ch)
    int             ch;
{
    static const int      aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
	20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Make 1st level prediction averaged over specified weight sets */
void
predict(int argc, char **argv)
{
    int             aa, i, j, k, n, winpos, winpos2, midpos, ws, seqsep;
    float **predsum;
    double prob, x;

    predsum = allocmat(target.length, target.length, sizeof(float));

    for (ws=8; ws<argc; ws++)
    {
	load_wts(argv[ws]);
	
	for (winpos = 0; winpos < target.length; winpos++)
	    for (winpos2 = winpos+5; winpos2 < target.length; winpos2++)
	    {
		for (j = 0; j < NUM_IN; j++)
		    activation[j] = 0.0;
		
		for (j = WINL; j <= WINR; j++)
		    if (j + winpos >= 0 && j + winpos < target.length)
		    {
			for (aa = 0; aa < 21; aa++)
			    activation[(j - WINL) * IPERGRP + aa] = target.profile[j + winpos][aa];
			activation[(j - WINL) * IPERGRP + 21] = target.cprob[j + winpos];
			activation[(j - WINL) * IPERGRP + 22] = target.hprob[j + winpos];
			activation[(j - WINL) * IPERGRP + 23] = target.eprob[j + winpos];
			activation[(j - WINL) * IPERGRP + 24] = target.psolv[j + winpos];
			activation[(j - WINL) * IPERGRP + 25] = target.entropy[j + winpos];
		    }
		    else
			activation[(j - WINL) * IPERGRP + 26] = 1.0;
		
		for (j = WINL; j <= WINR; j++)
		    if (j + winpos2 >= 0 && j + winpos2 < target.length)
		    {
			for (aa = 0; aa < 21; aa++)
			    activation[(WINR-WINL+1) * IPERGRP + (j - WINL) * IPERGRP + aa] = target.profile[j + winpos2][aa];
			activation[(WINR-WINL+1) * IPERGRP + (j - WINL) * IPERGRP + 21] = target.cprob[j + winpos2];
			activation[(WINR-WINL+1) * IPERGRP + (j - WINL) * IPERGRP + 22] = target.hprob[j + winpos2];
			activation[(WINR-WINL+1) * IPERGRP + (j - WINL) * IPERGRP + 23] = target.eprob[j + winpos2];
			activation[(WINR-WINL+1) * IPERGRP + (j - WINL) * IPERGRP + 24] = target.psolv[j + winpos2];
			activation[(WINR-WINL+1) * IPERGRP + (j - WINL) * IPERGRP + 25] = target.entropy[j + winpos2];
		    }
		    else
			activation[(WINR-WINL+1) * IPERGRP + (j - WINL) * IPERGRP + 26] = 1.0;
		
		midpos = (winpos + winpos2) / 2;
		
		for (j = CWINL; j <= CWINR; j++)
		    if (j + midpos >= 0 && j + midpos < target.length)
		    {
			for (aa = 0; aa < 21; aa++)
			    activation[2*(WINR-WINL+1) * IPERGRP + (j - CWINL) * IPERGRP + aa] = target.profile[j + midpos][aa];
			activation[2*(WINR-WINL+1) * IPERGRP + (j - CWINL) * IPERGRP + 21] = target.cprob[j + midpos];
			activation[2*(WINR-WINL+1) * IPERGRP + (j - CWINL) * IPERGRP + 22] = target.hprob[j + midpos];
			activation[2*(WINR-WINL+1) * IPERGRP + (j - CWINL) * IPERGRP + 23] = target.eprob[j + midpos];
			activation[2*(WINR-WINL+1) * IPERGRP + (j - CWINL) * IPERGRP + 24] = target.psolv[j + midpos];
			activation[2*(WINR-WINL+1) * IPERGRP + (j - CWINL) * IPERGRP + 25] = target.entropy[j + midpos];
		    }
		    else
			activation[2*(WINR-WINL+1) * IPERGRP + (j - CWINL) * IPERGRP + 26] = 1.0;
		
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP] = target.mimat[winpos][winpos2];
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP+1] = target.mipmat[winpos][winpos2];
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP+2] = target.potmat[winpos][winpos2];
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP+3] = target.psicov[winpos][winpos2];
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP+4] = target.evfold[winpos][winpos2];
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP+5] = target.ccmpred[winpos][winpos2];
		
		seqsep = winpos2 - winpos;
		
		if (seqsep < 5)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 6] = 1.0;
		else if (seqsep < 14)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 2 + seqsep] = 1.0;
		else if (seqsep < 18)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 16] = 1.0;
		else if (seqsep < 23)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 17] = 1.0;
		else if (seqsep <= 28)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 18] = 1.0;
		else if (seqsep <= 38)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 19] = 1.0;
		else if (seqsep <= 48)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 20] = 1.0;
		else
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 21] = 1.0;
		
		for (aa=0; aa<21; aa++)
		    activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + aa] = target.aacomp[aa];
		
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 21] = target.cmean;
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 22] = target.hmean;
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 23] = target.emean;
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 24] = target.smean;
		
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 25] = log((double)target.length);
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 26] = log((double)target.nseq);
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 27] = log((double)target.effnseq);
		activation[2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 28] = target.entropmean;
		
#if (2*(WINR-WINL+1) * IPERGRP + (CWINR-CWINL+1) * IPERGRP + 22 + 29) != NUM_IN
#error "MISMATCHING NUMBER OF INPUTS!"
#endif
		
		compute_output();
		
		predsum[winpos][winpos2] += activation[TOTAL - NUM_OUT];
	    }
    }

    for (winpos = 0; winpos < target.length; winpos++)
	for (winpos2 = winpos+5; winpos2 < target.length; winpos2++)
	{
	    predsum[winpos][winpos2] /= (float)(argc - 8);
	    
	    x = predsum[winpos][winpos2];

#if 0
	    /* Transform to posterior prob by scaled log transform */

#define A 3.6512190000514688E-01
#define B 1.5491157018510503E+01
#define C 9.8856650825604842E-01	    

	    prob = A * log(B * x + C);
	    
	    if (prob < 0.001)
		prob = 0.001;
	    if (prob > 0.999)
		prob = 0.999;

	    printf("%d %d 0 8 %f\n", winpos+1, winpos2+1, prob);
#else
	    printf("%d %d 0 8 %f\n", winpos+1, winpos2+1, x);
#endif
	}
}

/* Read PSI AA frequency data */
int             getmtx(FILE *lfil)
{
    int             aa, i, j, naa;
    char            buf[256], *p;

    if (fscanf(lfil, "%d", &naa) != 1)
      fail("Bad mtx file - no sequence length!");

    if (naa > MAXSEQLEN)
      fail("Input sequence too long!");

    if (fscanf(lfil, "%s", seq) != 1)
      fail("Bad mtx file - no sequence!");

    while (!feof(lfil))
      {
	if (!fgets(buf, 65536, lfil))
	  fail("Bad mtx file!");
	if (!strncmp(buf, "-32768 ", 7))
	  {
	    for (j=0; j<naa; j++)
	      {
		if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &profile[j][ALA],  &profile[j][CYS], &profile[j][ASP],  &profile[j][GLU],  &profile[j][PHE],  &profile[j][GLY],  &profile[j][HIS],  &profile[j][ILE],  &profile[j][LYS],  &profile[j][LEU],  &profile[j][MET],  &profile[j][ASN],  &profile[j][PRO],  &profile[j][GLN],  &profile[j][ARG],  &profile[j][SER],  &profile[j][THR],  &profile[j][VAL],  &profile[j][TRP],  &profile[j][TYR]) != 20)
		  fail("Bad mtx format!");
		aa = aanum(seq[j]);
		if (!fgets(buf, 65536, lfil))
		  break;
	      }
	  }
      }

    return naa;
}

/* Read in data */
void            read_dat(char *colname, char *pairname, char *psiname, char *evfname, char *ccmpfname, char *ssname, char *solvname)
{
    FILE           *ifp, *tfp, *cfp, *pfp, *sfp;
    char            buf[512];
    int             i, j, k, l, m, ncon, nncon, nres, npairs, ctoggle = 0;
    float val, consv[22], aventropy;

    if (!(cfp = fopen(colname, "r")))
	fail("Cannot open column stats file!");
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    sscanf(buf, "%d", &nres);
    
    target.length = nres;
    target.cprob = allocvec(nres, sizeof(float));
    target.hprob = allocvec(nres, sizeof(float));
    target.eprob = allocvec(nres, sizeof(float));
    target.entropy = allocvec(nres, sizeof(float));
    target.psolv = allocvec(nres, sizeof(float));
    target.aacomp = allocvec(21, sizeof(float));
    target.profile = allocmat(nres, 21, sizeof(float));

    npairs = nres * (nres - 1) / 2;
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    sscanf(buf, "%d", &target.nseq);
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    sscanf(buf, "%f", &target.effnseq);
    
    if (!fgets(buf, 512, cfp))
	fail("Bad column stats file!");
    
    if (sscanf(buf, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
	       &consv[0], &consv[1], &consv[2], &consv[3], &consv[4],
	       &consv[5], &consv[6], &consv[7], &consv[8], &consv[9],
	       &consv[10], &consv[11], &consv[12], &consv[13],
	       &consv[14], &consv[15], &consv[16], &consv[17],
	       &consv[18], &consv[19], &consv[20]) != 21)
	fail("Bad AAcomp records in column stats file!");
    
    for (j = 0; j < 21; j++)
	target.aacomp[j] = consv[j];
    
    for (i=0; i<nres; i++)
    {
	if (!fgets(buf, 512, cfp))
	    fail("Bad column stats file!");
	
	if (sscanf(buf, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
		   &consv[0], &consv[1], &consv[2], &consv[3], &consv[4],
		   &consv[5], &consv[6], &consv[7], &consv[8], &consv[9],
		   &consv[10], &consv[11], &consv[12], &consv[13],
		   &consv[14], &consv[15], &consv[16], &consv[17],
		   &consv[18], &consv[19], &consv[20], &consv[21]) != 22)
	{
	    puts(buf);
	    fail("Bad profile records in column stats file!");
	}
	
	for (aventropy = j = 0; j < 21; j++)
	    target.profile[i][j] = consv[j];
	aventropy += (target.entropy[i] = consv[21]);
    }
    
    fclose(cfp);

    aventropy /= nres;
    
    target.entropmean = aventropy;
    
    if (!(pfp = fopen(pairname, "r")))
	fail("Cannot open pair stats file!");
    
    target.potmat = allocmat(nres, nres, sizeof(float));
    target.mimat = allocmat(nres, nres, sizeof(float));
    target.mipmat = allocmat(nres, nres, sizeof(float));
    
    for (k=0; k<npairs; k++)
    {
	if (!fgets(buf, 512, pfp))
	    fail("Bad pair stats file!");
	
	if (sscanf(buf, "%d%d%f%f%f", &i, &j, &consv[0], &consv[1], &consv[2]) != 5)
	    fail("Bad pair stats file!");
	
	target.potmat[i-1][j-1] = target.potmat[j-1][i-1] = consv[0];
	target.mimat[i-1][j-1] = target.mimat[j-1][i-1] = consv[1];
	target.mipmat[i-1][j-1] = target.mipmat[j-1][i-1] = consv[2];
    }
    
    fclose(pfp);

    target.psicov = allocmat(nres, nres, sizeof(float));
    
    if ((pfp = fopen(psiname, "r")))
    {
	for (;;)
	{
	    if (!fgets(buf, 512, pfp))
		break;
	    
	    if (sscanf(buf, "%d%d%*s%*s%f", &i, &j, &consv[0]) != 3)
		fail("Bad PSICOV file!");
	    
	    target.psicov[i-1][j-1] = target.psicov[j-1][i-1] = consv[0];
	}
    
	fclose(pfp);
    }
    
    target.evfold = allocmat(nres, nres, sizeof(float));
    
    if ((pfp = fopen(evfname, "r")))
    {
	for (;;)
	{
	    if (!fgets(buf, 512, pfp))
		break;
	    
	    if (sscanf(buf, "%d%*s%d%*s%*s%f", &i, &j, &consv[0]) != 3)
		fail("Bad evfold file!");
	    
	    target.evfold[i-1][j-1] = target.evfold[j-1][i-1] = consv[0];
	}
    
	fclose(pfp);
    }
    
    target.ccmpred = allocmat(nres, nres, sizeof(float));

    if ((pfp = fopen(ccmpfname, "r")))
    {
	if (fscanf(pfp, "%f", &consv[0]) == 1)
	{
	    rewind(pfp);
	    
	    for (i=0; i<nres; i++)
		for (j=0; j<nres; j++)
		{
		    if (fscanf(pfp, "%f", &consv[0]) != 1)
			fail("Bad CCMPRED file!");
		    
		    /* NOTE: Reading full NxN matrix for CCMPRED! */
		    target.ccmpred[i][j] = consv[0];
		}
	}
    
	fclose(pfp);
    }
    
    if (!(sfp = fopen(ssname, "r")))
	fail("Cannot open SS2 file!");
    
    if (!fgets(buf, 512, sfp) || !fgets(buf, 512, sfp))
	fail("Bad SS2 file!");
    
    i = 0;
    while (!feof(sfp))
    {
	if (!fgets(buf, 512, sfp))
	    break;
	if (sscanf(buf, "%*s%*s%*s%f%f%f", target.cprob+i, target.hprob+i, target.eprob+i) != 3)
	    fail("Bad SS2 record!");
	target.cmean += target.cprob[i];
	target.hmean += target.hprob[i];
	target.emean += target.eprob[i];
	i++;
    }
    
    fclose(sfp);
    
    nres = i;
    
    if (!(sfp = fopen(solvname, "r")))
	fail("Cannot open solv file!");
    
    i = 0;
    while (!feof(sfp))
    {
	if (!fgets(buf, 512, sfp))
	    break;
	if (sscanf(buf, "%*s%*s%f", target.psolv+i) != 1)
	    fail("Bad solv record!");
	target.smean += target.psolv[i];
	i++;
    }
    
    fclose(sfp);
    
    if (i != nres)
	fail("Solv length mismatch!");
    
    target.cmean /= (float)nres;
    target.hmean /= (float)nres;
    target.emean /= (float)nres;
    target.smean /= (float)nres;
}

main(int argc, char **argv)
{
    int             i, j, niters;
    FILE *ifp;

    /* malloc_debug(3); */
    if (argc < 9)
	fail("usage : metapsicov colstats-file pairstats-file psicov-file evfold-file ccmpred-file ss2-file solv-file weight-file1 ... weight-filen");

    read_dat(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);

    init();

    predict(argc,argv);

    return 0;
}
