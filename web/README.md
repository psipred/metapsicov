# Web process for MetaPSICOV

These are the details of how the web version of metapsicov works. This is
a translation of the run_metapsicov.pl script and omits the CASP server
additions (consip4)

# Prerequisites

1. HH-suite
2. HH databases pdb70 and uniprot20
3. Blast+
4. Uniref90 blastdb
5. hmmer and uniref100
6. Psipred (with solv pred)
7. cblas library https://github.com/xianyi/OpenBLAS pass CFLAGS CPPFLAG AND LDFLAGS
8. CCMPRED git clone --recursive https://github.com/soedinglab/CCMpred
9. Freecontact wget ftp://rostlab.org/free/freecontact-1.0.21.tar.xz
     ./configure CFLAGS=-I/opt/OpenBLAS/include LIBS=-L/opt/OpenBLAS/lib LDFLAGS=-L/opt/OpenBLAS/lib CPPFLAGS=-I/opt/OpenBLAS/include

# MetaPSICOV Protocol

Adjust the number of threads/cpus appropriately to your environment

1. Run psiblast
`/scratch0/NOT_BACKED_UP/dbuchan/Applications/ncbi-blast-2.2.31+/bin/psiblast -evalue 0.001 -num_descriptions 2000 -num_alignments 0 -num_threads 1 -num_iterations 3 -inclusion_ethresh 0.001 -query example/test.fsa -db /scratch1/NOT_BACKED_UP/dbuchan/uniref/uniref90.fasta -out test.bls -out_pssm test.chk`
2. Run chkparse
`/opt/Code/psipred/bin/chkparse test.chk > test.mtx`

3. Run psipred
`/opt/Code/psipred/bin/psipred test.mtx /opt/Code/psipred/data/weights.dat /opt/Code/psipred/data/weights.dat2 /opt/Code/psipred/data/weights.dat3 > test.ss`
4. Run psipass
`/opt/Code/psipred/bin/psipass2 /opt/Code/psipred/data/weights_p2.dat 1 1.0 1.0 test.ss2 test.ss > test.horiz`
5. Run solvpred
`/opt/Code/metapsicov/bin/solvpred test.mtx /opt/Code/metapsicov/data/weights_solv.dat > test.solv`

6. Run hhblits over uniprot20
`/scratch0/NOT_BACKED_UP/dbuchan/Applications/hh-suite/bin/hhblits -i example/5ptpA.fasta -d /scratch1/NOT_BACKED_UP/dbuchan/hhblitsdb/uniprot20_2016_02/uniprot20_2016_02 -oa3m test.a3m -n 3 -maxfilt 500000 -diff inf -id 99 -cov 60 > test.hhblog`
7. Send hhblits alignment lines to an alignment file
`egrep -v '^>' test.a3m | sed 's/[a-z]//g' > test.aln`
8.  Test how many alignments lines we have
`cat test.aln | wc -l`;
9. IF there are more than 5 sequences we run the contact predictors, steps 10 to 12
10. Run PSICOV
`timeout 86400 /opt/Code/psicov/bin/psicov -o -d 0.03 test.aln > test.psicov`
11. Run Freecontact
`/opt/Applications/freecontact-1.0.21/src/freecontact < test.aln > test.evfold`
12. Run CCMpred
`/opt/Applications/CCMpred/bin/ccmpred -t 24 test.aln test.ccmpred`
13. Now run the metapsicov predictor over the contactfasta
14. Run alnstats
`/opt/Code/metapsicov/bin/alnstats test.aln test.colstats test.pairstats`
15. Run metapsicov stage1
`/opt/Code/metapsicov/bin/metapsicov test.colstats test.pairstats test.psicov test.evfold test.ccmpred test.ss2 test.solv /opt/Code/metapsicov/data/weights_6A.dat /opt/Code/metapsicov/data/weights_65A.dat /opt/Code/metapsicov/data/weights_7A.dat /opt/Code/metapsicov/data/weights_75A.dat /opt/Code/metapsicov/data/weights_8A.dat /opt/Code/metapsicov/data/weights_85A.dat /opt/Code/metapsicov/data/weights_9A.dat /opt/Code/metapsicov/data/weights_10A.dat /opt/Code/metapsicov/data/weights_811A.dat /opt/Code/metapsicov/data/weights_1012A.dat | sort -n -r -k 5 > test.metapsicov.stage1`
16. Run metapsicov stage 2
`/opt/Code/metapsicov/bin/metapsicovp2 test.colstats test.metapsicov.stage1 test.ss2 test.solv /opt/Code/metapsicov/data/weights_pass2.dat | sort -n -r -k 5 > test.metapsicov.stage2`
17. Run metapsicovhb
`/opt/Code/metapsicov/bin/metapsicovhb test.colstats test.metapsicov.stage1 test.ss2 test.solv /opt/Code/metapsicov/data/weights_hbpass2.dat | sort -n -r -k 5 > test.metapsicov.hb`
18. Compile the outputs
19. Plot contact map
`tail -n +2 example/5ptpA.fasta > test.contactnewfasta`
then
`/opt/Code/metapsicov/bin/plotmap test.contactnewfasta test.metapsicov.stage1  test.metapsicov.stage1`
20. copy intermediary files to final
`cp $tempdir/$jobid.metapsicov.stage1  $htmldir/$jobid.metapsicov.stage1.txt`
`cp $tempdir/$jobid.metapsicov.stage2  $htmldir/$jobid.metapsicov.stage2.txt`
`cp $tempdir/$jobid.metapsicov.hb  $htmldir/$jobid.metapsicov.hb`
`cp $tempdir/$jobid.psicov  $htmldir/$jobid.psicov.txt`
