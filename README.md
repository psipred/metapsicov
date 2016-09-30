# Prerequisites

1. HH-suite
2. HH databases pdb70 and uniprot20
3. Blast+
4. Uniref90 blastdb
5. hmmer and uniref100
6. Psipred (with solv pred)
7. cblas library https://github.com/xianyi/OpenBLAS pass CFLAGS CPPFLAG AND LDFLAGS
7. CCMPRED git clone --recursive https://github.com/soedinglab/CCMpred
8. Freecontact wget ftp://rostlab.org/free/freecontact-1.0.21.tar.xz


# MetaPSICOV Protocol

Adjust the number of threads/cpus appropriately to your environment

##. Run the initial HHBLits to find areas with hits
`hhblits -i example/test.fsa -n 3 -e 0.001 -d /scratch0/NOT_BACKED_UP/dbuchan/hhblitsdb/pdb70 -cpu 4 -o test.hhr > test.pdbhhblog`

##. Running the MetaPSICOV components

1. Run blast+ and build matrix
`/scratch0/NOT_BACKED_UP/dbuchan/Applications/ncbi-blast-2.2.31+/bin/psiblast -evalue 0.001 -num_descriptions 2000 -num_alignments 0 -num_threads 4 -num_iterations 3 -inclusion_ethresh 0.001 -query example/test.fsa -db /scratch0/NOT_BACKED_UP/dbuchan/uniref/uniref90.fasta -out test.bls -out_pssm test.pssm`

`./bin/chkparse test.pssm > test.mtx`

2. Run HHBlits vs unprot20 db - HH-Suite soeding lab
`hhblits -i example/test.fsa -n 3 -e 0.001 -d /scratch0/NOT_BACKED_UP/dbuchan/hhblitsdb/uniprot20_2016_02/uniprot20_2016_02 -cpu 10 -oa3m test.a3m -diff inf -cov 50 -id 99 > test.hhlog`

3. Get all hhblits alignments lines(!?) and count
Generates the PSICOV alignment file format, just every aligned sequence one on top of the other.
`grep -v '^>' test.a3m | sed 's/[a-z]//g' > test.hhbaln`
`cat test.hhbaln | wc -l`

5. If you don't have at least 3000 alignments lines then
run jack_hhblits to perform an hhblits like search using jackhmmer
http://hmmer.org/download.html
Requires openmpi
Also requires uniref100

`./bin/jack_hhblits test /scratch0/NOT_BACKED_UP/dbuchan/Applications/hmmer-3.1b2  /cs/research/bioinf/home1/green/dbuchan/Code/metapsicov /cs/research/bioinf/home1/green/dbuchan/Code/consip3/data/blast/uniref100 example/test.fsa /cs/research/bioinf/home1/green/dbuchan/Code/metapsicov > test.jacklog`

6. Check which ever of hhbaln or jackaln has the most hits and copy that to the PSICOV input format
`cat test.hhbaln | wc -l`
`cat test.jackaln | wc -l`

`cp test.? test.aln`

7. Run PSIPRED + SOLVPRED components
`../psipred/bin/psipred test.mtx ../psipred/data/weights.dat ../psipred/data/weights.dat2 ../psipred/data/weights.dat3 > test.ss`

`../psipred/bin/psipass2 ../psipred/data/weights_p2.dat 1 1.0 1.0 test.ss2 test.ss`

`./bin/solvpred test.mtx data/weights_solv.dat > test.solv`

8. Run PSICOV pre-processor components

`ulimit -s unlimited`
`./bin/alnstats test.aln test.colstats test.pairstats`

9. Run Contact predictors
if `cat $tempdir/$prefix.aln | wc -l` is greater than 10 do these steps

`timeout 86400 bin/psicov -z 6 -o -d 0.03 test.aln > test.psicov 2>&1`
`timeout 86400 /cs/research/bioinf/home1/green/dbuchan/Code/CCMpred/bin/ccmpred -t 6 test.aln test.ccmpred > /dev/null 2>&1`
`/cs/research/bioinf/home1/green/dbuchan/Code/freecontact-1.0.21/src/freecontact -a 8 < test.aln > test.evfold`

10. Run metapsicov elements
`touch test.psicov test.evfold test.ccmpred`

`bin/metapsicov test.colstats test.pairstats test.psicov test.evfold test.ccmpred test.ss2 test.solv data/weights_6A.dat data/weights_65A.dat data/weights_7A.dat data/weights_75A.dat data/weights_8A.dat data/weights_85A.dat data/weights_9A.dat data/weights_10A.dat data/weights_811A.dat data/weights_1012A.dat > $test.metapsicov.stage1`
`bin/metapsicovp2 test.colstats test.metapsicov.stage1 test.ss2 test.solv data/weights_pass2.dat | sort -n -r -k 5 | head -5000 > test.metapsicov.stage2`
