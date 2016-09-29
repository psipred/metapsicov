# Prerequisites

1. HH-suite
2. HH databases pdb70 and uniprot20
3. Blast+
4. Uniref90 blastdb
5. hmmer and uniref100
6. Psipred (with solv pred)
7. CCMPRED
8.


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

`./bin/jack_hhblits test /scratch0/NOT_BACKED_UP/dbuchan/Applications/hmmer-3.1b2 /scratch0/NOT_BACKED_UP/dbuchan/Applications/hh-suite /cs/research/bioinf/home1/green/dbuchan/Code/metapsicov /cs/research/bioinf/home1/green/dbuchan/Code/consip3/data/blast/uniref100 example/test.fsa /cs/research/bioinf/home1/green/dbuchan/Code/metapsicov > test.jacklog`



6. Check which ever of hhbaln or jackaln has the most hits and copy that to the PSICOV input format

`cat test.hhbaln | wc -l`
`cat test.jackaln | wc -l`

`cp test.? test.aln`

7. Run PSIPRED + SOLVPRED components
`../psipred/bin/psipred test.mtx ../psipred/data/weights.dat ../psipred/data/weights.dat2 ../psipred/data/weights.dat3 > test.ss`

`../psipred/bin/psipass2 ../psipred/data/weights_p2.dat 1 1.0 1.0 test.ss2 test.ss`

`./bin/solvpred test.mtx data/weights_solv.dat > test.solv`

8. Run PSICOV components

`ulimit -s unlimited`
`./bin/alnstats test.aln test.colstats test.pairstats`
