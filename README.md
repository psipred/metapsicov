# MetaPSICOV Protocol

Adjust the number of threads/cpus appropriately to your environment

1. Run blast+
`/scratch0/NOT_BACKED_UP/dbuchan/Applications/ncbi-blast-2.2.31+/bin/psiblast -evalue 0.001 -num_descriptions 2000 -num_alignments 0 -num_threads 4 -num_iterations 3 -inclusion_ethresh 0.001 -query example/test.fsa -db /scratch0/NOT_BACKED_UP/dbuchan/uniref/uniref90.fasta -out test.bls -out_pssm test.pssm`

2. Build matrix
`./bin/chkparse test.pssm > test.mtx`

3. Run HHBlits vs unprot20 db - HH-Suite soeding lab
`hhblits -i example/test.fsa -n 3 -e 0.001 -d /scratch0/NOT_BACKED_UP/dbuchan/hhblitsdb/uniprot20_2016_02/uniprot20_2016_02 -cpu 10 -oa3m test.a3m -diff inf -cov 80 -id 99 > test.hhlog`

4. Get all hhblits alignments lines(!?) and count

`grep -v '^>' test.a3m | sed 's/[a-z]//g' > test.hhbaln`
`cat test.hhbaln | wc -l`

5. run jack_hhblits to perform an hhblits like search using jackhmmer
http://hmmer.org/download.html

`./bin/jack_hhblits test /scratch0/NOT_BACKED_UP/dbuchan/Applications/hmmer-3.1b2 /scratch0/NOT_BACKED_UP/dbuchan/Applications/hh-suite ./ /cs/research/bioinf/home1/green/dbuchan/Code/metapsicov_web/consip/data/blast/uniref100 > test.jacklog`
