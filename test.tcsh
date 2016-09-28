#!/usr/bin/tcsh
set echo
set tmpdir = "/cs/research/bioinf/home1/green/dbuchan/Code/metapsicov"
set HHLIB = "/scratch0/NOT_BACKED_UP/dbuchan/Applications/hh-suite"
set jobid = test
python $HHLIB/scripts/hhsuitedb.py --cpu 4 -o $tmpdir/$jobid-mydb/mydb --ia3m "${tmpdir}/${jobid}-mya3m/*" --force
