#!/bin/tcsh -f
set jobid = $1
echo $jobid
foreach x ($jobid-mya3m/*.hhm)
    tail -n +2 $x > tmpfile
    mv tmpfile $x
    # set filename = `echo $x | sed -e 's/\.a3m//'`
    # echo $filename
    # /scratch0/NOT_BACKED_UP/dbuchan/Applications//hh-suite/build/src/hhmake -i $x -o $filename.hhm
    # $HHLIB/build/src/hhmake -i $x -o $jobid-mya3m/$x.hhm
end
