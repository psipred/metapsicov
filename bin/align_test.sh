#!/bin/sh
value=`cat $1 | wc -l`
if [ $value -gt 5 ]
then
    echo "$1 has sufficient alignments"
    exit 0
else
    echo "$1 has too few alignments"
    exit 127
fi
