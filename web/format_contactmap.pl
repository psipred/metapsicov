#!/usr/bin/perl
use strict;
use warnings;

#my $input_dir = addSlash($ARGV[0]);


open(INPUT,"<", "$ARGV[0]") or die "cant open Large_Database_file: $!";

#open(OUTPUT,">", $ARGV[1]) or die "cant open the output file $!"; #open output file
while (defined (my $line = <INPUT>))

{

next if ($line=~/^>/);
#print OUTPUT $1 if ($line =~/(\w+)/);
print  $1 if ($line =~/(\w+)/);
}

#close (OUTPUT);
close (INPUT);


