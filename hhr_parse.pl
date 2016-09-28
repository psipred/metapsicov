#!/usr/bin/perl -w
#
# Work out the alignmenbt range
#

use strict;
use Data::Dumper;

my @recs;
my @range;

open(SEQFILE, $ARGV[0]) or die;
my $seq = do { local $/; <SEQFILE> };
$seq =~ s/>.*\n//;
$seq =~ s/[^A-za-z]//g;
close SEQFILE;

my $nres = length($seq);

my $masked = $seq;
#print Dumper $masked;

#open HHR results and find all the results with a prob>=99. Then get the HHM
#range and build a masked sequence that length of that range
open(HHRFILE, $ARGV[1]) or die;
while (<HHRFILE>)
{
    if (/^..[0-9] /)
    {
      #print $_;
	    @recs = split(' ', substr($_, 34));
      #print Dumper @recs;
	    if ($recs[0] >= 99.0)
	    {
	      @range = split('-', $recs[6]);

  	    # print $range[0]-1, "-", $range[1]-1, "\n";
  	    for (my $i = $range[0]-1; $i <= $range[1]-1; $i++)
  	    {
	  	    substr($masked, $i, 1) = ' ';
          #print $i." ".$masked."\n";
	      }
	    #print $range[0], " ", $range[1], "\n";
	    }
    }

    if (/^>/)
    {
	    last;
    }
}

while ($masked =~ /([A-Za-z]+)/g)
{
  my $roffset = pos($masked) - length($1);
  print($roffset);
  my $domseq = $1;
  print $domseq;
  
  if (length($domseq) >= 30)
  {
    open(DOMSOUT, ">$ARGV[2]/$ARGV[3].$roffset.fasta") or die;
    print DOMSOUT ">Domain$roffset\n$domseq\n";
    close(DOMSOUT);
  }
}
# print Dumper @range;
# print Dumper $masked;