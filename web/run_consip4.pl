#!/usr/bin/perl

# Script to execute CONSIP3 CASP prediction
# Original written by David Jones in 1998
# Current version written by David Jones in April 2016

use Digest::MD5 qw(md5_hex);
use BSD::Resource;
use IO::Handle;
use Text::Wrap;
use File::Basename;

setpriority 0, 0, 1;

setrlimit(RLIMIT_AS, RLIM_INFINITY, RLIM_INFINITY);
setrlimit(RLIMIT_DATA, RLIM_INFINITY, RLIM_INFINITY);
setrlimit(RLIMIT_STACK, RLIM_INFINITY, RLIM_INFINITY);

$maindir = dirname(__FILE__);
$bindir = "$maindir/bin";
$tempdir = "/var/www/tmp/consip";

# sub MailError {
#   my($msg,$email) = @_;
#   $subjectline = $subject . " CONSIP3 *ERROR*";
# #  open(MAILOUT, "| /usr/bin/Mail -s '$subjectline' $email") or die "Cannot open to send mail";
#   open(MAILOUT, "| /usr/bin/Mail -s '$subjectline' d.t.jones\@ucl.ac.uk") or die "Cannot open to send mail";
#   print MAILOUT $msg;
#   close(MAILOUT);
#   if (open(LOGF, ">>/var/www/tmp/consip.log")) {
#       flock(LOGF, 2) or die;
#       $date = localtime;
#       print LOGF $date," : ERROR $jobid\n";
#       close(LOGF);
#   }
# }

sub RunMetaPSICOV
{
    my($prefix) = @_;

    if (system("$bindir/blastpgp -a 4 -b 0 -v 2000 -j 3 -h 0.001 -e 0.001 -d $maindir/data/blast/unifilt -i $tempdir/$prefix.fasta -C $tempdir/$jobid.chk > $tempdir/$jobid.blast") != 0)
    {
	&MailError("CONSIP3 ERROR 01 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    open(MMOUT, ">$tempdir/$jobid.pn") or die;
    print MMOUT "$jobid.chk";
    close(MMOUT);
    open(MMOUT, ">$tempdir/$jobid.sn") or die;
    print MMOUT "$prefix.fasta";
    close(MMOUT);

    if (system("$bindir/makemat -P $tempdir/$jobid > /dev/null") != 0)
    {
	&MailError("CONSIP3 ERROR 02 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    if (system("$bindir/hhblits -i $tempdir/$prefix.fasta -n 3 -e 0.001 -d $maindir/data/hhsuite/uniprot20_2016_02 -cpu 4 -oa3m $tempdir/$prefix.a3m -diff inf -cov 50 -id 99 > $tempdir/$prefix.hhblog 2>&1") != 0)
    {
	&MailError("CONSIP3 ERROR 03 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    if (system("grep -v '^>' $tempdir/$prefix.a3m | sed 's/[a-z]//g' > $tempdir/$prefix.hhbaln") != 0)
    {
	&MailError("CONSIP3 ERROR 04 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    $naln_hhblits = `cat $tempdir/$prefix.hhbaln | wc -l`;

    if ($naln_hhblits < 3000)
    {
	# Try to find more sequences via jackhmmer/UNIREF100
	if (system("$bindir/jack_hhblits $prefix $bindir $tempdir $maindir/data/blast/uniref100 > $tempdir/$prefix.jacklog") != 0)
	{
	    &MailError("CONSIP3 ERROR 05 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	    exit 0;
	}

	$naln_jack = `cat $tempdir/$prefix.jackaln | wc -l`;
    }
    else
    {
	$naln_jack = 0;
    }

    if ($naln_jack > $naln_hhblits)
    {
	system("/bin/cp -f $tempdir/$prefix.jackaln $tempdir/$prefix.aln")
    }
    else
    {
	system("/bin/cp -f $tempdir/$prefix.hhbaln $tempdir/$prefix.aln")
    }

    if (system("$bindir/psipred $tempdir/$jobid.mtx $maindir/data/psipred/weights.dat $maindir/data/psipred/weights.dat2 $maindir/data/psipred/weights.dat3 > $tempdir/$jobid.ss") != 0)
    {
	&MailError("CONSIP3 ERROR 06 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    if (system("$bindir/psipass2 $maindir/data/psipred/weights_p2.dat 1 1.0 1.0 $tempdir/$prefix.ss2 $tempdir/$jobid.ss > /dev/null") != 0)
    {
	&MailError("CONSIP3 ERROR 07 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    if (system("$bindir/solvpred $tempdir/$jobid.mtx $maindir/data/weights_solv.dat > $tempdir/$prefix.solv") != 0)
    {
	&MailError("CONSIP3 ERROR 08 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    if (system("$bindir/alnstats $tempdir/$prefix.aln $tempdir/$prefix.colstats $tempdir/$prefix.pairstats > /dev/null") != 0)
    {
	&MailError("CONSIP3 ERROR 09 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    $naln = `cat $tempdir/$prefix.aln | wc -l`;

    if ($naln >= 10)
    {
	$psicovret = system("timeout 86400 $bindir/psicov -z 6 -o -d 0.03 $tempdir/$prefix.aln > $tempdir/$prefix.psicov 2>&1");

	if ($psicovret != 0 && $psicovret >> 8 != 124)
	{
	    &MailError("CONSIP3 ERROR 10 ($psicovret) - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	    exit 0;
	}

	$ccmret = system("timeout 86400 $bindir/ccmpred -t 6 $tempdir/$prefix.aln $tempdir/$prefix.ccmpred > /dev/null 2>&1");

	if ($ccmret != 0 && $ccmret >> 8 != 124)
	{
	    &MailError("CONSIP3 ERROR 11 ($ccmret) - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	    exit 0;
	}

	if (system("$bindir/freecontact -a 8 < $tempdir/$prefix.aln > $tempdir/$prefix.evfold") != 0)
	{
	    &MailError("CONSIP3 ERROR 12 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	    exit 0;
	}
    }

    if (system("touch $tempdir/$prefix.psicov $tempdir/$prefix.evfold $tempdir/$prefix.ccmpred") != 0)
    {
	&MailError("CONSIP3 ERROR 13 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    if (system("$bindir/metapsicov $tempdir/$prefix.colstats $tempdir/$prefix.pairstats $tempdir/$prefix.psicov $tempdir/$prefix.evfold $tempdir/$prefix.ccmpred $tempdir/$prefix.ss2 $tempdir/$prefix.solv $maindir/data/weights_6A.dat $maindir/data/weights_65A.dat $maindir/data/weights_7A.dat $maindir/data/weights_75A.dat $maindir/data/weights_8A.dat $maindir/data/weights_85A.dat $maindir/data/weights_9A.dat $maindir/data/weights_10A.dat $maindir/data/weights_811A.dat $maindir/data/weights_1012A.dat > $tempdir/$prefix.metapsicov.stage1") != 0)
    {
	&MailError("CONSIP3 ERROR 14 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }

    if (system("$bindir/metapsicovp2 $tempdir/$prefix.colstats $tempdir/$prefix.metapsicov.stage1 $tempdir/$prefix.ss2 $tempdir/$prefix.solv $maindir/data/weights_pass2.dat | sort -n -r -k 5 | head -5000 > $tempdir/$prefix.metapsicov.stage2") != 0)
    {
	&MailError("CONSIP3 ERROR 15 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
	exit 0;
    }
}


$ENV{'BLASTMAT'} = "$maindir/data/blast";
$ENV{'BLASTDB'} = "$maindir/data/blast";
$ENV{'HHLIB'} = "$maindir/data/hhsuite/lib/hh";

$jobid = $ARGV[0];
$email = $ARGV[1];
$subject = $ARGV[2];

if (open(LOGF, ">>/var/www/tmp/consip.log")) {
    flock(LOGF, 2) or die;
    $date = localtime;
    print LOGF $date," : RUNNING $jobid $email $subject\n";
    close(LOGF);
}

chdir '/var/www/cgi-bin/consip3';

# Check sequence still exists

if (! -e "$tempdir/$jobid.fasta") {
    &MailError("CONSIP3\n\nInternal queue error. Please resubmit request.\n", $email);
    exit 0;
}

open(SEQFILE, "$tempdir/$jobid.fasta") or die;
$seq = do { local $/; <SEQFILE> };
$seq =~ s/>.*\n//;
$seq =~ s/[^A-za-z]//g;
close SEQFILE;

$nres = length($seq);

$masked = $seq;

if (system("$bindir/hhblits -i $tempdir/$jobid.fasta -n 3 -e 0.001 -d $maindir/data/hhsuite/pdb70 -cpu 4 -o $tempdir/$jobid.hhr > $tempdir/$jobid.pdbhhblog 2>&1") != 0)

{
    &MailError("CONSIP3 ERROR PH1 - please report error to psipred\@cs.ucl.ac.uk\n", $email);
    exit 0;
}

open(HHRFILE, "$tempdir/$jobid.hhr") or die;

while (<HHRFILE>)
{
    if (/^..[0-9] /)
    {
	@recs = split(' ', substr($_, 34));

	if ($recs[0] >= 98.0)
	{
	    @range = split('-', $recs[6]);

	    #print $range[0]-1, "-", $range[1]-1, "\n";
	    for ($i = $range[0]-1; $i <= $range[1]-1; $i++)
	    {
		substr($masked, $i, 1) = ' ';
	    }
	    #print $range[0], " ", $range[1], "\n";
	}
    }

    if (/^>/)
    {
	last;
    }
}

# First run complete sequence
&RunMetaPSICOV($jobid);

open(S2FILE, "<$tempdir/$jobid.metapsicov.stage2") or die;
while (<S2FILE>)
{
    @recs = split(' ');
    $contacts[$recs[0]][$recs[1]] = $recs[4];
}
close S2FILE;

# Now run detected domains (if any) that could not be matched to PDB templates
while ($masked =~ /([A-Za-z]+)/g)
{
    $domlen = length($1);
    $roffset = pos($masked) - $domlen;
    $domseq = $1;

    if (length($domseq) >= 30 && length($domseq) < $nres)
    {
	open(DOMSOUT, ">$tempdir/$jobid.$roffset.fasta") or die;
	print DOMSOUT ">Domain$roffset\n$domseq\n";
	close(DOMSOUT);

	&RunMetaPSICOV("$jobid.$roffset");

	# Erase any previously predicted contacts within this domain
	for ($i=$roffset+1; $i<=$roffset+$domlen; $i++)
	{
	    for ($j=$i+1; $j<=$roffset+$domlen; $j++)
	    {
		$contacts[$i][$j] = 0.0;
	    }
	}

	# Copy in newly predicted contacts for this domain
	open(S2FILE, "<$tempdir/$jobid.$roffset.metapsicov.stage2") or die;
	while (<S2FILE>)
	{
	    @recs = split(' ');
	    $contacts[$recs[0] + $roffset][$recs[1] + $roffset] = $recs[4];
	}
	close S2FILE;
    }
}

open(S3FILE, ">$tempdir/$jobid.metapsicov.stage3") or die;
for ($i=1; $i<=$nres; $i++)
{
    for ($j=$i+5; $j<=$nres; $j++)
    {
	$prob = $contacts[$i][$j];
	if ($prob > 0.0)
	{
	    print S3FILE "$i $j 0 8 $prob\n";
	}
    }
}
close S3FILE;

$message = "PFRMAT RR\nTARGET " . $subject . "\nAUTHOR MetaPSICOV\nMODEL 1\n";
$subjectline = $subject . " MetaPSICOV Results";
open(MAILOUT, "| /usr/bin/Mail -s '$subjectline' $email") or die;
print MAILOUT $message;

$Text::Wrap::columns = 50;
open(SEQFILE, "$tempdir/$jobid.fasta");
<SEQFILE>;
print MAILOUT Text::Wrap::fill( '', '', join '', <SEQFILE>);
close SEQFILE;
print MAILOUT "\n";

$results = `cat $tempdir/$jobid.metapsicov.stage3 | sort -n -r -k 5 | head -5000`;
print MAILOUT $results, "END\n";

close(MAILOUT);

if (open(LOGF, ">>/var/www/tmp/consip.log")) {
    flock(LOGF, 2) or die;
    $date = localtime;
    print LOGF $date," : COMPLETED $jobid\n";
    close(LOGF);
}

exit 0;
