#!/usr/bin/perl

#reads erpin output files and filters solutions with 3Ns or more, or below
#a given score cutoff.

sub usage{print "erpin output | $0 [-fasta] [-nogap] [-c <cutoff>] [-e <Evalue>] [-t <Training set>]\n";exit;}

$cutoff=-99999999;
$evalue=-$cutoff;
$badsol=0;
$goodsol=0;

while ($_ = shift){
    /^-[\s]*[\?h]$/ && usage();
    /^-fasta$/ && ($fasta=1);
    /^-nogap$/ && ($nogap=1);
    /^-c$/ && ($_=shift) && ($cutoff=$_);
    /^-t$/ && ($_=shift) && ($trset=$_);
    /^-e$/ && ($_=shift) && ($evalue=$_);
}

usage() if eof();

if (defined $trset) {
    print STDERR "training set = $trset\n";
    open(TR, $trset) || die ("Can't open OK $trset\n");
    while ($line = <TR>) {
	if ($line !~ /^>/) {
	    chomp $line;
	    $line =~ s/[-\.\s]//gi;
	    $training .= $line." ";
	}
    }
   close TR;
}

if($cutoff != -99999999){print STDERR "Score cutoff = $cutoff\n";}
if($evalue != 99999999){print STDERR "E-value cutoff = $evalue\n";}

$line=<>;
while($line) {
    if ($line =~ /^>/) {
	$name=$line;
	chomp $name;
	$line=<>;
	readsol();
	while ($line && ($line!~ /^>/) && ($line !~ /^\s*$/)) {
	    readsol();
	}
    }else{
	$line=<>;
    }
}


print STDERR "$goodsol alignments retained ($badsol rejected).\n";

sub readsol {
    $line =~ /(\S+)\s+\S+\s+(\d+)\.\.(\d+)\s+(\S+)\s+(\S+)/;
    $sens=$1;
    $pos1=$2;
    $pos2=$3;
    $sco=$4;
    $eva=$5;

    $seq=<>;
    chomp $seq;
    $Ncnt = $seq =~ tr/N/N/;
    if (($Ncnt < 3) && ($sco > $cutoff)&&($eva<=$evalue)) {
	$goodsol++;
	if ($fasta) {
	    $seq =~ s/\.//g;
	}
	if($nogap){
	    $seq =~ s/\-//g;
	}
	print substr($name,0,100), " ", $sens, " ", $pos1, "-", $pos2, " ", $sco,," ",$eva, "\n";
	print $seq, "\n";
    }else{
	$badsol++;
    }
    $line=<>;
}



