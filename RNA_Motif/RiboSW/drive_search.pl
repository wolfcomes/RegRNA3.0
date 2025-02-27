#!/usr/bin/perl -w

use strict;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Long;
use Cwd;

use constant BULGE      => 2;
use constant DEBUG      => 0;
use constant SHOW_RESULT=> 0;
use constant SHOW_SEQ   => 0;

# get current program directory
my $program_dir     = $0;
$program_dir        =~ s|(.+)/drive_search\.pl|$1|;
if (scalar @ARGV < 2) {
    print "\n";
    print "Usage: drive_search.pl [-a] fasta_file descriptor_file\n";
    print "\tArgument:\t-a: search the anti-sense sequence\n";
    print "\n";
}
(scalar @ARGV) >= 2 or exit 1;# "need more than 1 argument";

# for antisense search, check for "-a" parameter
my $antisense       = '';
GetOptions ('antisense|a' => \$antisense );
my $cwd             = getcwd;

my $inputFastaFile  = $ARGV[0];
my $descr_file      = $ARGV[1];

# check descriptor file exists
if ( !( -e $descr_file ) ) {
    print "NO DESCRIPTOR FILE PROVIDED!\n";
    exit 1;
}

open (DESCR, $descr_file) or die "$descr_file cannot be open! $!";

# get descriptor content and remove comments
my @descr_content = <DESCR>;
for (@descr_content) { chomp; if(/^#/) {$_=''} }
my $descr_content = join ('', @descr_content);


close (DESCR);

open (FASTA, $inputFastaFile) or die "$inputFastaFile cannot be open! $!";

# variables
my @raw_data    = <FASTA>;	# fasta file content
my @sequence    = ();		# array of multiple line sequence
my @AoSeqDivided = ();		# divide long sequence into shorter to process
my $fasta_name  = '';		# fasta name
my $seq_length  = 0;		# total sequences length
my $startPoint  = 0;		# predict start position
my $procSize    = 100000;	# length to divide long length into shorter

close (FASTA);

my $seq = '';
my $rev_seq = '';
my ($fh, $tmpfile);

# processing input sequence
for (@raw_data) {
    chomp;

    # if process line is a fasta name
    if ($_ =~ />/) {
        $seq = join ('', @sequence);
        ($fh, $tmpfile) = tempfile();
        open ($fh, '>', $tmpfile);
        $seq_length = length $seq;
        #  divide seq into reasonable length
        @AoSeqDivided = &divideSeq ($seq);
	$startPoint = 0;
        
	# fold the sense strand, each AoSeqDivided size is $procSize, divide long sequence to process
        for (@AoSeqDivided) {
            print $fh "$_\n";
            system ("$cwd/$program_dir/a.out", $fasta_name, $tmpfile, $descr_content, $startPoint) if ($seq ne '');
            $startPoint = $startPoint + $procSize;
        }
        
	# fold the antisense strand, same as sense strand process
        if ($antisense) {
            $startPoint = 0;
            $rev_seq = &reverse_seq($seq);
            @AoSeqDivided = &divideSeq ($rev_seq);
            
            for (@AoSeqDivided) {
                print $fh "$_\n";
                system ("$cwd/$program_dir/a.out", $fasta_name, $tmpfile, $descr_content, $startPoint, $seq_length) if ($seq ne '');
                $startPoint = $startPoint + $procSize;
            }
        }
        
        close ($fh);
        unlink $tmpfile;
        s/>//g;
        $fasta_name = $_;
        @sequence = ();
        $seq = '';
    }
    # if process line is sequence data, remove not nucleotide characters and translate all to lowercase
    elsif ($_ !~ />/) {
        print "$_\n" if SHOW_SEQ;
        s/[^ATUCGatucg]//g;
        tr/A-Z/a-z/;
        s/t/u/g;
        print "$_\n\n" if SHOW_SEQ;
        push @sequence, $_;
    }
}
## last line data will be processed here... same as other line process above
$seq = join ('', @sequence);
($fh, $tmpfile) = tempfile();
open ($fh, '>', $tmpfile);
$seq_length = length $seq;
#  divide seq into reasonable length
@AoSeqDivided = &divideSeq ($seq);
$startPoint = 0;

for (@AoSeqDivided) {
    print $fh "$_\n";

    system ("$cwd/$program_dir/a.out", $fasta_name, $tmpfile, $descr_content, $startPoint) if ($seq ne '');
    $startPoint = $startPoint + $procSize;
}
$startPoint = 0;

if ($antisense) {
    $startPoint = 0;
    $rev_seq = &reverse_seq($seq);
    @AoSeqDivided = &divideSeq ($rev_seq);
    
    for (@AoSeqDivided) {
        print $fh "$_\n";
        system ("$cwd/$program_dir/a.out", $fasta_name, $tmpfile, $descr_content, $startPoint, $seq_length) if ($seq ne '');
        $startPoint = $startPoint + $procSize;
    }
}

close ($fh);
unlink $tmpfile;


# transform the sense strand sequence into antisense strand sequence
sub reverse_seq {
    my $seq = shift(@_);
    my $rev_one = reverse $seq;
    $rev_one =~ tr/aucg/uagc/;

    return $rev_one;
}

# divide long sequence into shorter sequence to process, size defined in $procSize
sub divideSeq {
    my $seq = shift(@_);
    my @AoSeq = ();

    my $index = 0;
    my $seqLeng = length $seq;
    while ($index < $seqLeng) {
        push @AoSeq, (substr $seq, $index, $procSize + 1000);
        $index = $index + $procSize;
    }

    if ( ((scalar @AoSeq) > 1) && (length $AoSeq[-1]) < 1000) {
        pop @AoSeq;
    }

    return @AoSeq;
}
