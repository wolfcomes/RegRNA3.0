#!/usr/bin/perl -w

use strict;

if ( (scalar @ARGV) < 3 ) {
    print "Incorrect number of arguments!\n";
    print "Usage: HMM_proc.pl fasta_file hmm_file e-value\n";
    exit 1;
}

# fasta file processing
my $input       = $ARGV[0];
# location of HMM model file
my $hmm_file    = $ARGV[1];
# e-value
my $e_value     = $ARGV[2];

# grep sequences and structures for HMMER
if ( !(-e "$input.myResult.grep") ) {
    `pcregrep -v "^# " $input.myResult > $input.myResult.grep`;
} else {
    print "$input.myResult.grep already EXIST!\n";
}

# run HMMER
if ( !(-e "$input.E$e_value.hmmResult") ) {
    `/usr/local/bin/hmmsearch -E $e_value $hmm_file $input.myResult.grep > $input.E$e_value.hmmResult`;
} else {
    print "$input.E$e_value.hmmResult already EXIST!\n";
}

print "Done! $input.E$e_value.hmmResult is ready!\n";
