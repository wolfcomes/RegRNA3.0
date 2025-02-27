#!/usr/bin/perl
my $input_file=$ARGV[0];
my $renamed_file=$ARGV[1];
my $index_file=$ARGV[2];

my $seq_no=0;

open(IN, $input_file) or die "$input_file: $!";
open(RENAMED_FILE, ">$renamed_file");
open(INDEX_FILE,">$index_file");

while (defined(my $line = <IN>))
{
  chomp($line);
  if((index($line,">")==0))
  {
        $seq_no++;
	print RENAMED_FILE ">seq$seq_no\n";
	print INDEX_FILE "$line!!!>seq$seq_no\n";
  }
  else
  {
	print RENAMED_FILE  "$line\n";
  }
}
close(IN);
close(RENAMED_FILE);
close(INDEX_FILE);
