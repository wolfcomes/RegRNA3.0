#!/usr/bin/perl
my $input_file=$ARGV[0];
my $index_file=$ARGV[1];
 
my $seq_no=0;

my $index_no=0;
my @indexname;
my @indexseqno;
open(INDEX_FILE,"$index_file");	
while (defined(my $line2 = <INDEX_FILE>))
{
	$index_no++;
	chomp($line2);
	my @temp=split("!!!",$line2);
	$indexname[$index_no]=$temp[0];
	$indexseqno[$index_no]=$temp[1];
	
}
close(INDEX_FILE);

#for(my $i=1;$i<=$index_no;$i++)
#{
#	print "$indexname[$i] | $indexseqno[$i]\n";
#}


open(IN, $input_file) or die "$input_file: $!";

while (defined(my $line = <IN>))
{
  chomp($line);
  if((index($line,">")==0))
  {
    $seq_no++;
		
		recall_seqname($line);
		
  }
  else
  {
		print "$line\n";
  }
}
close(IN);



sub recall_seqname
{
	my $current_seqno=$_[0]; # e.g. >seq1|Purine|Start:186|Free energy: -13.40|HMM E-value:1.2e-05
	
	my @temp=split('\|',$current_seqno);
	my $real_seqno=$temp[0];
	my $real_seqno_len=length($real_seqno);
	
	for(my $i=1;$i<=$index_no;$i++)                   
	{		
		#print "----i=$i--\n";		
		#print "$real_seqno == $indexseqno[$i]\n";
		if(index($real_seqno,$indexseqno[$i])==0 && length($real_seqno)==length($indexseqno[$i]))
		{
			print "$indexname[$i]";
			my $seqname_suffix= substr $current_seqno, $real_seqno_len; length(current_seqno)-$real_seqno_len;
			
			#print "<font color=red><b>";
			print "$seqname_suffix\n";
			#print "</b></font>";
			
			
			$i=99999;
		}
	}      
}
