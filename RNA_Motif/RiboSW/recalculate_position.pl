#!/usr/bin/perl 

###讀取Seq file###
my $input_file=$ARGV[0];

open(IN, $input_file) or die "$input_file: $!";
while (defined(my $line = <IN>)) 
{	
	if(index($line,">")==0) #碰到>開頭，處理seqname
	{
		# >X83878|B.subtilis###168###|Purine|Start:18|Free energy: -15.00|HMM E-value:5.5e-06
		
		
		chomp($line);
		my @split_line=split('###',$line);
		my $cmhmm_result_startp=$split_line[1];	#取出168
		
		my @split_line2=split('Start:',$split_line[2]);
		my @split_line3=split('\|',$split_line2[1]);		
		my $RiboSW_startp=$split_line3[0];
		
		#print "$cmhmm_result_startp | $RiboSW_startp\n";
		my $true_startp=$cmhmm_result_startp + $RiboSW_startp;
		
		print $split_line[0].$split_line2[0]."Start:".$true_startp."|".$split_line3[1]."|".$split_line3[2]."\n";
		
	}	
	else	
	{
		print "$line";
	}
}

close(IN);
