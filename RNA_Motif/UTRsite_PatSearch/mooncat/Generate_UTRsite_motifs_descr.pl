#!/usr/bin/perl 

#usage: ./Generate_UTRsite_motifs_descr.pl UTRsite_motifs_detail.txt

###讀取Seq file###
my $input_file=$ARGV[0]; #UTRsite_motifs_detail.txt

my $ID="";

open(IN, $input_file) or die "$input_file: $!";
while (defined(my $line = <IN>)) 
{	
	if(index($line,"ID")==0) #ID:
	{		
		# ID: U0001; date: 30-07-1997 
		chomp($line);
		my @split_line=split(' ',$line); # ID: U0001; date: 30-07-1997 
		my @split_line2=split(';',$split_line[1]); #U0001;
		$ID=$split_line2[0]; #U0001		
	}
	
	if(index($line,"Pattern")==0) #Pattern:
	{		
		
		#Pattern:     r1={au,ua,gc,cg,gu,ug}
    #             0...1 mmmm p1=ggyyy u hhuh a r1~p1 mm 0...3

		open(OUT,">$ID.descr") or die "$ID.descr: $!";
		
		#處理Pattern的第一行
		chomp($line);
		my @split_line=split(' ',$line); #Pattern:     r1={au,ua,gc,cg,gu,ug}
		
		for(my $i=1;$i<@split_line;$i++)
		{			
			print OUT "$split_line[$i] ";
		}
		print OUT "\n";
		
		#處理Pattern的第2行&...
		while (1) 
		{
			$line = <IN>;		
			if(index($line,"Taxon_Range")==0)
			{
				last;
			}
			else
			{
				print OUT $line;
			}
		}
		
		close(OUT);
	}
}	
close(IN);
