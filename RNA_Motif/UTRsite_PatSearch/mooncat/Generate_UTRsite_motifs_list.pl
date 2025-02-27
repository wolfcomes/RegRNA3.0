#!/usr/bin/perl 

#usage: ./Generate_UTRsite_motifs_list.pl UTRsite_motifs_detail.txt > UTRsite_motifs_list.txt

###Åª¨úSeq file###
my $input_file=$ARGV[0]; #UTRsite_motifs_detail.txt



open(IN, $input_file) or die "$input_file: $!";
while (defined(my $line = <IN>)) 
{	
	if(index($line,"ID")==0) #ID:
	{
		my $ID="";
		# ID: U0001; date: 30-07-1997 
		chomp($line);
		my @split_line=split(' ',$line); # ID: U0001; date: 30-07-1997 
		my @split_line2=split(';',$split_line[1]); #U0001;
		$ID=$split_line2[0]; #U0001
		print "$ID|";
	}
	if(index($line,"Name")==0) #Name:
	{
		my $Name="";
		#Name:            Histone 3'UTR stem-loop structure (HSL3)
		chomp($line);
		my @split_line=split(' ',$line); #Name:            Histone 3'UTR stem-loop structure (HSL3)
		
		for(my $i=1;$i<@split_line;$i++)
		{			
			$Name=$Name.$split_line[$i]." ";
		}
		#print "$Name\n";
		print "$Name|";
	}
	if(index($line,"Taxon_Range")==0) #Taxon_Range:
	{
		my $Taxon_Range="";
		
		#Taxon_Range: 3HUM, 3'UTRs from human mRNAs
    #             3INV, 3'UTRs from invertebrate mRNAs
    #             3MAM, 3'UTRs from mammalian (not human/not rodent) mRNAs
    #             3MUS, 3'UTRs from mouse mRNAs
    #             3ROD, 3'UTRs from rodent (not mouse) mRNAs
    #             3VRT, 3'UTRs from vertebrate (not mammal) mRNAs	
    #Description:
		
		chomp($line);
		my @split_line=split(' ',$line); #Taxon_Range: 3HUM, 3'UTRs from human mRNAs
		my @split_line2=split(',',$split_line[1]); #3HUM,
		$Taxon_Range=$Taxon_Range.$split_line2[0];
				
		while (1) 
		{
			$line = <IN>;		
			if(index($line,"Description")==0)
			{
				last;
			}
			else
			{
				chomp($line);
				my @split_line3=split(' ',$line); #             3INV, 3'UTRs from invertebrate mRNAs
				my @split_line4=split(',',$split_line3[0]); #3INV,
				$Taxon_Range=$Taxon_Range.",".$split_line4[0];
			}
		}
		
		print "$Taxon_Range\n";
	}
	
}	
close(IN);
