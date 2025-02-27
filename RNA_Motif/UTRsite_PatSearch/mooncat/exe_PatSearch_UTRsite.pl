#!/usr/bin/perl 

#usage: ./exe_PatSearch_UTRsite.pl UTRsite_motifs_list.txt mooncat.input.fas mooncat.UTRsite.temp mooncat.UTRsite.result
#e.g. ./exe_PatSearch_UTRsite.pl UTRsite_motifs_list.txt input.fas temp.txt Result.txt

###讀取Seq file###
my $UTRsite_file=$ARGV[0]; #UTRsite_motifs_list.txt

my $input_file=$ARGV[1]; #mooncat.input.fas
my $temp_file=$ARGV[2]; #mooncat.UTRsite.temp
my $out_file=$ARGV[3]; #mooncat.UTRsite.result

my $first_motif_flag=0;

#system("echo \"\" > $out_file"); #先清乾淨$out_file(好讓一個個結果能夠附加>>進來)

open(IN, $UTRsite_file) or die "$UTRsite_file: $!";
while (defined(my $line = <IN>)) 
{		
	#U0001|Histone 3'UTR stem-loop structure (HSL3) |UTR
	chomp($line);
	my @split_line=split('\|',$line);	#UTRsite_ID
	my $UTRsite_ID=$split_line[0];
	my $UTRsite_Name=$split_line[1];
	my $UTRsite_Taxon_Range=$split_line[2];		
	
	if(!($UTRsite_ID eq "U0033")) #UTRsite U0033是在計算Open readin frame (ORF)這個fRNAfinder已經有算了，所以不要
	{	
		my $descr_file="$UTRsite_ID.descr";	#U0001.descr
		#print "$descr_file\n";
			
		#print "./PatSearch -inputfile $input_file -outputfile $temp_file -commandfile $descr_file -noreport -nooutcommand > /dev/null\n";
		system("./PatSearch -inputfile $input_file -outputfile $temp_file -commandfile $descr_file -noreport -nooutcommand  > /dev/null");
		
		open(IN2, $temp_file) or die "$temp_file: $!";
		my $line_no=0;
		while (defined(my $line2 = <IN2>)) 
		{
			$line_no++;
			
			if($line_no>=5) #第五行開始才會是結果 
			{
				#NM_000032 : [99,121] : guu c guccu cagugc agggc aac				
				chomp($line2);
				my @split_line2=split('\[',$line2);
				my @split_line3=split('\]',$split_line2[1]);
				my @split_line4=split('\,',$split_line3[0]);
				my $motif_start=$split_line4[0];
				my $motif_end=$split_line4[1];
				
				my @split_line5=split('\:',$line2);
				my $motif_seq=$split_line5[2];
				
				if($first_motif_flag==0)
				{
					system("echo \"$motif_start $motif_end | $UTRsite_ID | $UTRsite_Name| $UTRsite_Taxon_Range |$motif_seq\" > $out_file");
					$first_motif_flag=1;
				}
				else
				{
					system("echo \"$motif_start $motif_end | $UTRsite_ID | $UTRsite_Name| $UTRsite_Taxon_Range |$motif_seq\" >> $out_file");
				}
			}		
		}
		close(IN2);	
	}
}	
close(IN);