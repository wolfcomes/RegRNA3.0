#!/usr/bin/env perl

my $now_time = $ARGV[1];
#print "id=$now_time\n";
use warnings;
use strict;

my $DEBUG = 0;

use Data::Dumper;

# get top number of result
my $getResultNo = 3;
# hash to store the top number results
my %seqHoA = ();
# e-value of each different structure sequences
my %e_value = ();

#my @in_queue_seqname;
#my $in_queue_seqno=1;

my $in_queue_seq_no=1;
my @in_queue_seqname;
my @in_queue_start_position;
my @in_queue_free_energy;
my @in_queue_fasta_locus;
my @in_queue_shape;

# help messages
if ( (scalar @ARGV) < 1 ) {
    print "\nIncorrect number of arguments!\n";
    print "Usage: rbprune.pl hmmResult_file\n\n";
    exit 1;
}
            

# input hmm file to process
my $infile = $ARGV[0];
# change name to be grep file
my $grepfile = $infile;
# store the riboswitch name
my $which_riboswitch = $infile;
$which_riboswitch =~ s/^.*input\.(.+)\.fas\..+\.hmmResult/$1/;
$grepfile =~ s/(.+)\.E\d+\.hmmResult/$1\.myResult/;




while ( <> ) {
    # Get result in HMMER file and parse it
    if ( /(.+=\d+)(s\d+):.+\s+(.+)$/) {
        if ( !exists($seqHoA{$1}) ) { @{$seqHoA{$1}} = (); }
        next if scalar @{ $seqHoA{$1} } >= $getResultNo;
#print "$1 = @{$seqHoA{$1}}\n";
        $e_value{$1.$2} = $3 if scalar @{ $seqHoA{$1} } <= $getResultNo-1;
#print "$1 + $2 = $e_value{$1.$2}\n";
        push @{ $seqHoA{$1} }, $2 if scalar @{ $seqHoA{$1} } <= $getResultNo-1;
#print "$1 = @{$seqHoA{$1}}\n";
    }
}

# print for debug mode
for my $d ( keys %seqHoA ) {
    print "$d: @{$seqHoA{$d}}\n" if $DEBUG != 0;
}

# get the original fasta name to trace back to grep file and get sequence and structure data
my $grep_data = "";
for my $d ( keys %seqHoA ) {
    for ( @{ $seqHoA{$d} } ) {
        $grep_data .= "$d$_\$|";
    }
}
$grep_data =~ s/\|$//;
exit if $grep_data eq "";
$grep_data = "\^\$" if $grep_data eq "";
#print "\n--------test0-------\n$grep_data\n";

my @split_grep_data=split('\|',$grep_data);
#my @split_grep_data=split(/\$\|/,$grep_data);

#my @split_grep_dataxxx=split(/\$\|/,$grep_data);
#print "$split_grep_dataxxx[0]----- $split_grep_dataxxx[1]---- $split_grep_dataxxx[2]\n";

my @result;

#print "==============================================\n";
open (RESULT, ">$infile.RESULT");
for(my $i=0;$i<@split_grep_data;$i++)
{
  #to avoid the problem of pcregrep searching while  '|' in fasta seq_name
  if(index($split_grep_data[$i],'$')==length($split_grep_data[$i])-1)
  {
	#print "\n---------test1-------\n$split_grep_data[$i]\n";
	my @temp_array=();
	@temp_array = `/usr/local/bin/pcregrep -A 3 \"$split_grep_data[$i]\" $grepfile`;

	#print "\n---------test2-------\n$result[$i]\n";
        
        push @result, @temp_array;	
  }
#}
	
	
	
	# count data for RNAeval
	my $shapeCount = 0;
	# data for RNAeval
	my $RNAeval = "";
	# result of RNAeval
	my $RNAevalResult = "";
	# fasta name of other data for printing result
	my $fasta_locus = "";
	# sequence and structure result for printing result
	my $shape = "";
	# e-value for printing result
	my $HMM_e_value = "";
	
	# result file
	#open (RESULT, ">$infile.RESULT");
	
	# process each line of results
	for ( @result ) 
	{
    # get fasta nam
    if ( /^>/ ) {
        chomp;
        s/^>//;
        $fasta_locus = $_;
    }
    # process sequence and structure
    if ( /#/ ) 
    {
        $shapeCount++;
        chomp;
        s/# //;
        $RNAeval = "$_\n".$RNAeval;
	
        # send data to RNAeval to evaluate if sequences and structure is ready
        if ( $shapeCount == 2 ) 
        {
            $shape = $RNAeval;
            $RNAeval =~ s/^\n//;
            $RNAevalResult = `echo -e \"$RNAeval\" | /usr/local/bin/RNAeval`;
            $RNAeval = "";
            $shapeCount = 0;
        }
    }
    # print result if RNAeval is not empty
    if ( $RNAevalResult ne "" ) 
    {
        # get free energy
        $RNAevalResult =~ s/.+\(\s*(\-*\d+\.\d+)\).+/$1/s;
        # get e-value
        $HMM_e_value = $e_value{$fasta_locus};
        # add free energy and e-value into fasta name line
        $fasta_locus =~ s/^(.+)$/>$1/;
        $fasta_locus =~ s/=(\d+)s\d+/|$which_riboswitch|Start:$1|/;
        $fasta_locus .= "Free energy: $RNAevalResult|HMM E-value:$HMM_e_value";
        # print results



        #print "\n\n---------New Round-----\n\n";


        my @temp_subname=split("Free energy: ",$fasta_locus);
        #print "$temp_subname[0] $temp_subname[1]\n";
        my @temp_free_energy=split('\|',$temp_subname[1]);
        #print "$temp_free_energy[0]\n";

        my @temp_subname2=split("Start:",$fasta_locus);
        #print "$temp_subname2[0] $temp_subname2[1]\n";
        my @temp_start_position=split('\|',$temp_subname2[1]);
        #print "$temp_start_position[0]\n";

        if($in_queue_seq_no==1)
        {
				  #print "------first results seq--------\n";
          $in_queue_seqname[$in_queue_seq_no]=$temp_subname2[0];
          $in_queue_start_position[$in_queue_seq_no]=$temp_start_position[0];
          $in_queue_free_energy[$in_queue_seq_no]=$temp_free_energy[0];

          $in_queue_fasta_locus[$in_queue_seq_no]=$fasta_locus;
          $in_queue_shape[$in_queue_seq_no]=$shape;
          
          $in_queue_seq_no++;
        }
        else
        {
        	 #print "####test1######\n";
           #my $temp_dis=abs($in_queue_start_position-$temp_start_position[0]);
           #print "########$temp_dis##########\n";
           #print "^^^^^^^^".index($temp_subname2[0],$in_queue_seqname)."^^^^^^^^^^\n";
           my $same_seq_flag=0;
				   for(my $j=1;$j<$in_queue_seq_no;$j++)
				   {
				   		my $temp_dis=abs($in_queue_start_position[$j]-$temp_start_position[0]);
				   	
					   	if(index($temp_subname2[0],$in_queue_seqname[$j])==0 && $temp_dis<=10) #same region result
				      {
				      	$same_seq_flag=1;
				        if( $temp_free_energy[0]<$in_queue_free_energy[$j]) #if has lower free energy, replace with new result in queue
				        {
				        	#print "==========test replace===============\n";
				        	$in_queue_seqname[$j]=$temp_subname2[0];
						      $in_queue_start_position[$j]=$temp_start_position[0];
						      $in_queue_free_energy[$j]=$temp_free_energy[0];
				       
						      $in_queue_fasta_locus[$j]=$fasta_locus;
						      $in_queue_shape[$j]=$shape;
						    }
						    else
						    {
						    	#print "******test abandon*******\n";
						    }
								$j=999999;
				      }
			     }
			     #else #different region, then print result seq in queue, and replace with new result in queue
			     if($same_seq_flag==0)
			     {
							#print "---------test new result-----------\n";
							#print "$in_queue_fasta_locus\n";
			      	#print "$in_queue_shape\n";
							#print RESULT "$in_queue_fasta_locus\n";
			        #print RESULT "$in_queue_shape\n";
			  
			        $in_queue_seqname[$in_queue_seq_no]=$temp_subname2[0];
			        $in_queue_start_position[$in_queue_seq_no]=$temp_start_position[0];
			        $in_queue_free_energy[$in_queue_seq_no]=$temp_free_energy[0];
			  
			        $in_queue_fasta_locus[$in_queue_seq_no]=$fasta_locus;
			        $in_queue_shape[$in_queue_seq_no]=$shape;
			  
			        $in_queue_seq_no++;
			     }          
			  }

        #if($in_queue_seqno==1 || index($fasta_locus, $in_queue_seqname[$in_queue_seqno-1])==-1)
        #{
        #  $in_queue_seqname[$in_queue_seqno]=$fasta_locus;
        #  print "$fasta_locus\n";
        #  print RESULT "$fasta_locus\n";
        #  print "$shape\n";
        #  print RESULT "$shape\n";
        #  $in_queue_seqno++;
        #}

				#顯示所有結果在網頁上
        #print "$fasta_locus\n";        
        #print "$shape\n";
        
        #print RESULT "$fasta_locus\n";
        #print RESULT "$shape\n";

        $RNAevalResult = "";
        $fasta_locus = "";
        $shape = "";
    }
  }
}

#print "---------test2-----------\n";
for(my $j=1;$j<$in_queue_seq_no;$j++)
{
	#print "($in_queue_seqname[$j] | $in_queue_start_position[$j] | $in_queue_free_energy[$j] )\n";
	
	#for package #直接輸出在螢幕上
	#print "$in_queue_fasta_locus[$j]\n";
	#print "$in_queue_shape[$j]\n";
	
	#for web version #需傳遞變數給figure.php (START)
	my $type;
	if(index($infile,"Purine")!=-1)	
	{
		$type="Purine";
	}
	if(index($infile,"PreQ1")!=-1)	
	{
		$type="PreQ1";
	}
	if(index($infile,"Lysine")!=-1)	
	{
		$type="Lysine";
	}
	if(index($infile,"Cobalamin")!=-1)	
	{
		$type="Cobalamin";
	}
	if(index($infile,"TPP")!=-1)	
	{
		$type="TPP";
	}
	if(index($infile,"Glycine")!=-1)	
	{
		$type="Glycine";
	}
	if(index($infile,"glmS")!=-1)	
	{
		$type="glmS";
	}
	if(index($infile,"FMN")!=-1)	
	{
		$type="FMN";
	}
	if(index($infile,"SAM")!=-1)	
	{
		$type="SAM";
	}
	if(index($infile,"SAM_alpha")!=-1)	
	{
		$type="SAM_alpha";
	}
	if(index($infile,"yybP-ykoY")!=-1)	
	{
		$type="yybP-ykoY";
	}
	if(index($infile,"ykkC-yxkD")!=-1)	
	{
		$type="ykkC-yxkD";
	}	
	#print "<a href=./figure.php?type=$type&id=$now_time&file=$infile.RESULT&Seqno=$j style=text-decoration:none target=_blank><font color=blue>";

	##print out in > $now_time.temp
	#print "<a href=./figure.php?type=$type&id=$now_time&Seqno=$j style=text-decoration:none target=_blank><font color=blue>\n";
	#print "$in_queue_fasta_locus[$j]\n";
	#print "$in_queue_shape[$j]\n";
	#print "</font></a>";
	#for web version #需傳遞變數給figure.php (END)	
	
	print RESULT "$in_queue_fasta_locus[$j]\n";
	print RESULT "$in_queue_shape[$j]\n";
}

close RESULT;
