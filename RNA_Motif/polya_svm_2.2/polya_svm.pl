#! /usr/bin/perl -w
# Created by Yiming Cheng on 2006/02/22
# Function: Predict the poly(A) site of the input file containing FASTA format sequences.
# 	    The output format depends on the option [-m mode]
#	    If the [-l location] option is given, it will only do one hit model
# Comments: 2006/05/04: modified by Yiming for simplification.
# 	    2006/05/22: last modifed by Yiming for integration.
#           2006/07/19: add the loading of train.model.range, need to give the train.model file.
#     	    		add one more option for the selection of min score of the cis-element or -15. 
#           2006/08/19: change the default value to inf, hpr:32 and cutoff value 6.
#			add the leave one out option
#	    2006/08/20: change the output of cutoff value and add some code for the case that if the sequence is less than 232.
#			use 0.8781261 which is average value of the default setup: 2^(-6/32).
#Changes    2006/11/28:	added load_range function to load range (mean/sd)
#by			added load_settings function to load settings.inf which contains various default parameters such as paths
#Michael Tsai		print functions for consistant 3 decimal formatted output 
#			added parameter for a common in/out-path which will be concated to various in/out-files
#	    2006/11/29	added new mode 'b' which will generate e and s outputs and concat into a single out file
#	    2006/12/04	added range functionality taking parameters x and y from command line.
#	buggy		new functionality limits the size of the score matrix by not calculating values more than
#			100nt around x and y (denoted as $x1 $x2 in code)
#			added positional information to mode s printout and added spaces every 10 nt
#	    2006/12/07	added paslookup function and yespas variable, command line parameter is -p yes or no, no by default
#			added usage of it in scanseq
#			added function negseq(seq) which reverses/translates a sequence to it's negative form.
#				note that range is applied AFTER the sequence is reversed
#			added to matchmat(), command line parameter is -n yes, no by default
#			positional information for negative sequences = Total length - position
#	    2006/12/20	fixed a bug with match mode and ranges - output was shifted
#			fixed a bug with match mode and ranges - neg sequences now displays correctly			
#			added location now uses range for optimizing speed
#	    2006/12/26	changed min-score appended to 5' of mat-matrix to -6.28
#			changed sub match_mat so appended min-score is not visited by hpr function (appended min scores are used for
#			calculations of latter scores but not saved to score mat file. 
# 	    2006/12/27	added -g option, for usage in gen_traindata.pl scripts, only generates score mat file
#			changed is_pos_hpr: now looks for the max probability of a whole region where all probabilities are >cutoff
#				ie. if an hpr size is 32 and a region of 64 is found, rather than 2 predictions, a single prediction is made
#				on the maximum value of the whole 64nt region.	
#	    2007/01/04	changed neg sequence option. Range is applied after the sequence is reversed, with respect to the reversed sequence. 
#			thus if seq is 0-400, reversed with range 100-200, it is applied to region 200-300 opposite of the original sequence
#			and 100-200 of the reversed sequence 
#	    2007/01/09	changed is_pos_hpr to use logs for product of vec thereby increasing the max limit that cutoff can be before overflow.
#			changed is_pos_m to use logs for product (same as above)
#			changed when range<hpr dynamically calc avg instead of using hard coded .878 value
#	    2007/01/17  fixed a misplaced bracket for -g option that interfered with mode s.
#	    2007/01/25	added function to translate all U to T
#	    2007/07/13	fixed a bug regarding non-ATCG sequences
#			added filtering criteria based on support vector data in human models
#				scaled mat must have at least 2 positive scores for AUE..ADE && the average must be >-1
#	    2007/07/18  fixed a bug in read_num_file, caused problems for single predictions (3 element array ended up as a 2 element array)

use strict;
use Getopt::Std;
my %args;
getopt("iocmlthsdrxypngI",\%args);
my $switch=1; #for s mode printout, 0 for single +|:.-, 1 for extending the symbols to the length of the cis element
my ($input,$output,$cutoff,$mode,$location,$train_model,$hpr,$min_score,$element_out);
##load settings variables
my @track_id=();
my ($x1,$x2,$rangeoffset);		#12/04/06 variables for range function
my $inname="settings.inf";		#settings file
my $FOS=60;				#number of min-scores to add to the Front Of Sequence
my $yespas=0;				#create a file containing pas hexamer locations
my $yesneg=0;				#generate a negative sequence
my %hlen=();				#hash length
my %phash;				#pas hash
my @array;
my @break;
my $modelpath="";
my $inpath="";
my $outpath="";
my $svmpath="";
my $out="junk";
my $scutoff=6;
my $in;
my $smode="e";
my $slocation=0;
my $range_file;
my $scale_ave_score_limit=-1;
my $scale_pos_count_limit=2;
&load_settings($inname);
&ini_options();
my $p_cutoff = 2**(-$cutoff);
##VARIABLES for cis elements
my @file_ext = ("AUE1","AUE2","AUE3","AUE4","CUE1","CUE2","CDE1","CDE2","CDE3","CDE4","ADE1","ADE2","ADE3","ADE4","ADE5");
my $num_ext=@file_ext;
my $op=0;
my @reg_start = (0,0,0,0,60-$op,60-$op,100,100,100,100,140-$op,140-$op,140-$op,140-$op,140-$op);
my @reg_end = (60+$op,60+$op,60+$op,60+$op,100,100,140+$op,140+$op,140+$op,140+$op,200,200,200,200,200);

my $num_file = $input.".num";  				#output file "<input>.num" 
my $score_file = $input.".score";  			#output file "<input>.score" 
my $output_mat_score = $input.".mat.score"; 		# "<input>.mat.score" pwd  	the "score mat file"
my $score_2_svm = $input.".svm.score";  		# "<input>.svm.score" 		
my $predict_junk = $input.".pre.junk";  		# "<input>.pre.junk":

open(NUM_FILE,">$num_file");				#open write only file
open(SCORE,">$score_file");				#open write only file

&update_element_out();					
my (%ele_mat,%ele_mat_size,%score_mat,%ele_inf);	#hash tables
&load_element();					#load element file ->pssm (located at bottom of page

#if (defined($args{t})){
if (defined($range_file)){
&load_range();						#load the range file, training model file sent to svm predict
}
my $inf_file = "element.inf";				#human element files
my $mat_file = "element.mat";				#currently not used, they're hard coded to __DATA__
#&load_element_inf();   #initiate %ele_inf
#&load_element_mat();   #initiate %ele_mat,%ele_mat_size

if ($yespas>0){open(PAS,">$output.pas");}
&scan_seq();	       #initiate %score_mat
		       #if mode's', generate $input.score
		       #if mode'e', generate $input.mat.score

if ($yespas>0){close(PAS);}
close(NUM_FILE);
close(SCORE);

if (defined($args{g})){if ($args{g} =~/y/) {exit;}}	#for gen_traindata.pl script, exits prematurely to save on comp time. gen_traindata.pl only needs mat score file.

if ($mode eq 's'){
	system("mv $score_file $output");
	system("rm $num_file -f");
	system("rm $output_mat_score -f");
	exit;
}


my @posmap=@{&scale($output_mat_score,$score_2_svm)};		#normalize sequence scores
##keep in mind ./svm-predict not WIN compatible in this instance
system("svm-predict -b 1 $score_2_svm $train_model $predict_junk");
#system($svmpath."svm-predict -b 1 $score_2_svm $train_model $predict_junk");
my @pro_data;						#used in is_pos_hpr
my $num_line = 0;
my $pos=0;
&read_input($predict_junk);  				#initiate @pro_data
&read_num_file(); 					#generate the temporary output file and do the optimization method
open(OUT,">$output");
&print_output();
close(OUT);

if ($mode eq 'c'){
	&mode_c();
}
if ($mode eq 'b'){
system("cat $output $score_file>$output.b");
system("mv $output.b $output");
}
system("rm $score_file -f");
system("rm $num_file -f");
system("rm $output_mat_score -f");
system("rm $score_2_svm -f");
system("rm $predict_junk -f");

sub ini_options{
	$input = $args{i};				#input file
	if (!defined($input)){				#make sure there's an input
		&exit_help();
	}
	$input = $inpath.$args{i};			#path + input file
	if (!-e $input){print "no such input file: $input\n";exit;}	#make sure file even exists
	$output = $args{o};
	if (!$output){					#make sure there's an output name, otherwise make one.
		$output = $outpath."junk.out";
	}
	else{$output=$outpath.$args{o};}
	$yespas=$args{p};
	if (defined($yespas)){				#yes i want a pas location file
		if ($yespas=~/^y/){$yespas=1;}
		else{$yespas=0;}
	}else{$yespas=0;}
	$yesneg=$args{n};
	if (defined($yesneg)){				#yes i want to polyasvm the negative sequence
		if ($yesneg=~/^y/){$yesneg=1;}
		else{$yesneg=0;}
	}else{$yesneg=0;}
	$cutoff = $args{c};				#polya SVM cutoff value
	if (!defined($cutoff)){
		$cutoff = 6;				#default value is optimal at 6
		if ($scutoff>0){$cutoff=$scutoff;}
	} elsif ($cutoff <0){
		&exit_help();
	}
	$mode = $args{m};
	if (!defined($mode)){
		#$mode = 'e';
		$mode=$smode;
	}
	$x1=$args{x};					#the start range
	$x2=$args{y};					#the end range
	$location = $args{l};
	if (!defined($location)){
		$location=$slocation;
		if ((!defined($location))||!($location>=0)){$location=0;} #make sure a location exists or 0
	}
	if ($location!=0){$x1=$location;$x2=$location;$FOS=0;}	#range speed optimization, only calc in the area, and there's no need for FOS
	if (defined($x1)&&defined($x2)){		#range sanity check
		if ($x1<100){print "Range start position needs to be greater than 100\n";exit;}
		if ($x1>$x2){print "Range start needs to be less than Range end\n";exit;}
		if ((!($x2-$x1)>0)&&($location==0)){print "Range end needs to be bigger than Range start\n";exit;}
		$rangeoffset=$x1;
	}
	if (!defined($train_model)){
	$train_model = $args{t};
	}
	if (!defined($train_model)){
		$train_model = $modelpath."human.scaled.model";	#get a train model
	}else{$train_model=$modelpath.$train_model;}
	if (!-e $train_model){print "PolyA_SVM can't find any training models!\n";exit;}
	if (defined($args{r})){$range_file=$args{r};}
	if (!defined($range_file)){
		$range_file="$train_model.range";	#get a range of the model
	}else{$range_file=$modelpath.$range_file;}
	if (!-e $train_model&&!defined($args{g})){print "Can't find train model: $train_model\n";exit;}
	if (!-e $range_file&&!defined($args{g})){print "Can't find range file: $range_file\n";exit;}
	$hpr=$args{h};
	if (!defined($hpr)){
		$hpr = 32;				#default HPR size is optimal at 32
	}
	$min_score = $args{s};
	if (!defined($min_score)){
		$min_score = "inf";
	}
	$element_out = $args{d};
	if (!defined($element_out)){
		$element_out = 0;
	}
}

sub exit_help{
	print "\nUsage: perl polya_svm.pl -i input [-o output] [-c cutoff] [-m mode] [-t train_model] [-l location] [-h hpr] [-s min_score] [-d element_out] \n";
	print "Options:\n";
	print "	-i input --- file containing the FASTA format sequences. \n";
	print "	-o output --- the output file, the output format is based on the mode. The default output is junk.out\n";
	print "	-c cutoff value --- the default evalue is 6, the lower the evalue, the more accuracy. Need to be positive \n";
	print "	-m mode --- if 'e', prediction mode, which returns prediction results. This is the default mode.\n";
	print "		    if 's', matching element mode, which returns matching results for 15 cis-elements.\n";
	print "	-t train_model --- the train model file. Always with the train_model.range file.\n";
	print "	-l location --- if this is given, it will do one hit model, just predict if the given location is a poly(A) site or not.\n";
	print "	-h hpr --- the HPR size. The default size is 32. \n";
	print "	-s min_score --- the selection of the min score for Inf. Default is inf, it can also be m15,min,zero.\n";
	print "	-d element_out --- the cis-element out, 1 means AUE1, 2 means AUE2 ... 15 means ADE4.\n";
	print "	-n yes/no --- generate negative sequence then polya_svm it. **the range is applied after the reversal\n";	
	print "	-x <num> --- start of the range. When -x and -y are specified, the range option is initiated\n";
	print "	-y <num --- end of the range\n";
	print "\n";
	exit;
}

sub update_element_out{
	if ($element_out ne '0'){
		my (@temp_file_ext,@temp_reg_start,@temp_reg_end);
		my @out_arr = split(/ /,$element_out);
		my $temp_index=0;
		for (my $j=0;$j<$num_ext;$j++){
			if (($temp_index <= $#out_arr)&&($out_arr[$temp_index] == ($j+1))){
				$temp_index++;
				next;
			}
			push @temp_file_ext,$file_ext[$j];
			push @temp_reg_start,$reg_start[$j];
			push @temp_reg_end,$reg_end[$j];
		}
		@file_ext=();
		@file_ext = @temp_file_ext;
		@reg_start=();
		@reg_start = @temp_reg_start;
		@reg_end = ();
		@reg_end = @temp_reg_end;
		$num_ext = @file_ext;
	}
}

sub load_element{		#reads __DATA__ at bottom of file, populates %ele_mat, %ele_inf, @nts
	my ($cur_element,@nts);
        while (<DATA>){
        	chomp;
		if ($_ =~ /^[AC]/){		##traps lines AUE1 AUE2 .. CUE1..ADE5
							#AUE1	8.07	5.36	2.69	1389	1.031025	4.69755043500543
                	my @line =split(/\t/,$_);
                	$ele_inf{$line[0]}{'25'}=$line[3];
                	$ele_inf{$line[0]}{'50'}=$line[2];
                	$ele_inf{$line[0]}{'75'}=$line[1];
			$ele_inf{$line[0]}{'mean'} = $line[5];
			$ele_inf{$line[0]}{'sd'} = $line[6];
		} elsif ($_ =~ /([\w\d]+).MAT/){	#traps line #AUE1.MAT
			$cur_element = $1;		#cur_element=AUE1 or whatever was matched
			@nts =();
		} elsif ($_=~/pos/i){			#traps line #pos	A	C	G	U
                	$_=~s/U/T/g;    			#replace U with T
                        @nts = split /\t/, $_;  # @nts = "pos A C G T"
                        shift @nts;       # @nts = "A C G T"
                } elsif ($_=~/^\d/){		#traps digits
                        my @vals=split /\t/, $_;	#split on tabs; 4	-2.00	Inf	2.22	Inf
                        my $pos=shift @vals;		#STORE index number 1 2 3 ...
                        for (my $i=0; $i<@vals;$i++){	
                        	next if($vals[$i] =~/Inf/i);	#skip block below if "Inf"
                                $ele_mat{$cur_element}{$pos}{$nts[$i]}=$vals[$i];	#populates ele_mat without Infs
						# $ele_mat{cis element}{index in elements file}{ACGT} ie. $ele_mat{AUE1}{1}{A}
                        }
                        $ele_mat_size{$cur_element}++;
                } else {
		}
        }
#	foreach my $ele (sort keys %ele_inf){
#		print "\n$ele\n";
#		foreach my $i (sort keys %{$ele_inf{$ele}}){
#			print "$i\[$ele_inf{$ele}{$i}]";
#		}
#	}

#print element matrix that was loaded into memory	
#	foreach my $ele (sort keys %ele_mat){
#		print "$ele\t";
#		foreach my $ind (sort keys %{$ele_mat{$ele}}){
#			print "$ind ";
#			foreach my $nt (sort keys %{$ele_mat{$ele}{$ind}}){
#				print "$nt\[$ele_mat{$ele}{$ind}{$nt}]";
#			}print "\n";
#		}
#	}
	
}

#DEFUNCT
sub load_element_mat{

	open(MAT,"<$mat_file")||die "Can't find the file $mat_file\n";
	my ($cur_element,@nts);
        while (<MAT>){
        	chomp;
		if ($_ =~ /([\w\d]+).MAT/){
			$cur_element = $1;
			@nts =();
		} elsif ($_=~/Pos/i){
                	$_=~s/U/T/g;    
                        @nts = split /\s+/, $_;  # @nts = "pos A C G T"
                        shift @nts;       # @nts = "A C G T"
                } elsif ($_=~/^\d/){
                        my @vals=split /\s+/, $_;
                        my $pos=shift @vals;
                        for (my $i=0; $i<@vals;$i++){
                        	next if($vals[$i] =~/Inf/i);
                                $ele_mat{$cur_element}{$pos}{$nts[$i]}=$vals[$i];
                        }
                        $ele_mat_size{$cur_element}++;
                } else {
			die "Never reach in element.mat";
		}
        }
	close(MAT);
}

##DEFUNCT
sub load_element_inf{

        open(INF,$inf_file)||die "Can't find the file $inf_file.\n";
        my @temp_arr;
        while(<INF>){
                chomp;
                if ($_  =~ /^#/){
                        @temp_arr = split(/\t/,$_);
                        shift @temp_arr;
                } elsif ($_ =~ /^\w/) {
                        my @line =split(/\s+/,$_);
                        #for (my $i=0;$i<=$#temp_arr;$i++){
                        #     $ele_inf{$line[0]}{$temp_arr[$i]}=$line[$i+1];
			$ele_inf{$line[0]}{'25'}=$line[3];
                	$ele_inf{$line[0]}{'50'}=$line[2];
                	$ele_inf{$line[0]}{'75'}=$line[1];
			$ele_inf{$line[0]}{'mean'} = $line[5];
			$ele_inf{$line[0]}{'sd'} = $line[6];
                        #}
                }
        }
        close(INF);
#	if (-e "$train_model.range"){
#		open(RANGE,"<$train_model.range");
#		my $k=0;
#		while(<RANGE>){
#			chomp;
#			my @line = split(/\t/,$_);
#			$ele_inf{$file_ext[$k]}{'mean'} = $line[1];
#			$ele_inf{$file_ext[$k]}{'sd'} = $line[2];
#			$k++;
#		}
#		close(RANGE);
#	} else {
#		print "$train_model.range\n";
#		die "can't find the range file";
#	}

}

sub load_range{
        if (-e $range_file){
                open(RANGE,"<$range_file");
		my $k=0;
		while(<RANGE>){
			chomp;
			my @line = split(/\t/,$_);
			$ele_inf{$file_ext[$k]}{'mean'} = $line[1]; ##assume cis element order is preserved
			$ele_inf{$file_ext[$k]}{'sd'} = $line[2];
			$k++;
		}
		close(RANGE);
	} else {
                print "$range_file\n";
		die "can't find the range file $train_model.range\n";
	}
}

sub scan_seq{
        open (IN,"<$input") or die "Input file not found";	#read from the sequence input file
        open (MATOUT,">$output_mat_score") or die "couldn't open $output_mat_score\n";			#write out scores
push @track_id,[0,0];
	my $cnt_seq=0;  
        my ($id,$seq,$anno);
        while (<IN>){
                chomp;
                if ($_=~/^>/){
			#$_=~s/\t/ /g;
                        $cnt_seq++;
                        if ($cnt_seq >1){
			if ($seq !~ /[ATCG]/){print "Seq $id doesn't contain ATCG\n";}
			else{
push(@track_id,[$track_id[-1][0]+length($seq),$id]);
                                &match_mat($id,$seq);
			}
                        }
                        ($id,$anno)=split /\t/, $_,2;  #for example, id = ">p.2987.1"
                                                     #anno = "LLID:2987|CONTIG:NT_004559|ST:1|LOC:4512837|CLEVC:6|ACCC:103"
                        if ($id=~/^>lcl\|([\w\.]+)/){    # lcl
                                $id=$1;
                        }elsif($id=~/^>([\w\. ]+)/){    
                                $id=$1;
                        }else{
                                die("problem");
                        }
			$id=substr($_,1);
                        $seq='';
                }else{
			$_=~s/\W//g;
                        $seq .= $_;
                }
        }
	if ($seq !~ /[ATCG]/){print "Seq $id doesn't contain ATCG\n";}else{
push(@track_id,[$track_id[-1][0]+length($seq),$id]);
        &match_mat($id,$seq);
	}
        close(IN);
  	close(MATOUT);
}

#############
sub match_mat{
	my $matsize;
        my $temp_id=$_[0];      #id ="p.2987.1";
        my $seq=$_[1];          #seq ="ATGCATGCCCCCCCCCAGGGTTC..."
        $seq =~ s/\s//g;   #delete empty char in the sequence
        my @seq_nts;
	$hlen{$temp_id}=length($seq);  
	if ($yesneg>0){
		$seq=negseq($seq);
	}
	if ($yespas>0){
		%phash=&paslookup($seq);
		&print_pas($temp_id);
		%phash=();
	}
        $seq=uc($seq);			#make case insensitive
	$seq=~s/U/T/g;
        @seq_nts=split //, $seq;	#load sequence as array of characters
        my $seq_len= @seq_nts;	#length of sequence
        foreach my $mat_name (@file_ext){	#@file_ext = AUE1 AUE2 etc.., initializes array
                @{$score_mat{$mat_name}}=();	
        }
	if (defined($x1)&&defined($x2)){
		if ($x1>=length($seq)){print "Start Range for $temp_id falls outside the specified range\n";exit;}
		if ((length($seq)-$x2)<100){print "End Range for $temp_id falls outside the specified range\n";exit;}
	}
	&getscore($seq);		#populate score matrix for each seq
		
	if ($mode eq 's'||$mode eq 'b'){
	############## Start printing score;
	my %temp_score_print=();
	for (my $i =0;$i<$num_ext;$i++){			#number of cis elements
		my $temp_name = $file_ext[$i];		#copy AUE1..ADE5
		$matsize=$ele_mat_size{$temp_name};
		my $prev_emission="-";
		my $prev_value=-1;
		my @temp_arr_score = @{$score_mat{$temp_name}}; #copy score_mat[AUE][]..[]
		my @new_arr_score;
		for (my $j=0;$j<=$#temp_arr_score;$j++){
			if (defined($rangeoffset)){next if ($j<$FOS);}
				if ($temp_arr_score[$j]>$ele_inf{$temp_name}{'75'}){
					$prev_emission="+";
					$matsize=$ele_mat_size{$temp_name};
					$prev_value=75;
					push @new_arr_score,"+";
						
				} elsif ($temp_arr_score[$j]>$ele_inf{$temp_name}{'50'}){
					$matsize=$ele_mat_size{$temp_name};
					if ($switch>0){
						if ($prev_value<50){
							push @new_arr_score,"|";
							$prev_emission="|";
							$matsize=$ele_mat_size{$temp_name};
							$prev_value=50;
						}else{push @new_arr_score,$prev_emission;}
					}else{push @new_arr_score,"|";}
				} elsif ($temp_arr_score[$j]>$ele_inf{$temp_name}{'25'}){
					$matsize=$ele_mat_size{$temp_name};
					if ($switch>0){
						if ($prev_value<25){
							push @new_arr_score,":";
							$prev_emission=":";
							$matsize=$ele_mat_size{$temp_name};
							$prev_value=25;
						}else{push @new_arr_score,$prev_emission;}
					}else{push @new_arr_score,":";}
				} elsif ($temp_arr_score[$j]>0){
					$matsize=$ele_mat_size{$temp_name};
					if ($switch>0){
						if ($prev_value<0){
							push @new_arr_score,".";
							$prev_emission=".";
							$matsize=$ele_mat_size{$temp_name};
							$prev_value=0;
						}else{push @new_arr_score,$prev_emission;}
					}else{push @new_arr_score,".";}
				} else {
					if ($switch>0){					
						if ($matsize>0){push @new_arr_score, $prev_emission;}
						else{
							push @new_arr_score,"-";
							$prev_emission="-";
							$prev_value=-1;
						}
					}else{push @new_arr_score,"-";}
				}
				$matsize--;
		}#end for
		for (my $j=($#new_arr_score+1);$j<$seq_len;$j++){
			if ($matsize>0){
				push @new_arr_score, $prev_emission;$matsize--;
			}else{                	
				push @new_arr_score,"-";
			}
		}
		push @{$temp_score_print{$temp_name}},@new_arr_score;
	}#end outter for
	my ($start,$end);
	my $line_len = 80;
	my $position=1;
	my $posend=0;
	#if ($yesneg>0){
	#	$position=$hlen{$temp_id};
	#	$posend=$position+1;
	#}
	my $tempstring; 
	my @temparray;
	my $num = int($seq_len/$line_len);
	my $rnum;
	if (defined($rangeoffset)){
		$rnum=int(($x2-$x1+201)/$line_len);
	#if ($yesneg==0){$position+=$x1-100;$posend+=$x1-100;}
		$position+=$x1-100;$posend+=$x1-100;
	}
	print SCORE ">$temp_id\n";
	for (my $i=0; $i<($num+1); $i++){
#new addition
		if (defined($rnum)){
			last if ($i>=$rnum);
		}
		$start=$line_len*$i;
		if ($start+$line_len>$seq_len){
			$end=($seq_len-1);
		}else{
			$end=$start+$line_len-1; #CHANGED to len-1 instead of just len 12/05/06
		}
		last if ($start>=$end);
		#print SCORE "\n    \t",(join "",@seq_nts[$start..$end]);
		if (defined($rangeoffset)){
			#if ($yesneg==0){$tempstring=join "",@seq_nts[$start+$rangeoffset-100..$end+$rangeoffset-100];}
			#else{
				$tempstring=join "",@seq_nts[$start..$end];
			#}
		}else
		{$tempstring=join "",@seq_nts[$start..$end];}
		#if($yesneg>0){
		#	$posend-=length($tempstring);
		#}else{
			$posend+=length($tempstring);
		#}
		@temparray=split /(..........)/, $tempstring;
		print SCORE "\npos:$position\t@temparray ",($posend);
			for (my $j=0;$j<$num_ext;$j++){
				#print SCORE "\n",$file_ext[$j],"\t",(join "",@{$temp_score_print{$file_ext[$j]}}[$start..$end]);
				$tempstring=join "",@{$temp_score_print{$file_ext[$j]}}[$start..$end];
				@temparray=split /(..........)/, $tempstring;
				print SCORE "\n$file_ext[$j]\t@temparray ",($posend);
			}
		#if ($yesneg>0){$position-=$line_len;}else{
		$position+=$line_len;
		#}
	}#end for
	print SCORE "\n\n";
	################### End printing score
	}#end if mode s 
	if ($mode eq 'e'||$mode eq 'b'){	#default mode e
		my $index = 1; 		 #for the score matrix, 
if (defined($args{I})){$index=$args{I};}
		################# Start generate matrix for svm
		my %mat_max_score=();				#init empty hash
		my ($start_site,$end_site);			#new vars start/end
		if ($seq_len<120){				#seq too small 
			print NUM_FILE $temp_id,"\t",0,"\n";
			return(0);
		} elsif ($seq_len<=200){
			print MATOUT $index,"\t";
			my ($start_site,$end_site);
			for (my $i=0;$i<$num_ext;$i++){	#for i=0< #cis elements
				my @one_score = @{$score_mat{$file_ext[$i]}};
				my $mid_num = int(scalar(@one_score)/2);
				if ($file_ext[$i]=~/AUE/){
					$start_site = 0;
					$end_site = 60-$ele_mat_size{$file_ext[$i]};
				} elsif ($file_ext[$i]=~/CUE/){
					$start_site = $mid_num-40;
					$end_site = $mid_num-$ele_mat_size{$file_ext[$i]};
				} elsif ($file_ext[$i]=~/CDE/){
					$start_site = $mid_num;
					$end_site = $mid_num+40-$ele_mat_size{$file_ext[$i]};
				} elsif ($file_ext[$i]=~/ADE/){
					$start_site = $seq_len-60;
					$end_site = $seq_len-$ele_mat_size{$file_ext[$i]};
				} else {
					die "Not reachable\n";
				}
				print MATOUT &get_max(join(":",@one_score[$start_site..$end_site])),"\t";
			}
			print MATOUT "\n";
			print NUM_FILE $temp_id,"\t",1,"\t",0,"\n";
			return(0);
		}
		### sequence length will be greater than 200 following
		## location==0, location is anywhere
		if ($location == 0){
		my $stopcondition=$seq_len-199;
		if (defined($rangeoffset)){$stopcondition+=$FOS;}
			for (my $i=0;$i<$num_ext;$i++){ 			#for # of AUE1..ADE5
				my @one_score = @{$score_mat{$file_ext[$i]}}; #load stripe of entire AUE1..
				$start_site = $reg_start[$i];
				#size of cis block-length of cis logos (the # of positions)
				$end_site = ($reg_end[$i]-$ele_mat_size{$file_ext[$i]}); #Determine Boundries of AUE1..
				my @first_arr = @one_score[$start_site..$end_site];	#load cis sub-regions AUE block etc..
				my $curr_max = &get_max(join(":",@first_arr));	#gets max score in the block
				push @{$mat_max_score{0}},$curr_max;		#push the max score
				my $first_ele;
$end_site++;
				#for (my $j=1;$j<($seq_len-199);$j++) 
				#account for cis Blocks 1..60,61..100,101..140,141..200
				for (my $j=1;$j<$stopcondition;$j++){ #account for cis Blocks 1..60,61..100,101..140,141..200
					if (defined($rangeoffset)){if ($j>=($x2-$x1+$FOS+1)){last;}}
					$first_ele = shift @first_arr; 		#dump first element
					push @first_arr,$one_score[$end_site];	#fill gap from the shift earlier
					$end_site++;					#first ele= [0], [start..end]=[1..61]
					$start_site++;					#move cis sub-region block forward 1 nt.
					#if first_ele>=max, push a new max of the entire block onto mat_max_score
					#otherwise update curr_max if end_site is bigger, but always push 
					#current max onto mat_max_score. 198 times
next if (!defined($one_score[$end_site]));
					if ($curr_max<=$first_ele){			#populate mat_max_score
						#get a new max for a new block (when max of first cis block is found->overlap blocks)
						$curr_max = &get_max(join(":",@one_score[$start_site..$end_site]));
						push @{$mat_max_score{$j}},$curr_max;
					} else {
						 if ($one_score[$end_site]>$curr_max){
							$curr_max = $one_score[$end_site];
						}					
							push @{$mat_max_score{$j}},$curr_max;
					}
				}#end for j			
			}#end for i
			for (my $i=0;$i<$stopcondition;$i++){
				if (defined($rangeoffset)){if ($i>=($x2-$x1+$FOS+1)){last;}}
				next if ($i<$FOS&&defined($rangeoffset));
				print MATOUT $index,"\t",join("\t",@{$mat_max_score{$i}}),"\n";
			}
			if (!defined($rangeoffset)){
				print NUM_FILE $temp_id,"\t",($seq_len-199),"\n";
			}else{
				print NUM_FILE $temp_id,"\t",($x2-$x1+1),"\n";
			}
		} else { #else location !=0
			if (($location <100)||($location >($seq_len-100))){
				die "The location you specified is beyond the range of the sequence $temp_id.\n";
			}
			my $i=$location-100;
			if (defined($rangeoffset)){$i=$location-$rangeoffset+$FOS;}
			print MATOUT "$index\t";
			for (my $j=0;$j<$num_ext;$j++){
				my @one_score = @{$score_mat{$file_ext[$j]}};
				$start_site = $i+$reg_start[$j];
				$end_site = $i+$reg_end[$j]-$ele_mat_size{$file_ext[$j]}; 
				print MATOUT &get_max(join(":",@one_score[$start_site..$end_site])),"\t";
			} 
			print MATOUT "\n";
			print NUM_FILE $temp_id,"\t",1,"\t",$i+100,"\n";
		}
	}#end if mode e
} #end match_mat

sub get_max{
	my @temp_score = split(/:/,$_[0]);
	my $max_score = -99;
	for (my $i=0;$i<=$#temp_score;$i++){
		if ($temp_score[$i]>$max_score){
			$max_score = $temp_score[$i];
		}
	}
	if ($max_score == -99){
		die "Never reach in get_max!.\n";
	}
	return($max_score);
}

        sub get_max2{ #parameters array reference, start site, end site
                my $max_score = -99;
		my $max_pos=0;
                for (my $i=$_[1];$i<=$_[2];$i++){
                        if (${$_[0]}[$i]>$max_score){
                                $max_score = ${$_[0]}[$i];
				$max_pos=$i;
                        }
                }
                if ($max_score == -99){
                        die "Never reach in get_max!.\n";
                }
                #return($max_score);
                return($max_score,$max_pos);
        }



sub getscore{ #Scores the sequence
        my $sequence = $_[0];
	my $ms=-6.28;
	my @edge=();
	if (defined($rangeoffset)){	#if there is a range, redefine range to include 100nt zone around the range to predict.
		for (my $i=0;$i<$FOS;$i++){push @edge,$ms;}
		if ((length($sequence)-$x2)<100){
			#this condition should never be met
			print "range end position needs to be more than 100nt from the end of the sequence, will predict as much as possible.\n";
			$x2=length($sequence)-100;
			$sequence=substr($sequence,($x1-99));
		}else{
			$sequence=substr($sequence, ($x1-100),($x2-$x1+201)); #area from x1-100 to x2+100 inclusive
		}
	}
        my @sequence_nts = split //,$sequence;	#character array for sequence, also constructed in &match_mat->maybe can optimize here
        my $seq_length = scalar(@sequence_nts);	#get length of the array
        my $score = 0;					#init score
#ie. for AUE1, for every nt of seq, get total score with PSSM, push into score_mat
	foreach my $mat_name (keys %ele_mat_size){	#for each AUE, ADE etc.. there're 15 of them
                for (my $j=0;$j<=($seq_length-$ele_mat_size{$mat_name});$j++){ #length - #cis element positions in say AUE1
		#basically iterate through sequence
                        $score=0;
			if (substr($sequence,$j,$ele_mat_size{$mat_name}) !~ /[ATCG]/){
				#print "$j $mat_name\n";
				$score=-98;
			}else{
                        for (my $k=$j; $k<($j+$ele_mat_size{$mat_name}); $k++){ #loops through #cis ele hit times
				#scoring procedure, generate score matrix
				#if no hit...
				#k-j+1=index in ele_mat
                                if (!exists($ele_mat{$mat_name}{($k-$j+1)}{$sequence_nts[$k]})){
                                        if ($min_score eq "m15"){
						$score = -15;
						last;
					} elsif ($min_score eq "min"){
						$score = $ele_inf{$mat_name}{'min'};
						last;
					} elsif ($min_score eq "zero"){
						$score = 0;
						last;
					} elsif ($min_score eq "inf"){
						$score += -6.28;
					}
                                } else {		# ele_mat[AUE1..][k-j+1][ACGT] nucleotide IS in TABLE
	                                $score += $ele_mat{$mat_name}{($k-$j+1)}{$sequence_nts[$k]};
                        	}
				}#end for
			}
			if ($min_score eq "zero"){
				if ($score <0){
					$score = 0;
				}
			}	#assign output to 3 decimal places and push into array
                        push @{$score_mat{$mat_name}},sprintf("%.3f",$score);
                }#end for length
		if (defined($rangeoffset)){
			my @temp=();
			push @temp,@edge;
			push @temp, @{$score_mat{$mat_name}};
			$score_mat{$mat_name}=[@temp];
		}
	}#end for cis elements
}

sub scale{
	my $scale_in =$_[0];	#$input.mat.score
	my $scale_out = $_[1];	#$input.svm.score
my (@pmap,@tempout,$runave,$countpos);
	open(SIN,"<$scale_in");
        open(SOUT,">$scale_out");
	my $current_pos=0;
        while(<SIN>){
                chomp;  
                my @items = split(/\t/,$_);
                my $num_items = @items; 
		@tempout=();
		$runave=0;
		$countpos=0;
                push @tempout, "$items[0] ";
                for (my $i=1;$i<$num_items;$i++){ #for AUE1..ADE5
                        my $item_i = ($items[$i]-$ele_inf{$file_ext[$i-1]}{"mean"})/$ele_inf{$file_ext[$i-1]}{"sd"};
                        push @tempout, "$i:",sprintf("%.3f",$item_i)," ";
			$runave+=$item_i;
			if ($item_i>0){$countpos++;}
                }       
		$runave/=$num_items;
		if (($runave<$scale_ave_score_limit)&&($countpos<$scale_pos_count_limit)){
			#my $id;
			#for (my $i=1;$i<@track_id;$i++){if ($.>$track_id[$i-1][0] && $.<$track_id[$i][0]){$id=$track_id[$i][1];last;}}
			#print "$id $. ",join('',@tempout),"\n";
			push @pmap,-1;
			next;
		}else{push @pmap,$current_pos++;}
		print SOUT join("",@tempout),"\n";
        }       
        close(SIN);
        close(SOUT);
return \@pmap;
}

sub mode_c{
	### Complex output
	my (%all_loc,@all_pre,@all_score);
	open(COUT,"<$output");
	my $temp_id;
	while(<COUT>){
		chomp;
		if ($_ =~/^>([\d\w\.]+)/){
			$temp_id = $1;
			$all_loc{$temp_id}='';
		} elsif ($_ =~/^\+/){
			my @line = split(/\t/,$_);
			if ($line[1]=~/M/){
				$all_loc{$temp_id} = 1;
			} else {
				my $junk_junk = $all_loc{$temp_id};
				if ($junk_junk eq ''){
					$all_loc{$temp_id} =($line[3]-100);
				} else {
					$all_loc{$temp_id} = $junk_junk."+".($line[3]-100);
				}
			}
		} elsif ($_ =~/^\-/){
			$all_loc{$temp_id}=0;
		}
	}
	close(COUT);
	open(TEMP_OUT,"<$predict_junk");
	while(<TEMP_OUT>){
		chomp;
		next if ($.==1);
		my @line = split(/ /,$_);
		push @all_pre,sprintf("%.4f",$line[1]);
	}
	close(TEMP_OUT);
	open(RES,"<$output_mat_score");
	while(<RES>){
		chomp;
		my @line =split(/\t/,$_);
		push @all_score,join("\t",@line[1..$#line]);
	}
	close(RES);
	open(NUM,"<$num_file")||die "Can't open num file.\n";
	open(COM,">complex.junk");
	print COM "#+/-\tp-pro\tLocation\tAUE1\tAUE2\tAUE3\tAUE4\tCUE1\tCUE2\tCDE1\tCDE2\tCDE3\tCDE4\tADE1\tADE2\tADE3\tADE4\tADE5\n";
	my $j=0;
	while(<NUM>){
		chomp;
		my @line = split(/\t/,$_);
		print COM ">$line[0]\n";
		if ($line[1]==0){
			print COM "Sequence too short. No prediction\n";
			next;
		} elsif ($line[1] == 1){
			if ($all_pre[$j]>0.5){
				if ($line[2] ==0){
					print COM "+++++","\t",$all_pre[$j],"\tMiddle\t",$all_score[$j],"\n";
				} else {
					print COM "+++++","\t",$all_pre[$j],"\t$line[2]\t",$all_score[$j],"\n";
				}
			} else {
				if ($line[2] ==0){
					print COM "-","\t",$all_pre[$j],"\tMiddle\t",$all_score[$j],"\n";
				} else {
					print COM "-","\t",$all_pre[$j],"\t$line[2]\t",$all_score[$j],"\n";
				}
			}
		} else {
			my $temp_loc = $all_loc{$line[0]};
			my @temp_loc_arr;
			my $k= -1;
			if (($temp_loc ne 0)&&($temp_loc ne 1)){
				@temp_loc_arr = split(/\+/,$temp_loc);
				$k=0;
			}
			if ($k<0){
				for (my $i=0;$i<$line[1];$i++){
					if ($all_pre[$i+$j]>=0.5){
						print COM "+";
					} else {
						print COM "-";
					}
					print COM "\t",$all_pre[$i+$j],"\t",$i+100,"\t",$all_score[$i+$j],"\n";
				}
			} else {
				for (my $i=0;$i<$line[1];$i++){
					if (($k<=$#temp_loc_arr)&&($i==$temp_loc_arr[$k])){
						print COM "+++++\t",$all_pre[$i+$j],"\t",$i+100,"\t",$all_score[$i+$j],"\n";
						$k++;
					} else {
						if ($all_pre[$i+$j]>=0.5){
							print COM "+","\t",$all_pre[$i+$j],"\t",$i+100,"\t",$all_score[$i+$j],"\n";
						} else {
							print COM "-","\t",$all_pre[$i+$j],"\t",$i+100,"\t",$all_score[$i+$j],"\n";
						}
					}
				}
			}
		}	
		$j = $j+ $line[1];
	}
	close(NUM);
	close(COM);
	system("mv complex.junk $output");
}


sub read_input{
	my $r_input = $_[0];
	open(PIN,"<$r_input");
	my $pos_col;
	my $bit_pos=0;
	while(<PIN>){
		chomp;
		my @items = split(/ /,$_);
		if ($. == 1){
			if ($items[1] == 1){
				$pos_col = 1;
			} elsif ($items[1] == -1){
				$pos_col = 2;
			} else {
				die "Data has some problem.\n";
			}
			next;
		}
		if ($bit_pos<$.-2){$bit_pos=$.-2;}	#-2 (1st line = 'labels 1 -1', $. is 1 based, array are 0 based.)
		while($posmap[$bit_pos++]<0&&$bit_pos<@posmap){
#print "$.\n";
			push @pro_data,0;
		}
		push @pro_data,$items[$pos_col];
	}
	if (@pro_data<$num_ext){for (my $i=@pro_data;$i<=$num_ext;$i++){push @pro_data,0;}}
	close(PIN);
}

sub read_num_file{
	open(NUM,"<$num_file")||die "Can't open the number file.\n";
	my $cur_num = 0;
	my @temp_return;
	my $output_junk = $output.".cym";
	open(OUT_JUNK,">$output_junk");
	while(<NUM>){
		chomp;
		$num_line++;
		my @temp = split(/\t([\w]+)$/,$_);  # id, num_hit,[location]		#MT FIX 07/18/07, makes sure there're 2-3 items in num_tmp
		my @num_tmp = split(/\t([\w]+)$/,$temp[0]);  # id, num_hit,[location]	#MT FIX 07/18/07
		push @num_tmp,$temp[1];							#MT FIX 07/18/07
		undef(@temp);								#MT FIX 07/18/07
		#my @num_tmp = split(/\t/,$_);  # id, num_hit,[location]
		print OUT_JUNK ">",$num_tmp[0],"\n";
		if ($num_tmp[1] ==0){
			print OUT_JUNK "Seq too short.\n";
			next;
		} elsif ($num_tmp[1] ==1){ #location probably specified
			if ($pro_data[$cur_num]>=0.5){
				if ($num_tmp[2] != 0){
					#print OUT_JUNK ($num_tmp[2]-100),":",($num_tmp[2]-100),":",($num_tmp[2]-100),":",$pro_data[$cur_num],"\n";
					print OUT_JUNK ($num_tmp[2]-100),":",($num_tmp[2]-100),":",($num_tmp[2]-100),":",-log($pro_data[$cur_num])/log(2),"\n";
				} else {
					#print OUT_JUNK "M:M:M:",$pro_data[$cur_num],"\n";
					print OUT_JUNK "M:M:M:",-log($pro_data[$cur_num])/log(2),"\n";
				}
				$pos++;
			} else {
				if ($location !=0){
					print OUT_JUNK "N\t",$pro_data[$cur_num],"\n";
				} else {
					print OUT_JUNK "N\n";
				}
			}
		} else {
			if ((defined($hpr))&&($hpr >10)){
				@temp_return = &is_pos_hpr($cur_num,$num_tmp[1]);
			} else {
				@temp_return = &is_pos_m($cur_num,$num_tmp[1]);
			}
			if ($temp_return[$#temp_return]>0){
				$pos++;
				print OUT_JUNK join("\t",@temp_return[0 .. ($#temp_return-1)]),"\n";
			} else {
				print OUT_JUNK "N\n";
			}
		}
		$cur_num = $cur_num + $num_tmp[1];
	}
	close(NUM);
	close(OUT_JUNK);
}

sub print_output{
	my $counter=0;
        my $tempid=undef;
        my $shift=undef;
        my $output_junk = $output.".cym";
        open(OUT_JUNK,"<$output_junk")||die "Can not open file.cym\n";
        while(<OUT_JUNK>){
                chomp;
                if ($_ =~ /^>/){
                        print OUT $_,"\n";
                        $tempid=substr($_,1);
                        print OUT "#+/-\tHPR_Fr\tHPR_To\tSite\tPolyA_SVM Score\n";
                } elsif ($_ =~ /^S/){
                        print OUT "Seq too short, need to be greater than 120 nts","\n\n";
                } elsif ($_ =~ /^N/){   # Negative
                        my @junk_arr = split(/\t/,$_);
                        if ($location !=0) {
                                print OUT "-\t",$location,"\t",$location,"\t",$location,"\t",sprintf("%.3f",12*(1-$junk_arr[1])),"\n\n";
                        } else {
                                print OUT "-\n\n";
                        }
                } elsif ($_ =~ /^M/){
                        my @line_junk = split(/:/,$_);
                        #print OUT "+\tM\tM\tM\t",sprintf("%.3f",(-log($line_junk[3])/log(2))),"\n";
                        print OUT "+\tM\tM\tM\t",sprintf("%.3f",$line_junk[3]),"\n";
                } elsif ($_ =~ /^\d+/){
                        my @line_junk = split(/\t/,$_);
                        for (my $i=0;$i<($#line_junk+1);$i++){
                                if ($line_junk[$i] =~ /^(\d+):(\d+):(\d+):([\d\.]+)/){
					$counter++;
                                        if (!defined($rangeoffset)){
                                                #print OUT "+ \t",($1+100),"\t",($2+100),"\t",($3+100),"\t\t",sprintf("%.3f",(-log($4)/log(2))),"\n";
                                                print OUT "+ \t",($1+100),"\t",($2+100),"\t",($3+100),"\t\t",sprintf("%.3f",$4),"\n";
                                        }
                                        else{
                                                #print OUT "+ \t",($1+$rangeoffset),"\t",($2+$rangeoffset),"\t",($3+$rangeoffset),"\t\t",sprintf("%.3f",(-log($4)/log(2))),"\n";
                                                print OUT "+ \t",($1+$rangeoffset),"\t",($2+$rangeoffset),"\t",($3+$rangeoffset),"\t\t",sprintf("%.3f",$4),"\n";
                                        }
				}
                        }
                        print OUT "\n";
                } else {
                        die "Never reach in print_output.\n";
                }
        }

        print OUT '#########';
        print OUT "\n\#Summary\#\n";
        print OUT '#########';
        print OUT "\n";
        print OUT "#Number of input sequences: $num_line.\n";
        print OUT "#Number of positive sequences: $pos.\n";
        print OUT "#Number of positive hits: $counter\n";
        #print OUT "#The PolyA_SVM score cutoff value is $cutoff.\n\n";
        print OUT "\n";
        #print OUT "#Total number of sequences: $num_line.\n";
        #print OUT "#Positive sequences: $pos.\n";
        #print OUT "#The polyasvm score cutoff value is $cutoff.\n\n";

        close(OUT_JUNK);
        system("rm $output_junk -f");
}


sub print_output2{
	my $tempid=undef;
	my $shift=undef;
	my $output_junk = $output.".cym";
	open(OUT_JUNK,"<$output_junk")||die "Can not open file.cym\n";
	while(<OUT_JUNK>){
		chomp;
		if ($_ =~ /^>/){
			print OUT $_,"\n";
			$tempid=substr($_,1);
			print OUT "#+/-\tFrom\tTo\tMid-Location\tPolyaSVM cutoff score\n";
		} elsif ($_ =~ /^S/){
			print OUT "Seq too short, need to be greater than 120 nts","\n\n";
		} elsif ($_ =~ /^N/){   # Negative
			my @junk_arr = split(/\t/,$_);
			if ($location !=0) {
				print OUT "-\t",$location,"\t",$location,"\t",$location,"\t",sprintf("%.3f",12*(1-$junk_arr[1])),"\n\n";
			} else {
				print OUT "-\n\n";
			}
		} elsif ($_ =~ /^M/){
			my @line_junk = split(/:/,$_);
			#print OUT "+\tM\tM\tM\t",sprintf("%.3f",(-log($line_junk[3])/log(2))),"\n";
			print OUT "+\tM\tM\tM\t",sprintf("%.3f",$line_junk[3]),"\n";
		} elsif ($_ =~ /^\d+/){
			my @line_junk = split(/\t/,$_);
			for (my $i=0;$i<($#line_junk+1);$i++){
				if ($line_junk[$i] =~ /^(\d+):(\d+):(\d+):([\d\.]+)/){
					if (!defined($rangeoffset)){
						#print OUT "+ \t",($1+100),"\t",($2+100),"\t",($3+100),"\t\t",sprintf("%.3f",(-log($4)/log(2))),"\n";
						print OUT "+ \t",($1+100),"\t",($2+100),"\t",($3+100),"\t\t",sprintf("%.3f",$4),"\n";
					}
					else{
						#print OUT "+ \t",($1+$rangeoffset),"\t",($2+$rangeoffset),"\t",($3+$rangeoffset),"\t\t",sprintf("%.3f",(-log($4)/log(2))),"\n";
						print OUT "+ \t",($1+$rangeoffset),"\t",($2+$rangeoffset),"\t",($3+$rangeoffset),"\t\t",sprintf("%.3f",$4),"\n";
					}
				}
			}
			print OUT "\n";
		} else {
			die "Never reach in print_output.\n";
		}
	}	
	print OUT "#Total number of sequences: $num_line.\n";
	print OUT "#Positive sequences: $pos.\n";
	print OUT "#The polyasvm score cutoff value is $cutoff.\n\n";
	
	close(OUT_JUNK);
	system("rm $output_junk -f");
}

sub is_pos_m{
	my $start = $_[0];
	my $len = $_[1];
	my $pos_win_size;
	my $time_up_size;
	if ($len > 50){
		$pos_win_size = 30;
		$time_up_size = 10;
	} elsif ($len >1){
		$pos_win_size = int(0.592*$len+0.408);
		$time_up_size = int(0.183*$len+0.817);
	}

	my @vec = @pro_data[$start .. ($start+$len-1)];
for (my $i=0;$i<=$#vec;$i++){$vec[$i]=log($vec[$i])/log(2);}
	my $vec_len = @vec;
my $prod=0;
	my $product=1;
	my $max_product=$p_cutoff;
$max_product=-$cutoff;
	my ($k,$j,$m);
	my @return_arr;
	my $return_pos=0;
	my $pos_size=0;
	my $region_win=0;
	my $region_win_start =0;
	for ($k=0;$k<$vec_len;$k++){
		#if ($vec[$k]>=0.5){
		if ($vec[$k]>=-1){
			if ($region_win == 0){
				$region_win_start = $k;
			}
			$region_win++;
		}
		#if (($k==($vec_len-1))||($vec[$k]<0.5)){
		if (($k==($vec_len-1))||($vec[$k]<-1)){
			if ($region_win>=$pos_win_size){
				#$product = 1;
				$prod = 0;
				my @region_pro=@vec[$region_win_start .. $k];
				my $max_j =0;
				for ($j=0;$j<($region_win-$time_up_size+1);$j++){
					for ($m=0;$m<$time_up_size;$m++){
						#$product = $product*$region_pro[$j+$m];
						$prod = $prod + $region_pro[$j+$m];
					}
					if ($prod>$max_product){
						#$max_product = $product;
						$max_product = $prod;
						$max_j = $j+int($time_up_size/2);
					}
					$prod=0;
					#$product=1;	
				}
				#if ($max_product >$p_cutoff){
				if ($max_product >-$cutoff){
					#my $temp_ret = $region_win_start.":".($region_win_start+$region_win).":".($region_win_start+$max_j).":".sprintf("%.3f",-$max_product);
					my $temp_ret = $region_win_start.":".($region_win_start+$region_win-1).":".($region_win_start+$max_j).":".sprintf("%.3f",-$max_product);
					push(@return_arr,$temp_ret);
					$return_pos++;
				}
				$max_product = $p_cutoff;
$max_product=-$cutoff;
				$region_win = 0;
			} else {
				$region_win=0;
			}
		}

	}
	push(@return_arr,$return_pos);
	return(@return_arr);
}

sub is_pos_hpr{
	my $start = $_[0];
	my $len = $_[1];
	#if (defined($rangeoffset)){$len=$x2-$x1+$FOS+1;}
	if (defined($rangeoffset)){$len=$x2-$x1+1;}
	my @index=($start .. ($start+$len-1));
	my @rev_posmap;
	my $rpc=0;
	for (my $i=0;$i<@posmap;$i++){
		if ($posmap[$i]>=0){
			push @rev_posmap,$i;
			if ($rpc!=$posmap[$i]){print "rpc $rpc != posmap $posmap[$i]\n";exit;}
			$rpc++;
		}
	}
#for (my $i=0;$i<@rev_posmap;$i++){print "$i $rev_posmap[$i]\n";}
#print "size = ",scalar(@index),"\n";
	my $time_up_size = $hpr; #### input from the program.
	my @vec = @pro_data[@index];
#for(my $i=0;$i<@vec;$i++){print "$i $vec[$i]\n";}
	my $temp=log(2);
	for (my $i=0;$i<=$#vec;$i++){
		if (!defined($vec[$i])){$vec[$i]=-100;next;}
		if ($vec[$i]==0){$vec[$i]=-100;next;}
		$vec[$i]=log($vec[$i])/$temp;
	}
	if ($len<$hpr){
		$temp=-$cutoff/$hpr;
		for (my $i=$#vec+1;$i<$hpr;$i++){
			#$vec[$i] = 0.8781261;
			$vec[$i]=$temp;
		}
	}
	my ($k,$j);
	my @return_arr;
	my $return_pos=0;
	my $pos_size=0;
	$k=0;
	my $prod=0;
	for ($j =0;$j<$time_up_size;$j++){
		$prod+=$vec[$j];
	}
	if ($len<$hpr){
		if ($prod>-$cutoff){
			my $temp_ret ="$k:".($k+$len-1).":".($k+int($len/2)).":".sprintf("%.3f",-$prod);
			$return_pos++;
			push @return_arr,$temp_ret;
			push @return_arr,$return_pos;
			return(@return_arr);
		}
	}
	my @parray=();
	my @ppos=();
	my $state=0;
	my $laststate=-1;
	my $temp_ret;
	while ($k<($len-$time_up_size)){
		push @parray, $prod;
		if ($prod>=-$cutoff){$state=1;}else{$state=-1;}
		if ($state!=$laststate){
			push @ppos, $k;	
		}
		$laststate=$state;
if (!defined($vec[$k+$time_up_size])){print "vec end not defined! end=".($k+$time_up_size)."\n";exit;}
if (!defined($vec[$k])){print "vec start not defined!\n";exit;}
        	$prod = $prod+$vec[$k+$time_up_size] -$vec[$k];
        	$k++;
        }      
	if ($state==1&&$laststate==1){push @ppos,$k-1;} 
	$k=$#ppos;
	my ($maxscore,$maxpos);
	while($k>0){
		$laststate=shift @ppos; #front of hpr window
		$state=(shift @ppos);	#end of hpr window
		##merge regions that overlap
		($maxscore,$maxpos)=&get_max2(\@parray,$laststate,$state);
		$temp_ret = "$maxpos:".($maxpos+$time_up_size-1).":".($maxpos+int($time_up_size/2)).":".sprintf("%.3f",-$maxscore);
		$return_pos++;
		push @return_arr,$temp_ret;
		$k-=2;
	}
	push(@return_arr,$return_pos);
	return(@return_arr);
}

sub load_settings{
open(READ,"<$_[0]");
while(<READ>){
        chomp;
        next if ($_ =~ /^\#/ || length($_)<1);
        push @array,$_;
}
for my $line (@array){
chomp $line;
@break=split /=/,$line;
if ($break[0] eq "input path"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$inpath=$break[1];}   
}
if ($break[0] eq "model path"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$modelpath=$break[1];}   
}
if ($break[0] eq "output path"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$outpath=$break[1];}
}
if ($break[0] eq "svm predict path"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$svmpath=$break[1];}
}
if ($break[0] eq "output file"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$out=$break[1];}
}
if ($break[0] eq "input file"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])<1){die "no input file";}else{$in=$break[1];}
}
if ($break[0] eq "mode"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$mode=$break[1];}
}
if ($break[0] eq "cutoff"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$cutoff=$break[1];}
}
if ($break[0] eq "location"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$location=$break[1];}
}
if ($break[0] eq "model"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$train_model=$break[1];}
}
if ($break[0] eq "range"){
	next if (!defined($break[1]));
        $break[1]=~ s/\n//g;
        if (length($break[1])>=1){$range_file=$break[1];}
}
}
close(READ);
}
sub get_scoremat{
my $pos=$_[0];
my @stripe=("pos:$pos");
my @returnArray=();

if (defined($_[1])){push @returnArray,join "\t", ("cis ele",@file_ext);}
foreach my $var (@file_ext){
	push @stripe, $score_mat{$var}[$pos];
}
push @returnArray, join "\t",@stripe;
return join "\n",@returnArray;

}


sub paslookup{
my $temp="";
my @pas=("AATAAA","ATTAAA","TATAAA","AGTAAA","AAGAAA","CATAAA","GATAAA","TTTAAA","ACTAAA","AATATA","AATACA","AATGAA","AATAGA");
my %pashash;
my $pos=0;
	foreach my $hex (@pas){
	while(1){
		$temp=index($_[0],$hex,$pos);
		$pos=$temp+1;
		if ($temp<0){last;}
		if (!defined($pashash{$hex})){$pashash{$hex}=[$temp];}else{push @{$pashash{$hex}}, $temp;}
	}
	$pos=0;
	}
	return %pashash;
}

sub print_pas{
	#assume global %phash is populated and PAS file handle is open
	if (defined($_[0])){print PAS "\#$_[0]\n";}
	foreach my $pas (keys %phash){
		print PAS ">$pas\n";
		foreach my $position (@{$phash{$pas}}){
			print PAS "$position\t";
		}
		print PAS "\n";
	}	
}

sub negseq{
$_[0]=~s/\s//g;
$_[0]=uc($_[0]);
$_[0]=~tr/ATCG/TAGC/;
$_[0]=reverse($_[0]);
return $_[0];
}



###
__DATA__
AUE1	8.07	5.36	2.69	1389	1.031025	4.69755043500543
AUE2	6.75	4.50	2.26	2086	3.591545	2.62519776295469
AUE3	5.91	4.49	3.21	92	1.27952	4.17734050041748
AUE4	6.92	4.62	2.31	1871	2.74084500000001	3.26155992165348
CUE1	6.31	4.21	2.11	2662	2.04321	3.32550269583261
CUE2	5.44	3.63	1.82	1380	2.38214	4.12573749471991
CDE1	7.59	5.06	2.54	2614	2.707035	3.01568103831248
CDE2	5.03	3.36	1.68	1826	2.52439499999999	3.06693043043164
CDE3	6.58	4.39	2.20	1420	1.05726	4.14213301637727
CDE4	6.75	4.48	2.25	1265	1.81868	3.87491316465502
ADE1	8.98	7.69	5.83	65	-0.134514999999999	5.23918225299664
ADE2	8.66	4.05	2.30	414	-0.956804999999999	5.08174206283353
ADE3	8.50	5.67	2.84	2377	2.67058	3.86396103254555
ADE4	8.33	6.44	4.39	402	0.47064	5.00757801828435
ADE5	9.05	6.04	3.02	3168	3.08603	3.28865995220556
#AUE1.MAT
#pos	A	C	G	U
1	-0.73	0.41	0.80	-0.54
2	-0.33	-1.05	-0.98	0.84
3	-2.59	-2.70	2.09	-2.16
4	-2.00	Inf	2.22	Inf
5	Inf	Inf	2.32	Inf
6	Inf	Inf	2.32	Inf
7	0.10	-0.29	0.12	-0.01
8	-0.97	-0.94	1.39	-0.49
9	-0.18	0.14	0.45	-0.29
#AUE2.MAT
#pos	A	C	G	U
1	-0.02	-0.00	0.50	-0.40
2	-0.64	-0.84	-0.45	0.79
3	-0.36	-0.78	0.97	-0.24
4	-5.30	-4.27	-4.05	1.57
5	-1.97	-0.87	-1.44	1.19
6	Inf	Inf	Inf	1.61
7	-1.97	-2.24	1.80	-0.77
8	-1.60	-1.31	-1.48	1.21
9	-0.01	-0.42	-0.10	0.26
10	-0.03	-0.00	-0.16	0.12
#AUE3.MAT
#pos	A	C	G	U
1	0.00	-1.54	-1.30	0.81
2	Inf	Inf	Inf	1.61
3	1.19	Inf	Inf	0.17
4	-1.07	-1.79	Inf	1.31
5	1.85	Inf	Inf	Inf
6	Inf	Inf	Inf	1.61
7	0.69	-0.90	0.09	-0.49
#AUE4.MAT
#pos	A	C	G	U
1	0.34	-0.05	-0.08	-0.28
2	0.55	-0.27	-0.49	-0.17
3	0.64	-0.07	-1.44	-0.08
4	-5.10	-4.27	-4.63	1.57
5	Inf	Inf	2.32	Inf
6	Inf	Inf	Inf	1.61
7	1.72	-3.07	-2.35	-3.64
8	0.28	0.19	-1.88	0.24
9	0.46	-0.36	-0.81	0.11
10	-0.10	-0.14	0.03	0.13
#CUE1.MAT
#pos	A	C	G	U
1	0.01	-0.22	-0.10	0.12
2	-0.05	-0.28	-0.12	0.21
3	-0.29	-0.27	0.24	0.28
4	-0.71	-0.91	-1.02	0.94
5	0.43	-1.34	0.49	-0.56
6	Inf	-1.95	Inf	1.59
7	-2.65	-5.20	-3.88	1.54
8	-3.72	-3.88	-3.21	1.57
9	-2.41	-3.11	-0.71	1.36
10	0.29	-0.42	-0.32	-0.09
11	0.27	-0.50	-0.55	0.05
#CUE2.MAT
#pos	A	C	G	U
1	0.10	-0.24	-0.04	-0.00
2	0.19	0.20	-0.47	-0.17
3	1.04	-1.77	-1.34	-1.56
4	1.26	-3.19	-3.03	-2.91
5	Inf	Inf	Inf	1.65
6	1.38	Inf	Inf	Inf
7	1.14	-2.29	-1.77	-2.05
8	0.91	-1.38	-1.20	-1.05
9	0.02	0.06	0.36	-0.27
10	-0.03	-0.10	-0.07	0.11
#CDE1.MAT
#pos	A	C	G	U
1	-0.26	-0.06	0.01	0.17
2	-0.29	0.98	-0.39	-0.37
3	-1.36	1.33	0.67	-1.54
4	-4.75	-1.94	-1.76	1.24
5	Inf	-1.34	2.22	Inf
6	-1.30	0.42	-0.07	0.32
7	-3.52	2.30	-2.98	-2.68
8	-2.11	-1.36	-2.34	1.15
9	-0.52	0.65	0.44	-0.45
10	-0.01	-0.18	0.15	0.01
#CDE2.MAT
#pos	A	C	G	U
1	-0.52	-0.56	-0.54	0.60
2	-2.60	-1.11	-0.15	0.89
3	-4.46	-4.69	-4.17	1.37
4	-0.57	1.15	0.73	-1.88
5	Inf	Inf	Inf	1.42
6	Inf	-1.07	0.17	0.88
7	-1.65	-1.77	-1.56	1.10
8	-0.21	0.06	-0.27	0.22
9	-0.03	-0.18	-0.11	0.15
#CDE3.MAT
#pos	A	C	G	U
1	-0.07	0.24	0.14	-0.17
2	-0.44	0.73	-0.03	-0.21
3	-2.45	-3.79	-1.35	1.21
4	Inf	-1.51	2.23	Inf
5	Inf	Inf	Inf	1.42
6	Inf	Inf	2.33	Inf
7	-2.59	-2.17	-0.75	1.09
8	-1.11	0.12	0.90	-0.24
9	0.06	-0.10	0.02	-0.00
#CDE4.MAT
#pos	A	C	G	U
1	-0.38	0.37	-0.26	0.13
2	-1.44	0.76	0.49	-0.21
3	-3.39	-1.38	-3.23	1.24
4	Inf	2.46	Inf	Inf
5	Inf	Inf	Inf	1.42
6	Inf	-0.75	2.16	Inf
7	-1.72	-0.15	0.19	0.51
8	-0.46	-0.03	-0.18	0.32
9	-0.14	-0.11	0.54	-0.23
#ADE1.MAT
#pos	A	C	G	U
1	-1.57	1.40	-1.51	0.02
2	Inf	2.36	Inf	Inf
3	-0.79	Inf	Inf	1.43
4	Inf	2.36	Inf	Inf
5	Inf	2.36	Inf	Inf
6	-0.58	2.06	Inf	Inf
7	-0.79	1.35	-0.71	-0.57
#ADE2.MAT
#pos	A	C	G	U
1	-0.07	0.17	-0.02	-0.04
2	-0.46	0.63	0.30	-0.42
3	Inf	2.36	Inf	Inf
4	Inf	2.36	Inf	Inf
5	0.83	1.38	Inf	Inf
6	Inf	Inf	2.21	Inf
7	-4.11	2.31	-4.17	-6.28
8	-1.96	0.24	-0.90	0.90
#ADE3.MAT
#pos	A	C	G	U
1	-0.00	-0.21	0.21	-0.03
2	-0.59	-0.84	0.89	-0.02
3	-2.11	-1.60	1.63	-0.64
4	-3.86	1.22	1.23	-3.95
5	Inf	Inf	-1.52	1.57
6	Inf	Inf	2.21	Inf
7	-3.65	-4.14	2.13	-4.37
8	-1.13	-0.67	1.47	-1.04
9	0.03	-0.51	0.71	-0.46
10	0.12	-0.08	0.01	-0.07
#ADE4.MAT
#pos	A	C	G	U
1	-0.35	-0.28	0.67	-0.16
2	-0.88	-0.83	1.46	-1.11
3	Inf	Inf	2.21	Inf
4	Inf	Inf	2.21	Inf
5	Inf	-0.24	1.94	Inf
6	Inf	2.36	Inf	Inf
7	1.01	-0.09	-3.30	-0.39
8	-0.93	0.02	1.22	-1.02
#ADE5.MAT
#pos	A	C	G	U
1	0.03	0.05	-0.04	-0.04
2	-0.49	0.49	0.19	-0.15
3	0.35	-0.33	0.02	-0.20
4	-1.53	-1.59	1.79	-1.83
5	-1.58	-5.27	2.04	-4.95
6	-2.17	-0.03	1.79	Inf
7	1.71	-1.59	-4.29	-3.95
8	-2.78	-2.69	2.00	-2.34
9	0.16	-1.71	1.18	-1.15
10	-0.55	-0.15	1.00	-0.63
11	-0.03	0.26	-0.07	-0.10
###


