#!/usr/bin/perl 

my $input_file = $ARGV[0];
my $reading_flag = 0;
my $line_no = 0;

my $seqname;
my $startp;
my $endp;
my $Evalue;
my $start_read_flag = 0;
my $score;
our $global_motif = "" ;
open(IN, $input_file) or die "$input_file: $!";

while (defined(my $line = <IN>)) {
    chomp($line);
    my @split_line = split(/\s+/, $line);
      # 使用 \s+ 来分割多个空格
    if ($split_line[0] eq "Query:"){
        #print "$line\n" ;
    }

    if ($split_line[0] eq "Description:") {
        my @motif = split(/\s+/, $line);
        $global_motif = join(" ", @motif[1..$#motif]);  # 直
    }
}
close(IN);

open(IN, $input_file) or die "$input_file: $!";
while (defined(my $line = <IN>)) {
    chomp($line);
    my @split_line = split(/\s+/, $line);
      # 使用 \s+ 来分割多个空格
    
    # 处理序列名
    if (index($line, ">>") == 0) {
        $seqname = $split_line[1];
        $line = <IN>;  # 读取一行 
        $line = <IN>;  # 读取一行
        $line = <IN>;  # 读取一行
		my @split_line = split(/\s+/, $line);
		if (@split_line >= 13) {
        # 获取start and end 
            my $Evalue= $split_line[3] ;
        	my $startp = $split_line[10];  # 第10个元素
        	my $endup = $split_line[11];    # 第11个元素
            my $score = $split_line[4] ;
            #print "E-value: $Evalue\n";
            #print "Score: $score\n";
            print "$startp $endup $global_motif $score\n" ;
    	
        }
    }
    else 
    {   

        if ($split_line[1] eq $seqname) 
        {
            #print "$split_line[3]";
   
        }
    }
}


close(IN);
