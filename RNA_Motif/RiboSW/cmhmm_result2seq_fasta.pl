#!/usr/bin/perl 

###Ū��Seq file###
my $input_file=$ARGV[0];

my $reading_flag=0;
my $line_no=0;

my $seqname;
my $startp;
my $endp;

my $start_read_flag=0;

open(IN, $input_file) or die "$input_file: $!";
while (defined(my $line = <IN>)) 
{	
	$line_no++;	
	chomp($line);
	my @split_line=split(' ',$line);		
	
	if(index($line,">")==0) #�I��>�}�Y�A��seqname���X��
	{
		$seqname=$line;
	}
	elsif(index($split_line[0],"Query")==0) #�I��Query�}�Y�A��᭱Target = 168 - 267����m���X��
	{
		$startp=$split_line[7];
		$endp=$split_line[9];
		print "$startp $endup";		
	}	
	elsif(index($split_line[0],"Score")==0)
	{
		$start_read_flag=1;		
	}
	else	
	{
		if($start_read_flag==1 && @split_line==0) #����bScore�U���A�B�u���@�Ӵ���Ÿ�(�h�U���|�檺�̫�@��h��result seq���@)
		{
			$line = <IN>;
			$line = <IN>;
			$line = <IN>;
			$line = <IN>;			
			my @split_line=split(' ',$line);					
			#print "$split_line[3]";
			
			if($split_line[2]==$endp) #����̫�@�q�����Ʀr��$endp�@�ˮɡA��result seq�~����Q����
			{
				print "\n";
				$start_read_flag=0;
			}
		}
	}
}



close(IN);
