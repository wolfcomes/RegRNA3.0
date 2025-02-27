#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use File::Spec;
use File::Path qw(make_path);  # 引入 File::Path 模块

my %opts;
my $scriptdir=$ENV{PWD};
getopts('f:s:o:h',\%opts);  # 增加 -o 选项

if ($opts{'h'}) {
    print "STD  output
      -f   fasta file.
      -s   secondary structure file (optional)(format like example/1ddy_A_dot.txt).
      -o   output directory (optional).
      -h   help.\n";
    exit;
}

exit() unless ($opts{'f'} || $opts{'s'});

# 设置输出目录，默认使用当前目录
my $out_dir = $opts{'o'} || $scriptdir;  # 需要声明
# 如果输出目录不存在，则创建目录
unless (-d $out_dir) {
    make_path($out_dir) or die "Failed to create directory $out_dir: $!\n";
}

if($opts{'f'}) {
    my ($file1, $base_dir, $fasta_name) = File::Spec->splitpath($opts{'f'});
    $fasta_name=~s/\.fas.upper//g;
    my $name=$fasta_name;

    `/home/RegRNA/.conda/envs/regrna/bin/RNAfold --noPS -i $base_dir/$name.fas.upper >$out_dir/$name\_dot.txt`; # This should be replaced with your RNAfold directory
    `$scriptdir/dot2motif.pl $out_dir/$name`;
    `$scriptdir/dot2motif2.pl $out_dir/$name`;
    `$scriptdir/process_motif.pl $out_dir/$name`;
    `$scriptdir/process_motif2.pl $out_dir/$name`;
    `$scriptdir/search_motif.pl $out_dir/$name`;
    `$scriptdir/search_R-BIND.pl $out_dir/$name`;
    `$scriptdir/search_miRBase.pl $out_dir/$name`;
}

if($opts{'s'}) {
    my ($file, $base_dir, $ss_name) = File::Spec->splitpath($opts{'s'});
    $ss_name=~s/\.(.+)//g;
    my $name=$ss_name;
    my $format=$1;

    `cp "$out_dir/$name$format" "$out_dir/$name\_dot.txt"`;
    `$scriptdir/dot2motif.pl $out_dir/$name`;
    `$scriptdir/dot2motif2.pl $out_dir/$name`;
    `$scriptdir/process_motif.pl $out_dir/$name`;
    `$scriptdir/process_motif2.pl $out_dir/$name`;
    `$scriptdir/search_motif.pl $out_dir/$name`;
    `$scriptdir/search_R-BIND.pl $out_dir/$name`;
    `$scriptdir/search_miRBase.pl $out_dir/$name`;
}

