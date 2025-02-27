#!/usr/bin/perl

# Converts parenthesized alignment to .epn format
# This program understands parenthesized alignemnt with pseudoknots denoted by curly or square brackets.
# Typical input format:
#
#  >secondary structure 
#  (((((-----((((((((----(((((((-------------))))--)))---))))))-))-(((((((--((((((((---))))))))--)))))))---)))))
#  >Ver.spinosum
#  GGUGACGAUAGCGAGAAGGUCACACCCGUUCCCAUACCGAACACGGAAGUGAAGCUUCUCAGCGCCGACGGUAGUUCGGGGCUGCCUCCCGUGAGAGUAGGACGUUGCC
#  >Emb.brevis
#  CGGCGGCCAUAGCGGCAGGGAAACGCCCGGUCCCAUUCCGAACCCGGAAGCUAAGCCUGCCAGCGCCGAUGGUACUGCACUCGUGGGGUGUGGGAGAGUAGGUCACCGC
#  >F.johnsoniae
#  GGUGGCGAUAGCAGAGAGGUCACACCCGUUCCCAUGCCGAACACGGAAGUUAAGUUCUCUAGCGCCGAUGGUAGUGUGGACUUUGUCCCUGUGAGAGUAGGACGUUGCC
#  >Bor.burgdorferi
#  GUGGUUAA-AGAAAAGAGGAAACACCUGUUAUCAUUCCGAACACAGAAGUUAAGCUCUUAUUCGCUGAUGGUACUGCG---AGU---CGCGGGAGAGUAGG-UUAUUGC
#  >Trp.pallidum
#  GUUGCCAU-GGUGGAGAGGUCAUACCCGUUCCCAUCCCGAACACGGAAGUCAAGCUCUCCUACGCCGAUGAUACUGCU---CCU---CGCGGGAAAGUAGG-UAGUAGC
#


$open{")"} = "(";
$open{"]"} = "[";
$open{"}"} = "{";


if ($#ARGV < 0) {
    print "$0 <parenthesized alignment>\n";
    print "   - Reads parenthesized alignment (fasta format) and translates it into .epn (ERPIN2) format. \n";
    print "   - Checks that helices have no gap. Otherwise change positions to single-stranded.\n";
    print "   - Deletes columns with gaps only\n";
    print "   See Perl source file for more information on parenthesized alignments\n";
    print "   Daniel Gautheret, 2002\n";
    exit;
}

if (!open(F, $ARGV[0])) {
    print "Can't open \"$ARGV[0]\"\n";
    exit;
}

while (($line = <F>) !~ /^>/) {}
$namestr = $line;
while (($line = <F>) && ($line !~ /^>/)) {
    chomp $line;
    $str .=  $line;
}

#===============================================================
# gap check

# constructs pair list
# p, q, u distinguish 5' helix, 3' helix and single strand

for ($i=0; $i<length($str); $i++) {
    $c=substr($str,$i,1);
    if ($c =~ /[\(\[\{]/ ) {  
	push @{ $pstack{$c} }, $i;
    }
    if ($c =~ /[\)\}\]]/) {  
	$j = pop @{ $pstack{$open{$c}} };
	$pair[$i]=$j;
	$pair[$j]=$i;
    }
}

# Now reads each sequence and looks for gap in helices. 
# Modifies parenthesis list accordingly

$seqnum=0;
while ($line) {
    if ($line =~ /^>/) {
	$seqname[$seqnum] = $line;
	$seq="";
	while (($line = <F>) && ($line !~ /^>/)) {
	    chomp $line;
	    $seq .=  $line;
	}

	$lali = length($seq);
	for ($i=0; $i<$lali; $i++) {	
	    if ((defined ($pair[$i])) && (substr($seq,$i,1) eq "-")) {
		$j=$pair[$i];
		print STDERR "Warning: Gap in helix. Seq $seqnum (setting $i:$j as single stranded)\n";
		substr($str,$i,1,"-");
		substr($str,$j,1,"-");
	    }
	}
	$seqs[$seqnum]=$seq;
	$seqnum++;
#	if ($seqnum%50==0){print STDERR $seqnum,"\n";}
    } else {
	$line = <F>;
    }
}

close F;

# ===========================================
# Now checks for fully gapped columns

foreach $seq (@seqs) {
    for ($i=0; $i<$lali; $i++) {	
	if (substr($seq,$i,1) ne "-") {
	    $nogap[$i]=1;
	}
    }
}

$removed = 0;
for ($i=0; $i<$lali; $i++) {	
    if ($nogap[$i] != 1) {
	$removed++;
	print STDERR "Removing column $i (gaps only) \n";
	substr ($str,$i,1,"Z");
	foreach $seq (@seqs) {
	    substr ($seq, $i, 1, "Z");
	}
    }
}

$str =~ s/Z//g;
foreach $seq (@seqs) {
    $seq =~ s/Z//g;
}


$lali = length ($str);

print STDERR "$lali columns\n";

print STDERR $str, "\n";

#===============================================================
# REconstructs pair list
# p, q, u distinguish 5' helix, 3' helix and single strand

for ($i=0; $i<length($str); $i++) {
    $c=substr($str,$i,1);
    if ($c =~ /[\(\[\{]/ ) {  
	push @{ $pstack{$c} }, $i;
	$status[$i]="p";
    }
    if ($c =~ /[\)\}\]]/) {  
	$j = pop @{ $pstack{$open{$c}} };
	$pair[$i]=$j;
	$pair[$j]=$i;
	$status[$i]="q";
    }
    if ($c !~ /[\(\[\{\)\}\]]/) {  
	$status[$i]="u";
    }
}

#===============================================================

#construct helix/strand numbering from parenthesis list, for printing

$b="X";
$lasthelix=0;
$laststrand=-1;
$curnum=0;

for ($i=0; $i<length($str); $i++) {
    #change element
    if (($i==0) || ($status[$i] ne $status[$i-1])){
	if ($status[$i] eq "u"){
	    $curnum = $laststrand+2;
	    $laststrand=$curnum;
	    $erpinstruc[$i]=$curnum;
	}else{
	    $j=$pair[$i];
	    if ($i<$j) {
		$curnum = $lasthelix+2;
		$lasthelix=$curnum;
		$erpinstruc[$i]=$curnum;
		$erpinstruc[$j]=$curnum;
	    }
	}
    #same element (except if helix with bulge on 3')
    }else{
	if ($status[$i] eq "u"){
	    $erpinstruc[$i]=$curnum;
	}else{
	    $j=$pair[$i];
	    if ($i<$j) {
		if ($status[$j] eq $status[$j+1]){		
		    $erpinstruc[$i]=$curnum;
		    $erpinstruc[$j]=$curnum;
		}else{
		    $curnum = $lasthelix+2;
		    $lasthelix=$curnum;
		    $erpinstruc[$i]=$curnum;
		    $erpinstruc[$j]=$curnum;
		}
	    }
	}
    }
}

# Descriptor formating in N raws
undef(@epn_desc);

for ($i=0; $i<length($str); $i++) {
    undef(@epn_desc_pos);
    @epn_desc_pos = split("",sprintf '%0*d',length($curnum),$erpinstruc[$i]);
    for($j=0;$j<length($curnum);$j++){
	$epn_desc[$j][$i]=$epn_desc_pos[$j];
    }
}

print "$namestr";
for($i=0;$i<length($curnum);$i++){
    for($j=0;$j<length($str);$j++){
	print $epn_desc[$i][$j];
    }
    print "\n";
}

$k=0;
foreach $seq (@seqs) {
    print $seqname[$k++];
    print $seq, "\n";
}







