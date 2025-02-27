#!/usr/bin/perl 

#usage: ./Generate_UTRsite_motifs_descr.pl UTRsite_motifs_detail.txt

sub trim($)
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

@files = glob "./UTRSite/*";

open(OUT1,">UTRsite_motifs_list.txt") or die "UTRsite_motifs_list.txt: $!";

foreach $file(@files)
{
	if(index($file,"./UTRSite/UX")!=0)
	{
		@pageID = glob $file."/pageID/*";

		open(IN ,"$pageID[0]/index.html");
		
		while(chomp($line=<IN>))
		{
			if($line=~ m|<td class=\"fieldValue\" colspan=\"2\" style=\"font-weight:bold;\">(.*)</td>|s)
			{
				$MotifName = $1;

				open(OUT,">$MotifName.descr") or die "$MotifName.descr: $!";
				print OUT1 $MotifName;
			}
			if(index($line,"                        <td class=\"fieldValue\" colspan=\"2\"><pre>")!=-1 || index($line,"                        <td class=\"fieldValue\"><pre>")!=-1)
			{
				if(index($line,"</pre></td>")!=-1)
				{
					$line=~ m|<pre>(.*)</pre>|s;

					$pattern = $1;
					
					$pattern=~ s/&gt;/>/g;

					#print $MotifName."\n".$pattern."\n";
					print OUT $pattern;
					close(OUT);
				}
				else
				{
					$line=~ m|<pre>(.*)|s;

					$pattern = $1;
					
					while(chomp($line=<IN>))
					{
						if(index($line,"</pre></td>")!=-1)
						{
							if(index($line,"</pre></td>")==0)
							{
								last;
							}
							else
							{
								$line=~ m|(.*)</pre>|s;

								$pattern = $pattern."\n".$1;

								last;
							}
						}
						else
						{
							$pattern = $pattern."\n".$line;
						}
					}
					
					$pattern=~ s/&gt;/>/g;

					#print $MotifName."\n".$pattern."\n";
					print OUT $pattern;
					close(OUT);
				}
			}
			if(index($line,"Standard Name")!=-1)
			{
				chomp($line=<IN>);

				$line=~ m|>(.*)</td>|s;

				print OUT1 "|".$1;
			}
			if(index($line,"UTR Region")!=-1)
			{
				chomp($line=<IN>);

				$line=~ m|>(.*)</td>|s;

				print OUT1 "|".trim($1);
			}
			if(index($line,"Taxon Range")!=-1)
			{
				chomp($line=<IN>);
				chomp($line=<IN>);

				$line=~ m|(.*)</td>|s;

				print OUT1 " ".trim($1)."\n";
			}
		}
		close(IN);
	}
}
close(OUT1);