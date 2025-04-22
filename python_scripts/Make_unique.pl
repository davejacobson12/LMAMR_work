open(fh,"$ARGV[0]");
%count;

while($line=<fh>.<fh>)
{
	chomp $line;
	@data=split("\n",$line);
	chomp $data[0];
	chomp $data[1];
	$count{$data[1]}++;
}

$seq=1;

foreach $int (keys %count)
{
	chomp $int;
	chomp $count{$int};
	print ">Seq",$seq,";size=",$count{$int},";\n";
	print "$int\n";
	$seq++;
}
close fh;

