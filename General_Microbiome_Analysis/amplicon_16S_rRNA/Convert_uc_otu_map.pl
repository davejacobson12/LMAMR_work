open(fh,"$ARGV[0]");
my %Array;
open(fail,">$ARGV[2]");
open(ftable,">$ARGV[1]");
while($line=<fh>)
{
	chomp $line;
	@data=split("\t",$line);
	chomp $data[0];
	if($data[0] eq "H")
	{
		chomp $data[8];
		chomp $data[9];
		push (@{ $Array {$data[9]} }, $data[8]);
		
	}
	else
	{
		chomp $data[8];
		print fail "$data[8]\n";
	}
	
}
close fh;

foreach $key(keys %Array)
{
	print ftable "$key\t";
	foreach (@{$Array{$key}})
	{
		chomp $_;
		print ftable "$_\t";
	}
	print ftable "\n";
}

close fail, ftable;
