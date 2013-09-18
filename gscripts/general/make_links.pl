use strict;
use warnings;

my $files_path = shift @ARGV;

open (DATAPATHS, "<".$files_path);
while( my $line  = <DATAPATHS>){
	chomp($line);
	my @values = split(/\t/, $line);

	if($values[1]){
		`ln -s $values[1] ./$values[0]`
	}else{

		`ln -s $values[0]`;

	}


}
