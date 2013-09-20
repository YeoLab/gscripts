###
#
#   Provide a tab delimited file of link filename and path to 
#   original file. Will go through file and make links for each
#   pair of link filename and path. One per line. 
#
###

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
