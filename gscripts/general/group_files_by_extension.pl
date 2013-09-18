####
#
#   Run in directory to group files by their extension (string after last period).
#   Will ignore dashes and things following it (err and out output of array jobs)
#
####

use strict;
use warnings;


my @files = <*>;

foreach (@files){


	my @filename_values = split(/\./, $_);
	
	if($filename_values[-1] and -f $_ and $filename_values[-1] ne $_){
		my $extension = $filename_values[-1];
        my @ext_values = split(/\-/, $extension);
        if (scalar(@ext_values) > 1){
            $extension = $ext_values[0];
        }

		my $folder_name = $extension."_files";

		print $extension."\n";
		unless(-e $folder_name){
			`mkdir $folder_name`;

		}

		`mv $_ $folder_name`;
	}
}
		
