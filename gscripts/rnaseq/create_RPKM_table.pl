####
#
#	concatenate multiple single-lane RPKM files
#	filter flagged and lowly-expressed genes
# 
#
#	ppliu
#
####

# filter and replacement values for flagged and low RPKMs
my $RPKM_CUTOFF = 0.5;
my $RPKM_REPLACEMENT = 0.1;

# global variables
my %all_rpkms;			# hash to hold all RPKM values from each input file
my @file_handler_list; 		# list of all file handlers created from inputted file paths
my @file_paths;			# list of filename paths of RPKM files

# initalize variables 
# grab list of file paths passed from command line
my @file_paths = @ARGV;

# populate RPKM values hash with values from each file taking into 
# account cutoffs and value replacements
foreach my $file (@file_paths){

	# file existence check 
	unless (-e $file){
		die("$file does not exist!\n");
	}

	# open file handler;
	open(RPKM_FILE, "<".$file);

	# skip RPKM header
	<RPKM_FILE>;

	while (<RPKM_FILE>){
		
		# remove whitespace
		chomp $_;

		# split line
		my ($cluster_id, $flag, $rpkm) = split("\t", $_);

		# skip IDs with flag 2, expected to be Ensembl predictive genes
		if ($flag eq 2){

			next;
		}

		# store RPKM value in hash by cluster id and filename keys
		$all_rpkms{$cluster_id}{$file} = $rpkm;

	}

}

# print header
print("gene_ID\t");
print join("\t", @file_paths);
print "\n";

# filter then print values
foreach my $cluster_id (keys %all_rpkms){

	# flag if at least one RPKM value is above the cutoff
	my $above_cutoff = 0;

	foreach $file (@file_paths){
	
		if ($all_rpkms{$cluster_id}{$file} < $RPKM_CUTOFF){

			$all_rpkms{$cluster_id}{$file} = $RPKM_REPLACEMENT;

		}else{

			$above_cutoff = 1;

		}

	}

	if ($above_cutoff){
		
		print $cluster_id;

		foreach (@file_paths){

			print "\t";
			print $all_rpkms{$cluster_id}{$_};

		}
		print "\n";


	}else{

		print STDERR "$cluster_id did not meet the cutoff RPKM value of $RPKM_CUTOFF\n";


	}


}
