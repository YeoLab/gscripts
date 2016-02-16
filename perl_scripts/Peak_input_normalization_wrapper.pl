use warnings;
use strict;
use Statistics::Basic qw(:all);


my %rbp_list;
my @peak_filelist;
my %clipshort_to_peakshort;
my @clip_bam_fi_shorts;
my @peak_filelist_shorts;

my $filist_fi = $ARGV[0];
my $output_folder = $ARGV[1];

unless ($filist_fi && $output_folder) {
    print STDERR "Usage Error\n";
    print STDERR "Usage: perl Peak_input_normalization_wrapper.pl filist_fi output_folder\n";
    exit;
}

## Reads in mapped_read_num file if it already exists - useful for re-running analysis more quickly on the same manifest file
my %mapped_num;
my $filist_fi_mappedreadnums = $filist_fi.".mapped_read_num";
if (-e $filist_fi_mappedreadnums) {
    open(READNUM,"$filist_fi_mappedreadnums");
    for my $line (<READNUM>) {
	chomp($line);
	my ($fi,$num) = split(/\t/,$line);
	$mapped_num{$fi} = $num;
    }
    close(READNUM);
}

unless ($output_folder =~ /\/$/) {
    $output_folder = $output_folder."/";
}

my @rep_listing;
my $type_flag;

my %peakfi2uid;

## First, read in manifest file & do pre-processing to pair datasets with paired inputs, get usable read #s if not done yet, and define replicate structure
open(F,$filist_fi);
for my $line (<F>) {
    chomp($line);
    $line =~ s/\r//g;
    my @tmp = split(/\s+/,$line);
    next unless ($tmp[0]);
    next if ($tmp[0] eq "uID");

## IMPORTANT: this script is written for two manifest file structures: 
## 1) a 5 column structure where column 4 is CLIP and column 5 is SMInput
## 2) a 6 column structure where column 4 is CLIP rep1, column 5 is CLIP rep2, and column 6 is a single SMInput used for both
## This will not run properly if the manifest is structured differently

    if (scalar(@tmp) == 6) {
	$type_flag = "two_replicate_ENCODEstyle";
    } elsif (scalar(@tmp) == 5) {
	$type_flag = "one_replicate";
    }

    my $uid = shift(@tmp);
    my $rbp = shift(@tmp);
    my $cellline = shift(@tmp);
    my %CLIP;
    $CLIP{"_01"} = shift(@tmp);
    $CLIP{"_02"} = shift(@tmp) if ($type_flag eq "two_replicate_ENCODEstyle");
    my $input = shift(@tmp);

    $CLIP{"_01"} =~ s/\.bam$//;
    $CLIP{"_02"} =~ s/\.bam$// if ($type_flag eq "two_replicate_ENCODEstyle");
    $input =~ s/\.bam$//;

    if ($type_flag eq "two_replicate_ENCODEstyle") {
	@rep_listing = ("_01","_02");
    } elsif ($type_flag eq "one_replicate") {
	@rep_listing = ("_01");
    } else {
	print STDERR "TYPE flag is not set properly!!!!\n";
    }

    next if ($CLIP{"_01"} eq "NA" || $input eq "NA");
    if ($type_flag eq "two_replicate_ENCODEstyle") {
	next if ($CLIP{"_02"} eq "NA");
    }
    for my $rep (@rep_listing) {

	my @clip_fi_split = split(/\//,$CLIP{$rep});
	my $clip_fi_short = $clip_fi_split[$#clip_fi_split];

	my $clip_bam_fi_short = $clip_fi_short.".bam";
	my $clip_bam_fi_softlink = $output_folder.$clip_bam_fi_short;
	my $clip_bam_fi = $CLIP{$rep}.".bam";
	push @clip_bam_fi_shorts,$clip_bam_fi_short;
	
	my @input_fi_split = split(/\//,$input);
	my $input_fi_short = $input_fi_split[$#input_fi_split];

	my $input_bam_fi_short = $input_fi_short.".bam";
	my $input_bam_fi_softlink = $output_folder.$input_bam_fi_short;
	my $input_bam_fi = $input.".bam";
	
	my $peak_fi_short = $clip_fi_short.".peaks.bed";
	my $peak_fi_softlink = $output_folder.$peak_fi_short;
	my $peak_fi = $CLIP{$rep}.".peaks.bed";

	$peakfi2uid{$peak_fi} = $uid.$rep;
	push @peak_filelist,$peak_fi;
	push @peak_filelist_shorts,$peak_fi_short;
	$clipshort_to_peakshort{$clip_bam_fi_short} = $peak_fi_short;
	
	my $output = $output_folder.$uid.$rep.".basedon_".$uid.$rep.".peaks.l2inputnormnew.bed";
	unless (-e $clip_bam_fi && -e $input_bam_fi) {
	    print STDERR "ERROR ERROR one of these files doesn't exist $clip_bam_fi $input_bam_fi\n";
	    next;
	}
	system("ln -s $clip_bam_fi $clip_bam_fi_softlink") unless (-e $clip_bam_fi_softlink);
	system("ln -s $input_bam_fi $input_bam_fi_softlink") unless (-e $input_bam_fi_softlink);
	system("ln -s $peak_fi $peak_fi_softlink") unless (-e $peak_fi_softlink);

	unless (exists $mapped_num{$clip_bam_fi_short}) {
	    my $mapped_read_num = `samtools view -c -F 4 $clip_bam_fi_softlink`;
	    $mapped_num{$clip_bam_fi_short} = $mapped_read_num;
	}
	unless (exists $mapped_num{$input_bam_fi_short}) {
	    my $input_read_num = `samtools view -c -F 4 $input_bam_fi_softlink`;
	    $mapped_num{$input_bam_fi_short} = $input_read_num;
	}

	if (-e $output) {
	    print STDERR "skipping $output alreadydone\n";
	    next;
	}
    }
}
close(F);


open(READNUMOUT,">$filist_fi_mappedreadnums");
for my $fi (keys %mapped_num) {
    chomp($mapped_num{$fi});
    print READNUMOUT "$fi\t$mapped_num{$fi}\n";
}
close(READNUMOUT);

my %foldenriched;
open(F,$filist_fi);
for my $line (<F>) {
    chomp($line);
    $line =~ s/\r//g;
    my @tmp = split(/\s+/,$line);
    next if ($tmp[0] eq "uID");


    my $uid = shift(@tmp);
    my $rbp = shift(@tmp);
    my $cellline = shift(@tmp);
    my %CLIP;
    $CLIP{"_01"} = shift(@tmp);
    $CLIP{"_02"} = shift(@tmp) if ($type_flag eq "two_replicate_ENCODEstyle");
    my $input = shift(@tmp);

    $CLIP{"_01"} =~ s/\.bam$//;
    $CLIP{"_02"} =~ s/\.bam$// if ($type_flag eq "two_replicate_ENCODEstyle");
    $input =~ s/\.bam$//;

    
    my $next_flag = 0;
    for my $rep (@rep_listing) {
	$next_flag = 1 if ($CLIP{$rep} eq "NA");
    } 
    $next_flag = 1 if ($input eq "NA");
    next if ($next_flag == 1);
    
    for my $rep (@rep_listing) {

        my @clip_fi_split = split(/\//,$CLIP{$rep});
        my $clip_fi_short = $clip_fi_split[$#clip_fi_split];

        my $clip_bam_fi_short = $clip_fi_short.".bam";
        my $clip_bam_fi_softlink = $output_folder.$clip_bam_fi_short;
        my $clip_bam_fi = $CLIP{$rep}.".bam";

	my @input_fi_split = split(/\//,$input);
        my $input_fi_short = $input_fi_split[$#input_fi_split];

        my $input_bam_fi_short = $input_fi_short.".bam";
        my $input_bam_fi_softlink = $output_folder.$input_bam_fi_short;
	my $input_bam_fi = $input.".bam";

        my $peak_fi_short = $clip_fi_short.".peaks.bed";
        my $peak_fi_softlink = $output_folder.$peak_fi_short;
	my $peak_fi = $CLIP{$rep}.".peaks.bed";

	unless (-e $peak_fi && -e $input_bam_fi && -e $clip_bam_fi) {
	    print STDERR "ERROR ERROR ERROR one of these doesn't exist $peak_fi $clip_bam_fi $input_bam_fi\n";
	    next;
	}

	my @desired_peak_filelist;
	$desired_peak_filelist[0] = $peak_fi;


	for my $compare_peakfi (@desired_peak_filelist) {
	    
	    my @compare_peakfii = split(/\//,$compare_peakfi);
	    my $compare_peakfi_short = $compare_peakfii[$#compare_peakfii];
	    my $compare_peakfi_softlink = $output_folder.$compare_peakfi_short;
	    my $output = $output_folder.$uid.$rep.".basedon_".$peakfi2uid{$compare_peakfi}.".peaks.l2inputnormnew.bed";
	    my $output_compressed = $output.".compressed.bed";
	    
	    ## For faster re-running - will pick up where it left off and not overwrite files. Output to a new output_directory if full re-run is desired
	    if (-e $output) {
		print STDERR "skipping $output alreadydone\n";
		if (-e $output_compressed) {
		} else {
		    system("perl compress_l2foldenrpeakfi.pl $output");
		}
	    } else {
		system("perl overlap_peakfi_with_bam_PE.pl $clip_bam_fi_softlink $input_bam_fi_softlink $compare_peakfi_softlink $filist_fi_mappedreadnums $output");
		system("perl compress_l2foldenrpeakfi.pl $output");
	    }

	}
    }
}
close(F);

