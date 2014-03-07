use strict;
use warnings;

#read GrOUPED, drops reads into genes and convert coordinates into BED
my %GENES;
my $spec = $ARGV[0]; #hg18, mm8, ce6
my $file = $ARGV[1]; #"AGO2_293T_comb";
#open(GRP,"/nas3/yeolab/Genome/ucsc_output/hg18data4/grouped") || die " can't open file: $!";
open(GRP,"/nas3/yeolab/Genome/ensembl/AS_STRUCTURE/hg19data4/grouped") || die "can't open file: $!";
while (<GRP>)
{
   chomp;
   my ($chr,$genes) = split(/\t/,$_);
   my @genes = split(/\,/,$genes);
   my $first = $genes[0];
   my $last = $genes[scalar(@genes)-1];
   my ($gene_f,$x_f,$y_f,$strand_f) = split(/\|/,$first);
   my ($gene_l,$x_l,$y_l,$strand_l) = split(/\|/,$last);
   $GENES{$chr}{$x_f."-".$y_l} = $genes; #gene start and stop
}
close(GRP);

## convert to bowtie output but relative to genes to  be compatible with rest of scripts

my $counts_sense_antisense = $file.".anti";

open(OUT,">".$file.".notingenes.BED") || die;

printf OUT "track name\=".$file."_genomic description\=".$file."_genomic visibility\=2 itemRgb\=\"On\" useScore\=1\n";

open(BED,">".$file.".ingenes.BED") || die;

printf BED "track name\=".$file."_genic description\=".$file."_genic visibility\=2 itemRgb\=\"On\" useScore\=1\n";

my %target_counts;
open(FI,"<".$file.".unique.rmsk") || die;
while (<FI>)
{
   chomp;
   my $line = $_;
   my ($qname,$qsign,$chr,$x,$y,$seq) = split(/\t/,$line);
   my $qstrand = -1;
   if ($qsign eq '+')
   {
      $qstrand = 1;
   }
   $chr =~ s/chr//g;
   my $flag = -1;
   foreach my $grouped_loc (keys %{$GENES{$chr}})
   {
      my ($x_f,$y_l)= split(/\-/,$grouped_loc);
      if ($x_f <= $x && $y <= $y_l) 
      {
	 foreach my $gene_locs ($GENES{$chr}{$grouped_loc}){
	    my @gene_locs = split(/\,/,$gene_locs);
	    foreach my $gene_loc (@gene_locs)
	    {
	       my ($geneid,$genex,$geney,$genestrand) = split(/\|/,$gene_loc);
	       if ($genex <= $x && $y <= $geney )
	       {
		  if ($genestrand eq $qstrand) {
		     $target_counts{$geneid}{'+'}++; 
		     printf BED "chr".$chr."\t".$x."\t".$y."\t".$qname.",".$gene_loc."\t1\t".$qsign."\t".$x."\t".$y."\t0,255,0\n";
		     $flag = 1;
		     last;
		  }		
		  else
		  {
		     $target_counts{$geneid}{'-'}++;
		  }
	       }
	    }	
	 }	
      }
   }
   if ($flag eq -1)
   {
      printf OUT "chr".$chr."\t".$x."\t".$y."\t".$qname.",NA\t1\t".$qsign."\t".$x."\t".$y."\t0,0,255\n";
   }
}
open(CTS,">".$counts_sense_antisense) || die;
foreach my $id (keys %target_counts)
{
   if (!$target_counts{$id}{'+'})
   {
      $target_counts{$id}{'+'} = 0;
   }
   if (!$target_counts{$id}{'-'})
   {
      $target_counts{$id}{'-'} = 0;
   }
   printf CTS $id."\t".$target_counts{$id}{'+'}."\t".$target_counts{$id}{'-'}."\n";
}
close(CTS);
undef %target_counts;
