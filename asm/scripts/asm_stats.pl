#!/usr/bin/perl
#
use File::Spec;
use strict;
use warnings;

my %stats;

my $dir = shift || 'genomes';
my $outfile = shift || "assembly_stats.tsv";
my $lookup_names = shift || "strain_samples.tsv";

my %lookups;
if ( -f $lookup_names )  {
    open(my $lfh => $lookup_names) || die "$lookup_names: $!";
    my $lhdrline = <$lfh>;
    my $n = 0;
    chomp($lhdrline);
    my %lhdr = map { $_ => $n++ } split(/\t/,$lhdrline);
    
    while(<$lfh>) {
	chomp;
	my @row = split(/\t/,$_);
	$lookups{ $row[$lhdr{RunAcc}] } = $row[ $lhdr{Strain} ];
    }
}
my $outfh;
if ( $outfile eq '-') { 
    $outfh = \*STDOUT;
} elsif ( -w $outfile ) {
    open($outfh => ">$outfile") || die $!;
} else {
    die "Cannot write to $outfile";
    
}

my %cols = map { $_ => 1 }  qw(z1_BUSCO_Complete z2_BUSCO_Single
    z3_BUSCO_Duplicate z4_BUSCO_Fragmented
    z5_BUSCO_Missing z6_BUSCO_NumGenes);

opendir(DIR,$dir) || die $!;

foreach my $file ( readdir(DIR) ) {
 next unless ( $file =~ /(\S+)\.stats.txt$/);
 my $stem = $1;
 $stem =~ s/\.sorted//;
 open(my $fh => "$dir/$file") || die $!;
 while(<$fh>) {
  next if /^\s+$/;
  s/^\s+//;
  chomp;
  if ( /(.+)\s+=\s+(\d+(\.\d+)?)/ ) {
      $stats{$stem}->{$1} = $2;
#      warn($1," ", $2,"\n");
      $cols{$1}++;
  }
 }

 my $busco_file = File::Spec->catfile("BUSCO",sprintf("run_%s",$stem),
				      sprintf("short_summary_%s.txt",$stem));
				      
 if ( -f $busco_file ) {
     
     open(my $fh => $busco_file) || die $!;
     while(<$fh>) {	 
	 if (/^\s+C:(\d+\.\d+)\%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)/ ) {
	     $stats{$stem}->{"z1_BUSCO_Complete"} = $1;
	     $stats{$stem}->{"z2_BUSCO_Single"} = $2;
	     $stats{$stem}->{"z3_BUSCO_Duplicate"} = $3;
	     $stats{$stem}->{"z4_BUSCO_Fragmented"} = $4;
	     $stats{$stem}->{"z5_BUSCO_Missing"} = $5;
	     $stats{$stem}->{"z6_BUSCO_NumGenes"} = $6;
	 } 
     }

 } else {
    warn("Cannot find $busco_file");
 }
}


my @cols = sort keys %cols;
my @cols2 = sort keys %cols;
for my $m ( @cols2 ) {
	$m =~ s/^z\d+_//;
}
print $outfh join("\t", qw(SampleID Strain), @cols2), "\n";
foreach my $sp ( sort keys %stats ) {    
    print $outfh join("\t", $sp, $lookups{$sp} || $sp,
			 map { $stats{$sp}->{$_} || '-' } @cols), "\n";
}
