#!/usr/bin/env perl
#
use strict;
use warnings;

my %samps;
my $in = shift || '../variantcall/SRA_samples.csv';
open(my $fh => $in) || die "cannot open $in: $!";
while(<$fh>){ 
	chomp;
	my ($sra,$strain,$biosample,@rest) = split(/,/,$_);
	push @{$samps{$biosample}->{$strain}}, $sra;
}

for my $samp ( keys %samps ) {
    while( my ($strain, $sralst) = each %{$samps{$samp}} ) {
	my @sra = (sort @$sralst);
	next if ( @sra <= 1 ) ;
	my $strainnospace = $strain;
	$strainnospace =~ s/\s+/_/g;

	foreach my $strand ( qw(R1 R2) ) {
	    if (  -f sprintf("input/%s_%s.fq.gz",$sra[0],$strand) ) {

		printf("if [ ! -f %s_%s.fq.gz ]; then\n",$strainnospace,$strand);
		printf("   zcat %s | gzip -c > %s_%s.fq.gz\n",
		       join(" ", map { sprintf("%s_%s.fq.gz", $_, $strand) } @sra),
		       $strainnospace, $strand);
		print "fi\n";
		printf "rm %s\n",join(" ", 
				      map { sprintf("%s_%s.fq.gz",$_,$strand) } @sra);
	    }
	}
    }
}
