#!/usr/bin/env perl
#
use strict;
use warnings;

my %samps;

while(<>){ 
	chomp;
	my ($sra,$strain,$biosample,@rest) = split(/,/,$_);
	push @{$samps{$biosample}->{$strain}}, $sra;
}

for my $samp ( keys %samps ) {
	my $strain_sample = $samps{$samp};
	for my $strain ( keys %{$strain_sample}) {
		my @sra = sort values %{$strain_sample{$strain}};
		for $strand ( qw(R1 R2) ) {
			printf("zcat %s | gzip -c > %s_%s.fq.gz\n",
				join(" ", map { sprintf("%s_%s.fq.gz", $_, $strand) } @sra));
		}
	}
}
