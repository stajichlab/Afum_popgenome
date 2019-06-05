#!/usr/bin/env perl
use strict;
use warnings;

my $file = shift || 'strain_samples.tsv';

open(my $fh => $file ) || die $!;
my %acc2strain;
while(<$fh>) {
    next if /^RunAcc/;
    chomp;
    my ($srr_acc, $strain) = split(/\t/,$_);
    $strain =~ s/[ \/]/_/g; 
    push @{$acc2strain{$strain}}, $srr_acc;
}

for my $str ( sort { scalar @{$acc2strain{$b}} <=> 
			 scalar @{$acc2strain{$a}} } keys %acc2strain ) {
    my $accs = $acc2strain{$str};
    if ( scalar @$accs > 1 ) {
	print join("\t", $str, $str,'Ascomycota'), "\n";
    } else {
	print join("\t", $accs->[0], $str,'Ascomycota'), "\n";
    }
}
