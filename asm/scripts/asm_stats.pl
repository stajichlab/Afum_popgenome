#!/usr/bin/perl
#
use File::Spec;
use strict;
use warnings;

my %stats;

my $readlen = 150;    # assume reads are 150bp?

my $read_map_stat = 'mapping_report';
my $dir           = shift || 'genomes';
my %cols;
my @header;
my %header_seen;

opendir( DIR, $dir ) || die $!;
my $first = 1;
foreach my $file ( readdir(DIR) ) {
    next unless ( $file =~ /(\S+)\.stats.txt$/ );
    my $stem = $1;
    $stem =~ s/\.sorted//;
    open( my $fh => "$dir/$file" ) || die $!;
    while (<$fh>) {
        next if /^\s+$/;
        s/^\s+//;
        chomp;
        if (/(.+\S+)\s+=\s+(\d+(\.\d+)?)/) {
            $stats{$stem}->{$1} = $2;
            #warn("'$1' '$2'\n");
            $cols{$1}++;
            if ( !exists $header_seen{$1} ) {
                push @header, $1;
                $header_seen{$1} = 1;
            }
        }
   }

    if ($first) {
        push @header, qw(BUSCO_Complete BUSCO_Single BUSCO_Single BUSCO_Single
          BUSCO_Fragmented BUSCO_Missing BUSCO_NumGenes
        );
    }

    my $busco_file = File::Spec->catfile(
        "BUSCO", $stem,
        sprintf( "run_%s",               $stem ),
        sprintf( "short_summary_%s.txt", $stem )
    );

    if ( -f $busco_file ) {

        open( my $fh => $busco_file ) || die $!;
        while (<$fh>) {
            if (
/^\s+C:(\d+\.\d+)\%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)/
              )
            {
                $stats{$stem}->{"BUSCO_Complete"}   = $1;
                $stats{$stem}->{"BUSCO_Single"}     = $2;
                $stats{$stem}->{"BUSCO_Duplicate"}  = $3;
                $stats{$stem}->{"BUSCO_Fragmented"} = $4;
                $stats{$stem}->{"BUSCO_Missing"}    = $5;
                $stats{$stem}->{"BUSCO_NumGenes"}   = $6;
            }
        }

    }
    else {
        warn("Cannot find $busco_file");
    }

    my $sumstatfile = File::Spec->catfile( $read_map_stat,
        sprintf( "%s.bbmap_summary.txt", $stem ) );
    if ( -f $sumstatfile ) {
        open( my $fh => $sumstatfile ) || die "Cannot open $sumstatfile: $!";
        my $read_dir   = 0;
        my $base_count = 0;
        $stats{$stem}->{'Mapped reads'} = 0;
        while (<$fh>) {
            if (/Read (\d+) data:/) {
                $read_dir = $1;
            }
            elsif ( $read_dir && /^mapped:\s+(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/ ) {
                $base_count += $4;
                $stats{$stem}->{'Mapped_reads'} += $2;
            }
            elsif (/^Reads Used:\s+(\S+)/) {
                $stats{$stem}->{'Reads'} = $1;
            }

        }
        if ( ! exists $stats{$stem}->{'TOTAL LENGTH'} ) {
            warn("no total length for $stem");
        }
        $stats{$stem}->{'Average_Coverage'} =
          sprintf( "%.1f", $base_count / $stats{$stem}->{'TOTAL LENGTH'} );
          $stats{$stem}->{'Mapping_Percent'} =
            sprintf( "%.2f", 100 * $stats{$stem}->{'Mapped_reads'}/ $stats{$stem}->{'Reads'} );
    }
    if ($first) {
        push @header, ( 'Reads', 'Mapped_reads', 'Mapping_Percent','Average_Coverage' );
    }
    $first = 0;
}

print join( "\t", qw(SampleID), @header ), "\n";
foreach my $sp ( sort keys %stats ) {
    print join( "\t", $sp, map { $stats{$sp}->{$_} || 'NA' } @header ), "\n";
}
