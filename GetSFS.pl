#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_v $opt_o $opt_u $opt_w $opt_l $opt_e $opt_m);

# Usage
my $usage = "
GetSFS.pl

Calculates site frequency spectrum

Copyright (C) 2025 by Jacob A Tennessen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl GetSFS.pl options
 required:
  -v	(path to) allele counts table
  -o  (path to) output file
 optional:
  -u maximum ratio between het expected and het observed [default = 1.5]
  -w maximum difference between het expected and het observed [default = 5]
  -l comma-delimited list of samples to include [samples start at 1]
  -e [comma-delimited list of] regions to exclude [format CHROM:SITE-SITE,CHROM:SITE-SITE]
  -m max missing samples [default = no limit]
";

#############

getopts('v:o:u:w:l:e:m:');
die $usage unless ($opt_v);
die $usage unless ($opt_o);

my ($alleles, $outfile, $maxratio, $maxdif, @list, @exclude, $maxmiss);

$alleles = $opt_v if $opt_v;

$outfile = $opt_o if $opt_o;

if (defined $opt_l) {
  @list = split ",", $opt_l;
}

if (defined $opt_u) {
  $maxratio = $opt_u;
} else {
  $maxratio = 1.5;
}

if (defined $opt_m) {
  $maxmiss = $opt_m;
}

if (defined $opt_e) {
  @exclude = split ",", $opt_e;
}

if (defined $opt_w) {
  $maxdif = $opt_w;
} else {
  $maxdif = 5;
}

my %HoHexcludestart;

my %HoHexcludeend;

foreach my $e (@exclude) {
  my @chromsite = split ":", $e;
  my @sites = split "-", $chromsite[1];
  $HoHexcludestart{$chromsite[0]}{$chromsite[1]} = $sites[0];
  $HoHexcludeend{$chromsite[0]}{$chromsite[1]} = $sites[1];
}

open(IN, $alleles) || die "can't open $alleles\n";

my $gametes;

my %sfs;

my $allmin = 1;

while (<IN>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @origdata = split "\t", $line;
  my @realdata;
  if (defined $list[0]) {
    push @realdata, $origdata[0];
    push @realdata, @origdata[@list];
  } else {
    push @realdata, @origdata;
  }
  if ($realdata[0] =~ /Site/) {
    $gametes = 2*(scalar(@realdata) - 1);
    if (defined $maxmiss) {
      $allmin = $gametes/2 - $maxmiss;
    }
  } else {
    my @complexsite = split "_", $realdata[0];
    my $site = pop @complexsite;
    my $chrom = join "_", @complexsite;
    my $bannedregion;
    foreach my $excluderegion ( keys %{ $HoHexcludestart{$chrom} } ) {
      if (($site >= $HoHexcludestart{$chrom}{$excluderegion})&&($site <= $HoHexcludeend{$chrom}{$excluderegion})) {
        $bannedregion = 1;
      }
    }
    my $altcount = 0;
    my $gametecount = 0;
    my $hetcount = 0;
    for (my $i = 1; $i < (scalar(@realdata)); $i++) {
      if ($realdata[$i] =~ /\d/) {
        if ($realdata[$i] == 0) {
          $gametecount +=2;
        } elsif ($realdata[$i] == 1) {
          $altcount +=1;
          $gametecount +=2;
          $hetcount += 1;
        } elsif ($realdata[$i] == 2) {
          $altcount +=2;
          $gametecount +=2;
        }
      }
    }
    if ($gametecount >= 2*$allmin) {
      my $freq = $altcount/$gametecount;
      my $het = 2*($freq)*(1-$freq);
      my $expectedhet = $het*($gametecount/2);
      unless (($freq == 0)||($freq == 1)) {
        my $ratio = $hetcount/$expectedhet;
        my $difference = abs($hetcount - $expectedhet);
        if ((($ratio < $maxratio)&&($ratio > (1/$maxratio)))||($difference <= $maxdif)) {
          my $folded = $freq;
          if ($folded > 0.5) {
            $folded = 1 - $freq;
          }
          my $estimatedcount = sprintf "%.0f", $folded*$gametes;
          if (defined $sfs{$estimatedcount}) {
            $sfs{$estimatedcount} += 1;
          } else {
            $sfs{$estimatedcount} = 1;
          }
        }
      }
    }
  }
}

close (IN);

my @out;

for (my $c = 1; $c <= $gametes/2; $c++) {
  unless (defined $sfs{$c}) {
    $sfs{$c} = 0;
  }
  push @out, "$c\t$sfs{$c}";
}

my $result = join "\n", @out;

unless ( open(META, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print META "Freq\tCount\n$result";
close (META);

###################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}
 