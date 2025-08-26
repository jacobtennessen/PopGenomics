#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_a $opt_b $opt_o $opt_t $opt_n $opt_y $opt_w $opt_m $opt_r $opt_g $opt_v $opt_x $opt_d $opt_s);

# Usage
my $usage = "
HudsonFstAverageWindows.pl

Calculates pairwise Fst (Hudson) in windows

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

Usage: perl HudsonFstAverageWindows.pl options
 required:
  -f  (path to) an allele table (numerical genotypes 0, 1, or 2)
  -a  comma-delimited list of sample freq positions [samples start at 1], first group
  -b  comma-delimited list of sample freq positions [samples start at 1], second group
  -o  (path to) output file
 optional:
  -g  comma-delimited list of chromosome(s) to examine
  -r  range of positions to examine on chromosomes START-END
  -m  minimum joint frequency [default = none]
  -t  minimum Fst to report [default = none]
  -n  minimum number of samples per group [default = 5]
  -y  comma-delimited list of male positions [samples start at 1]
  -w  window size [default = 100000]
  -v  comma-delimited list of regions to avoid [CHROM:Site-Site,CHROM:Site-Site]
  -x  exit once target range is seen
  -d  only consider sites segregating in both populations
  -s  minimal count of SNPs per window [default = 100]

";

#############

# command line processing.
getopts('f:a:b:o:t:n:y:w:m:g:r:v:s:xd');
die $usage unless ($opt_f);
die $usage unless ($opt_a);
die $usage unless ($opt_b);
die $usage unless ($opt_o);

my ($genotypes, $grouppos1, $grouppos2, $outfile, $fstthresh, $freqthresh, $minsamp, @males, $window, @chromstouse, @range, @avoid, $exit, $double, $minsnps);

$genotypes	= $opt_f if $opt_f;
$grouppos1 = $opt_a if $opt_a;
$grouppos2 = $opt_b if $opt_b;
$outfile = $opt_o if $opt_o;

if (defined $opt_t) {
  $fstthresh = $opt_t;
}

if (defined $opt_s) {
  $minsnps = $opt_s;
} else {
  $minsnps = 100;
}

if (defined $opt_m) {
  $freqthresh = $opt_m;
}

if (defined $opt_n) {
  $minsamp = $opt_n;
} else {
  $minsamp = 5;
}

if (defined $opt_y) {
  @males = split ",", $opt_y;
}

if (defined $opt_g) {
  @chromstouse = split ",", $opt_g;
}

if (defined $opt_v) {
  @avoid = split ",", $opt_v;
}

if (defined $opt_x) {
  $exit = 1;
}

if (defined $opt_d) {
  $double = 1;
}

if (defined $opt_r) {
  @range = split "-", $opt_r;
} else {
  @range = (0,1000000000000);
}

if (defined $opt_w) {
  $window = $opt_w;
} else {
  $window = 100000;
}

my %males;

my %chromstouse;

my %HoHavoidstart;

my %HoHavoidend;

foreach my $c (@chromstouse) {
  $chromstouse{$c} = 1;
}

foreach my $av (@avoid) {
  my @chromdata = split ":", $av;
  my @sites = split "-", $chromdata[1];
  $HoHavoidstart{$chromdata[0]}{$chromdata[1]} = $sites[0];
  $HoHavoidend{$chromdata[0]}{$chromdata[1]} = $sites[1];
}

foreach my $m (@males) {
  $males{$m} = 1;
}

my @g1f = split ",", $grouppos1;

my @g2f = split ",", $grouppos2;

my @out;

print "\n\nExamining Populations\n$grouppos1\n$grouppos2\n";

if (defined $chromstouse[0]) {
  my $chromlist = join " ", @chromstouse;
  print "Chroms $chromlist\n\n";
} else {
  print "Full genome\n\n";
}

my $lastchrom = "N";

my $lastsite = 0;

unless ( open(OUT, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OUT "Chrom\tWindow\tNum\tDen\tFst\tSNPs";
close (OUT);

my $sitecount = 0;
my $numtotal = 0;
my $dentotal = 0;

my @fsts;

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[0] =~ /Site/) {
      next;
    }
    my @allchromdata = split "_", $data[0];
    my $site = pop @allchromdata;
    my $chrom = join "_", @allchromdata;
    my $forbidden;
    unless ((!defined $chromstouse[0])||((defined $chromstouse{$chrom})&&($site >= $range[0])&&($site <= $range[1]))) {
      if ((defined $exit)&&($sitecount > 0)) {
        last;
      }
      next;
    }
    foreach my $avr (keys %{ $HoHavoidstart{$chrom} } ) {
      if (($site >= $HoHavoidstart{$chrom}{$avr})&&($site <= $HoHavoidend{$chrom}{$avr})) {
        $forbidden = 1;
      }
    }
    if (defined $forbidden) {
      next;
    }
    my $roundsite = $window*int($site/$window);
    unless (($chrom =~ /$lastchrom/)&&($roundsite == $lastsite)) {
      if ($sitecount >= $minsnps) {
        my $numav = $numtotal/$sitecount;
        my $denav = $dentotal/$sitecount;
        my $fstav = 0;
        unless ($denav == 0) {
          $fstav = $numav/$denav;
        }
        push @fsts, $fstav;
        push @out, "$lastchrom\t$lastsite\t$numav\t$denav\t$fstav\t$sitecount";
      }
      $lastchrom = $chrom;
      $lastsite = $site;
      $sitecount = 0;
    }
    my $group1alleletotal = 0;
    my $good1samples = 0;
    foreach my $g1 (@g1f) {
      if ($data[$g1] =~ /\d/) {
        if ((defined $males[0])&&(defined $males{$g1})&&($data[0] =~ /X/)) {
          if ($data[$g1] > 0) {
            $group1alleletotal += 1;
          }
          $good1samples += 0.5;
        } else {
          $group1alleletotal += $data[$g1];
          $good1samples +=1;
        }
      }
    }
    my $group2alleletotal = 0;
    my $good2samples = 0;
    foreach my $g2 (@g2f) {
      if ($data[$g2] =~ /\d/) {
        if ((defined $males[0])&&(defined $males{$g2})&&($data[0] =~ /X/)) {
          if ($data[$g2] > 0) {
            $group2alleletotal += 1;
          }
          $good2samples += 0.5;
        } else {
          $group2alleletotal += $data[$g2];
          $good2samples +=1;
        }
      }
    }
    if (($good1samples >= $minsamp)&&($good2samples >= $minsamp)) {
      my $g1size = $good1samples*2;
      my $g2size = $good2samples*2;
      my $group1freq = $group1alleletotal/$g1size;
      my $group2freq = $group2alleletotal/$g2size;
      my $meanfreq = ($group1freq + $group2freq)/2;
      if ($meanfreq > 0.5) {
        $meanfreq = 1 - $meanfreq;
      }
      if ($meanfreq > 0) {
        unless ((defined $double)&&(($group1freq == 0)||($group1freq == 1)||($group2freq == 2)||($group2freq == 1))) {
          if ($data[0] =~ /Mt_/) {
            $g1size = $good1samples;
            $g2size = $good2samples;
          }
          my $fstnum = ($group1freq-$group2freq)*($group1freq-$group2freq) - $group1freq*(1-$group1freq)/($g1size-1) - $group2freq*(1-$group2freq)/($g2size-1);
          my $fstden = $group1freq*(1-$group2freq) + $group2freq*(1-$group1freq);
          $sitecount +=1;
          $numtotal += $fstnum;
          $dentotal += $fstden;
        }
      }
    }
}

close (IN);

if ($sitecount >= $minsnps) {
  my $numav = $numtotal/$sitecount;
  my $denav = $dentotal/$sitecount;
  my $fstav = 0;
  unless ($denav == 0) {
    $fstav = $numav/$denav;
  }
  push @fsts, $fstav;
  push @out, "$lastchrom\t$lastsite\t$numav\t$denav\t$fstav\t$sitecount";
}

my $windowcount = scalar(@fsts);

my $mean = get_avg(\@fsts);
my $stddev = get_stddev(\@fsts);

print "In $windowcount windows, Fst mean = $mean and standard deviation = $stddev\n";

my $result = join "\n", @out;

unless ( open(OUT, ">>$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OUT "\n$result";
close (OUT);

#######################

sub get_stddev {
  return sqrt(get_disp(@_));
}

sub get_disp {
  my ($array_ref) = @_;
  my $mean = get_avg($array_ref);
  my $count = @$array_ref;
  my $sum = 0;
  
  for my $num (@$array_ref) {
      $sum += (($num - $mean) ** 2);
  }
  return $sum / $count;
}

sub get_avg {
  my ($array_ref) = @_;
  my $count = @$array_ref;
  my $sum = 0;
  for my $num (@$array_ref) {
      $sum += $num;
  }
  return $sum / $count;
}