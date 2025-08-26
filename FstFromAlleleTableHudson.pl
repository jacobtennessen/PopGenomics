#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_a $opt_b $opt_o $opt_t $opt_n $opt_y $opt_s $opt_m $opt_c $opt_p $opt_x);

# Usage
my $usage = "
FstFromAlleleTableHudson.pl

Calculates pairwise Fst (Hudson)

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

Usage: perl FstFromAlleleTableHudson.pl options
 required:
  -f  (path to) an allele table (numerical genotypes 0, 1, or 2)
  -a  comma-delimited list of sample freq positions [samples start at 1], first group
  -b  comma-delimited list of sample freq positions [samples start at 1], second group
  -o  (path to) output file
 optional:
  -c  table with chrom and position for each contig
  -m  minimum joint frequency [default = none]
  -t  minimum Fst to report [default = none]
  -n  minimum number of samples per group [default = 5]
  -y  comma-delimited list of male positions [samples start at 1]
  -s  spacer size [default = 100000]
  -p  name of haploid chromosome
  -x  name of X chromosome

";

#############

# command line processing.
getopts('f:a:b:o:t:n:y:s:m:c:p:x:');
die $usage unless ($opt_f);
die $usage unless ($opt_a);
die $usage unless ($opt_b);
die $usage unless ($opt_o);

my ($genotypes, $grouppos1, $grouppos2, $outfile, $fstthresh, $freqthresh, $minsamp, @males, $spacer, $postable, $haploid, $xname);

$genotypes	= $opt_f if $opt_f;
$grouppos1 = $opt_a if $opt_a;
$grouppos2 = $opt_b if $opt_b;
$outfile = $opt_o if $opt_o;

if (defined $opt_t) {
  $fstthresh = $opt_t;
}

if (defined $opt_c) {
  $postable = $opt_c;
}

if (defined $opt_p) {
  $haploid = $opt_p;
} else {
  $haploid = "ignore";
}

if (defined $opt_x) {
  $xname= $opt_x;
} else {
  $xname = "ignore";
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

if (defined $opt_s) {
  $spacer = $opt_s;
} else {
  $spacer = 100000;
}

my %males;

foreach my $m (@males) {
  $males{$m} = 1;
}

my $cumul = 0 - $spacer;

my @g1f = split ",", $grouppos1;

my @g2f = split ",", $grouppos2;

my @out;

my %contigpos;

if (defined $postable) {

  open(PT, $postable) || die "can't open $postable\n";
  
  while (<PT>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[2] =~ /\d/) {
      $contigpos{$data[0]} = $data[1];
    }
  }
  
  close (PT);

}

my $lastchrom = "N";

my $lastsite = 0;

unless ( open(OUT, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OUT "Marker\tChrom\tSite\tPos\tSize1\tSize2\tFreq1\tFreq2\tFst";
close (OUT);

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[0] =~ /Site/) {
      next;
    }
    my $group1alleletotal = 0;
    my $good1samples = 0;
    foreach my $g1 (@g1f) {
      if ($data[$g1] =~ /\d/) {
        if (($data[0] =~ /$haploid/)||((defined $males[0])&&(defined $males{$g1})&&($data[0] =~ /$xname/))) {
          if ($data[$g1] > 0) {
            $group1alleletotal += 1;
          }
          $good1samples += 0.5;
        } else {
          if ($data[$g1] =~ /l/) {
            if ($data[$g1] =~ /2l/) {
              $group1alleletotal += 1;
            }
            $good1samples +=0.5;
          } else {
            $group1alleletotal += $data[$g1];
            $good1samples +=1;
          }
        }
      }
    }
    my $group2alleletotal = 0;
    my $good2samples = 0;
    foreach my $g2 (@g2f) {
      if ($data[$g2] =~ /\d/) {
        if (($data[0] =~ /$haploid/)||((defined $males[0])&&(defined $males{$g2})&&($data[0] =~ /$xname/))) {
          if ($data[$g2] > 0) {
            $group2alleletotal += 1;
          }
          $good2samples += 0.5;
        } else {
          if ($data[$g2] =~ /l/) {
            if ($data[$g2] =~ /2l/) {
              $group1alleletotal += 1;
            }
            $good2samples +=0.5;
          } else {
            $group2alleletotal += $data[$g2];
            $good2samples +=1;
          }
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
      unless ((defined $freqthresh)&&($meanfreq < $freqthresh)) {
        my $fstnum = ($group1freq-$group2freq)*($group1freq-$group2freq) - $group1freq*(1-$group1freq)/($g1size-1) - $group2freq*(1-$group2freq)/($g2size-1);
        my $fstden = $group1freq*(1-$group2freq) + $group2freq*(1-$group1freq);
        my $fstat = sprintf "%.6f", $fstnum/$fstden;
        unless ((defined $fstthresh)&&($fstat < $fstthresh)) {
            my @allchromdata = split "_", $data[0];
            my $site = pop @allchromdata;
            my $chrom = join "_", @allchromdata;
            $group1freq = sprintf "%.6f", $group1freq;
            $group2freq = sprintf "%.6f", $group2freq;
            unless ($chrom =~ /^$lastchrom$/) {
              $cumul += $lastsite + $spacer;
            }
            my $pos = $cumul+$site;
            if (defined $postable) {
              $pos = $contigpos{$chrom}+$site;
            }
            push @out, "$data[0]\t$chrom\t$site\t$pos\t$g1size\t$g2size\t$group1freq\t$group2freq\t$fstat";
            $lastchrom = $chrom;
            $lastsite = $site;
        }
      }
    }
    if ((scalar(@out)) >= 1000) {
        my $result = join "\n", @out;
        unless ( open(OUT, ">>$outfile") ) {
            print "Cannot open file \"$outfile\" to write to!!\n\n";
            exit;
        }
        print OUT "\n$result";
        close (OUT);
        @out = ();
    }
}

close (IN);

if (defined $out[0]) {
    my $result = join "\n", @out;
    unless ( open(OUT, ">>$outfile") ) {
        print "Cannot open file \"$outfile\" to write to!!\n\n";
        exit;
    }
    print OUT "\n$result";
    close (OUT);
}