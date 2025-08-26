#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_v $opt_o $opt_r $opt_a $opt_b $opt_c $opt_u $opt_w $opt_p $opt_l );

# Usage
my $usage = "
GetPi.pl

Calculates pi and heterozygosity

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

Usage: perl GetPi.pl options
 required:
  -v	(path to) allele counts table
  -o  (path to) outout file
optional:
  -r window size [default = 1000000]
  -a comma-delimted list of samples to be examned as a group [samples start at 1]
  -b comma-delimted list of samples to be examned as a group [samples start at 1]
  -c comma-delimted list of samples to be examned as a group [samples start at 1]
  -u maximum ratio between het expected and het observed [default = 1.5]
  -w maximum difference between het expected and het observed [default = 5]
  -p chrom to examine [defualt = all]
  -l range [START-END]
";

#############

getopts('v:o:r:a:b:c:u:w:p:l:');
die $usage unless ($opt_v);
die $usage unless ($opt_o);

my ($alleles, $outfile, $window, @group1, @group2, @group3, $maxratio, $maxdif, $chromtouse, @range);

$alleles = $opt_v if $opt_v;

$outfile = $opt_o if $opt_o;

if (defined $opt_r) {
  $window = $opt_r;
} else {
  $window = 1000000;
}

if (defined $opt_u) {
  $maxratio = $opt_u;
} else {
  $maxratio = 1.5;
}

if (defined $opt_w) {
  $maxdif = $opt_w;
} else {
  $maxdif = 5;
}

if (defined $opt_p) {
  $chromtouse = $opt_p;
}

if (defined $opt_l) {
  @range = split "-", $opt_l;
}

if (defined $opt_a) {
  @group1 = split ",", $opt_a;
}

if (defined $opt_b) {
  @group2 = split ",", $opt_b;
}

if (defined $opt_c) {
  @group3 = split ",", $opt_c;
}

my %HoHbakery;

my %HoHsnps;

my %HoHhw;

my %HoHoHbakery;

my %highest;

my %hets;

my %names;

my $sitecount = 0;

open(IN, $alleles) || die "can't open $alleles\n";

while (<IN>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ($data[0] =~ /Site/) {
    for (my $i = 1; $i < (scalar(@data)); $i++) {
      $names{$i} = $data[$i];
    }
  } else {
    my @chromsitedata = split "_", $data[0];    
    my $site = pop @chromsitedata;
    my $chrom = join "_", @chromsitedata;
    unless ((!defined $chromtouse)||($chrom =~ /^$chromtouse$/)) {
      next;
    }
    unless ((!defined $range[0])||(($range[0] <= $site)&&($range[1] >= $site))) {
      next;
    }
    my $roundedsite = (int($site/$window))*$window;
    unless ((defined $highest{$chrom})&&($highest{$chrom} > $roundedsite)) {
      $highest{$chrom} = $roundedsite;
    }
    my $altcount = 0;
    my $goodcount = 0;
    my $hetcount = 0;
    for (my $i = 1; $i < (scalar(@data)); $i++) {
      if ($data[$i] =~ /\d/) {
        if ($data[$i] == 0) {
          $goodcount +=1;
        } elsif ($data[$i] == 1) {
          $altcount +=1;
          $goodcount +=1;
          $hetcount += 1;
        } elsif ($data[$i] == 2) {
          $altcount +=2;
          $goodcount +=1;
        }
      }
    }
    if ($goodcount > 0) {
      my $freq = $altcount/(2*$goodcount);
      my $het = 2*($freq)*(1-$freq);
      my $expectedhet = $het*$goodcount;
      unless (($freq == 0)||($freq == 1)) {
        my $ratio = $hetcount/$expectedhet;
        my $difference = abs($hetcount - $expectedhet);
        if ((($ratio < $maxratio)&&($ratio > (1/$maxratio)))||($difference <= $maxdif)) {
          for (my $i = 1; $i < (scalar(@data)); $i++) {
            if ($data[$i] =~ /\d/) {
              if ($data[$i] == 1) {
                if (defined $hets{$i}) {
                  $hets{$i} += 1;
                } else {
                  $hets{$i} = 1;
                }
                if (defined $HoHoHbakery{$chrom}{$roundedsite}{$i}) {
                  $HoHoHbakery{$chrom}{$roundedsite}{$i} += 1;
                } else {
                  $HoHoHbakery{$chrom}{$roundedsite}{$i} = 1;
                }
              }
            }
          }
          if (defined $HoHsnps{$chrom}{$roundedsite}) {
            $HoHsnps{$chrom}{$roundedsite} += 1;
          } else {
            $HoHsnps{$chrom}{$roundedsite} = 1;
          }
          if (defined $HoHbakery{$chrom}{$roundedsite}) {
            $HoHbakery{$chrom}{$roundedsite} += $het;
          } else {
            $HoHbakery{$chrom}{$roundedsite} = $het;
          }
          if (defined $group1[0]) {
            my $goodcount1 = 0;
            my $altcount1 = 0;
            foreach my $g1 (@group1) {
              if ($data[$g1] =~ /\d/) {
                if ($data[$g1] == 0) {
                  $goodcount1 +=1;
                } if ($data[$g1] == 1) {
                  $altcount1 +=1;
                  $goodcount1 +=1;
                } elsif ($data[$g1] == 2) {
                  $altcount1 +=2;
                  $goodcount1 +=1;
                }
              }
            }
            if ($goodcount1 > 0) {
              my $freq1 = $altcount1/(2*$goodcount1);
              my $het1 = 2*($freq1)*(1-$freq1);
              if (defined $HoHoHbakery{$chrom}{$roundedsite}{"g1"}) {
                $HoHoHbakery{$chrom}{$roundedsite}{"g1"} += $het1;
              } else {
                $HoHoHbakery{$chrom}{$roundedsite}{"g1"} = $het1;
              }
            }
          }
          if (defined $group2[0]) {
            my $goodcount2 = 0;
            my $altcount2 = 0;
            foreach my $g2 (@group2) {
              if ($data[$g2] =~ /\d/) {
                if ($data[$g2] == 0) {
                  $goodcount2 +=1;
                } if ($data[$g2] == 1) {
                  $altcount2 +=1;
                  $goodcount2 +=1;
                } elsif ($data[$g2] == 2) {
                  $altcount2 +=2;
                  $goodcount2 +=1;
                }
              }
            }
            if ($goodcount2 > 0) {
              my $freq2 = $altcount2/(2*$goodcount2);
              my $het2 = 2*($freq2)*(1-$freq2);
              if (defined $HoHoHbakery{$chrom}{$roundedsite}{"g2"}) {
                $HoHoHbakery{$chrom}{$roundedsite}{"g2"} += $het2;
              } else {
                $HoHoHbakery{$chrom}{$roundedsite}{"g2"} = $het2;
              }
            }
          }
          if (defined $group3[0]) {
            my $goodcount3 = 0;
            my $altcount3 = 0;
            foreach my $g3 (@group3) {
              if ($data[$g3] =~ /\d/) {
                if ($data[$g3] == 0) {
                  $goodcount3 +=1;
                } if ($data[$g3] == 1) {
                  $altcount3 +=1;
                  $goodcount3 +=1;
                } elsif ($data[$g3] == 2) {
                  $altcount3 +=2;
                  $goodcount3 +=1;
                }
              }
            }
            if ($goodcount3 > 0) {
              my $freq3 = $altcount3/(2*$goodcount3);
              my $het3 = 2*($freq3)*(1-$freq3);
              if (defined $HoHoHbakery{$chrom}{$roundedsite}{"g3"}) {
                $HoHoHbakery{$chrom}{$roundedsite}{"g3"} += $het3;
              } else {
                $HoHoHbakery{$chrom}{$roundedsite}{"g3"} = $het3;
              }
            }
          }
          $sitecount +=1;
        } else {
          if (defined $HoHhw{$chrom}{$roundedsite}) {
            $HoHhw{$chrom}{$roundedsite} += 1;
          } else {
            $HoHhw{$chrom}{$roundedsite} = 1;
          }
        }
      }
    }
  }
}

close (IN);

my @names = sort by_number (keys (%names));

if ($sitecount > 0) {
  foreach my $i (@names) {
    unless (defined $hets{$i}) {
      $hets{$i} = 0;
    }
    my $het =  sprintf "%.8f", $hets{$i}/$sitecount;
    print "$names{$i}\t$het\n";
  }
}

my @bakeryout;

my @title = ("Chrom","Site","Pi");

foreach my $i (@names) {
  push @title, $names{$i};
}

if (defined $group1[0]) {
  push @title, "Group1";
}

if (defined $group2[0]) {
  push @title, "Group2";
}

if (defined $group3[0]) {
  push @title, "Group3";
}

push @title, "SNPs";

push @title, "HW";

my $title = join "\t", @title;

push @bakeryout, $title;

foreach my $chrom (sort (keys %highest)) {
  my $startw = 0;
  my $endw = $highest{$chrom};
  if (defined $range[0]) {
    $startw = (int($range[0]/$window))*$window;
    $endw = (int($range[1]/$window))*$window;
  }
  for (my $x = $startw; $x <= $endw; $x+= $window) {
    my @bakerytemp;
    if (defined $HoHbakery{$chrom}{$x}) {
      my $pi = sprintf "%.8f", $HoHbakery{$chrom}{$x}/$window;
      push @bakerytemp, "$chrom\t$x\t$pi";
    } else {
      push @bakerytemp, "$chrom\t$x\t0";
    }
    foreach my $i (@names) {
      if (defined $HoHoHbakery{$chrom}{$x}{$i}) {
        my $het = sprintf "%.8f", $HoHoHbakery{$chrom}{$x}{$i}/$window;
        push @bakerytemp, $het;
      } else {
        push @bakerytemp, 0;
      }
    }
    if (defined $group1[0]) {
      if (defined $HoHoHbakery{$chrom}{$x}{"g1"}) {
        my $het1 = sprintf "%.8f", $HoHoHbakery{$chrom}{$x}{"g1"}/$window;
        push @bakerytemp, $het1;
      } else {
        push @bakerytemp, 0;
      }
    }
    if (defined $group2[0]) {
      if (defined $HoHoHbakery{$chrom}{$x}{"g2"}) {
        my $het2 = sprintf "%.8f", $HoHoHbakery{$chrom}{$x}{"g2"}/$window;
        push @bakerytemp, $het2;
      } else {
        push @bakerytemp, 0;
      }
    }
    if (defined $group3[0]) {
      if (defined $HoHoHbakery{$chrom}{$x}{"g3"}) {
        my $het3 = sprintf "%.8f", $HoHoHbakery{$chrom}{$x}{"g3"}/$window;
        push @bakerytemp, $het3;
      } else {
        push @bakerytemp, 0;
      }
    }
    unless (defined $HoHsnps{$chrom}{$x}) {
      $HoHsnps{$chrom}{$x} = 0;
    }
    push @bakerytemp, $HoHsnps{$chrom}{$x};
    unless (defined $HoHhw{$chrom}{$x}) {
      $HoHhw{$chrom}{$x} = 0;
    }
    push @bakerytemp, $HoHhw{$chrom}{$x};
    my $bakerytemp = join "\t", @bakerytemp;
    push @bakeryout, $bakerytemp;
  }
}

my $result = join "\n", @bakeryout;

unless ( open(META, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print META $result;
close (META);

###################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}

