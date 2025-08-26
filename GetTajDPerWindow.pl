#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_v $opt_o $opt_a $opt_b $opt_c $opt_u $opt_w $opt_l $opt_t $opt_s $opt_m);

# Usage
my $usage = "
GetTajDPerWindow.pl

Calculates Tajima's D per window

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

Usage: perl GetTajDPerWindow.pl options
 required:
  -v	(path to) allele counts table
  -o  (path to) output file
optional:
  -a comma-delimted list of samples to be examined as a group [samples start at 1]
  -b comma-delimted list of samples to be examined as a group [samples start at 1]
  -c comma-delimted list of samples to be examnied as a group [samples start at 1]
  -u maximum ratio between het expected and het observed [default = 1.5]
  -w maximum difference between het expected and het observed [default = 5]
  -l comma-delimited list of samples to include [samples start at 1]
  -t translate chromosome
  -s window size [default = 10000]
  -m max missing samples [default = no limit]
";

#############

getopts('v:o:a:b:c:u:w:l:s:m:t');
die $usage unless ($opt_v);
die $usage unless ($opt_o);

my ($alleles, $outfile, @group1, @group2, @group3, $maxratio, $maxdif, @list, $translate, $window, $maxmiss);

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

if (defined $opt_w) {
  $window = $opt_w;
} else {
  $window = 10000;
}

if (defined $opt_w) {
  $maxdif = $opt_w;
} else {
  $maxdif = 5;
}

if (defined $opt_m) {
  $maxmiss = $opt_m;
}

my $g1size;

if (defined $opt_a) {
  @group1 = split ",", $opt_a;
  $g1size = scalar(@group1);
}

my $g1min = 1;
if (defined $maxmiss) {
  $g1min = $g1size - $maxmiss;
}

my $g2size;

if (defined $opt_b) {
  @group2 = split ",", $opt_b;
  $g2size = scalar(@group2);
}

my $g2min = 1;
if (defined $maxmiss) {
  $g2min = $g2size - $maxmiss;
}

my $g3size;

if (defined $opt_c) {
  @group3 = split ",", $opt_c;
  $g3size = scalar(@group3);
}

my $g3min = 1;
if (defined $maxmiss) {
  $g3min = $g3size - $maxmiss;
}

my %chromids = (
  "NC_064873.1" => "X",
  "NC_064874.1" => "2RL",
  "NC_064875.1" => "3RL",
);

my $allsize = 0;

my $allmin = 1;

my %bakery;

my %snps;

my %HoHsnps;

my %hw;

my %HoHmetabakery;

my %names;

my %sitecount;

my $linecount = 0;

my $anno;

my %windows;

open(IN, $alleles) || die "can't open $alleles\n";

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
    if ($realdata[-1] =~ /Annotation/) {
      $anno = pop @realdata;
    }
    for (my $i = 1; $i < (scalar(@realdata)); $i++) {
      $names{$i} = $realdata[$i];
      $allsize +=1;
    }
    if (defined $maxmiss) {
      $allmin = $allsize - $maxmiss;
    }
  } else {
    if (defined $anno) {
      my $junk = pop @realdata;
    }
    my @sitedata = split "_", $realdata[0];
    my $site = pop @sitedata;
    my $chromname = join "_", @sitedata;
    my $startwindow = $window*(int($site/$window));
    my $endwindow = $startwindow + $window;
    my $localwindow = "$chromname\t$startwindow\t$endwindow";
    $windows{$localwindow} = 1;
    $linecount +=1;
    my $altcount = 0;
    my $goodcount = 0;
    my $hetcount = 0;
    for (my $i = 1; $i < (scalar(@realdata)); $i++) {
      if ($realdata[$i] =~ /\d/) {
        if ($realdata[$i] == 0) {
          $goodcount +=1;
        } elsif ($realdata[$i] == 1) {
          $altcount +=1;
          $goodcount +=1;
          $hetcount += 1;
        } elsif ($realdata[$i] == 2) {
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
          if ($goodcount >= $allmin) {
            for (my $i = 1; $i < (scalar(@realdata)); $i++) {
              if ($realdata[$i] =~ /\d/) {
                if ($realdata[$i] == 1) {
                  if (defined $HoHmetabakery{$localwindow}{$i}) {
                    $HoHmetabakery{$localwindow}{$i} += 1;
                  } else {
                    $HoHmetabakery{$localwindow}{$i} = 1;
                  }
                }
              }
            }
            if (defined $snps{$localwindow}) {
              $snps{$localwindow} += 1;
            } else {
              $snps{$localwindow} = 1;
            }
            if (defined $bakery{$localwindow}) {
              $bakery{$localwindow} += $het;
            } else {
              $bakery{$localwindow} = $het;
            }
          }
          if (defined $group1[0]) {
            my $goodcount1 = 0;
            my $altcount1 = 0;
            foreach my $g1 (@group1) {
              if ($realdata[$g1] =~ /\d/) {
                if ($realdata[$g1] == 0) {
                  $goodcount1 +=1;
                } if ($realdata[$g1] == 1) {
                  $altcount1 +=1;
                  $goodcount1 +=1;
                } elsif ($realdata[$g1] == 2) {
                  $altcount1 +=2;
                  $goodcount1 +=1;
                }
              }
            }
            if ($goodcount1 >= $g1min) {
              my $freq1 = $altcount1/(2*$goodcount1);
              unless (($altcount1 == 0)||($altcount1 == (2*$goodcount1))) {
                if (defined $HoHsnps{$localwindow}{"g1"}) {
                  $HoHsnps{$localwindow}{"g1"} += 1;
                } else {
                  $HoHsnps{$localwindow}{"g1"} = 1;
                }
                my $het1 = 2*($freq1)*(1-$freq1);
                if (defined $HoHmetabakery{$localwindow}{"g1"}) {
                  $HoHmetabakery{$localwindow}{"g1"} += $het1;
                } else {
                  $HoHmetabakery{$localwindow}{"g1"} = $het1;
                }
              }
            }
          }
          if (defined $group2[0]) {
            my $goodcount2 = 0;
            my $altcount2 = 0;
            foreach my $g2 (@group2) {
              if ($realdata[$g2] =~ /\d/) {
                if ($realdata[$g2] == 0) {
                  $goodcount2 +=1;
                } if ($realdata[$g2] == 1) {
                  $altcount2 +=1;
                  $goodcount2 +=1;
                } elsif ($realdata[$g2] == 2) {
                  $altcount2 +=2;
                  $goodcount2 +=1;
                }
              }
            }
            if ($goodcount2 >= $g2min) {
              my $freq2 = $altcount2/(2*$goodcount2);
              unless (($altcount2 == 0)||($altcount2 == (2*$goodcount2))) {
                if (defined $HoHsnps{$localwindow}{"g2"}) {
                  $HoHsnps{$localwindow}{"g2"} += 1;
                } else {
                  $HoHsnps{$localwindow}{"g2"} = 1;
                }
                my $het2 = 2*($freq2)*(1-$freq2);
                if (defined $HoHmetabakery{$localwindow}{"g2"}) {
                  $HoHmetabakery{$localwindow}{"g2"} += $het2;
                } else {
                  $HoHmetabakery{$localwindow}{"g2"} = $het2;
                }
              }
            }
          }
          if (defined $group3[0]) {
            my $goodcount3 = 0;
            my $altcount3 = 0;
            foreach my $g3 (@group3) {
              if ($realdata[$g3] =~ /\d/) {
                if ($realdata[$g3] == 0) {
                  $goodcount3 +=1;
                } if ($realdata[$g3] == 1) {
                  $altcount3 +=1;
                  $goodcount3 +=1;
                } elsif ($realdata[$g3] == 2) {
                  $altcount3 +=2;
                  $goodcount3 +=1;
                }
              }
            }
            if ($goodcount3 >= $g3min) {
              my $freq3 = $altcount3/(2*$goodcount3);
              unless (($altcount3 == 0)||($altcount3 == (2*$goodcount3))) {
                if (defined $HoHsnps{$localwindow}{"g3"}) {
                  $HoHsnps{$localwindow}{"g3"} += 1;
                } else {
                  $HoHsnps{$localwindow}{"g3"} = 1;
                }
                my $het3 = 2*($freq3)*(1-$freq3);
                if (defined $HoHmetabakery{$localwindow}{"g3"}) {
                  $HoHmetabakery{$localwindow}{"g3"} += $het3;
                } else {
                  $HoHmetabakery{$localwindow}{"g3"} = $het3;
                }
              }
            }
          }
          $sitecount{$realdata[0]} = 1;
        } else {
          if (defined $hw{$localwindow}) {
            $hw{$localwindow} += 1;
          } else {
            $hw{$localwindow} = 1;
          }
        }
      }
    }
  }
}

close (IN);

my $sitecount = scalar (keys %sitecount);

print "$sitecount SNPs assessed.\n";

my @names = sort by_number (keys (%names));

my @bakeryout;

foreach my $localwindow (keys %windows) {
  my @bakerytemp;
  unless (defined $snps{$localwindow}) {
    $snps{$localwindow} = 0;
  }
  if ((defined $bakery{$localwindow})&&($snps{$localwindow} > 0)) {
    my $pi = sprintf "%.6f", $bakery{$localwindow}/$window;
    push @bakerytemp, "$localwindow\t$pi";
    push @bakerytemp, sprintf "%.3f", tajima_D_counts(2*$allsize,$snps{$localwindow},$bakery{$localwindow});
  } else {
    push @bakerytemp, "$localwindow\t0";
    push @bakerytemp, "NA";
  }
  if (defined $group1[0]) {
    if ((defined $HoHmetabakery{$localwindow}{"g1"})&&($HoHsnps{$localwindow}{"g1"} > 0)) {
      my $het1 = sprintf "%.8f", $HoHmetabakery{$localwindow}{"g1"}/$window;
      push @bakerytemp, $het1;
      push @bakerytemp, sprintf "%.3f", tajima_D_counts(2*$g1size,$HoHsnps{$localwindow}{"g1"},$HoHmetabakery{$localwindow}{"g1"});
    } else {
      push @bakerytemp, 0;
      push @bakerytemp, "NA";
    }
  }
  if (defined $group2[0]) {
    if ((defined $HoHmetabakery{$localwindow}{"g2"})&&($HoHsnps{$localwindow}{"g2"} > 0)) {
      my $het2 = sprintf "%.8f", $HoHmetabakery{$localwindow}{"g2"}/$window;
      push @bakerytemp, $het2;
      push @bakerytemp, sprintf "%.3f", tajima_D_counts(2*$g2size,$HoHsnps{$localwindow}{"g2"},$HoHmetabakery{$localwindow}{"g2"});
    } else {
      push @bakerytemp, 0;
      push @bakerytemp, "NA";
    }
  }
  if (defined $group3[0]) {
    if ((defined $HoHmetabakery{$localwindow}{"g3"})&&($HoHsnps{$localwindow}{"g3"} > 0)) {
      my $het3 = sprintf "%.8f", $HoHmetabakery{$localwindow}{"g3"}/$window;
      push @bakerytemp, $het3;
      push @bakerytemp, sprintf "%.3f", tajima_D_counts(2*$g3size,$HoHsnps{$localwindow}{"g3"},$HoHmetabakery{$localwindow}{"g3"});
    } else {
      push @bakerytemp, 0;
      push @bakerytemp, "NA";
    }
  }
  push @bakerytemp, $snps{$localwindow};
  unless (defined $hw{$localwindow}) {
    $hw{$localwindow} = 0;
  }
  push @bakerytemp, $hw{$localwindow};
  my $bakerytemp = join "\t", @bakerytemp;
  push @bakeryout, $bakerytemp;
}

my @title = ("Chrom","Start","Stop","Pi","D");

if (defined $group1[0]) {
  push @title, "Group1Pi";
  push @title, "Group1D";
}

if (defined $group2[0]) {
  push @title, "Group2Pi";
  push @title, "Group2D";
}

if (defined $group3[0]) {
  push @title, "Group3Pi";
  push @title, "Group3D";
}

push @title, "SNPs";

push @title, "HW";

my $title = join "\t", @title;

unshift @bakeryout, $title;

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
 
######################
 
sub tajima_D_counts {
    my ($n,$seg_sites,$pi) = @_;
    my $a1 = 0; 
    for(my $k= 1; $k < $n; $k++ ) {
    $a1 += ( 1 / $k );
    }
 
     my $a2 = 0;
     for(my $k= 1; $k < $n; $k++ ) {
     $a2 += ( 1 / $k**2 );
     }
     
    my $b1 = ( $n + 1 ) / ( 3* ( $n - 1) );
    my $b2 = ( 2 * ( $n ** 2 + $n + 3) ) / 
         ( ( 9 * $n) * ( $n - 1) );
    my $c1 = $b1 - ( 1 / $a1 );
    my $c2 = $b2 - ( ( $n + 2 ) /
             ( $a1 * $n))+( $a2 / $a1 ** 2);
    my $e1 = $c1 / $a1;
    my $e2 = $c2 / ( $a1**2 + $a2 );
     
    my $denom = sqrt ( ($e1 * $seg_sites) + (( $e2 * $seg_sites) * ( $seg_sites - 1)));
    return if $denom == 0;
    my $D = ( $pi - ( $seg_sites / $a1 ) ) / $denom;
    return $D;
}