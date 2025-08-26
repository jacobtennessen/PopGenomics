#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_v $opt_o $opt_a $opt_b $opt_c $opt_u $opt_w $opt_g $opt_l $opt_s);

# Usage
my $usage = "
GetPiAnnotatedPerGene.pl

Calculate pi per gene at site types

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

Usage: perl GetPiAnnotatedPerGene.pl options
 required:
  -v	(path to) allele counts table (annotated)
  -o  (path to) output file
  -g  (path to) gff
optional:
  -a comma-delimted list of samples to be examined as a group [samples start at 1]
  -b comma-delimted list of samples to be examined as a group [samples start at 1]
  -c comma-delimted list of samples to be examnied as a group [samples start at 1]
  -u maximum ratio between het expected and het observed [default = 1.5]
  -w maximum difference between het expected and het observed [default = 5]
  -l comma-delimited list of samples to include [samples start at 1]
  -s (path to) table of synonymous site proportion per gene (product of GetSiteTypesPerGene.pl)
";

#############

getopts('v:o:a:b:c:u:w:g:l:s:');
die $usage unless ($opt_v);
die $usage unless ($opt_o);
die $usage unless ($opt_g);

my ($alleles, $outfile, @group1, @group2, @group3, $maxratio, $maxdif, $gff, @list, $syntable);

$alleles = $opt_v if $opt_v;

$outfile = $opt_o if $opt_o;

$gff = $opt_g if $opt_g;

if (defined $opt_s) {
  $syntable = $opt_s;
}

if (defined $opt_l) {
  @list = split ",", $opt_l;
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

if (defined $opt_a) {
  @group1 = split ",", $opt_a;
}

if (defined $opt_b) {
  @group2 = split ",", $opt_b;
}

if (defined $opt_c) {
  @group3 = split ",", $opt_c;
}

my %chromids = (
  "NC_064873.1" => "X",
  "NC_064874.1" => "2RL",
  "NC_064875.1" => "3RL",
);

my %synscale;

if (defined $syntable) {

  open(SYN, $syntable) || die "can't open $syntable\n";

  while (<SYN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    if ($data[4] =~ /\d/) {
      $synscale{$data[0]} = $data[4];
    }
  }

  close (SYN);

}

if ($gff =~ /gz$/) {
  open(GFF, "gunzip -c $gff |") || die "can’t open pipe to $gff";
} else {
  open(GFF, $gff) || die "can't open $gff\n";
}

my %products;

my %HoHoHoHcds;

my %locations;

my %HoHcds;

while (<GFF>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if (defined $data[2]) {
    if ((defined $data[8])&&($data[3] =~ /\d/)&&($data[4] =~ /\d/)&&(defined $chromids{$data[0]})) {
      my $chrom = $data[0];
      if (defined $chromids{$data[0]}) {
        $chrom = $chromids{$data[0]};
      }
      if ($data[2] =~ /^CDS$/) {
        my $shortsitestart = int($data[3]/100000);
        my $shortsiteend = int($data[4]/100000);
        my @genedata1 = split "gene=", $data[8];
        my @genedata2 = split ";", $genedata1[1];
        my $gene = $genedata2[0];
        my $span = "$data[3]\t$data[4]";
        for (my $ss = $shortsitestart; $ss <= $shortsiteend; $ss++) {
          $HoHoHoHcds{$chrom}{$ss}{$gene}{$span} = 1;
        }
        $HoHcds{$gene}{$span} = 1;
      } elsif ($data[2] =~ /^gene$/) {
        my @genedata1 = split "gene=", $data[8];
        my @genedata2 = split ";", $genedata1[1];
        my $gene = $genedata2[0];
        my @productdata1 = split "description=", $data[8];
        my @productdata2 = split ";", $productdata1[1];
        $products{$gene} = $productdata2[0];
        $locations{$gene} = "$chrom\t$data[3]\t$data[4]";
      }
    }
  }
}
close (GFF);

print "GFF data uploaded.\n";

my %lengths;

foreach my $gene (keys %products) {
  my %coding;
  foreach my $cds ( keys %{ $HoHcds{$gene} } ) {
    my @span = split "\t", $cds;
    for (my $s = $span[0]; $s <= $span[1]; $s++) {
      $coding{$s} = 1;
    }
  }
  $lengths{$gene} = scalar(keys %coding);
  unless (defined $synscale{$gene}) {
    $synscale{$gene} = 0.25;
  }
}

print "Gene lengths calculated.\n";

my %HoHbakery;

my %HoHsnps;

my %hw;

my %HoHoHmetabakery;

my %sitecount;

my $linecount = 0;

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
  unless ($realdata[0] =~ /Site/) {
    my @sitedata = split "_", $realdata[0];
    my $siteexact = pop @sitedata;
    my $shortsite = int($siteexact/100000);
    my $formalchrom = join "_", @sitedata;
    my %genes;
    if (defined $chromids{$formalchrom}) {
      my $chrom = $chromids{$formalchrom}; 
      foreach my $gene ( keys %{ $HoHoHoHcds{$chrom}{$shortsite} } ) {
        foreach my $cds ( keys %{ $HoHoHoHcds{$chrom}{$shortsite}{$gene} } ) {
          my @span = split "\t", $cds;
          if (($span[0] <= $siteexact)&&($span[1] >= $siteexact)) {
            $genes{$gene} = 1;
          }
        }
      }
    }
    $linecount +=1;
    foreach my $gene (keys %genes) {
      my $altcount = 0;
      my $goodcount = 0;
      my $hetcount = 0;
      for (my $i = 1; $i < ((scalar(@realdata))-1); $i++) {
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
      my $freq = $altcount/(2*$goodcount);
      my $het = 2*($freq)*(1-$freq);
      my $expectedhet = $het*$goodcount;
      unless (($freq == 0)||($freq == 1)) {
        my $ratio = $hetcount/$expectedhet;
        my $difference = abs($hetcount - $expectedhet);
        if ((($ratio < $maxratio)&&($ratio > (1/$maxratio)))||($difference <= $maxdif)) {
          if (defined $HoHsnps{$gene}{$origdata[-1]}) {
            $HoHsnps{$gene}{$origdata[-1]} += 1;
          } else {
            $HoHsnps{$gene}{$origdata[-1]} = 1;
          }
          if (defined $HoHbakery{$gene}{$origdata[-1]}) {
            $HoHbakery{$gene}{$origdata[-1]} += $het;
          } else {
            $HoHbakery{$gene}{$origdata[-1]} = $het;
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
            my $freq1 = $altcount1/(2*$goodcount1);
            my $het1 = 2*($freq1)*(1-$freq1);
            if (defined $HoHoHmetabakery{$gene}{"g1"}{$origdata[-1]}) {
              $HoHoHmetabakery{$gene}{"g1"}{$origdata[-1]} += $het1;
            } else {
              $HoHoHmetabakery{$gene}{"g1"}{$origdata[-1]} = $het1;
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
            my $freq2 = $altcount2/(2*$goodcount2);
            my $het2 = 2*($freq2)*(1-$freq2);
            if (defined $HoHoHmetabakery{$gene}{"g2"}{$origdata[-1]}) {
              $HoHoHmetabakery{$gene}{"g2"}{$origdata[-1]} += $het2;
            } else {
              $HoHoHmetabakery{$gene}{"g2"}{$origdata[-1]} = $het2;
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
            my $freq3 = $altcount3/(2*$goodcount3);
            my $het3 = 2*($freq3)*(1-$freq3);
            if (defined $HoHoHmetabakery{$gene}{"g3"}{$origdata[-1]}) {
              $HoHoHmetabakery{$gene}{"g3"}{$origdata[-1]} += $het3;
            } else {
              $HoHoHmetabakery{$gene}{"g3"}{$origdata[-1]} = $het3;
            }
          }
          $sitecount{$realdata[0]} = 1;
        } else {
          if (defined $hw{$gene}) {
            $hw{$gene} += 1;
          } else {
            $hw{$gene} = 1;
          }
        }
      }
    }
  }
}

close (IN);

my $sitecount = scalar (keys %sitecount);

print "$sitecount SNPs assessed.\n";

my @bakeryout;

foreach my $gene (sort keys %products) {
  my @bakerytemp;
  push @bakerytemp, $gene;
  if (defined $locations{$gene}) {
    push @bakerytemp, $locations{$gene};
  } else {
    push @bakerytemp, "NA";
  }
  push @bakerytemp, $lengths{$gene};
  push @bakerytemp, $products{$gene};
  my $piN = 0;
  my $piS = 0;
  if (defined $HoHbakery{$gene}{"nonsyn"}) {
    $piN += sprintf "%.6f", $HoHbakery{$gene}{"nonsyn"}/($lengths{$gene}*(1-$synscale{$gene}));
  }
  if (defined $HoHbakery{$gene}{"splice"}) {
    $piN += sprintf "%.6f", $HoHbakery{$gene}{"splice"}/($lengths{$gene}*(1-$synscale{$gene}));
  }
  if (defined $HoHbakery{$gene}{"syn"}) {
    $piS += sprintf "%.6f", $HoHbakery{$gene}{"syn"}/($lengths{$gene}*$synscale{$gene});
  }
  push @bakerytemp, $piN;
  push @bakerytemp, $piS;
  if (defined $group1[0]) {
    my $piN1 = 0;
    my $piS1 = 0;
    if (defined $HoHoHmetabakery{$gene}{"g1"}{"nonsyn"}) {
      $piN1 += sprintf "%.6f",$HoHoHmetabakery{$gene}{"g1"}{"nonsyn"}/($lengths{$gene}*(1-$synscale{$gene}));
    }
    if (defined $HoHoHmetabakery{$gene}{"g1"}{"splice"}) {
      $piN1 += sprintf "%.6f", $HoHoHmetabakery{$gene}{"g1"}{"splice"}/($lengths{$gene}*(1-$synscale{$gene}));
    }
    if (defined $HoHoHmetabakery{$gene}{"g1"}{"syn"}) {
      $piS1 += sprintf "%.6f", $HoHoHmetabakery{$gene}{"g1"}{"syn"}/($lengths{$gene}*$synscale{$gene});
    }
    push @bakerytemp, $piN1;
    push @bakerytemp, $piS1;
  }
  if (defined $group2[0]) {
    my $piN2 = 0;
    my $piS2 = 0;
    if (defined $HoHoHmetabakery{$gene}{"g2"}{"nonsyn"}) {
      $piN2 += sprintf "%.6f",$HoHoHmetabakery{$gene}{"g2"}{"nonsyn"}/($lengths{$gene}*(1-$synscale{$gene}));
    }
    if (defined $HoHoHmetabakery{$gene}{"g2"}{"splice"}) {
      $piN2 += sprintf "%.6f", $HoHoHmetabakery{$gene}{"g2"}{"splice"}/($lengths{$gene}*(1-$synscale{$gene}));
    }
    if (defined $HoHoHmetabakery{$gene}{"g2"}{"syn"}) {
      $piS2 += sprintf "%.6f", $HoHoHmetabakery{$gene}{"g2"}{"syn"}/($lengths{$gene}*$synscale{$gene});
    }
    push @bakerytemp, $piN2;
    push @bakerytemp, $piS2;
  }
  if (defined $group3[0]) {
    my $piN3 = 0;
    my $piS3 = 0;
    if (defined $HoHoHmetabakery{$gene}{"g3"}{"nonsyn"}) {
      $piN3 += sprintf "%.6f",$HoHoHmetabakery{$gene}{"g3"}{"nonsyn"}/($lengths{$gene}*(1-$synscale{$gene}));
    }
    if (defined $HoHoHmetabakery{$gene}{"g3"}{"splice"}) {
      $piN3 += sprintf "%.6f", $HoHoHmetabakery{$gene}{"g3"}{"splice"}/($lengths{$gene}*(1-$synscale{$gene}));
    }
    if (defined $HoHoHmetabakery{$gene}{"g3"}{"syn"}) {
      $piS3 += sprintf "%.6f", $HoHoHmetabakery{$gene}{"g3"}{"syn"}/($lengths{$gene}*$synscale{$gene});
    }
    push @bakerytemp, $piN3;
    push @bakerytemp, $piS3;
  }
  my $snpsN = 0;
  my $snpsS = 0;
  if (defined $HoHsnps{$gene}{"nonsyn"}) {
    $snpsN += $HoHsnps{$gene}{"nonsyn"};
  }
  if (defined $HoHsnps{$gene}{"splice"}) {
    $snpsN += $HoHsnps{$gene}{"splice"};
  }
  if (defined $HoHsnps{$gene}{"syn"}) {
    $snpsS += $HoHsnps{$gene}{"syn"};
  }
  push @bakerytemp, $snpsN;
  push @bakerytemp, $snpsS;
  unless (defined $hw{$gene}) {
    $hw{$gene} = 0;
  }
  push @bakerytemp, $hw{$gene};
  my $bakerytemp = join "\t", @bakerytemp;
  push @bakeryout, $bakerytemp;
}

my @title = ("Gene","Chrom","Start","Stop","Length","Description","PiN","PiS");

if (defined $group1[0]) {
  push @title, "Group1PiN";
  push @title, "Group1PiS";
}

if (defined $group2[0]) {
  push @title, "Group2PiN";
  push @title, "Group2PiS";
}

if (defined $group3[0]) {
  push @title, "Group3PiN";
  push @title, "Group3PiS";
}

push @title, "SNPsN";

push @title, "SNPsS";

push @title, "HW";

my $title = join "\t", @title;

unshift @bakeryout, $title;

my $result = join "\n", @bakeryout;

$result =~ s/ /_/g;

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