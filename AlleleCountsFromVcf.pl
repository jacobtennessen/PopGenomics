#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_v $opt_o $opt_d $opt_e $opt_l $opt_c $opt_m $opt_f $opt_z $opt_h $opt_n $opt_s $opt_q);

# Usage
my $usage = "
AlleleCountsFromVcf.pl

Converts a vcf into allele counts for R

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

Usage: perl AlleleCountsFromVcf.pl options
 required:
  -v	(path to) vcf file
  -o  (path to) output file baseline name
 optional:
  -d  minimum depth (default = 5)
  -e  maximun depth (default = 50)
  -l  comma-delimited list of sample positions to use (0-based, first sample column = 0)
  -c  comma-delimited list of chromosomes to use
  -m  number of missing genotypes allowed per marker [default = 0]
  -f  mimimum count of each allele [default = 1]
  -z  file is zipped
  -h  phred-scaled threshold for excess heterozygosity [default = none]
  -n  only normal SNPs (two alleles, not indels)
  -s  only SNPs (multiallelic ok; if indel is called in multiallelic, ignore it)
  -q  genotype quality threshold [default = none]
";

#############

getopts('v:o:d:e:l:c:m:f:z:h:q:ns');
die $usage unless ($opt_v);
die $usage unless ($opt_o);

my ($genotypes, $out, $mindepth, $maxdepth, @samplelist, @chromlist, $missingthresh, $mincount, $hetthresh, $genoqual, $normalsnps, $allsnps);

$genotypes = $opt_v if $opt_v;

$out = $opt_o if $opt_o;

if (defined $opt_d) {
  $mindepth = $opt_d;
} else {
  $mindepth = 5;
}

if (defined $opt_q) {
  $genoqual = $opt_q;
}

if (defined $opt_e) {
  $maxdepth = $opt_e;
} else {
  $maxdepth = 50;
}

my $outfile = "$out"."_$mindepth"."_$maxdepth.txt";

if (defined $opt_l) {
  @samplelist = split ",", $opt_l;
}

if (defined $opt_c) {
  @chromlist = split ",", $opt_c;
}

if (defined $opt_m) {
  $missingthresh = $opt_m;
} else {
  $missingthresh = 0;
}

if (defined $opt_f) {
  $mincount = $opt_f;
} else {
  $mincount = 1;
}

if (defined $opt_h) {
  $hetthresh = $opt_h;
}

if (defined $opt_n) {
  $normalsnps = 1;
}

if (defined $opt_s) {
  $allsnps = 1;
}

if (defined $opt_z) {
  open(IN, "gunzip -c $genotypes |") || die "can’t open pipe to $genotypes";
} else {
  open(IN, $genotypes) || die "can't open $genotypes";
}

my %chroms;

foreach my $c (@chromlist) {
  $chroms{$c} = 1;
}

my $outputsize = 1000;

my @allnames;

my @allsamples;

my $namelist;

my $linecount = 0;

my @out;

my $cleared;

while (<IN>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ($data[0] =~ /CHROM/) {
    my $samplecount = 0;
    my @allsamplelist;
    for (my $fnx = 9; $fnx < (scalar(@data)); $fnx ++) {
        push @allnames, $data[$fnx];
        push @allsamples, $fnx;
        push @allsamplelist, $samplecount;
        $samplecount +=1;
    }
    unless (defined $samplelist[0]) {
      push @samplelist, @allsamplelist;
    }
    $namelist = join "\t", @allnames[@samplelist];
    push @out, "Site\t$namelist";
    next;
  } elsif ($line =~ /^#/) {
      next;
  }
  if (defined $chromlist[0]) {
    unless (defined $chroms{$data[0]}) {
      next;
    }
  }
  if (defined $hetthresh) {
    my @hetdata = split "ExcessHet=", $data[7];
    if (defined $hetdata[1]) {
      my @hetdata2 = split ";", $hetdata[1];
      if ($hetdata2[0] > $hetthresh) {
        next;
      }
    }
  }
  my %indelallele;
  if (defined $normalsnps) {
    unless ((length($data[3]) == 1)&&(length($data[4]) == 1)&&($data[3] =~ /[a-zA-Z]/)&&($data[4] =~ /[a-zA-Z]/)) {
      next;
    }
  } elsif (defined $allsnps) {
    my $justindel;
    if ($data[4] =~ /,/) {
      if ((length($data[3]) == 1)&&($data[3] =~ /[a-zA-Z]/)) {
        my $snp;
        my @alts = split ",", $data[4];
        my $acount = 1;
        foreach my $a (@alts) {
          if ((length($a) == 1)&&($a =~ /[a-zA-Z]/)) {
            $snp = 1;
          } else {
            $indelallele{$acount} = 1;
          }
          $acount +=1;
        }
        unless (defined $snp) {
          $justindel = 1;
        }
      } else {
        $justindel = 1;
      }
    } else {
      unless ((length($data[3]) == 1)&&(length($data[4]) == 1)&&($data[3] =~ /[a-zA-Z]/)&&($data[4] =~ /[a-zA-Z]/)) {
        $justindel = 1;
      }
    }
    if (defined $justindel) {
      next;
    }
  }
  my @inds = @data[@allsamples[@samplelist]];
  my $genoflag = 0;
  my @genos;
  my $seenref = 0;
  my $seenalt = 0;
  foreach my $i (@inds) {
    my @idata = split ":", $i;
    unless ((defined $idata[2])&&($idata[2] =~ /\d/)) {
      $idata[2] = 0;
    }
    unless ((defined $idata[3])&&($idata[3] =~ /\d/)) {
      $idata[3] = 0;
    }
    if ((defined $genoqual)&&($idata[3] < $genoqual)) {
      push @genos, "NA";
      $genoflag += 1;
    } else {
      my $excludeindel;
      foreach my $indel (keys %indelallele) {
        if ($idata[0] =~ /$indel/) {
          $excludeindel = 1;
        }
      }
      if (($idata[2] > $maxdepth)||($idata[2] < $mindepth)||(defined $excludeindel)) {
        push @genos, "NA";
        $genoflag += 1;
      } else {
        if (($idata[0] =~ /0\|0/)||($idata[0] =~ /0\/0/)) {
          push @genos, 0;
          $seenref += 2;
        } elsif (($idata[0] =~ /0\|1/)||($idata[0] =~ /1\|0/)||($idata[0] =~ /2\|0/)||($idata[0] =~ /0\|2/)||($idata[0] =~ /0\/1/)||($idata[0] =~ /1\/0/)||($idata[0] =~ /2\/0/)||($idata[0] =~ /0\/2/)) {
          push @genos, 1;
          $seenref += 1;
          $seenalt += 1;
        } elsif (($idata[0] =~ /1\|1/)||($idata[0] =~ /1\|2/)||($idata[0] =~ /2\|1/)||($idata[0] =~ /2\|2/)||($idata[0] =~ /1\/1/)||($idata[0] =~ /1\/2/)||($idata[0] =~ /2\/1/)||($idata[0] =~ /2\/2/)) {
          push @genos, 2;
          $seenalt += 2;
        } else {
          push @genos, "NA";
          $genoflag += 1;
        }
      }
    }
  }
  if (($genoflag <= $missingthresh)&&($seenref >= $mincount)&&($seenalt >= $mincount)) {
    my $genolist = join "\t", @genos;
    push @out, "$data[0]_$data[1]\t$genolist";
  }
  if (scalar (@out) >= $outputsize) {
      my $result = join "\n", @out;
      if (defined $cleared) {
          unless ( open(META, ">>$outfile") ) {
              print "Cannot open file \"$outfile\" to write to!!\n\n";
              exit;
          }
          print META "$result\n";
          close (META);
      } else {
          unless ( open(META, ">$outfile") ) {
              print "Cannot open file \"$outfile\" to write to!!\n\n";
              exit;
          }
          print META "$result\n";
          close (META);
          $cleared = 1;
      }
      @out = ();
  }
  $linecount +=1;
}

close (IN);

if (scalar (@out) >= 0) {
    my $result = join "\n", @out;
    if (defined $cleared) {
        unless ( open(META, ">>$outfile") ) {
            print "Cannot open file \"$outfile\" to write to!!\n\n";
            exit;
        }
        print META "$result\n";
        close (META);
    } else {
        unless ( open(META, ">$outfile") ) {
            print "Cannot open file \"$outfile\" to write to!!\n\n";
            exit;
        }
        print META "$result\n";
        close (META);
        $cleared = 1;
    }
    @out = ();
}