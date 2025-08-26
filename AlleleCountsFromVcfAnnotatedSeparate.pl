#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_v $opt_o $opt_d $opt_e $opt_l $opt_c $opt_m $opt_f $opt_z $opt_h $opt_a $opt_r $opt_y $opt_s $opt_q);

# Usage
my $usage = "
AlleleCountsFromVcfAnnotatedSeparate.pl

Converts a vcf into allele counts for R, using a separate vcf for annotation information

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

Usage: perl AlleleCountsFromVcfAnnotatedSeparate.pl options
 required:
  -v	(path to) vcf file
  -a	(path to) annotated vcf file
  -o  (path to) outout file baseline name
 optional:
  -d  minimum depth (default = 5)
  -e  maximun depth (default = 50)
  -l  comma-delimited list of sample positions to use (0-based, first sample column = 0)
  -c  comma-delimited list of chromosomes to use
  -m  number of missing genotypes allowed per marker [default = 0]
  -f  mimimum count of each allele [default = 1]
  -z  main vcf file is zipped
  -y  annotation vcf is zipped
  -h  phred-scaled threshold for excess heterozygosity [default = none]
  -r  site range to use [START-END]
  -s  only SNPs (multiallelic ok; if indel is called in multiallelic, ignore it)
  -q  genotype quality threshold [default = none]
";

#############

getopts('v:o:d:e:l:c:m:f:r:z:h:a:y:q:s');
die $usage unless ($opt_v);
die $usage unless ($opt_o);
die $usage unless ($opt_a);

my ($genotypes, $annotations, $out, $mindepth, $maxdepth, @samplelist, @chromlist, $missingthresh, $mincount, @range, $hetthresh, $genoqual, $allsnps);

$genotypes = $opt_v if $opt_v;

$annotations = $opt_a if $opt_a;

$out = $opt_o if $opt_o;

if (defined $opt_d) {
  $mindepth = $opt_d;
} else {
  $mindepth = 5;
}

if (defined $opt_e) {
  $maxdepth = $opt_e;
} else {
  $maxdepth = 50;
}

if (defined $opt_q) {
  $genoqual = $opt_q;
}

if (defined $opt_r) {
  @range = split "-", $opt_r;
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

if (defined $opt_s) {
  $allsnps = 1;
}

if (defined $opt_y) {
  open(ANN, "gunzip -c $annotations |") || die "can’t open pipe to $annotations";
} else {
  open(ANN, $annotations) || die "can't open $annotations";
}

my %HoHoHannotation;

my %chroms;

foreach my $c (@chromlist) {
  $chroms{$c} = 1;
}

my $outputsize = 1000;

while (<ANN>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  my @data = split "\t", $line;
  if ($line =~ /^#/) {
      next;
  }
  if (defined $chromlist[0]) {
    unless (defined $chroms{$data[0]}) {
      next;
    }
  }
  if (defined $range[0]) {
    if ($data[1] > $range[1]) {
      last;
    }
    unless ($data[1] >= $range[0]) {
      next;
    }
  }
  my @annot;
  my @allele;
  my @annodata1 = split "ANN=", $data[7];
  my @annodataalt = split /\|,/, $annodata1[1];
  foreach my $aalt (@annodataalt) {
    my @annodata2 = split /\|/, $aalt;
    if (($annodata2[1] =~ /frame/)||($annodata2[1] =~ /gained/)||($annodata2[1] =~ /lost/)||($annodata2[1] =~ /missense/)||($annodata2[1] =~ /initiator_codon_variant/)) {
      push @annot, "nonsyn";
      push @allele, $annodata2[0];
    } elsif (($annodata2[1] =~ /splice_donor/)||($annodata2[1] =~ /splice_acceptor/)||($annodata2[1] =~ /gene_fusion/)||($annodata2[1] =~ /gene_fusion/)||($annodata2[1] =~ /exon_loss/)||($annodata2[1] =~ /premature_start_codon_gain/)) {
      push @annot, "splice";
      push @allele, $annodata2[0];
    } elsif (($annodata2[1] =~ /synonymous/)||($annodata2[1] =~ /retained/)) {
      push @annot, "syn";
      push @allele, $annodata2[0];
    }
  }
  if (defined $annot[0]) {
    for (my $a = 0; $a < scalar(@annot); $a++) {
      $HoHoHannotation{$data[0]}{$data[1]}{$allele[$a]} = $annot[$a];
    }
  }
}

close (ANN);

print "Annotations uploaded.\n";

if (defined $opt_z) {
  open(IN, "gunzip -c $genotypes |") || die "can’t open pipe to $genotypes";
} else {
  open(IN, $genotypes) || die "can't open $genotypes";
}

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
    push @out, "Site\t$namelist\tAnnotation";
    next;
  } elsif ($line =~ /^#/) {
      next;
  }
  if (defined $chromlist[0]) {
    unless (defined $chroms{$data[0]}) {
      next;
    }
  }
  if (defined $range[0]) {
    if ($data[1] > $range[1]) {
      last;
    }
    unless ($data[1] >= $range[0]) {
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
  unless (defined $HoHoHannotation{$data[0]}{$data[1]}) {
    next;
  }
  my @alleles = split ",", $data[4];
  unless (defined $HoHoHannotation{$data[0]}{$data[1]}{$alleles[0]}) {
    next;
  }
  my $annot = $HoHoHannotation{$data[0]}{$data[1]}{$alleles[0]};
  my %indelallele;
  if (defined $allsnps) {
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
        if (($idata[0] =~ /0\|0/)||($idata[0] =~ /0\/0/)||($idata[0] =~ /2\|2/)||($idata[0] =~ /2\/2/)||($idata[0] =~ /2\|0/)||($idata[0] =~ /0\|2/)||($idata[0] =~ /2\/0/)||($idata[0] =~ /0\/2/)) {
          push @genos, 0;
          $seenref += 2;
        } elsif (($idata[0] =~ /0\|1/)||($idata[0] =~ /1\|0/)||($idata[0] =~ /0\/1/)||($idata[0] =~ /1\/0/)||($idata[0] =~ /1\|2/)||($idata[0] =~ /2\|1/)||($idata[0] =~ /1\/2/)||($idata[0] =~ /2\/1/)) {
          push @genos, 1;
          $seenref += 1;
          $seenalt += 1;
        } elsif (($idata[0] =~ /1\|1/)||($idata[0] =~ /1\/1/)) {
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
    push @out, "$data[0]_$data[1]\t$genolist\t$annot";
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