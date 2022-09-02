#!/usr/bin/perl -w
#------
# Name: xpl2xyz.pl
# Author:J. Michael Word
# Date Written:  7/ 5/2002
# Purpose: reformat an Xplor ascii map file for input to kin3Dcont

# usage: xpl2xyz.pl XPLORmap.txt >mapxyz.txt
#    or: xpl2xyz.pl XPLORmap.txt | kin3Dcont [opt...] -  >map.kin
# where the options to kin3Dcont should include -samp and -first
# and usually -kin
#
# Modifications:
#  7/ 5/02 - JM Word - v1.0 - first cut
#  7/22/02 - JM Word - v1.1 - got non-orthogonal axes working

use strict;              # complain if sloppy programming found

my $version = "$0 v1.2.030212";

if (! defined($ARGV[0])) {
   warn "Reformat an Xplor ascii map file for input to kin3Dcont.\n";
   warn "syntax: $0 XPLORmap.txt >mapxyz.txt\n";
   warn "${version}\n";
   die "command line parameter error, stopped";
}

#----------------------------------------------------------
<>; # first line is going to be blank so ignore
#----------------------------------------------------------
if (! defined($_ = <>)) { # read a line
   die "out of data reading header, stopped";
}
chomp; my ($ntitle) = split " "; # number of title lines
#----------------------------------------------------------
for (my $nt = 0; $nt < $ntitle; $nt++) {
   <>; # ignore title lines
}
#----------------------------------------------------------
if (! defined($_ = <>)) { # read a line
   die "out of data reading header, stopped";
}
chomp;
my @bnds = split " ";          # fill array with bounds

die "error reading bounds, stopped" if scalar(@bnds) < 9;

my ($na, $amin, $amax,
    $nb, $bmin, $bmax,
    $nc, $cmin, $cmax) = @bnds;   # break into individual variables
#----------------------------------------------------------
if (! defined($_ = <>)) { # read a line
   die "out of data reading header, stopped";
}
chomp;
my @cell = unpack "a12a12a12a12a12a12", $_; # fill unit cell array

my ($a_len, $b_len, $c_len, $alpha, $beta, $gamma) = @cell;

my $deg2rad  = 0.017453292519943295769;
my $cosalpha = cos($alpha * $deg2rad);
my $cosbeta  = cos($beta  * $deg2rad);
my $cosgamma = cos($gamma * $deg2rad);
my $singamma = sin($gamma * $deg2rad);

my $scsq = 1.0 - $cosalpha * $cosalpha
               - $cosbeta  * $cosbeta
               - $cosgamma * $cosgamma
               + 2.0 * $cosalpha * $cosbeta * $cosgamma;
if (abs($scsq) < 1.0E-9) { $scsq = 0.0; } # chop small (neg?) numbers

my $vol = $a_len * $b_len * $c_len * sqrt($scsq);

# unit cell grid spacing
my $a_gap = $a_len / $na;
my $b_gap = $b_len / $nb;
my $c_gap = $c_len / $nc;

# unit cell to cartesian multiplication factors
my $cp = $a_gap;
my $cq = $b_gap * $cosgamma;
my $cr = $c_gap * $cosbeta;
my $cs = $b_gap * $singamma;
my $ct = $c_gap * ($cosalpha - $cosbeta * $cosgamma) / $singamma;
my $cu = $c_gap * $vol / ($a_len * $b_len * $c_len * $singamma);
      
my $orthog = ((abs($cq) + abs($cr) + abs($ct)) < 1.0E-9);

#----------------------------------------------------------
if (! defined($_ = <>)) { # read a line
   die "out of data reading header, stopped";
}
chomp;
my ($mode) = split " ";           # data order
die "error in map mode- expected ZYX found $mode, stopped"
	    unless /^\s*ZYX/i;

#==========================================================
# start reading data

my $ksect = -1; # a ksect line will come after header
my $nslice = ($amax - $amin + 1) * ($bmax - $bmin + 1);
my $per = 6; # XPLOR rows are six wide
my $numfull   = int($nslice / $per);
my $remainder = $nslice % $per;

#(ic is the slowest dimension)

for (my $ic = $cmin; $ic <= $cmax; $ic++) {
   my @datarec  = (); # one lines worth of data
   my @datavals = (); # a whole section worth of data
   my ($i, $j);
   my $entry = 0.0;

   # consume data records
   #----------------------------------------------------------
   $_ = <>; chomp; ($ksect) = split " ";  # read the k-section number
   for ($i = 0; $i < $numfull; $i++) {
      #----------------------------------------------------------
      if (! defined($_ = <>)) { # read a line
         die "out of data reading ksect $ksect, stopped";
      }
      @datarec = unpack "a12a12a12a12a12a12", $_; # fill datarec
      for ($j = 0; $j < $per; $j++) {
         $entry = $datarec[$j] * 1.0;
         push @datavals, $entry;
      }
   }
   if ($remainder > 0) {
      #----------------------------------------------------------
      if (! defined($_ = <>)) { # read a line
         die "out of data reading ksect $ksect, stopped";
      }
      @datarec = unpack "a12a12a12a12a12a12", $_; # fill datarec
      for ($j = 0; $j < $remainder; $j++) {
         $entry = $datarec[$j] * 1.0;
         push @datavals, $entry;
      }
   }

   # produce output formatted for kin3Dcont: val, x, y, z

   my ($ia, $ib);

   for ($ib = $bmin; $ib <= $bmax; $ib++) {
      for ($ia = $amin; $ia <= $amax; $ia++) {
         $entry = shift @datavals
	          || die "too few data values in ksect $ksect, stopped";

         # a slight efficiency (saves 3 multiplies and 3 adds)
         if ($orthog) {
            printf("%.4g  %.6g %.6g %.6g\n", $entry,
	              $ia * $cp,
                      $ib * $cs,
                      $ic * $cu);
         }
         else {
            printf("%.4g  %.6g %.6g %.6g\n", $entry,
	              $ia * $cp + $ib * $cq + $ic * $cr,
                      $ib * $cs + $ic * $ct,
                      $ic * $cu);
         }
      }
   }
   if (scalar(@datavals) > 0) {
      my $extra = scalar(@datavals);
      die "$extra too many data values in ksect $ksect, stopped";
   }
}
