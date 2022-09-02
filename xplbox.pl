#!/usr/bin/perl -w
#------
# Name: xplbox.pl
# Author:J. Michael Word
# Date Written:  7/ 5/2002
# Purpose: read the header from an Xplor ascii map file
#          and format a box for the unit cell

# usage: xplbox.pl [-flags] XPLORmap.txt >map.kin
#
# Modifications:
#  7/22/02 - JM Word - v1.0 - first cut
#  7/23/02 - JM Word - v1.1 - added gap, mean, sd data

use strict; # complain if sloppy programming found

my $version = "$0 v1.2.030212";

my $def_contgrid  = 1.0; # default gaps per contour factor
my $def_smoothing = 1.5; # default smoothing factor (times grid size)
my $contgrid  = $def_contgrid;
my $smoothing = $def_smoothing;
my $calcstats = 0;
my $estimate  = 1;

while(defined($ARGV[0]) && $ARGV[0] =~ /^\-\S+/) {
   $_ = shift;
   if (/^-G(\.?\d+\.?\d*)/i) {
      $contgrid = $1 * 1.0;
   }
   elsif (/^-S(\.?\d+\.?\d*)/i) {
      $smoothing = $1 * 1.0;
   }
   elsif (/^-vs$/i) {
      $calcstats = 1;
      $estimate  = 0;
   }
   elsif (/^-q$/i) {
      $calcstats = 1;
      $estimate  = 1;
   }
   else {
      die "unknown flag $_, stopped";
   }
}
if (! defined($ARGV[0])) {
   warn "Read the header from an Xplor ascii map file\n";
   warn "and format a box for the unit cell.\n";
   warn "\n";
   warn "syntax: $0 [-flag] XPLORmap.txt > map.kin\n";
   warn "where flag can be\n";
   warn "   -q     one pass calculation of the mean and sd (quicker than -vs)\n";
   warn "   -vs    two pass calculation of the mean and sd (very slow)\n";
   warn "   -g#.#  scale recommended contour grid by # times map grid\n";
   warn "          (default $def_contgrid)\n";
   warn "   -s#.#  scale recommended contour smoothing (times contour grid)\n";
   warn "          (default $def_smoothing)\n";
   warn "${version}\n";
   die "command line parameter error, stopped";
}

my $infname = $ARGV[0]; # the name of the input file

my @remember; # stuff for the caption

#----------------------------------------------------------
# first line is going to be blank so ignore
if (! defined($_ = <>)) { # read a line
   die "out of data reading header, stopped";
}
#----------------------------------------------------------
if (! defined($_ = <>)) { # read a line
   die "out of data reading number of title lines, stopped";
}
chomp; my ($ntitle) = split " "; # number of title lines
#----------------------------------------------------------
for (my $nt = 0; $nt < $ntitle; $nt++) {
   if (! defined($_ = <>)) { # read a line
      die "out of data reading title lines, stopped";
   }
   chomp;
   # save title lines
   push @remember, $_;
}
#----------------------------------------------------------
if (! defined($_ = <>)) { # read a line
   die "out of data reading bounds, stopped";
}
chomp;
my @bnds = split " ";          # fill array with bounds

die "error reading bounds, stopped" if scalar(@bnds) < 9;

# save bounds line
push @remember, $_;

my ($na, $amin, $amax,
    $nb, $bmin, $bmax,
    $nc, $cmin, $cmax) = @bnds;   # break into individual variables
#----------------------------------------------------------
if (! defined($_ = <>)) { # read a line
   die "out of data reading unit cell, stopped";
}
chomp;
my @cell = unpack "a12a12a12a12a12a12", $_; # fill unit cell array

# save unit cell line
push @remember, $_;

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
   die "out of data reading mode, stopped";
}
chomp;
my ($mode) = split " ";           # data order
die "error in map mode- expected ZYX found $mode, stopped"
	    unless /^\s*ZYX/i;
# save mode line
push @remember, $_;

#==========================================================
# start reading data

my $ksect = -1; # a ksect line will come after header
my $nslice = ($amax - $amin + 1) * ($bmax - $bmin + 1);
my $per = 6; # XPLOR rows are six wide
my $numfull   = int($nslice / $per);
my $remainder = $nslice % $per;

my $min_val = undef;
my $max_val = undef;
my $val_cnt = 0;
my $avg_val = undef;
my $sigma   = 0.0;
my @nums    = ();
my $sum_val = 0.0;
my $ssq_val = 0.0;

if ($calcstats) {
   for (my $ic = $cmin; $ic <= $cmax; $ic++) {
      my @datarec  = (); # one lines worth of data
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
            &consume_ent($entry);
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
            &consume_ent($entry);
         }
      }
   }
   &calc_stats();
}
#==========================================================
# format output

print "\@kinemage 1\n";
print "\@caption\n";
print " $version\n";
foreach my $line (@remember) {
   print " $line\n";
}

#-----------------------------------------------
# WARNING: need to figure out if we should 
#          use -sampled if contgrid variable is not 1.0
#-----------------------------------------------
printf(" %s [-flags] %s > cont.kin\n", $0, $infname);
printf(" xpl2xyz.pl %s > mapxyz.txt\n", $infname);
printf(" kin3Dcont -group -noaxis -lens -first -align \\\n");
printf(" -gxyz %.5g %.5g %.5g",
         $contgrid*$a_gap,
	 $contgrid*$b_gap,
	 $contgrid*$c_gap);
printf(" -nosmooth \\\n");
#printf(" -sxyz %.2g %.2g %.2g \\\n",
#         $contgrid*$smoothing*$a_gap,
#	 $contgrid*$smoothing*$b_gap,
#	 $contgrid*$smoothing*$c_gap);
printf(" -slevel 2 sky 3 orange mapxyz.txt >> cont.kin\n");

if ($val_cnt > 1) {
   printf(" Levels at 1 and 2 sigma: -level %6g sky %6g orange\n",
            $avg_val+$sigma,
	    $avg_val+2.0*$sigma);
   printf(" (%d samples, min=%6g, max=%6g, avg=%6g, sigma=%6g)\n",
            $val_cnt, $min_val, $max_val, $avg_val, $sigma);
}
#-----------------------------------------------

print "\@perspective\n";
print "\@group {box}\n";
print "\@vectorlist {unit_cell} color=gray width=1\n";

&drawbox(0, $na, 0, $nb, 0, $nc);

print "\@vectorlist {data_range} color=yellow width=1\n";
&drawbox($amin, $amax, $bmin, $bmax, $cmin, $cmax);

sub drawbox() {
   my ($a_low,$a_high,$b_low,$b_high,$c_low,$c_high) = @_;

&writepoint("000", "P", $a_low,  $b_low,  $c_low);
&writepoint("100", "L", $a_high, $b_low,  $c_low);

&writepoint("001", "P", $a_low,  $b_low,  $c_high);
&writepoint("101", "L", $a_high, $b_low,  $c_high);

&writepoint("010", "P", $a_low,  $b_high, $c_low);
&writepoint("110", "L", $a_high, $b_high, $c_low);

&writepoint("011", "P", $a_low,  $b_high, $c_high);
&writepoint("111", "L", $a_high, $b_high, $c_high);

&writepoint("000", "P", $a_low,  $b_low,  $c_low);
&writepoint("001", "L", $a_low,  $b_low,  $c_high);
&writepoint("011", "L", $a_low,  $b_high, $c_high);
&writepoint("010", "L", $a_low,  $b_high, $c_low);
&writepoint("000", "L", $a_low,  $b_low,  $c_low);

&writepoint("100", "P", $a_high, $b_low,  $c_low);
&writepoint("101", "L", $a_high, $b_low,  $c_high);
&writepoint("111", "L", $a_high, $b_high, $c_high);
&writepoint("110", "L", $a_high, $b_high, $c_low);
&writepoint("100", "L", $a_high, $b_low,  $c_low);
}

sub writepoint() {
   my ($lbl, $type, $ia, $ib, $ic) = @_;

   my ($px, $py, $pz) = convertI2F($ia, $ib, $ic);

   printf("{%s}%s %.6g %.6g %.6g\n", $lbl, $type, $px, $py, $pz);
}

sub convertI2F() {
   my ($ia, $ib, $ic) = @_;

   return ($ia * $cp + $ib * $cq + $ic * $cr,
           $ib * $cs + $ic * $ct,
           $ic * $cu);
}
sub consume_ent() {
   my ($e) = @_;
   if ($val_cnt < 1) {
      $min_val = $e;
      $max_val = $e;
   }
   $min_val = $e if $e < $min_val;
   $max_val = $e if $e > $max_val;
   if ($estimate) {
      $sum_val += $e;
      $ssq_val += $e*$e;
   }
   else {
      push(@nums, $e);
   }
   $val_cnt++;
}
sub calc_stats() {
   if ($val_cnt > 1) {
      if ($estimate) {
         $avg_val = $sum_val / $val_cnt;
         my $var_est = ($ssq_val - ($sum_val*$sum_val/$val_cnt))/$val_cnt;
         $sigma = sqrt($var_est);
      }
      else {
#         @nums = sort absnums @nums;
         my $e       = 0.0;
         my $d       = 0.0;
         my $sum_e   = 0.0;
         foreach $e (@nums) {
           $sum_e += $e;
         }
         $avg_val = $sum_e / $val_cnt;
         my $sumdev_val = 0.0;
         my $ssqdev_val = 0.0;
         foreach $e (@nums) {
            $d = $e - $avg_val;
            $sumdev_val += $d;
            $ssqdev_val += $d*$d;
         }
         my $variance = ($ssqdev_val
	       - ($sumdev_val*$sumdev_val/$val_cnt))/($val_cnt - 1);
         $sigma = sqrt($variance);
      }
   }
}
# lets us work with smaller values first to lower roundoff error
sub absnums { abs($a) <=> abs($b) }
