#! /usr/bin/perl

use strict;

my ($dim, $finalt, $deltat, $order);
my (%ids, %fns, $id, $filename, $fh);
my ($gh, $first, $using);

$dim = 0;
$finalt = 0.;
$deltat = 0.;
$order = 0;
while (<>) {
  if (m/^\[p4est\]\s+Dimension is (\d+)/) {
    $dim = $1;
  }
  if (m/^\[p4est\]\s+finaltime\s+([\d\.eE+-]+)/) {
    $finalt = $1;
  }
  if (m/^\[p4est\]\s+deltat\s+([\d\.eE+-]+)/) {
    $deltat = $1;
  }
  if (m/^\[p4est\]\s+rkorder\s+(\d+)/) {
    $order = $1;
  }

  if (m/^\[p4est (\d+)\] T (.+) I (\d+) X (.+) V (.+)$/) {
    $id = $3;

    # check for existing file handle
    $fh = $ids{$id};
    if (!$fh) {

      # this is the first occurrence of this particle
      $filename = sprintf "I%06d.txt", $id;
      open $fh, ">", $filename or die "Open file for $id";
      print $fh "# P T I X3 V3\n" or die "Print first line for $id";

      $ids{$id} = $fh;
      $fns{$id} = $filename;
    }

    # print the line into appropriate file
    printf $fh "%d %g %d %s %s\n", $1, $2, $3, $4, $5 or die "Print for $id";
  }
}

# begin writing gnuplot script
open $gh, ">", "G.gnuplot";
print $gh "set term x11 size 500,500\n";
print $gh "set size ratio -1\n";
print $gh "set xtics .2,.2 offset -2\n";
print $gh "set ytics .2,.2 offset 1.5\n";
print $gh "set zrange [0:1]\n";
#print $gh "set xrange [0:1]\nset yrange [0:1]\n";
#print $gh "set xrange [0:1]\nset yrange [0:1]\nset zrange [0:1]\n";
print $gh "set view 65,30\n";
print $gh "set xyplane at .3\n";
print $gh "set key off\n";
print $gh "set title \"D=$dim T=$finalt dt=$deltat rk=$order\"\n";
if ($dim == 2) {
  print $gh "plot \\\n";
  $using = "4:5";
}
else {
  print $gh "splot \\\n";
  $using = "4:5:6";
}

# close data files
$first = "";
foreach $id (sort keys %ids) {
  close $ids{$id};

  $filename = $fns{$id};
  print $gh "$first  \"$filename\" using $using with lines";

  $first = ",\\\n";
}
print $gh "\n";

print $gh "set term tikz color solid tightboundingbox\n";
print $gh "set output \"G.tikz\"\n";
print $gh "set title\n";
print $gh "replot\n";

# finish writing gnuplot script
print $gh "pause -1\n";
close $gh;
