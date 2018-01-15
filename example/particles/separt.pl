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
  if (m/Into p([48])est_new/) {
    if ($1 == "4") {
      $dim = 2;
    }
    elsif ($1 == "8") {
      $dim = 3;
    }
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
      open $fh, ">", $filename;
      print $fh "# P T I X3 V3\n";

      $ids{$id} = $fh;
      $fns{$id} = $filename;
    }

    # print the line into appropriate file
    printf $fh "%d %g %d %s %s\n", $1, $2, $3, $4, $5;
  }
}

# begin writing gnuplot script
open $gh, ">", "G.gnuplot";
print $gh "set term x11 size 500,500\n";
print $gh "set size ratio -1\n";
#print $gh "set xrange [0:1]\nset yrange [0:1]\n";
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

# finish writing gnuplot script
print $gh "pause -1\n";
close $gh;
