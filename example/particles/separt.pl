#! /usr/bin/perl

use strict;

my (%ids, %fns, $id, $filename, $fh);
my ($gh, $first);

while (<>) {
  if (m/^\[p[48]est (\d+)\] T (.+) I (\d+) X (.+) V (.+)$/) {
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
#print $gh "set xrange [0:1]\nset yrange [0:1]\n";
print $gh "set size ratio -1\n";
print $gh "plot \\\n";

# close data files
$first = "";
foreach $id (sort keys %ids) {
  close $ids{$id};

  $filename = $fns{$id};
  print $gh "$first  \"$filename\" using 4:5 with lines";

  $first = ",\\\n";
}
print $gh "\n";

# finish writing gnuplot script
print $gh "pause -1\n";
close $gh;
