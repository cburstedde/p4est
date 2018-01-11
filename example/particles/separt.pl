#! /usr/bin/perl

use strict;

my (%ids, $id, $filename, $fh);

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
    }

    # print the line into appropriate file
    printf $fh "%d %g %d %s %s\n", $1, $2, $3, $4, $5;
  }
}

# close files
foreach $fh (keys %ids) {
  close $fh;
}
