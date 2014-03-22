#!/bin/sh

command -v perf >/dev/null 2>&1 || { echo >&2 "Could not find perf.  Aborting."; exit 1; }

for test in L1-dcache-load-misses L1-dcache-store-misses LLC-load-misses LLC-store-misses LLC-prefetch-misses dTLB-load-misses dTLB-store-misses dTLB-prefetch-misses iTLB-load-misses branch-load-misses; do
  for type in nodes lnodes base; do
    if [ $type = "nodes" ]; then
      othertype="lnodes"
    elif [ $type = "base" ]; then
      othertype="lnodes --skip-nodes"
    else
      othertype="nodes"
    fi
    for i in 0 1 2 3 4 5 6 7; do
      perf stat --event=$test ./p8est_timings -c unit --skip-$othertype --level $i 2> p8.$test.$type.$i > /dev/null
    done
  done
done

