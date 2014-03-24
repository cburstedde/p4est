#!/bin/bash

command -v perf >/dev/null 2>&1 || { echo >&2 "Could not find perf.  Aborting."; exit 1; }

STATS="instructions,cache-misses,branch-misses"
REPCOUNT=30

SPACESTATS=`echo $STATS | sed 's/,/ /g'`
echo "type level time $SPACESTATS" > p8.perf

for i in 5 6 7 8; do
  for type in "base" "nodes" "lnodes"; do
    if [ $type = "nodes" ]; then
      othertype="lnodes"
    elif [ $type = "base" ]; then
      othertype="lnodes --skip-nodes"
    else
      othertype="nodes"
    fi
    perf stat -o p8.$type.$i -e "$STATS" -r $REPCOUNT ./p8est_timings \
      -c unit --skip-$othertype --level $i > p8.$type.$i.stdout
    runtime=`egrep "[0-9.]+ seconds time elapsed" p8.$type.$i | awk '{print $1;}'`
    accum=$runtime
    for stat in $SPACESTATS; do
      count=`egrep "[0-9,]+ $stat" p8.$type.$i | awk '{print $1;}' | sed 's/,//g'`
      accum="$accum $count"
    done
    if [ $type = "base" ]; then
      baseaccum=$accum
    else
      read -a basearray <<< "$baseaccum"
      read -a array <<< "$accum"
      for j in "${!array[@]}"; do
        val=`echo "${array[j]}-${basearray[j]}" | bc -l`
        array[j]=$val
      done
      diff="${array[*]}"
      echo "$type $i $diff" >> p8.perf
    fi
  done
done

