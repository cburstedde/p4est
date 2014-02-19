#! /bin/sh

# Pipe in the contents of one or more log files written by p{4,8}est_timings

# Find path for the content filter script
TIMANAWK=`echo $0 | sed -e 's/\.sh$/.awk/'`
if test ! -x "$TIMANAWK" ; then
	echo "No executable $TIMANAWK"
	exit 1
fi

TMPF=`tempfile --prefix=timana`
#echo "Temp file: $TMPF"

$TIMANAWK > $TMPF

gnuplot <<EOF

set logscale xy
set xlabel "Elements"
set ylabel "Runtime / \#cores [1s]"
set key left

set output "timana.eps"
set term postscript color

# Fields 4..9: Refine Balance Partition Ghost Nodes Lnodes
plot "$TMPF" using (\$3 / \$1):(\$5) title "Balance", \
	1e-5 * x title "ideal scaling" with lines
plot "$TMPF" using (\$3 / \$1):(\$9) title "Lnodes", \
	2e-5 * x title "ideal scaling" with lines
# plot "$TMPF" using (\$3):(\$5 / (\$3 * 1e-6 / \$1)) title "Balance"
# plot "$TMPF" using (\$3):(\$9 / (\$3 * 1e-6 / \$1)) title "Lnodes"

EOF
