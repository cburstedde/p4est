#! /bin/sh

# Pipe in the contents of one or more log files written by p{4,8}est_timings
# Optionally let the first and second argument be --title "Title of plot"

# Find path for the content filter script
TIMANAWK=`echo $0 | sed -e 's/\.sh$/.awk/'`
if test ! -x "$TIMANAWK" ; then
	echo "No executable $TIMANAWK"
	exit 1
fi

TITLECOMM=
if test "x$1" = "x--title" ; then
	TITLECOMM="set title \"$2\""
	shift
	shift
fi

TMPF=`tempfile --prefix=timana`
#echo "Temp file: $TMPF"

$TIMANAWK $@ > $TMPF

gnuplot <<EOF

set logscale xy
set key left
set style data points

#set output "timana.eps"
set term postscript color
$TITLECOMM

# Fields 4..9: Refine Balance Partition Ghost Nodes Lnodes

set xlabel "\#elements / core"
set ylabel "runtime [1s]"

set output "timana-balance.eps"
plot "$TMPF" using (\$3 / \$1):(\$5) title "Balance" pt 5, \
	5e-6 * x title "linear scaling" with lines, \
	8e-4 * x**(2./3.) title "scaling with x**(2/3)" with lines, \
	1 title "1 second"

set output "timana-ghost.eps"
plot "$TMPF" using (\$3 / \$1):(\$7) title "Ghost" pt 5, \
	2e-4 * x**(2./3.) title "scaling with x**(2/3)" with lines, \
	1 title "1 second"

set output "timana-nodes.eps"
plot "$TMPF" using (\$3 / \$1):(\$8) title "Nodes" pt 5, \
	2e-5 * x title "linear scaling" with lines, \
	1 title "1 second"

set output "timana-lnodes.eps"
plot "$TMPF" using (\$3 / \$1):(\$9) title "Lnodes" pt 5, \
	2e-5 * x title "linear scaling" with lines, \
	1 title "1 second"

# set xlabel "Elements"
# set ylabel "Runtime / \#cores [1s]"
# plot "$TMPF" using (\$3):(\$5 / (\$3 * 1e-6 / \$1)) title "Balance"
# plot "$TMPF" using (\$3):(\$7 / (\$3 * 1e-6 / \$1)) title "Ghost"
# plot "$TMPF" using (\$3):(\$8 / (\$3 * 1e-6 / \$1)) title "Nodes"
# plot "$TMPF" using (\$3):(\$9 / (\$3 * 1e-6 / \$1)) title "Lnodes"

EOF
