#! /bin/sh

# Pipe in the contents of one or more log files written by p8est_tsearch

# Find path for the content filter script
TSRANAWK=`echo $0 | sed -e 's/\.sh$/.awk/'`
if test ! -x "$TSRANAWK" ; then
	echo "No executable $TSRANAWK"
	exit 1
fi

#TMPF=tsrana.tmp
TMPF=`tempfile --prefix=tsrana`
#echo "Temp file: $TMPF"

$TSRANAWK $@ > $TMPF

gnuplot <<EOF

set logscale xy
set key left

set output "tsrana.eps"
set term postscript color

set xlabel "points / core"
set ylabel "runtime [1s]"

# Fields: procs level shift refinelocal elements points Search_1 Search_N

set macros
#ignore = "(\$4 != 0 || \$6 / \$1 < 100)"
#ignore = "0"
# ignore the runs with --refine-local
ignore = "(\$4 != 0)"
ign262 = "(\$4 != 0 || \$6 == 262144)"
onl262 = "(\$4 != 0 || \$6 != 262144)"

plot	"$TMPF" using (@ign262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1", \
	"$TMPF" using (@onl262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1, P=262144", \
	"$TMPF" using (@ignore ? 1/0 : \$6 / \$1):(\$8) title "Search_N", \
	1e-3 * x title "linear scaling" with lines, \
	1 title "1 second"

EOF
