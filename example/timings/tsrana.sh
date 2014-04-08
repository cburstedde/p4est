#! /bin/sh

# Pipe in the contents of one or more log files written by p8est_tsearch
# We create three plots: one each for uniform, fractal, and localized ref.

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

set term postscript color

set logscale xy
set key left

set xlabel "points / core"
set ylabel "runtime [1s]"

# Fields: procs level shift refinelocal elements points Search_1 Search_N

set macros

# the @ignore variable is set below (using lazy evaluation)"
# except the replot command does not pick up the variable changes
ign262 = "(@ignore || \$6 == 262144)"
onl262 = "(@ignore || \$6 != 262144)"

# ignore the runs with --refine-local altogether
set output "tsrana-uniform.eps"
set title "Uniform refinement"
ignore = "(\$3 != 0 || \$4 != 0)"
plot	"$TMPF" using (@ign262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1", \
	"$TMPF" using (@onl262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1, P=262144", \
	"$TMPF" using (@ignore ? 1/0 : \$6 / \$1):(\$8) title "Search_N", \
	1e-3 * x title "linear scaling" with lines, \
	1 title "1 second"

set output "tsrana-fractal.eps"
set title "Fractal refinement"
ignore = "(\$3 != 4 || \$4 != 0)"
plot	"$TMPF" using (@ign262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1", \
	"$TMPF" using (@onl262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1, P=262144", \
	"$TMPF" using (@ignore ? 1/0 : \$6 / \$1):(\$8) title "Search_N", \
	1e-3 * x title "linear scaling" with lines, \
	1 title "1 second"

set output "tsrana-local.eps"
set title "Localized refinement"
ignore = "(\$3 != 4 || \$4 != 1)"
plot	"$TMPF" using (@ign262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1", \
	"$TMPF" using (@onl262 ? 1/0 : \$6 / \$1):(\$7) title "Search_1, P=262144", \
	"$TMPF" using (@ignore ? 1/0 : \$6 / \$1):(\$8) title "Search_N", \
	1e-3 * x title "linear scaling" with lines, \
	1 title "1 second"

EOF
