#! /usr/bin/awk -f

BEGIN {
	print "# procs level shift refinelocal elements points Search_1 Search_N"
}
/refine-local:/ {
	if ($3 == "true") { refinelocal = 1 } else {refinelocal = 0 };
}
/New p(4|8)est/ {
	if ($3 == "p4est") { dim = 2 }
	else if ($3 == "p8est") { dim = 3 }
	#print "We are in dimension " dim
}
/^\[p4est\] Processors [[:digit:]]+ level/ {
	procs = $3
	level = $5
	shift = $7
	points = $9
	elems = 0
	#printf "Procs %d level %d shift %d points %d\n", procs, level, shift, points
}
/^\[p4est\] .* balanced to/ {
	elems = $9
}
/^\[p4est\] Summary =/ {
	# Output runtimes for: Search_1 Search_N
	printf "%d %d %d %d %d %d %g %g\n",
	       procs, level, shift, refinelocal, elems, points, $9, $10
}
