#! /usr/bin/awk -f

/New p(4|8)est/ {
	if ($3 == "p4est") { dim = 2 }
	else if ($3 == "p8est") { dim = 3 }
	#print "We are in dimension " dim
}
/^\[p4est\] Processors [[:digit:]]+ level/ {
	procs = $3
	level = $5
	elems = 0
	#printf "Procs %d level %d\n", procs, level
}
/^\[p4est\] .* balanced to/ {
	elems = $9
}
/^\[p4est\] Summary =/ {
	# Output runtimes for: Refine Balance Partition Ghost Nodes Lnodes
	printf "%d %d %d %g %g %g %g %g %g\n", procs, level, elems,
		$5, $6, $33, $34, $35, $38
}
