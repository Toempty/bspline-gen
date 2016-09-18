
"""
Generate the "omega" matrix required to find the smoothest basis-spline
curve through a set of intervals.
"""

from sys import argv
from tabular import print_matrix
from bspline import BSpline

K = sorted([ float(v) for v in argv[1].split(',') ])
spline = BSpline( K, int(argv[2]) )

DS    = [ [ spline(si,ii).deriv() for ii in range(len(K)-1) ] for si in range(len(spline)) ]
DDS   = [ [             p.deriv() for p  in s               ] for s  in DS ]

# Following calculates
# \lbrace \Omega \rbrace_{jk} = \int B_j^{\prime\prime}(t) B_k^{\prime\prime}(t) dt
OMEGA = [
	[ sum([ (p[i]*q[i]).integrate( K[i], K[i+1] ) for i in range(len(K)-1) ])
	for p in DDS ]
	for q in DDS ]

print_matrix( OMEGA, sep="\t" )

