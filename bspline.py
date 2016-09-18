
from poly import Polynomial

class BSpline(object):
	"""
	Calculates the (matrix of) piecewise polynomials of a basis spline on an
	augmented set of intervals and provides indexed access to them.
	It precomputes the matrix of polynomials.
	Each spline basis function is a piecewise polynomial that covers one
	more interval that its degree. (The empty intervals at the ends formed
	by repeating the boundary knots count as intervals in this sense.)
	Below is a schematic of a cubic basis spline for 4 (real) intervals.

	Ki:       0    1    2    3    4 
	Ui: 0 1 2 3    4    5    6    7 8 9 10
	    |D|D|D|____|____|____|____|D|D|D|
	     0|0|0|0000|    |    |    |
	       1|1|1111|1111|    |    |
	         2|2222|2222|2222|    |
	          |3333|3333|3333|3333|
	          |     4444|4444|4444|4
	          |          5555|5555|5|5
	          |_______________6666|6|6|6
	"""
	def __init__( self, knot_seq, degree=3 ):
		# Preserve the distinction between the actual and...
		self.K = sorted([float(k) for k in knot_seq])
		assert all( [ self.K[i-1] < self.K[i] for i in range(1,len(self.K))] ),\
			"knots must be strictly increasing (in this code)"
		# ...augmented knot sequence which has dummy knots at both ends.
		self.U = [self.K[0],]*degree + self.K + [self.K[-1],]*degree
		self.deg = degree
		# Precompute the piecewise-polynomials
		self.polys = [
			[ self._poly_r( si, degree, (self.K[ii-1]+self.K[ii])/2.0 ) for ii in range(1,len(self.K)) ]
			for si in range(len(self)) ]
		# polys should be a diagonal banded matrix of polynomials with zero
		# polynomials in the upper-right and lower-left corners.

	def __len__( self ):
		return self.deg + (len(self.K)-1)

	def __call__( self, spline_index, interval_index ):
		"""
		Return the polynomial for the given spline in the given interval
		(both 0-based indices).
		"""
		return self.polys[ spline_index ][ interval_index ]


	def _poly_r( self, i:"spline_index", d, x ):
		"""
		Generate the ith B-spline basis function via the recurrence
		relation.
		The ith spline is that whose support begins in the ith interval
		This is method NOT the most efficient; it is intended to be the
		clearest by virtue of directly implementing the Cox-de Boor
		recurrence relation.
		"""
		u = self.U
		if d == 0:
			return Polynomial.haar( i, u, x, i >= len(u)-2 )
		else:
			l = 1.0/(u[i+d  ]-u[i  ]) if u[i+d  ] != u[i  ] else 0.0
			r = 1.0/(u[i+d+1]-u[i+1]) if u[i+d+1] != u[i+1] else 0.0
			return \
				Polynomial( [-u[i],     1], l) * self._poly_r( i,   d-1, x ) + \
				Polynomial( [ u[i+d+1],-1], r) * self._poly_r( i+1, d-1, x )

# Unit test...
if __name__=="__main__" and __debug__:
	"""
	The output of this Python script can be plotted in R using:

	A <- commandArgs( trailingOnly=T );
	dat <- read.table( A[[1]], header=F );
	spline.dat <- split( dat[,2:3], dat[,1]);
	png( A[[2]] );
	plot( range(dat[,2]), range(dat[,3]), xlab="X", ylab="Y" );
	for( i in 1:length(spline.dat) ) lines( spline.dat[[i]], col=rainbow(length(spline.dat))[i], lwd=2 )
	dev.off()
	"""
	from sys import argv,stdout,stdin
	K = sorted([ float(v) for v in argv[1].split(',') ])
	spline = BSpline( K, int(argv[2]) )
	DELTA = (max(K)-min(K))/100.0
	for si in range(len(spline)):
		# For each interval...
		for ival in range(len(K)-1):
			x = K[ ival ]
			ival_ub = K[ ival+1 ]
			# Precompute the polynomial that will be relevant on this interval
			# using any value of X in the interval.
			poly = spline( si, ival )
			print( "#{}/{}\t{}".format( si, ival, poly ) )
			while x < ival_ub:
				print( "{}\t{:.2f}\t{:.3f}".format( si, x, poly(x) ) )
				x += DELTA

