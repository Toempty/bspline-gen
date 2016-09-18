
from math import sqrt

def _prune( coef ):
	"""
	Remove 0.0's from the tail of a list.
	"""
	i = len(coef)
	while i > 0 and coef[i-1] == 0:
		i -= 1
	return coef[:i] if i > 0 else [0.0,]


class Polynomial(object):
	"""
	Polynomials represented as lists of coefficients in ASCENDING order of
	degree, so that degree(p) == len(self.C)-1. For example:
	      x-2: Polynomial([-2,1   ])
	ax^2+bx+c: Polynomial([ c,b,a ])
	"""

	@staticmethod
	def fromAscendingCSV( csv ):
		return Polynomial( [ float(v) for v in csv.split(',') ] )

	@staticmethod
	def haar( i, k:"interval boundaries", x:"", leftcont=False ):
		"""
		All Haar functions are continuous from right (closed on the left, open
		on the right), but...
		Make special provision for the rightmost spline since we are following
		de Boor's (arbitrary) convention of defining piecewise-polynomials that
		are continuous from the right. We need the rightmost polynomial to be
		evaluate-able AT it's upper bound, as well, so that the function is
		defined on the CLOSED basic interval.
		"""
		return Polynomial([1.0,]) if k[i] <= x and ( x < k[i+1] or (leftcont and x<=k[i+1]) ) else Polynomial([0.0,])

	def __init__( self, coeffs, scalar=None ):
		self.C = coeffs
		if scalar:
			self.scale( scalar )

	def __str__( self ):
		return ','.join([ "{:.3e}".format(coef) if coef != 0.0 else "0" for coef in self.C ])

	def __len__( self ):
		return len(self.C)

	def deg( self ):
		return len(self)-1

	def __add__( self, other ):
		p = self.C
		q = other.C if isinstance( other, Polynomial ) else [ float(other), ]
		if len(p) <= len(q):
			s = [ q[i] for i in range(len(q)) ]
			for i in range(len(p)):
				s[i] += p[i]
			return Polynomial( s )
		else:
			s = [ p[i] for i in range(len(p)) ]
			for i in range(len(q)):
				s[i] += q[i]
			return Polynomial( s )

	def __mul__( self, other ):
		p = self.C
		q = other.C if isinstance( other, Polynomial ) else [ float(other), ]
		c = [ 0 for i in range(len(p)+len(q)-1) ]
		for i in range(len(p)):
			for j in range(len(q)):
				c[i+j] += p[i]*q[j]
		return Polynomial( _prune(c) )


	def __call__( self, x ):
		p = 1.0 # = x^0
		s = 0.0
		for c in self.C:
			s += c*p
			p *= x
		return s


	def scale( self, a ):
		"""
		Mutates self (unlike the multiplication operator, __mul__).
		"""
		for i in range(len(self.C)):
			self.C[i] *= a
		return self


	def deriv( self ):
		"""
		Calculate the derivative of polynomal p.
		"""
		p = self.C
		return Polynomial( _prune( [ float(i)*p[i] for i in range(1,len(p)) ] ) if len(p) > 1 else [0,] )


	def integrate( self, a=0, b=0 ):
		"""
		Calculate the derivative of polynomal p.
		"""
		p = Polynomial( [0.0,] + [ self.C[i]/float(i+1) for i in range(len(self)) ] )
		if a == b:
			return p
		return p(b) - p(a)

	def roots( self, interval=None ):
		"""
		Only supports quadratic roots for now.
			 0  1  2
		p = [c, b, a, 0, 0, ... ]
		[ -b +/- sqrt(b*b-4ac) ] /2a
		"""
		p = self.C
		if not all([c==0 for c in p[3:]]):
			raise NotImplementedError( "expected quadratic poly, have ({0})".format(p) )
		if interval and (interval[1] <= interval[0]):
			raise RuntimeError( "intervals upper bound ({1}) <= lower bound ({0})".format(*interval) )
		RAD = p[1]*p[1] - 4.0*p[2]*p[0]
		DEN = (2.0*p[2])
		SQR = sqrt( abs(RAD) )
		if RAD >= 0.0:
			r1 = (-p[1]-SQR)/DEN
			r2 = (-p[1]+SQR)/DEN
			if interval:
				return tuple( filter( lambda r:(interval[0] <= r) and (r <= interval[1]), (r1,r2) ) )
			else:
				return ( complex(r1,0) , complex(r2,0) )
		else:
			if interval:
				return () # there are no real roots so they can't be in the interval.
			else:
				return ( complex( -p[1]/DEN, -SQR/DEN ), complex( -p[1]/DEN, +SQR/DEN ) )


############################################################################
# Unit testing

if __name__=="__main__" and __debug__:
	import sys
	p = Polynomial.fromAscendingCSV( sys.argv[1] )
	print( p.deriv().integrate() )
	print( p.integrate().deriv() )
	print( p.integrate( 0.0, 1.0 ) )

