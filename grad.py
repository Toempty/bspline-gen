
# Generate R code to calculate the gradient of the quadratic form:
# \Theta^T \Omega \Theta
# ...where \Omega is a (necessarily square) matrix and \Theta is a vector of
# parameters with respect to each of which the gradient is sought.
# All indices are 0-based.

def _grad_comp( x:"0-based parameter vector index", n:"dimension", mat="M", param="X" ):
	COLSUM = " + ".join(["{}[{}]*{}[{},{}]".format(param,i,mat,i,x)
			for i in range(n) if i != x ])
	ROWSUM = " + ".join(["{}[{}]*{}[{},{}]".format(param,j,mat,x,j)
			for j in range(n) if j != x ])
	QUADTERM = "2*{}[{}]*{}[{},{}]".format(param,x,mat,x,x)
	return "{} + {} + {}".format( COLSUM, ROWSUM, QUADTERM )

def grad( dim, vec:"parameter vector name"="X", mat:"matrix name"="M" ):
	return "grad <- function({},{}) {{ c(\n{}) }}".format( vec,mat,
		",\n".join([ _grad_comp( i, dim, mat, vec ) for i in range(dim)]) )

if __name__=="__main__":
	from sys import argv
	DIM = int(argv[1])
	VEC =     argv[2] if len(argv) > 2 else "X"
	MAT =     argv[3] if len(argv) > 3 else "M"
	print( grad( DIM, VEC, MAT ) )

