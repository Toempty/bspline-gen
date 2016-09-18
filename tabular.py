
from sys import stdout

def print_matrix( m, sep=" ", fp=stdout ):
	"""
	Print a matrix of cells convertible to string such that columns
	are aligned and fields are right-justified.
	"""
	# Assume m is row-major
	NC = len(m[0])
	column_widths = [ max([len(str(m[r][c])) for r in range(len(m)) ]) for c in range(NC) ]
	formats = [ "{{: >{}s}}".format(w) for w in column_widths ]
	for row in m:
		print( sep.join([formats[c].format(str(row[c])) for c in range(NC)]), file=fp )

