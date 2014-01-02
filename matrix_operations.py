''' 
James Ma


run program:
------------
python matrix_operations.py
'''
import math,random,sys,csv,operator
import numpy as np
import scipy.linalg
import scipy as sp



#-------------------------------------------
# Helper Functions
#-------------------------------------------
# Generate simulated matrix
def gen_matrix(size=100):
	'''
	The point of inflection of x**2 + 3x + 4 is:
	2x + 3 = 0
	x = -3/2
	
	input:
	------
	size : number of samples
	
	returns:
	--------
	x : x matrix
	y : y vector
	'''
	x = range(-50, -2)
	x.append(-3.0/2)
	x.extend(range(0,51))
	
	x = np.array(x)
	y = x**2 + 3*x + 4 + np.random.normal(0,.1,len(x))
	
	x_mat = np.matrix([x**2,x, np.ones(len(x))]).T
	return x_mat,np.matrix(y).T
	

# helper functions
# swap rows in matrix
def row_swap_matrix(mat,i,j):
	cur_val = mat[i,:].copy()
	mat[i,:] = mat[j,:] 
	mat[j,:] = cur_val

# swap items in list
def row_swap_list(list,i,j):
	cur_val = list[i]
	list[i] = list[j]
	list[j] = cur_val


# Convert the diagonals to 1
def convert_diag(A, val=1):
	# dimension nxm
	n = A.shape[0]
	
	# all diags to 1
	for i in xrange(n):
		A[i,i] = val
	return A


# L2 -norm
l2_norm = lambda x: math.sqrt(sum([i**2 for i in x]))


# create array of zeros
# x param : size of array
zeros = lambda x : np.array([0 for i in xrange(x)])





#------------------------------------------------
# Parts a-n
#------------------------------------------------

# (a) Scipy functions lu_factor and lu_solve to find B from your simulated data
sp_lu_solve = lambda x,y : sp.linalg.lu_solve(sp.linalg.lu_factor(x),y)





# (b) function for LU Decomposition
def  LUdecomp(A, z=0, r=0):
	'''
	LU Decomposition using design matrix X = L*U
	This method uses Partial Pivoting for better numerical accuracy.
	Additionally, the matrix A is replaced with the LU decomposition.
	
	input:
	------
	A : design matrix to be factored
	z : permutation vector
	r : -1 or +1 depending on number of row swaps
	
	returns:
	--------
	A : LU decomposition 
	z : permutation vector
	r : -1 or +1 depending on number of row swaps
	'''
	# Make sure it's a square matrix
	assert A.shape[0] == A.shape[1], "Need a Square Matrix"
	
	# permutation vector, check if need to create it
	if z:
		pass
	else:
		z = np.arange(A.shape[0])
	
	# dimension of mxm matrix
	m = A.shape[0]
	
	# main loop
	for j in xrange(m):
		
		# find uij's
		for i in xrange(j+1):
			A[i,j] = A[i,j] - A[i,:i].dot(A[:i,j])
		
		# assign diagonal
		ujj = A[j, j]
		# store each lij's - later to find max
		lij_star = [ujj]
		
		# find lij_stars 
		for i in xrange(j+1, m):
			# lij
			A[i,j] = abs(A[i,j] - A[i,:j].dot(A[:j, j]))
			lij_star.append(A[i,j])
	

		# get max index
		index = max(enumerate(lij_star),key=operator.itemgetter(1))[0] + j
		
		# swap if j+1 ...j+k is larger than j
		if (index > j):
			# row swap matrix rows
			row_swap_matrix(A, j, index)
		
			# row swap list rows
			row_swap_list(z, j, index)
			
			# increment r
			r += 1
			
		# set diagonal ujj
		ujj = A[j, j]
		
		# divide by largest element
		A[j+1:,j] = A[j+1:,j] / ujj
			
	# if r is even then r = 1
	# else r = -1
	if (r % 2 == 0) : r = 1
	else : r = -1
	
	return A,z,r




# c. LUSub
def LUSub(B,Y,z):
	'''
	LU Substitution taking an LU matrix B, and solves for X:
	Y = B*X
	
	It replaces the Y vector with the solution.
	
	input:
	------
	B : LU decomposition
	Y : Y
	z : permutation vector
	
	returns:
	--------
	Y : solution 
	'''
	
	assert (B.shape[0] == B.shape[1]), "B needs to be a Square Matrix"
	assert (Y.shape[0] == B.shape[0]), "Y needs to have the same number of rows"
	
	# index table - to keep track of permutations
	it = {j:i for i,j in enumerate(z)}
	
	# length of z
	n = len(z)
	
	# forward substitution - using L matrix
	for i in xrange(1,n):
		update = 0
		idx = it[i]
		for j in xrange(i):
			idx_y = it[j]
			update += B[idx,j] * Y[idx_y]
		# update betas
		Y[idx] = Y[idx] - update
	
	# backward substitution - using U matrix
	for i in xrange(n-1, -1, -1):
		update = 0
		idx = it[i]
		for j in xrange(n-1, i,-1):
			idx_y = it[j]
			update += B[idx,j] * Y[idx_y]
		# update betas
		Y[idx] = (Y[idx] - update) / float(B[idx,idx])
	
	return Y




# d. Show the LU Decoposition
def print_lu(x,y):
	'''
	Prints the LU decomposition -> LU sub solution.
	'''
	# store x.T * x
	x_lu = x.T*x
	
	print "The original A is:"
	print x_lu
	
	# LU decomposition
	a,z,r = LUdecomp(x_lu)
	print "\nThe LU decomposition is: "
	print x_lu
	
	
	# store x.T*y
	y1 = x.T*y
	
	print "\nThe original Y is:"
	print y1
	
	LUSub(x_lu, y1, z)
	print "\nThe Beta Estimates / Transformed Y: "
	print y1




# e. Matrix Inversion using LU Decomposition
def InvertMatrix(A):
	'''
	Method inverts matrix A using the LU decomposition.
	
	input:
	------
	A : design matrix
	
	returns:
	--------
	inv(A) : inverse of A
	'''
	assert (A.shape[0] == A.shape[1]), "Need a Square Matrix"
	
	# Initialize empty matrix to store inverse
	mat = []
	
	# n rows, cols
	n = A.shape[0]
	
	# Identity matrix for inverse
	id = np.matrix(np.identity(n))
	
	# LU decomposition
	x_lu,z,r = LUdecomp(A)
	
	# use LU substitution to find inverse
	for i in xrange(n):
		# solve for each column in identity matrix
		b = LUSub(x_lu, id[:,i], z)
		mat.append(b)
	return np.hstack((mat))





# f. Using scipy's Cholesky Factorization
sp_cho_solve = lambda x,y: sp.linalg.cho_solve(sp.linalg.cho_factor(x.T*x),x.T*y)





# g. Cholesky Decomposition - lower triangular
def CholeskyDecomp(A):
	'''
	Cholesky decomposition on design matrix A.  
	This method replaces matrix A with the lower triangular part of
	the decomposition.
	
	input:
	------
	A : design matrix
	
	returns:
	--------
	Chol(A) : Cholesky decomposition of A  
	'''
	assert (A.shape[0] == A.shape[1]), "Needs a Square Matrix"
	
	
	# number of rows, cols
	n = A.shape[0]
		
	for k in xrange(n):
		# lkk
		lkk = math.sqrt(A[k,k] - A[k,:k].dot(A[k,:k].T))
		A[k,k] = lkk
		
		# fill in non-diagonals / lower triangle
		for i in xrange(k+1, n):
			A[i,k] = 1.0 / lkk * (A[i,k] - A[i,:k].dot(A[k,:k].T))
	
	return A


# Cholesky substituion - similar to LUSub
# lower triangular matrix from Cholesky
def ChoSub(B,Y):
	'''
	Cholesky Substitution to solve for X in :
	B*X = Y
	
	It replaces the Y vector with the solution.
	
	input:
	------
	B : Cholesky decomposition
	Y : Y
	
	returns:
	--------
	Y : solution 
	'''
	assert (B.shape[0] == B.shape[1]), "B needs to be a Square Matrix"
	assert (Y.shape[0] == B.shape[0]), "Y needs to have the same number of rows"
	
	
	# length of z
	n = B.shape[0]
	
	# forward substitution - using L matrix
	for i in xrange(n):
		update = 0
		for j in xrange(i):
			update += B[i,j] * Y[j]
		# update betas
		Y[i] = (Y[i] - update) / float(B[i,i])
	
	# transpose for upper triangle
	B = B.T
	
	# backward substitution - using L.T matrix
	for i in xrange(n-1, -1, -1):
		update = 0
		for j in xrange(n-1, i,-1):
			update += B[i,j] * Y[j]
		
		# update betas
		Y[i] = (Y[i] - update) / float(B[i,i])
	
	return Y




# h. Find B using CholeskyDecomp and ChoSub
def b_cho_lu(x, y):
	'''
	Prints the solution.
	Cholesky(x) -> Cholesky substituion(x,y) 
	'''
	# store x.T*x
	A = x.T*x
	
	print "The Original Matrix A:"
	print A
	
	CholeskyDecomp(A)
	print "\nThe Intermediate step, Cholesky Decomp:"
	print A
	
	# store x.T*y
	b = x.T*y
	ChoSub(A, b)
	print "\nThe Beta Estimates:"
	print b






# i. Use Scipy's QR to find betas
# returns the q,r Decomposition using scipy
scipy_qr = lambda x : sp.linalg.qr(x)

# Solving for beta using scipy
scipy_qr_solve = lambda q,r,y : np.matrix(sp.linalg.inv(r.T.dot(r)))*r.T*q.T*y



# j. QR Decomposition 
def QRdecomp(A):
	'''
	QR decomposition on design matrix A.
	This method finds orthonormal vectors that span the column
	space of A - matrix Q.
	
	Then solves for matrix R - an upper triangular matrix.
	
	input:
	------
	A : design matrix
	
	returns:
	--------
	Q : orthonormal matrix
	R : upper triangular matrix
	'''
	# number of columns
	n = A.shape[1]

	# copy of matrix
	Q = np.array(A.copy())

	for i in xrange(1,n):
		# keep a copy of xi
		xi = Q[:,i].T.copy()
		
		# find the projections
		for j in xrange(i):
			vi = Q[:,j]
			Q[:,i] -= (xi.dot(vi) / vi.T.dot(vi)) * vi
	
	# norm each column vector
	for i in xrange(n):
		Q[:,i] = Q[:,i] / l2_norm(Q[:,i])
	
	# R part
	R = Q.T*A
	
	return Q,R






# k. QR Solve
def QRsolve(Q,R,Y):
	'''
	This method solves for X using matrices Q and R.
	Q*R*X = Y
	
	
	input:
	------
	Q : orthonormal matrix
	R : upper triangle
	Y : Y
	
	returns:
	--------
	X : solution 
	'''
	rb = Q.T*Y
	
	# index R (mxm)
	m = R.shape[1]
	
	# backward substition to solve betas
	for i in xrange(m-1, -1, -1):	
		update = 0
		
		for j in xrange(m-1, i, -1):
			update += rb[j] * R[i,j]
		# update betas
		rb[i] = (rb[i] - update) / R[i,i]
	
	return rb





# l. QRdecomp and QRsolve to find B
def qr_decomp_solve(x, y):
	'''
	Prints and solves for betas:
	QRdecomp(x) -> QRsolve(q,r,y)
	
	'''
	Q,R = QRdecomp(x)
	print "The QR Decomposition\nQ Matrix:"
	print np.array(Q)
	
	print "\nR Matrix:"
	print R
	
	b = QRsolve(Q,R,y)
	
	print "\nThe Beta Estimates:"
	print b



# m. QR decomposition using Givens rotations
def GivensRotation(A):
	'''
	This method uses the Givens Rotation to zero out elements
	that are below the upper-triangle, R, and implicitly finds
	the matrix Q.T
	
	
	input:
	------
	A : design matrix
	
	returns:
	--------
	Q : orthonormal matrix
	R : upper triangular matrix 
	'''
	# A is an nxm matrix
	n = A.shape[0]
	m = A.shape[1]
	
	# empty list of matrices to compute Q
	gs = []
	
	# intialize r
	R = A.copy()
	
	# initialize Q
	Q = np.identity(n)
	
	# loop to zero out non-upper-triangle entries
	for i in xrange(m):	
		for j in xrange(n-1, i, -1):
			
			# empty identity to store little matrix
			g = np.identity(n)
			
			# x1,x2 in little matrix
			x1 = R[j-1,i]
			x2 = R[j,i]
			
			# denominator
			denom = math.sqrt(x1**2 + x2**2)
			
			# the little matrix
			xs = np.array([[x1 / denom, x2 / denom],[-x2 / denom, x1 / denom]])
		
			# put the little matrix in identity
			g[j-1:j+1,j-1:j+1] = xs
			
			# Q.T = prod(gs)
			Q = Q.dot(g.T)
			
			# R matrix
			R = g.dot(R)
	return Q,R





# n. Use Givens rotation to find beta
def givens_solve(A,y):
	'''
	Using the givens rotation to find Q and R.  
	Then solves for betas using backward substitution.
	'''
	Q,R = GivensRotation(A)
	
	# rb matrix
	rb = Q.T*y
	
	# shape of betas
	m = A.shape[1]
	
	# remake R matrix
	R = R[:m,:m]
	
	# backward substition to solve betas
	for i in xrange(m-1, -1, -1):	
		update = 0
		
		for j in xrange(m-1, i, -1):
			update += rb[j] * R[i,j]
		# update coefficients
		rb[i] = (rb[i] - update) / R[i,i]
	
	return rb[:m]




### main function
# for printing results
def main():
	
	# generate x,y
	x,y = gen_matrix()
	
	
	# (a). Scipy lu_factor and lu_solve to find Betas
	print "(a). Using Scipy's lu_factor and lu_solve\nBeta Estimates:"
	print sp_lu_solve(x.T*x, x.T*y)
	
	
	# (b). LU Decomp
	# show it changes original matrix
	x_lu = x.T*x
	
	p,z,r = LUdecomp(x_lu)
	print "\n\n(b). LU decomposition:"
	print x_lu
	print "\nz vector:"
	print z
	print "\nr:"
	print r
	
	
	# (c). LU substitution
	# store x.T*y to show it saves in y
	y1 = x.T*y
	LUSub(x_lu,y1,z)
	print "\n\n(c) LU substition\nBetas"
	print y1
	
	
	# (d). Find beta from simulated data
	print "\n\n(d). Using LUdecomp and LUSub:"
	print_lu(x,y)
	
	
	# (e). Invert Matrix 
	print "\n\n(e) Write a function using LU to invert matrix\nWritten Inverted Matrix"
	print InvertMatrix(x.T*x)
	print "\nInvert from Numpy:"
	print sp.linalg.inv(x.T*x)
	
	
	# (f) Scipy's Cho_factor and Cho_solve
	print "\n\n(f) Scipy's Cho_factor and Cho_solve\nBeta Estimates:"
	print sp_cho_solve(x, y)
	
	
	# (g). Cholesky Decomposition
	print "\n\n(g) Cholesky Decomposition"
	x_cho = x.T*x
	CholeskyDecomp(x_cho)
	print "Cholesky decomposition for lower triangular part:"
	print x_cho
	

	# (h) Use ChokleskyDecomp and ChoSub to find Betas
	print "\n\n(h) CholeskyDecomp and ChoSub to find Betas"
	b_cho_lu(x,y)
	
	
	# (i). Use Scipy's qr decomposition to find Betas
	print "\n\n(i) Scipy's QR Decomposition to find betas\nScipy QR Beta Estimates"
	x,y = gen_matrix()
	q,r = scipy_qr(x)
	b = scipy_qr_solve(q,r,y)
	print b
	
	

	# (j). QR Decomposition  
	print "\n\n(j) Use QRdecomp"
	q,r = QRdecomp(x)
	print "Q Matrix:"
	print q
	print "\nR Matrix:"
	print r
	
	# (k). QR solve
	print "\n\n(k) Use QRsolve to find Betas"
	print "Beta Estamates"
	print QRsolve(q,r,y)
	
	
	
	# (l) Using functions QRdecomp and QRsolve to find Beta from simulated data
	print "\n\n(l) Use QRdecomp and QRsolve to find Betas"
	qr_decomp_solve(x,y)
	
	# (m). QR Decomposition using Givens Rotation
	print "\n\n(m) Givens rotation to find Q,R:\nQ Matrix"
	q,r = GivensRotation(x)
	print q
	print "\nR Matrix (truncated):"
	print r[:4]
	
	# (n). Using GivensRotation to find Beta
	print "\n\n(n) Use GivensRotation to find Beta\nBeta Estimates:"
	print givens_solve(x,y)



# run main
if __name__ == "__main__":
	main()




