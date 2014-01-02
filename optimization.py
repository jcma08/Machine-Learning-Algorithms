__author__ = "James Ma"
__email__ = "jcma08@gmail.com"


# test function, the derivative
f_x = lambda x : x**3 - 0.165*x**2 + 0.0003993
# true func : 1/4 * x**4 - 0.165/3 * x**3 + 0.0003993*x

# second derivative
f_prime = lambda x : 3*x**2 - 2*0.165*x



# test interval
# true = there exists a root, else false
def test_root(f_x, start, end):
	val = f_x(start) * f_x(end)
	
	if (val < 0): return True
	else : return False



# bisection method
# f_x is the derivative of f_x
def bisect(f_x, start, end, maxiter, es):
	
	if test_root(f_x, start, end): pass
	else: return
	
	# assume f_x(start) * f_x(end) > 0
	midpt = (start + end) / 2
	
	for iter in xrange(maxiter):
	
		# test roots
		func = f_x(start) * f_x(midpt)
		if (func < 0):
			end = midpt
		elif (func > 0):
			start = midpt
		else:
			return midpt
	
		midold = midpt
		midpt = (start + end) / 2
		
		# relative error
		ea = abs(midpt - midold) / midpt * 100
		
		if (ea < es):
			return midpt
	return midpt




# newton-raphson method
def newt_raph(f_x, f_prime, start, maxiter, es):
	
	for iter in xrange(maxiter):
	
		start_old = start
	
		# find root
		start = start - f_x(start) / f_prime(start)
	
		# relative error
		ea = abs(start - start_old) / start * 100
	
		if (ea < es):
			return start
	return start
	


# secant method
def secant(f_x, start, end, maxiter, es):
	
	for iter in xrange(maxiter):
	
		start_old = start
	
		# find root
		start = start - f_x(start)*(start - end) / (f_x(start) - f_x(end))
		
		# relative error
		ea = abs(start - start_old) / start * 100
		
		if (ea < es):
			return start
		
		end = start_old
	return start



import numpy as np

# conjugate gradient
def steepestdescent(f_x, x0, h, e , step):
	xold = x0
	err = e*1.1
	while(err > e):
		s = calcgradient(f_x,xold,h)
		xnew = xold - step*s
		err = abs(f_x(xnew) - f_x(xold))
		xold = xnew
	return xnew

def f1(v):
	return np.sum((v-10)**2+2)




def f1(v):
	return np.sum(v**2) + 2




def calcgradient(f_x, x0, h):
	n = x0.shape[0]
	grad = np.zeros(n)
	for i in xrange(n):
		x1 = np.copy(x0)
		x2 = np.copy(x0)
		x1[i] = x0[i] - h
		x2[i] = x0[i] - h
		grad[i] = (f(x2) - f(x1)) / (2*h)
	return grad


def f2(x):
	return x**3 - 3


def deriv(f,x0,h):
	return (f(x0 + h) - f(x0 - h)) / (2*h)



import numpy.random as npr

def NewPop(N,minD,maxD):
	r = len(bin(,ax([abs(mindD,abs(maxD)))))-1
	
	popmat = np.zeros([N,r])
	
	for popcount in range(N):
		val = npr.random_integers(minD, maxD, 1)[0]
		
		if(val > 0):
			popmat[popcount,0] = 1
			bval = str(bin(abs(val)))
			lbv = len(bval) - 2
			popmat[popcount, -lbv:] = list(bval)[-lbv:]
			popmat[popcount,range(1,r-lbv)] = 0
			
	return popmat



def f(x):
	return 300-(x-20)**2


def fitness1(z):
	return f(z)



def CalcFitness(gvec,w):
	L = len(gvec)
	b10 = 0
	for i in range(1,L):
		b10 = b10 + gvec[i]*2**(L-1-i)
	b10 = b10 * (-1) **(1 - gvec[0])
	return w(b10)


def ApplyRouletteFitness(popmat,w):
	r = popmat.shape[1]
	N = popmat.shape[0]
	newpop = np.empty([N,r])
	
	fitvec = [CalcFitness(popmat[i],w) for i in range(N)]
	
	fitsum = np.sum(fitvec)
	CSfitvec = np.cumsum(fitvec)
	
	print fitsum
	print CSfitvec
	
	for i in range(N):
		rv = np.random.uniform(0.0,fitsum)
		pick = sum(CSfitvec < rv)
		print rv,pick
		newpop[i,:] = popmat[pick,:]
	
	return newpop



minD = 5
maxD = 35


N = 20
endgen = 30
mu = 0.0001
pop = NewPop(N, minD, maxD)

maxvec = np.empty([endgen])
for gen in range(endgen):
	pop1 = ApplyRouletteFitness(pop, fitness1)
	pop2 = CrossOver(pop1)
	pop3 = Mutate(pop2, mu)
	maxvec[i] = np.max(pop3)
	




### main
def main():
	print "Bisection : %.4f\n" % bisect(f_x, 0, 0.11, 20, 1)

	print "Newton Raphson: %.4f\n" % newt_raph(f_x, f_prime, 0.1, 20, 1)
	
	print "Secant Method: %.4f\n" % secant(f_x, 0, 0.11, 20, 1)
	
	


### run if main
if __name__ == "__main__":
	main()
