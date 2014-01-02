''' 
James Ma

run program:
------------
python distributions.py
'''
import math,random,sys,csv,time





#-------------------------------------------
# Useful functions
#-------------------------------------------

### Functions used in integration
# (a) function
def func_a(x):
	return math.exp(-x**2 / 2.0) / (math.sqrt(2.0 * math.pi)) 

# (b) function
def func_b(x):
	return math.sin(2.0 * x)**2 / float(x**2)

# (c) function
def func_c(x):
	return x**(-1.0 / 2)





### Funcions borrowed from Assignment 1 - solutions
# Generate Standard Cauchy
def rcauchy(k,s):
	'''
	Draw random samples from Standard Cauchy distribution.
	
	### inputs:
	### k - number of random variables to generate
	### s - seed value
	'''
	### set seed
	random.seed(s)
	### creates random variable using transformation
	x = [math.tan(math.pi*(random.random() - .5)) for i in range(k)]
	### reset seed
	random.seed()
	return x



### Generate Normal
def rnorm(k,mu,sigma,s):
	'''
	Draw random samples from a Normal distribution.
	
	### inputs:
	### k - number of random variables to generate
	### mu - mean parameter
	### sigma - standard deviation parameter
	### s - seed value
	'''
	### set seed
	random.seed(s)
	
	### set size so that n/2 is whole number
	ni = (k + k%2)/2
	
	### random values
	x = [random.random() for i in range(2*ni)]
	
	### creates random variable using box meuller transformation
	x = ([mu + sigma*(math.sqrt(-2*math.log(x[i]))*math.cos(2*math.pi*x[i+ni/2])) for i in range(ni)] +
	[mu + sigma*(math.sqrt(-2*math.log(x[i]))*math.sin(2*math.pi*x[i+ni/2])) for i in range(ni)])
	
	### reset seed
	random.seed()
	return x[0:k]



## Generate chisquared with mean df
def rchisq(k,df,s):
	'''
	Draw random samples from chi-square distribution.
	
	### inputs:
	### k - number of random variables to generate
	### df - degrees of freedom
	### s - seed value
	'''
	random.seed(s)
	
	### creates k*df N(0,1)
	x = rnorm(k*df,0,1,s)
	### squares each value
	x = [x[i]**2 for i in range(df*k)]
	### creates random variable by adding each df size set
	x = [sum(x[df*i:df*(i+1)])for i in range(k)]
	
	# reset seed
	random.seed()
	return x



# Generate T-distribution
def rtdist(n,df,s):
	'''
	Draw random samples from T distrubtion
	
    Parameters:
    -----------
    n: number of samples
    df: degrees of freedom
    s: seed number
    '''
	# seed
	random.seed(s)
    
    # if n is one, return value 
	if (n == 1):
		x = rnorm(n, 0, 1, s)[0] * math.sqrt(float(df) / rchisq(n,df,s)[0])
		
		# reset seed
		random.seed()
		return x
	# else return list
	else:
		norms = rnorm(n, 0, 1, s)
		rchis = rchisq(n,df,s)
		x = [norms[i] * math.sqrt(float(df) / rchis[i]) for i in xrange(n)]
		
		# reset seed
		random.seed()
		return x



# Generate Uniform distribution
def runif(a, b, size,s):
    '''
    Draw samples from a Uniform distribution using python's 
    built-in function.
    
    Parameters:
    -----------
    a: min
    b: max
    size: number of samples
    s: seed number
    '''
    # set seed
    random.seed(s)
    
    # list of uniforms
    x = [random.uniform(a,b) for i in xrange(size)]
    
    # reset seed
    random.seed()
    return x


# Pdf for T-distribution 
def ptdist(t,v):
	'''
	Returns the pdf of T-distribution with v degrees of freedom.
	
    Parameters:
    -----------
    t: t-value
    v: degrees of freedom
    '''
	return gamma((v + 1) / 2.0) * (1 + t**2 / float(v)) ** (-(v + 1) / 2.0) / (math.sqrt(v*math.pi) * gamma(v/2.0))

# Gamma helper function
def gamma(x):
	'''
	Helper function for the pdf of T-distribution.
	
    Parameters:
    -----------
    x: value to apply gamma function to
    '''
	if ((x - math.floor(x)) == 0):
		return math.factorial(x - 1)
	elif ((x - math.floor(x)) == 0.5):
		x -= 0.5
		return math.sqrt(math.pi)*math.factorial(2*x) / ((4**x) * math.factorial(x))
	else:
		raise TypeError


# Pdf for Standard Cauchy distribution
def pcauchy(x):
	'''
    Parameters:
    -----------
    x: cauchy
    '''
	return 1.0 / (math.pi * (1 + x**2))


# Pdf for Uniform distribution 
def punif(a,b):
	'''
	Returns the pdf of the Uniform distribution in interval [a,b].
	
	
    Parameters:
    -----------
    a: min
    b: max
    '''
	return float(1) / (b - a)


# Trapezoidal rule for Rhomberg Integration
def trapint1(a, b, infunc,nmax):
	'''
	Returns the integral using the Trapezoidal rule.  
	
	Same as below without error checking and printing.  
	Used for the first step in Rhomberg integration.
	'''
	# I1, first iteration
	Iold = (b - a)/2.0 * (infunc(a) + infunc(b))
	Inew = Iold
	
	# Main loop while less than nmax
	for k in xrange(2,nmax + 1):
		f = 0
		# Calculate the sum of f's
		for i in xrange(1,2**(k - 2) + 1):
			f += infunc(a + ((2*i - 1) * (b - a)) / (2.0**(k - 1)))
			
		# calc Inew
		Inew = 1.0 / 2 * Iold + (b - a) / 2.0**(k - 1) * f 
		
		# assign previous value
		Iold = Inew
		
	return Inew



### Integral values used to calculate error
# www.wolframalpha.com
integrals = dict(func_a =  0.341345, func_b =  0.0517455, func_c =  2.82843)


# Integral values for extended midpoint
integrals2 = dict(func_a =  0.5, func_b = .15727, func_c =  2.82843)



### Export data
def outdat(x,fname):
	'''
	### creates a textfile called test.txt
	###
	and creates a writer object

	### inputs:
	### x - data
	### fname - the desired filename
	'''
	with open(fname,'w') as ofile:
		writer = csv.writer(ofile)
		for i in xrange(len(x)): writer.writerow([x[i]])



### Weights and nodes for GuassQuad using 32 points
# http://pomax.github.io/bezierinfo/legendre-gauss.html

weights = [0.0965400885147,0.0965400885147,0.0956387200793,0.0956387200793,
			0.0938443990808,0.0938443990808,0.0911738786958,0.0911738786958,
			0.0876520930044,0.0876520930044,0.0833119242269,0.0833119242269,
			0.0781938957871,0.0781938957871,0.0723457941088,0.0723457941088,
			0.0658222227764,0.0658222227764,0.0586840934785,0.0586840934785,
			0.0509980592624,0.0509980592624,0.0428358980222,0.0428358980222,
			0.034273862913,0.034273862913,0.0253920653093,0.0253920653093,
			0.0162743947309,0.0162743947309,0.00701861000947,0.00701861000947]

nodes = [-0.0483076656877,0.0483076656877,-0.144471961583,0.144471961583,
		-0.239287362252,0.239287362252,-0.331868602282,0.331868602282,
		-0.421351276131,0.421351276131,-0.506899908932,0.506899908932,
		-0.587715757241,0.587715757241,-0.66304426693,0.66304426693,
		-0.73218211874,0.73218211874,-0.794483795968,0.794483795968,
		-0.849367613733,0.849367613733,-0.896321155766,0.896321155766,
		-0.934906075938,0.934906075938,-0.964762255588,0.964762255588,
		-0.985611511545,0.985611511545,-0.997263861849,0.997263861849]







#------------------------------------------------
# Parts (1) - (8)
#------------------------------------------------

# 1. Extended Trapezoidal rule
def trapint(a, b, infunc, nmax, e):
	'''
	The Trapezoidal rule for calculating the integral over the range [a,b].
	
	I = h * (f(x0) + f(x1)) / 2  + ... + h * (f(x n-1) + f(xn)) / 2) + error
	
	The range [a,b] is divided into equally spaced mesh points with step size
	
	h = (b - a) / n
	
	And n trapezoids are calculated and added to estimate the integral.
	
	Parameters:
	-----------
	a: min value of integration
	b: max value of integration
	infunc: function to integrate
	nmax: max number of iterations
	e: error threshold
	'''
	# Assign integral value to calculate error
	if(infunc == func_a):
		integral_val = integrals['func_a']
	elif(infunc == func_b):
		integral_val = integrals['func_b']
	else:
		integral_val = integrals['func_c']
	
	
	# I1, first iteration
	Iold = (b - a) / 2.0 * (infunc(a) + infunc(b))
	print "iteration=%d, l value=%.4f, error=%.4e" % (1,Iold,1000)
	
	# Main loop 
	for k in xrange(2,nmax + 1):
	
		# Initialize f
		f = 0
		
		# Calculate the sum of f's, which are heights of the mini trapezoids
		for i in xrange(1,2**(k - 2) + 1):
			f += infunc(a + ((2*i - 1) * (b - a)) / (2.0**(k - 1)))
			
		# Calculate Inew, as a function of Iold.  
		# also multiply by h to get area of the trapezoids
		Inew = 1.0/2 * Iold + (b - a) / 2.0**(k - 1) * f 
		
		# Check for convergence
		if (abs(Iold - Inew) <= abs(Iold * e)):
			print ""
			return Inew
		
		# Calculate error
		error = abs(Inew - integral_val)
		print "iteration=%d, l value=%.4f, error=%.4e" % (k,Inew,error)
		
		# Store previous value
		Iold = Inew
	
	print ""
	return Inew






# 2. Romberg integration
def romint(a,b,infunc,m,e):
	'''
	The Rhomberg integration for calculating the integral over the range [a,b].
	
	R(n,1): Trapezoidal rule with n mesh points
	R(n,k): (4**(k - 1) * R(n,k-1) - R(n-1,k-1)) / (4**(k - 1) - 1)
	
	The first step of each row is the Trapezoidal rule.  Each step after that
	is an extension of Simpson's rule.
	
	
	Parameters:
	-----------
	a: min value of integration
	b: max value of integration
	infunc: function to integrate
	m: max number of iterations
	e: error threshold
	'''
	# Initialize lists to store values of rows (old row, current row)
	old_row = []
	curr_row = []
	Inew = 0
    
	# Main loop
	for i in xrange(1,m + 1):

		# Calculate R(n,1) using Trapezoidal rule
		first_row_val = trapint1(a,b,infunc,i)
		# Append R(n,1)
		curr_row.append(first_row_val)
	
		# Calculate values in row R(n,k)
		for j in xrange(1,i):
			Inew = (4**(j) * curr_row[j-1] - old_row[j-1]) / (4**(j) -1)
			# Append R(n,k)
			curr_row.append(Inew)
	
		# Check for convergence
		if(len(old_row)> 0):
			if(abs(curr_row[j-1] - curr_row[j]) <= e * curr_row[j-1]):
				print ""
				return curr_row[j]
		
		# print row
		print "iteration=%d,%s" % (i, ",".join(map(str,curr_row)))
		
		# Set old row = current row
		old_row = curr_row
		
		# Reset current row
		curr_row = []
		
	print ""
	return Inew




# 3. Gaussian Quadrature
def GaussQuad(a,b,infunc,weights,nodes):
    '''
	The 32-point Guassian Quadrature for calculating the 
	integral over the range [a,b].
	
	This method is different from the previous methods in that 
	it assumes that f(x) is some high-ordered polynomial and 
	does not use equally spaced points, but rather chosen in some 
	optimal way.
	
    
	Parameters:
	-----------
	a: min value of integration
	b: max value of integration
	infunc: function to integrate
	weights: weights used for 32-point GaussQuad
	nodes: nodes used for 32-point GuassQuad
    '''
    # Initialize values
    val = 0
    
    # Convert values a,b so that the function is in the interval [-1,1]
    m = (b - a) / 2.0
    c = (b + a) / 2.0
    
    for i in xrange(len(weights)):
    	val += weights[i]*infunc(m*nodes[i] + c)
    return val * m





# 4. Naive Monte Carlo
def mc1int(a,b,n,infunc,fmax,s):
    '''
	The Naive Monte Carlo for calculating the integral over the range [a,b].
	
	I = A * P(yi < f(xi)) + error
	
	The method samples an area and estimates the integral by 
	finding the proportion of sampled points below the integral
	and multiplying it by the sampled area.
	
    
	Parameters:
	-----------
	a: min value of integration
	b: max value of integration
	infunc: function to integrate
	fmax: upper value of samples, a number greater than max(f(x))
	s: seed number
    '''
    # Initialize sums
    sums = 0
    
    # set seed
    random.seed(s)
    
    # Calculate area
    A = fmax*(b - a)
    
    # Loop through n points
    for i in xrange(n):
    	# Generate xi,yi from uniform
    	xi = random.uniform(a,b)
    	yi = random.uniform(0,fmax)
    	
    	# Check condition to find P(yi < f(xi))
    	if(yi < infunc(xi)):
    		sums += 1
    # reset seed
    random.seed()
    
    return A * sums / float(n)




# 5. Better Monte Carlo
def bmc1int(a,b,n,infunc,fmax,s):
	'''
	The Better Monte Carlo for calculating the integral over the range [a,b].
	
	Ihat = (b - a) / n * ( f(x1) + f(x2) + ... + f(xn) )
	
	This method is similar to that above.  It is a bit better because it does
	not need to sample yi nor check conditions.
	

	Parameters:
	-----------
	a: min value of integration
	b: max value of integration
	infunc: function to integrate
	fmax: upper value of samples, a number greater than max(f(x))
	s: seed number
	'''
	# Initialize sums
	sums = 0
	
	# set seed	
	random.seed(s)
	
	# Generate n points from uniform
	for i in xrange(n):
		xi = random.uniform(a,b)
		sums += infunc(xi)
	
	# reset seed
	random.seed()
	
	return (b - a) / float(n) * sums





# 6. Better Monte Carlo - antithetic variables
def bmc2int(a,b,n,infunc,fmax,s):
	'''
	The  Monte Carlo using anthetic variables
	for calculating the integral over the range [a,b].
	
	This method is similar to that above.  It is a bit better because 
	it uses a variance reduction technique by sampling from negatively 
	correlated variables.
	

	Parameters:
	-----------
	a: min value of integration
	b: max value of integration
	infunc: function to integrate
	fmax: upper value of samples, a number greater than max(f(x))
	s: seed number
	'''
	# Initailize sums
	sums = 0
	
	# set seed
	random.seed(s)
	
	# Generate n/2 points from uniform and n/2 from negatively correlated uniform
	for i in xrange(n / 2):
		xi = random.uniform(a,b)
		x2 = b - xi
		sums += (infunc(xi) + infunc(x2))
	
	# reset seed
	random.seed()
	
	return float(b - a) / (n / 2) * 1 / 2.0 * sums




# 7. Extended midpoint rule
def xmidpt(a,b,infunc,nmax,e):
	'''
	The  Extended Midpoint rule
	for calculating the integral over the range [a,b].
	
	This method is similar to the Trapezoidal rule except it calculates
	the midpoints between mesh points.  This helps it avoid endpoints 
	that are problematic.  
	
	Similar to the Trapezoidal rule, it can reuse old values.
	
	
	Parameters:
	-----------
	a: min value of integration
	b: max value of integration
	infunc: function to integrate
	nmax: max number of iterations
	e: error threshold
	'''
	# Assign integral value to calculate error
	if(infunc == func_a):
		integral_val = integrals2['func_a']
	elif(infunc == func_b):
		integral_val = integrals2['func_b']
	else:
		integral_val = integrals2['func_c']
	
	
	# Check for infinity and assign a float value of 100
	if (b == float('inf')):
		b = 100.0
	
	
	# I1, first iteration
	Iold = (b - a) * (infunc((a + b) / 2.0))
	
	# print first iteration
	print "iteration=%d, l value=%.4f, error=%.4e" % (1,Iold,1000)
	
	
	# Main loop up to nmax
	for k in xrange(2,nmax + 1):
		f = 0
		j = 0
		
		# Calculate the sum of f's, which are heigts
		for i in xrange(0, 2*3**(k - 2)):
			# This ifelse allows the use of old Ivalues and skips
			# redundant calculations
			if((i % 2) == 0):
				j += 2
				f += infunc(a + ((j - 1) * (b - a)) / (2.0*3**(k - 1)))
			else:
				# skips previous values, increment j by 4
				j += 4
				f += infunc(a + ((j - 1) * (b - a)) / (2.0*3**(k - 1)))
		
		# Calculate Inew by using Iold
		Inew = 1.0 / 3 * Iold +  f * (b - a)  / (3.0**(k - 1))
	
		# Check for convergence
		if (abs(Iold - Inew) <= Iold * e):
			print ""
			return Inew
		
		# Calculate error
		error = abs(Inew - integral_val)
		
		# print iteration 
		print "iteration=%d, l value=%.4f, error=%.4e" % (k,Inew,error)
		
		# Assign previous value
		Iold = Inew
	
	print ""
	return Inew



# 8. Importance Sampling
# a. Using Uniform(-50,50)
def importance1(a,b,infunc,r_func,p_func,n,s):
	'''
	The Importance Sampling method using a Uniform distribution.
	
	Ihat = 1 / n * ( f(x1) / g(x1) + f(x2) / g(x2) + ... + f(xn) / g(xn))
	
	This method is identical to Better Monte Carlo except that 
	it will weight points f(xi) by the inverse pdf.

	
	Parameters:
	-----------
	a: min value for uniform
	b: max value for uniform
	infunc: function to integrate
	r_func: distribution to generate values
	p_func: calculate pdf 
	n: number of values to generate 
	s: seed number
	'''
	# Generate values from uniform
	x = r_func(a,b,n,s)
	
	# find f(x) for each x value and weight by the inverse pdf
	sums = [infunc(i) / p_func(a,b) for i in x]
	
	return float(sum(sums)) / n


# b. Using Standard Cauchy(-50,50)
def importance2(infunc,r_func,p_func,n,s):
	'''
	The Importance Sampling method using a Standard Cauchy distribution.
	
	Ihat = 1 / n * ( f(x1) / g(x1) + f(x2) / g(x2) + ... + f(xn) / g(xn))
	
	This method is identical to Better Monte Carlo except that 
	it will weight points f(xi) by the inverse pdf.
	
	
	Parameters:
	-----------
	infunc: function to integrate
	r_func: distribution to generate values
	p_func: calculate pdf 
	n: number of values to generate 
	s: seed number
	'''
	# generate from cauchy
	x = r_func(n,s)
	
	# find f(x) for each x value and weight by the inverse pdf, only if abs(x) < 50
	sums = [infunc(i) / p_func(i) for i in x if (abs(i) < 50)]
	return float(sum(sums)) / n


# c. Using T-distribution(30)
def importance3(df,infunc,r_func,p_func,n,s):
	'''
	The Importance Sampling method using a T-distribution.
	
	Ihat = 1 / n * ( f(x1) / g(x1) + f(x2) / g(x2) + ... + f(xn) / g(xn))
	
	This method is identical to Better Monte Carlo except that 
	it will weight points f(xi) by the inverse pdf.
	
	
	Parameters:
	-----------
	df: degrees of freedom for T-distribution
	infunc: function to integrate
	r_func: distribution to generate values
	p_func: calculate pdf 
	n: number of values to generate 
	s: seed number
	'''
	# generate from t-distribution(df)
	x = r_func(n,df,s)
	
	# find f(x) for each x value and weight by the inverse pdf, only if abs(x) < 50
	sums = [infunc(i) / p_func(i,df) for i in x if (abs(i) < 50)]
	return float(sum(sums)) / n




### main function
# for printing results
def main():

	# 1. Trapezoidal
	print "Trapezoidal Rule applied to integral a: "
	trapint(0,1,func_a,10,10**(-6))
	print "Trapezoidal Rule applied to integral b: "
	# calculate time for question 2
	start1 = time.time()
	trapint(math.pi,1.5*math.pi,func_b,10,10**(-6))
	end1 = time.time()
	
	
	# 2. Rhomberg
	print "Rhomberg Integration applied to integral a: "
	romint(0,1,func_a,10,10**(-6))
	print "Rhomberg Integration applied to integral b: "
	# calculate time for question 2
	start2 = time.time()
	romint(math.pi,1.5*math.pi,func_b,10,10**(-6))
	end2 = time.time()
	
	print "Running time for Trapezoid: %f\nRunning time for Rhomberg: %f" % (end1 - start1,end2 - start2)
	print "It seems that Rhomberg is faster and requires less computation\n"
	
	# 3. Gaussian Quadrature
	print "Gaussian Quadrature applied to integral a (32-points): %f\n" % GaussQuad(0,1,func_a,weights,nodes)
	
	
	# 4. Naive Monte Carlo
	print "Naive Monte Carlo applied to integral a: %f\n" % mc1int(0,1,1000,func_a,1,1)
	print "Naive Monte Carlo applied to integral b: %f\n" % mc1int(math.pi,1.5*math.pi,1000,func_b,1,1)
	print "Naive Monte Carlo applied to integral c: %f\n" % mc1int(0,2,1000,func_c,1,1)
	
	
	# 5. Better Monte Carlo
	print "Better Monte Carlo applied to integral a: %f\n" % bmc1int(0,1,1000,func_a,1,1)
	print "Better Monte Carlo applied to integral b: %f\n" % bmc1int(math.pi,1.5*math.pi,1000,func_b,1,1)
	print "Better Monte Carlo applied to integral c: %f\n" % bmc1int(0,2,1000,func_c,1,1)
	
	
	# 6. Better Monte Carlo - antithetic variables
	print "Better Monte Carlo (antithetic variables) applied to integral a: %f\n" % bmc2int(0,1,1000,func_a,1,1)
	print "Better Monte Carlo (antithetic variables) applied to integral b: %f\n" % bmc2int(math.pi,1.5*math.pi,1000,func_b,1,1)
	print "Better Monte Carlo (antithetic variables) applied to integral c: %f\n" % bmc2int(0,2,1000,func_c,1,1)

	# 7. Extended midpoint rule
	print "Extended midpoint rule applied to integral a: "
	xmidpt(0,float('inf'),func_a,15,10**(-6))
	print "Extended midpoint rule applied to integral b: "
	xmidpt(math.pi,float('inf'),func_b,10,10**(-6))
	print "Extended midpoint rule applied to integral c: "
	xmidpt(0,2,func_c,10,10**(-6))
	
	# 8. Importance sampling
	print "Importance sampling applied to integral a, using Uniform(-50,50): %f\n" % importance1(-50,50,func_a,runif,punif,1000,1)
	print "Importance sampling applied to integral a, using Standard Cauchy[-50,50]: %f\n" % importance2(func_a,rcauchy,pcauchy,1000,1)
	print "Importance sampling applied to integral a, using T-Distribution(30 df): %f\n" % importance3(30,func_a,rtdist,ptdist,1000,1)
	
	# Answering Question
	print "\nWhich importance function does best?"
	print "As expected, the T-Distribution does the best."
	print "It is very close to the distribution that we are approximating (N(0,1)). "
	print "The sampling distribution for the T-distribution is within +/- 0.01 of the true area. \n"
	print "The next best is the Standard Cauchy, which is not as close to the N(0,1)"
	print "but better than the Uniform.\n"
	print "The sampling distribution for Cauchy is within +/- 0.06 of the true area.\n"
	print "The Uniform gave the worst approximations with sampling distribution +/- 0.6\n\n"
	
	# 8. Histogram part
	# Uniform[-50,50]
	print "Getting values for histogram using Uniform(-50,50)"
	x = [importance1(-50,50,func_a,runif,punif,1000,i) for i in xrange(5000)]
	outdat(x,'unif.txt')
	
	# Standard Cauchy[-50,50]
	print "Getting values for histogram using Standard Cauchy"
	x = [importance2(func_a,rcauchy,pcauchy,1000,i) for i in xrange(5000)]
	outdat(x,'stdcauchy.txt')
	
	# T-distribution
	print "Getting values for histogram using T-Distribution"
	x = [importance3(30,func_a,rtdist,ptdist,1000,i) for i in xrange(5000)]
	outdat(x,'tdist.txt')
	


# run main
if __name__ == "__main__":
	main()



