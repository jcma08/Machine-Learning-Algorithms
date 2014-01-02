import numpy as np
import math


# f_prime
def f_prime(x):
	return 4 * x**3 - 9 * x**2



# gradient descent
def grad_desc(f_prime,eps=0.01,precision=0.0001):

	x_old = 0
	x_new = 6 # The algorithm starts at x=6
	eps = 0.01 # step size
	precision = 0.00001
		 
	while abs(x_new - x_old) > precision:
		x_old = x_new
		x_new = x_old - eps * f_prime(x_old)
	print "Local minimum occurs at ", x_new

# node
class tree(object):
	def __init__(self,left=None,right=None,key=None,val=None,item=None):
		# key is condition
		# val is value (float, int, str)
		# left if condition true, else right
		# item = [outcomes..] only at base
		self.left = left
		self.right = right
		self.key = key
		self.val = val
		self.item = item
	
	def __repr__(self):
		return "%s(%r)" % (self.__class__,self.__dict__)
	
	@staticmethod
	def get_val(x):
		node = self
		while (node.left):
			if (x < self.val):
				node = self.left
			else:
				node = self.right
		return node.item

# derivative of func
def deriv(f,x,step=0.0001):
	derivs = []
	for i in xrange(len(f)):
		if (isinstance(f,tree));
			derivs.append((f.get_val(x + step) - f.get_val(x)) / step)
		else:
			derivs.append(0)
	return derivs


# squared loss
def mse(x,y):
	assert (len(x)==len(y)), "lengths do not match"
	loss = 0
	for item in xrange(len(x)):
		loss += (x[item] - y[item])**2
	return loss / float(len(x))


# training stump
def train(x,y):
	# minimize square loss
	n = len(x)
	loss = []
	for i in xrange(n):
		thresh = x[i]
		loss1 = mse([np.mean(y[x<thresh])]*n,y[x<thresh])
		loss2 = mse([np.mean(y[x>=thresh])]*n,y[x>=thresh])
		loss.append((loss1 + loss2)/2)
	max_val = max(enumerate(loss),key=operator.itemgetter(0))
	
	left_item = np.mean(y[max_val<thresh])
	right_item =  np.mean(y[max_val>=thresh])
	h = tree(val = x[max_val])
	h.left.item = left_item
	h.right.item = right_item
	return h

# function f(x)
def func(f,x,gamma):
	val = 0
	for i in xrange(len(f)):
		if (i==0):
			val += f[i]
		else:
			val -= gamma[i-1]*gradient(f,x)*loss(y[i],f[i-1](x[i]))

# gbr
def gbr(x,y,m=10):
	gamma = np.mean(y)
	
	f = [gamma]
	stepsize=0.001
	n=len(x)
	
	residuals = np.zeros((n,1))
	
	for i in xrange(m):
		for j in xrange(n):
			residuals[i] = -deriv(f,x[j])
		
		h = train(x,residuals)






#----------------------------------------------------------------------------
import numpy as np

# sgd using linear regression
class linear_sgd(object):
	
	def __init__(self,alpha=0.0001):
		base.__init__(self)
		# parameter :  alpha
		self.alpha = alpha
		self.weights = 0.
	
	def __repr__(self):
		return "%s %r" % ("Linear Regression",self.weights)
		
	def fit(self,x,y):
		assert(x.shape[0] == y.shape[0]), "Dimension mismatch!"
		n = x.shape[0]
		m = x.shape[1]
		
		# may need to randomize x in [0..1]
		self.weights = np.matrix(np.ones(m+1))
		
		# add temporary ones to x matrix
		x = np.hstack((np.matrix(np.ones(n)).T,x))
		
		# update weights
		for idx in xrange(n):
			self.weights -= self.alpha * 2 * (self.weights.dot(x[idx,:].T) - y[idx]) *  x[idx,:]
		# done
		
	def predict(self,x):
		x = np.hstack((np.matrix(np.ones(n)).T,x))
		return x.dot(self.weights.T)
	
	def coef_(self):
		return self.weights


#----------------------------------------------------------------------------
# sgd using logistic regression
# only supports binary classification
class logistic_sgd(object):
	
	def __init__(self,alpha=0.0001):
		# parameter :  alpha
		self.alpha = alpha
		self.weights = np.array([])
		
	def fit(self,x,y):
		assert(x.shape[0] == y.shape[0]), "Dimension mismatch!"
		n = x.shape[0]
		m = x.shape[1]
		
		# may need to randomize x in [0..1]
		self.weights = np.matrix(np.ones(m+1))
		
		# add temporary ones to x matrix
		x = np.hstack((np.matrix(np.ones(n)).T,x))
		
		# update weights
		for idx in xrange(n):
			self.weights -= self.alpha *2*(self.__sigmoid__(self.weights,x[idx,:].T) - y[idx])*x[idx,:]
		# done
	
	def __sigmoid__(self,weights,x):
		hx = 1. / (1 + np.exp(-weights.dot(x)))
		return hx

	def predict(self,x):
		x = np.hstack((np.matrix(np.ones(n)).T,x))
		return self.__sigmoid__(x,self.weights.T)
	
	def coef_(self):
		return self.weights




#---------------------------------------------------------------------------
# decision tree
class node:
	def __init__(self,left=None,right=None,key=None,val=None):
		self.left = left
		self.right = right
		self.key = key
		self.val = val
	

class decision_tree:
	def __init__(self):
		pass
	
	def fit(self,x,y):
		
	def predict(self,x,y):
	







