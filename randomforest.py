import numpy as np
import operator

# node
class node:
	def __init__(self,left=None,right=None,key=None,val=None,dist=None):
		self.left = left
		self.right = right
		self.key = key
		self.val = val
		self.dist = dist

# decision tree, classification?
class decision_tree(object):
	'''
	CART implementation of decision tree.
	input   : X,Y
	returns : tree
	'''
	
	def __init__(self):
		self.node = node
	
	def fit(self,x,y):
		assert x.shape[0] == y.shape[0], "Dimension mismatch!"
		
		
		if len(y) <=1:
			return
		
		n = x.shape[0]
		m = x.shape[1]
		
		gini_dict = {}
		for feature in xrange(m):
			# partition y's
			# assuming binary classification
			unique_feature_vals = np.unique(x[:,feature])
			for val in xrange(len(unique_feature_vals)-1):
				# calc gini for feature
				left = self.gini(unique_feature_vals[unique_feature_vals<=val],y[unique_feature_vals<=val])
				right = self.gini(unique_feature_vals[unique_feature_vals>val],y[unique_feature_vals>val])
				nleft = len(y[unique_feature_vals<=val])
				nright = len(y[unique_feature_vals>val])
				nm = len(y)
				gini = nleft / float(nm) * left + nright / float(nm) * right
				# {feature:[val,gini]
				gini_dict[feature] = [val,gini]
		
		# find max to recurse on
		min_gini = 100000
		min_feature = 0
		min_val = 0
		for feature,val in gini_dict:
			if val[1] < min_gini:
				min_gini = val
				min_feature = key
				min_val = val[0]
		
		# recurse
		self.left = self.fit(x[x[:,min_feature] <=min_val,:],y[x[:,min_feature] <=min_val])
		self.right = self.fit(x[x[:,min_feature] >min_val,:],y[x[:,min_feature] >min_val])
		
		new_node = node(key=min_feature,val=min_val)
		 
		return new_node
	
	def gini(self,feature_vals,y):
		n = y.shape[0]
		f_i = {}
		for val in feature_vals:
			f_i[val] = sum(y[y==val])
		gini = 1 - sum([item / float(n) for item in f_i.values()])
		return gini
	
	def predict(self,x):
		n = x.shape[0]
		m = x.shape[1]
		
		for row in xrange(n):
			if x[row,self.node.key] < self.node.val & self.node.left != None:
				current_node = self.node.left
			else:
				pass # TODO
			


#--------------------------------------------------------------------------
# Random Forest

class random_forest(object):
	'''
	CART implementation of decision tree.
	input   : X,Y
	returns : tree
	'''
	
	def __init__(self,n_estimators=None):
		self.node = node
		self.n_estimators = n_estimators
		self.tree = [] # list of trees
		
	def fit(self,x,y):
		assert x.shape[0] == y.shape[0], "Dimension mismatch!"
		
		
		if len(y) <=1:
			return
		
		n = x.shape[0]
		m = x.shape[1]
		
		gini_dict = {}
		
		for iter in xrange(n_estimators):
			for feature in xrange(m):
				# partition y's
				# assuming binary classification
				unique_feature_vals = np.unique(x[:,feature])
				for val in xrange(len(unique_feature_vals)-1):
					# calc gini for feature
					left = self.gini(unique_feature_vals[unique_feature_vals<=val],y[unique_feature_vals<=val])
					right = self.gini(unique_feature_vals[unique_feature_vals>val],y[unique_feature_vals>val])
					nleft = len(y[unique_feature_vals<=val])
					nright = len(y[unique_feature_vals>val])
					nm = len(y)
					gini = nleft / float(nm) * left + nright / float(nm) * right
					# {feature:[val,gini]
					gini_dict[feature] = [val,gini]
			
			# find max to recurse on
			min_gini = 100000
			min_feature = 0
			min_val = 0
			for feature,val in gini_dict:
				if val[1] < min_gini:
					min_gini = val
					min_feature = key
					min_val = val[0]
			
			# recurse
			self.left = self.fit(x[x[:,min_feature] <=min_val,:],y[x[:,min_feature] <=min_val])
			self.right = self.fit(x[x[:,min_feature] >min_val,:],y[x[:,min_feature] >min_val])
			
			new_node = node(key=min_feature,val=min_val)
		 
		return new_node
	
	def gini(self,feature_vals,y):
		n = y.shape[0]
		f_i = {}
		for val in feature_vals:
			f_i[val] = sum(y[y==val])
		gini = 1 - sum([item / float(n) for item in f_i.values()])
		return gini
	
	def predict(self,x):
		n = x.shape[0]
		m = x.shape[1]
		
		for row in xrange(n):
			if x[row,self.node.key] < self.node.val & self.node.left != None:
				current_node = self.node.left
			else:
				pass # TODO
			

#--------------------------------------------------------------------------
# Naive Bayes
