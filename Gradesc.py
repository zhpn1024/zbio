#################################################   
# Regression: input   
# Author : zhangp   
# Date   : 2014-05-26   
# HomePage : http://   
# Email  :   
#################################################   
 
from numpy import *  
#import matplotlib.pyplot as plt  
import time  
from multiprocessing import Process, Queue

# calculate the sigmoid function   
#def sigmoid(inX):  
#    return 1.0 / (1 + exp(-inX))  

def lgstBinProb(x, theta, args):
	if 'p0' not in args: args['p0'] = 0.0005456521
	if 'p1' not in args: args['p1'] = 0.01
	#p0 = 0.0005456521, p1 = 0.01
	#print shape(x), shape(theta[1:])
	return args['p0'] + args['p1'] / (1 + exp(- x * theta))
	
def lgstBinPD(x, y, theta, args):
	if 'p0' not in args: args['p0'] = 0.0005456521
	if 'p1' not in args: args['p1'] = 0.01
	p0= args['p0']
	p1 = args['p1']
	numSamples, numFeatures = shape(x)
	#d=mat(theta.copy())
	#print 'p1=', p1
	p = lgstBinProb(x, theta, args)
	a = (p - p0) / p / (1 - p) / p1
	b = y[:,1] - multiply(p,y[:,0])
	c = multiply(a, b)
	s = sum(c) / numSamples
	#print 's=', s
	args['s'] = s
	d = x.transpose() * multiply(c, p1 + p0 - p) / numSamples
	return d

#	cox binomial distribution model
#	x is row vector or data matrix
#	theta is column vetor that represent weights of x
def coxBinProb(x, theta, args):
	if 'p0' not in args: args['p0'] = 0.0005456521
	return args['p0'] + exp(-x * theta)
	
#	Partial derivative of cox binomial distribution log likelyhood
def coxBinPD(x,y,theta,args):
	if 'p0' not in args: args['p0'] = 0.0005456521
	p0=args['p0']
	numSamples, numFeatures = shape(x)
	p = coxBinProb(x, theta, args)
	a = (p - p0) / p / (1 - p)
	b = multiply(p,y[:,0])-y[:,1]
	c = multiply(a, b)
	d = x.transpose() * c / numSamples
	return d
	
#global d
#d=1
def lognCk(n, k):
	#reduce(add, [log(n-i)-log(k-i) for i in range(k)])
	return sum([float(log(n-i)-log(k-i)) for i in range(k)])
def log0(p):
	numSamples, numFeatures = shape(p)
	out=p.copy()
	for i in range(numSamples):
		if p[i]>0: out[i] = log(p[i])
	return out
def coxBinLogL(x, y, theta, args):
	numSamples, numFeatures = shape(x)
	p = coxBinProb(x, theta, args)
	pFull = y[:,1] / y[:,0]
	pNull = sum(y[:,1]) / sum(y[:,0])
	print 'p', shape(p), p[0:3].transpose(),'pFull', shape(pFull), pFull[0:3].transpose(),'pNull', pNull
	nck=mat(map(lognCk, y[:, 0], y[:, 1])).transpose()
	print 'nck', shape(nck), nck[0:3]
	p1=multiply(y[:, 1], log(p))
	p1f = multiply(y[:, 1], log0(pFull))
	p1n = multiply(y[:, 1], log(pNull))
	print 'p1', shape(p1), p1[0:3]
	p2=multiply((y[:, 0]-y[:, 1]), log(1 - p))
	p2f=multiply((y[:, 0]-y[:, 1]), log(1 - pFull))
	p2n=multiply((y[:, 0]-y[:, 1]), log(1 - pNull))
	print 'p2', shape(p2), p2[0:3]
	logL=sum(nck + p1 + p2)
	logLf=sum(nck + p1f + p2f)
	logLn=sum(nck + p1n + p2n)
	#print 'logL', logL[0:3]
	return logL

def binLogL(x, y, theta, args, prob):
	numSamples, numFeatures = shape(x)
	p = prob(x, theta, args)
	pFull = y[:,1] / y[:,0]
	pNull = sum(y[:,1]) / sum(y[:,0])
	print 'p', shape(p), p[0:3].transpose(),'pFull', shape(pFull), pFull[0:3].transpose(),'pNull', pNull
	nck=mat(map(lognCk, y[:, 0], y[:, 1])).transpose()
	print 'nck', shape(nck), nck[0:3].transpose()
	p1=multiply(y[:, 1], log(p))
	o=log0(pFull)
	print 'o', o[0:3].transpose()
	p1f = multiply(y[:, 1], o)
	p1n = multiply(y[:, 1], log(pNull))
	print 'p1', shape(p1), p1[0:3].transpose()
	p2=multiply((y[:, 0]-y[:, 1]), log(1 - p))
	p2f=multiply((y[:, 0]-y[:, 1]), log0(1 - pFull))
	p2n=multiply((y[:, 0]-y[:, 1]), log(1 - pNull))
	print 'p2', shape(p2), p2[0:3].transpose()
	logL=sum(nck + p1 + p2)
	logLf=sum(nck + p1f + p2f)
	logLn=sum(nck + p1n + p2n)
	#print 'logL', logL[0:3]
	args['logL'] = logL
	args['logLf'] = logLf
	args['logLn'] = logLn
	return logL

lgstBin={'func':lgstBinProb,'drv':lgstBinPD}
coxBin={'func':coxBinProb,'drv':coxBinPD}

#def run(q, alpha, drv, x, y, theta, j):
#	dd=alpha * drv(x,y,theta,j)
#	q.put((dd,j))
	#print j,dd
# train a model using gradient descent algorithm   
# input: train_x is a mat datatype, each row stands for one sample   
#		 train_y is mat datatype too, each row is the corresponding result   
#		 theta is original input and trained step by step
#        opts is optimize option include step and maximum number of iterations   
def gradesc(train_x, train_y, theta, args={}, model=lgstBin, alpha= 0.05, maxIter= 20, highest=0.5,maxThr=1, brk=1e-7, like=True):  
	# calculate training time   
	startTime = time.time()  
	
	numSamples, numFeatures = shape(train_x)  
	print 'numSamples', numSamples, 'numFeatures', numFeatures
	k=0
	while k != maxIter:
		d = model['drv'](train_x, train_y, theta, args)
		m = sqrt((d.transpose() * d)[0,0])
		#print d.transpose()
		if m < brk: break
		if m * alpha > highest : alpha = highest / m 
		#print max(d),min(d),alpha,m/alpha,k
		#print alpha, k, m#, str(time.time() - startTime)
		if k%100==0: 
			print 'alpha=',alpha, 'k=',k, 'm=',m
		#print d.transpose()
		#print shape(d.transpose() * d)
		#print theta.transpose(),k
		ad = alpha * d
		rr=1.0
		if k>0 :
			r=(d.transpose()*lastd)[0,0]/lastm/lastm
			r2=m/lastm
			#print r,r2
			if 0.98<= r < 1 and r2 <= 1: 
				if lastrr < 1: rr = 1.5
				else: rr = 1.2
				#alpha*=1.1
				#lastd=d*1.1
			elif 0.95 <= r < 0.98 and r2 <= 1: 
				rr=1.1
			elif 1 <= r < 1.01 and r2 < 1.01:
				rr = 1.01
			elif 0.85 <= r< 0.95 and r2 > 0.95:
				rr/=1.2
				ad /= 2
			elif 0.8 <= r < 0.85:
				if lastrr >= 1: 
					rr /= 1.5
				else: rr /= 1.2
				#theta-=d/4
				ad /= 3
			elif  0< r < 0.8 : 
				if lastrr >= 1 : 
					rr /= 2
					ad /= 4
				else: 
					rr /= 1.5
					ad /= 3
				#theta-=d/3
				#ad /= 3
			#elif 0.6<= r < 0.8 : 
			#	rr /=2
			elif -1<= r < 0 : 
				rr /= 10
				#theta-=d/2
				ad /= 5
			elif r < -1 : 
				rr /= 20
				#theta-=d*3/4
				ad /= 10
			if k%100==0: print 'r',r, r2, rr#, '\n'
		theta += ad
		#print theta.transpose(), '\n'
		#time.sleep(2)
		alpha *= rr
		#if alpha>1: 
		#	rr*=1/alpha
		#	alpha=1
			
		lastd = d
		lastm = m
		lastrr = rr
		k+=1

	print 'Congratulations, training complete! Took %fs!' % (time.time() - startTime)  
	print theta.transpose()
	if like :
		print 'Log likelihood: %f' % binLogL(train_x, train_y, theta, args, model['func'])
	print 'args: ', args
	return theta  



def lgstBinProb1(x, theta, args):
	if 'p0' not in args: args['p0'] = 0.0005456521
	#if 'p1' not in args: args['p1'] = 0.01
	#p0 = 0.0005456521, p1 = 0.01
	#print shape(x), shape(theta[1:])
	return args['p0'] + theta[0] * 0.01 / (1 + exp(- x * theta[1:]))
	
def lgstBinPD1(x, y, theta, args):
	if 'p0' not in args: args['p0'] = 0.0005456521
	#if 'p1' not in args: args['p1'] = 0.01
	p0= args['p0']
	#p1 = args['p1']
	numSamples, numFeatures = shape(x)
	d=mat(theta.copy())
	#print 'p1=', p1
	p = lgstBinProb1(x, theta, args)
	a = (p - p0) / p / (1 - p) / theta[0]
	b = y[:,1] - multiply(p,y[:,0])
	c = multiply(a, b)
	d[0] = sum(c) * 0.01 / numSamples
	#print 's=', s
	d[1:] = x.transpose() * multiply(c, theta[0] + p0 - p) / numSamples
	return d

def lgstBinLogL1(x, y, theta, args):
	numSamples, numFeatures = shape(x)
	p = lgstBinProb1(x, theta, args)
	print 'p', shape(p), p[0:3]
	nck=mat(map(lognCk, y[:, 0], y[:, 1])).transpose()
	print 'nck', shape(nck), nck[0:3]
	p1=multiply(y[:, 1], log(p))
	print 'p1', shape(p1), p1[0:3]
	p2=multiply((y[:, 0]-y[:, 1]), log(1 - p))
	print 'p2', shape(p2), p2[0:3]
	logL=sum(nck + p1 + p2)
	#print 'logL', logL[0:3]
	args['logL'] = logL
	return logL
	
lgstBin1={'func':lgstBinProb1,'drv':lgstBinPD1}


# test your trained Logistic Regression model given test set   
#def testLogRegres(weights, test_x, test_y):  
#    numSamples, numFeatures = shape(test_x)  
#    matchCount = 0  
#    for i in xrange(numSamples):  
#        predict = sigmoid(test_x[i, :] * weights)[0, 0] > 0.5  
#        if predict == bool(test_y[i, 0]):  
#            matchCount += 1  
#    accuracy = float(matchCount) / numSamples  
#    return accuracy  


# show your trained logistic regression model only available with 2-D data   
'''def showLogRegres(weights, train_x, train_y, function):  
    # notice: train_x and train_y is mat datatype   
	numSamples, numFeatures = shape(train_x)  
	if numFeatures != 3:  
	print "Sorry! I can not draw because the dimension of your data is not 2!"  
	return 1  
	
   	# draw all samples   
	for i in xrange(numSamples):  
		if int(train_y[i, 0]) == 0:  
			plt.plot(train_x[i, 1], train_x[i, 2], 'or')  
		elif int(train_y[i, 0]) == 1:  
			plt.plot(train_x[i, 1], train_x[i, 2], 'ob')  
			
	# draw the classify line   
	min_x = min(train_x[:, 1])[0, 0]  
	max_x = max(train_x[:, 1])[0, 0]  
	weights = weights.getA()  # convert mat to array   
	y_min_x = float(-weights[0] - weights[1] * min_x) / weights[2]  
	y_max_x = float(-weights[0] - weights[1] * max_x) / weights[2]  
	plt.plot([min_x, max_x], [y_min_x, y_max_x], '-g')  
	plt.xlabel('X1'); plt.ylabel('X2')  
	plt.show()
'''
''''
d=mat(ones((numFeatures,1)))
q = Queue()

#d2={}
#alpha = opts['alpha']; maxIter = opts['maxIter']  
#weights = ones((numFeatures, 1))  
	t=[0]*maxThr
	""""
	ti=0
	for j in range(numFeatures):
		if j<maxThr:
			ti=j
		else:
			ti=j % maxThr
			t[ti].join()
			#if not t[ti].is_alive(): print t[ti].d
			#d[t[ti].j,0]=q.get()
			#print t[ti].exitcode
			#print t[ti].isrun(),t[ti].j,t[ti].d,d.transpose()
		#t[ti]=pro(alpha, coxBinPD, train_x, train_y, theta, j)
		t[ti]=Process(target=run,args=(q, alpha, coxBinPD, train_x, train_y, theta, j))
		t[ti].start()
		#print str(j)+" start.."+str(time.time() - startTime)
	"""
	j=0
	for ti in range(maxThr):
		t[ti]=Process(target=run,args=(q, alpha, coxBinPD, train_x, train_y, theta, j))
		t[ti].start()
		j+=1
	
	while j<numFeatures:
		for ti in range(maxThr):
			if not t[ti].is_alive():
				t[ti]=Process(target=run,args=(q, alpha, coxBinPD, train_x, train_y, theta, j))
				t[ti].start()
				j+=1
				if j>=numFeatures: break
		if j>=numFeatures: break
		time.sleep(0.1)
	for ti in range(maxThr):
		t[ti].join()
		#d[t[ti].j,0]=t[ti].d
		#print t[ti].j,t[ti].d,d.transpose()
		#d[j,0]=d2[j]
		#d[j,0]=alpha * coxBinPD(train_x,train_y,theta,j)
		#print theta1[j],d
		#theta1[j] += d
	#q.join()
	for j in range(numFeatures):
		(dd,jj)=q.get()
		d[jj,0]=dd
s=0.0
for i in range(len(x)):
	if x[i,j]==0 :continue
	p=coxBinProb(x[i],theta)
	#print i,j,y[i]
	s+=x[i,j]*(p-p0)*(y[i,0]*p-y[i,1])/p/(1-p)
return s[0,0]/len(x)
'''
"""		
33.        if opts['optimizeType'] == 'gradDescent': # gradient descent algorilthm   
34.            output = sigmoid(train_x * weights)  
35.            error = train_y - output  
36.            weights = weights + alpha * train_x.transpose() * error  
37.        elif opts['optimizeType'] == 'stocGradDescent': # stochastic gradient descent   
38.            for i in range(numSamples):  
39.                output = sigmoid(train_x[i, :] * weights)  
40.                error = train_y[i, 0] - output  
41.                weights = weights + alpha * train_x[i, :].transpose() * error  
42.        elif opts['optimizeType'] == 'smoothStocGradDescent': # smooth stochastic gradient descent   
43.            # randomly select samples to optimize for reducing cycle fluctuations    
44.            dataIndex = range(numSamples)  
45.            for i in range(numSamples):  
46.                alpha = 4.0 / (1.0 + k + i) + 0.01  
47.                randIndex = int(random.uniform(0, len(dataIndex)))  
48.                output = sigmoid(train_x[randIndex, :] * weights)  
49.                error = train_y[randIndex, 0] - output  
50.                weights = weights + alpha * train_x[randIndex, :].transpose() * error  
51.                del(dataIndex[randIndex]) # during one interation, delete the optimized sample   
52.        else:  
53.            raise NameError('Not support optimize method type!')  
54.      '''''"""
