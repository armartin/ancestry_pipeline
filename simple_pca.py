'''
Chris Gignoux

simple Marsland-esque PCA

ideal for post-processing pcamask data

input: a column of labels, columns from 2 on are features, each row is a single observation

'''

from numpy import *
from sys import argv

debugging = 0
if debugging:
	infile = 'tempcaitlin2_2.pca.txt'
else:
	infile = argv[-1]
	

print '%s\nRunning PCA post-PCAmask\n%s\n\n\n' % ('#'*50,'#'*50)

print 'reading in data from %s' % (infile)
try:
	data = [line.strip().split() for line in file(infile)]
except IOError:
	print 'sorry, input file %s not found' % (infile)
	# quit
	assert 1 == 0

labels = [line[0] for line in data]
features = array([line[1:] for line in data], dtype=float)


# filter out individuals with nan/no relevant ancestry
print 'finding non-missing values'
numberfilter = features[:,0] * 1 == features[:,0]

usedata = features[numberfilter,:]
print 'running PCA'
# covariance is teensy, can use eigenvalue decomposition
evals, evecs = linalg.eig(cov(usedata.T))
# sort output by PC
print 'post-processing PCA'
indices = argsort(evals)[::-1]
evals = evals[indices]
evecs = evecs[:,indices]

# normalize data
for i in range(evecs.shape[1]):
	evecs[:,i] /= linalg.norm(evecs[:,i]) * sqrt(evals[i])

x = dot(evecs.T,usedata.T).T
# revise output
outputvals = where(features,-999.,0)
outputvals[numberfilter,:] = x

outfile = file(infile.replace('.pca.txt','.pca.orthog.txt') , 'w')
print 'writing output to %s' % (outfile)
outfile.write('ID\t%s\n' % ('\t'.join(['PC%s' % (i+1) for i in range(outputvals.shape[1])])))

for i, line in enumerate(outputvals):
	outputdata = line
	if -999 in outputdata:
		outputdata = '\t'.join(['NA','NA'])
	else:
		outputdata = '\t'.join([str(val) for val in line])
	outfile.write('%s\t%s\n' % (labels[i], outputdata))

outfile.close()
print 'done'
