import bisect
import numpy



class feature:
	#start and end are inclusive!
	def __init__(self,lsp):
		#cp is the coding position of the first base in exon: 0,1,2. Former is the name of last gene
		#lsp=line.split()
		warnDir=False
		#print lsp
		#this handles cases where chromosomes have a "chrX" and those that have the form "X"
		self.chr=lsp[0].split('r')[-1]
		if len(lsp)>5:
			self.dir=lsp[5]
		
			if self.dir!='+' and self.dir!='-' and not warnDir:
				print "Warning: direction is not + or -: Verify bed file format!" 
				print "will ignore columns beyond 5"
				warnedDir=True
				self.dir=None
		else:
			self.dir=None
		self.begin=int(lsp[1])+1
		self.end=int(lsp[2])
		if len(lsp)>3:
			self.name=lsp[3]
		else:
			self.name=None
		if len(lsp)>4:
			self.rest=lsp[4:]
		else:
			self.rest=[]
		self.len=abs(self.begin-self.end)+1
		#self.phase=int(lsp[7])
		"""if(former==self.gen):
			self.codeposb=0
			self.codepose=(self.len-1)%3
		else:
			self.codeposb=cp
			self.codepose=(cp+self.len-1)%3"""
		#test if end >, < begin
		
		if(self.begin>self.end):
			print "Error: self.begin>self.end"
			print lsp
			print line
			raise ValueError
			
	# now define function	
	def dist(self, pos):
		#returns the distance to the feature
		if (self.begin<=pos and self.end>=pos):
			#print "begin ", self.begin
			return 0
			
		elif self.begin>pos:
			return self.begin-pos
		elif self.end<pos:
			return pos-self.end
	
	#def annot(self,pos):
	#	verdict=self.isin(pos)
	#	if ver
	#	return verdict
	



	
	
class bed(object):	
	
	def __init__(self,fileN):			
		self.chrdict={}				
		self.chrstarts={}			
		self.chrends={}				
		a=open(fileN,'r')
		lines=a.readlines()
		former="";
		genes=[]
		self.features=[]
		for lin in lines:
			lsp=lin.split()
			
			if len(lsp)<3 or lsp[0][0]=='#':
				continue
			"""from gtf files
			if len(lsp)<3 or lsp[0][0]=='#':
				continue
			if lsp[2]=="gene" or lsp[2]=="transcript":
				gen=gene(lsp)
				#print gen.chr
				try:
					self.chrdict[gen.chr].append(gen)
					self.chrstarts[gen.chr].append(gen.begin)
					self.chrends[gen.chr].append(gen.end)
				except KeyError:
					self.chrdict[gen.chr]=[gen]
					self.chrstarts[gen.chr]=[gen.begin]
					self.chrends[gen.chr]=[gen.end]
		
		
			elif lsp[2]=="CDS":		
				gen.CDS.append(exon(lsp))
				if gen.CDS[-1].end>gen.end or gen.CDS[-1].begin<gen.begin:
					print "Error: CD exceeds gene!"
			else:
				#print lsp"""
			feat=feature(lsp)
			try:
				self.chrdict[feat.chr].append(feat)
				self.chrstarts[feat.chr].append(feat.begin)
				self.chrends[feat.chr].append(feat.end)
			except KeyError:
				self.chrdict[feat.chr]=[feat]
				self.chrstarts[feat.chr]=[feat.begin]
				self.chrends[feat.chr]=[feat.end]

		
		
		for chr, lst in self.chrstarts.iteritems():
			
			
			if (sorted(lst)==lst):
				print "starts sorted: success!"
			else:
				print"starts are not sorted!\n\n"
				
				b=numpy.array(lst).argsort()
				self.chrstarts[chr]=list(numpy.array(lst)[b] )
				self.chrends[chr]=list(numpy.array(self.chrends[chr])[b] )
				self.chrdict[chr]=list(numpy.array(self.chrdict[chr])[b] )
				print "sorted according to starts"
				lst=self.chrstarts[chr]
				if (sorted(lst)==lst):
					print "starts sorted: success!"
				else:
					print "starts are not sorted!\n\n"
					sys.exit(2)
				
			lste=self.chrends[chr]
		
			if (sorted(lste)==lste):
				print "ends sorted!"
			#else:
			#	#print lste
				
		
		#Calculate orderdict, which contains the position of the intervals that contain the start of the interval at a given position. If no overlap, orderdict[chr][i]=[[i]] 
		self.orderdict={}
		for chr, starts in self.chrstarts.iteritems():
			
			self.orderdict[chr]=[[0]]
			
			for ii in range(1,len(starts)):
				#print "start is" , ii
				self.orderdict[chr].append([])
				
				#print "len is ", len(self.orderdict[chr][ii-1])
				
				for posprev in self.orderdict[chr][ii-1]:
					#print "compare ",self.chrends[chr][posprev],self.chrstarts[chr][ii-1]
					if self.chrends[chr][posprev]>=self.chrstarts[chr][ii]:
						self.orderdict[chr][ii].append(posprev)
				self.orderdict[chr][ii].append(ii)			
			self.orderdict[chr].append([])	
		print "parsed bed file"
		#print self.chrstarts
		#print self.chrends
		#print self.orderdict
	
#given intervals (perhaps overlapping) as two lists (begin, end), 
#find the position of intervals in which a fits. To do this:
#1-find maxpos, the position of the last start lower than a. We also want orderdict, which contains the position of the intervals that contain the start of maxpos
	
	def loc(self,a,beg,end):
		#print "position ", a#, beg
		
		#print "length ", len(beg)
		#print sorted(beg)==beg
		#print "bisection result ", bisect.bisect_right(beg,a)
		#foo=1110358
		#print "foo ",foo
		#print "bisection result ", bisect.bisect_right(beg,foo)
		#print beg
		
		#we cannot do a right bissection. Only consider the last 10 genes, assuming there is no overlap of more than 10 genes
		maxpos=bisect.bisect_right(beg,a)-1
		#print "maxpos", maxpos
		return maxpos
	def ann(self,a,chr):#annotate
		ls=[]
		#print self.chrstarts.keys(),self.chrends.keys()
		try:
			maxpos=self.loc(a,self.chrstarts[chr],self.chrends[chr])
		except KeyError:
			print "chrom not in bedfile"
			return None 
		if maxpos==-1: #then value before any field
			return None
		#print maxpos, len(self.orderdict[chr])
		for p in self.orderdict[chr][maxpos]:
			#print "annotation ", self.chrdict[chr][p].annot(a)
			try:
				#print self.chrdict[chr][p].annot(a)
				if self.chrdict[chr][p].begin>a:
					print "Error; should not be looking before start of segment"
					raise IndexError
				
				if self.chrdict[chr][p].end>=a:
					if self.chrdict[chr][p].name:
						return self.chrdict[chr][p].name
					else:
						return "noname"
					ls.append((annotation[0],annotation[1]))
				#print "begin: ",self.chrdict[chr][p].begin
			except IndexError:
				pass
		#print "returning ",ls
		return None
	def annfull(self,a,chr):
		ls=[]
		#print self.chrstarts.keys(),self.chrends.keys()
		try:
			maxpos=self.loc(a,self.chrstarts[chr],self.chrends[chr])
		except KeyError:
			print "chrom not in bedfile"
			return None 
		if maxpos==-1: #then value before any field
			return None
		#print maxpos, len(self.orderdict[chr])
		for p in self.orderdict[chr][maxpos]:
			#print "annotation ", self.chrdict[chr][p].annot(a)
			try:
				#print self.chrdict[chr][p].annot(a)
				if self.chrdict[chr][p].begin>a:
					print "Error; should not be looking before start of segment"
					raise IndexError
				
				if self.chrdict[chr][p].end>=a:
					return (self.chrdict[chr][p].name,self.chrdict[chr][p].rest)
				
					ls.append((annotation[0],annotation[1]))
				#print "begin: ",self.chrdict[chr][p].begin
			except IndexError:
				pass
		#print "returning ",ls
		return None
		
