import csv
import math

sent_dict = {}
relation_dict = {}
countryMap = {}
targetRelations = ['AG.LND.TOTL.K2','BN.KLT.DINV.CD','BX.GSR.MRCH.CD','EG.ELC.PROD.KH',
					'EN.ATM.CO2E.KT','EP.PMP.DESL.CD','FP.CPI.TOTL.ZG','IT.NET.USER.P2',
					'NY.GDP.MKTP.CD','SP.DYN.LE00.IN','SP.POP.TOTL']
relation_keywords_dict = {}
countryClassDict = {}
countryAliases = {}

def generateCountryAliases():
	global countryAliases
	with open('countryAliases.tsv') as tsvin:
		tsvin = csv.reader(tsvin, delimiter='\t')
		for row in tsvin:
			countryAliases[row[0]] = row[1]

#generate an indentitiy matrix of order nxn
def makeIdentity(n):
    result=[]
    for i in xrange(n):
    	result.append([])
    	for j in xrange(n):
    		if (i ==j):
        		result[i].append(1)
        	else:
        		result[i].append(0)
    return result

#transpose of a matrix
def transpose(a):
	tr = []
	col = len(a[0])
	for j in range(col):
		l = []
		for i in range(len(a)):
			l.append(a[i][j])
		tr.append(l)
	return tr

def multiplyMatrices(a, b):  # a is self.tr and b is self.index
    aRows = len(a)
    aCols = len(a[0])
    bRows = len(b)
    bCols = len(b[0])
    assert(aCols == bRows)
    rows = aRows
    cols = bCols
    #c = make2dList(rows, cols)
    c=[]
    for row in xrange(rows):
    	l = []
        for j in xrange(cols):
        	#c[row].append([])
        	dotProduct = 0.0
        	for k in xrange(aCols):
        		dotProduct += a[row][k]*b[k][j]
        	l.append(dotProduct)
        c.append(l)
    return c

def multiplyRowOfSquareMatrix(m, row, k):
    n = len(m)
    rowOperator = makeIdentity(n)
    rowOperator[row][row] = k
    return multiplyMatrices(rowOperator, m)

def addMultipleOfRowOfSquareMatrix(m, sourceRow, k, targetRow):
    n = len(m)
    rowOperator = makeIdentity(n)
    rowOperator[targetRow][sourceRow] = k
    return multiplyMatrices(rowOperator, m)

def invertMatrix(m):
    n = len(m)
    assert(len(m) == len(m[0]))
    inverse = makeIdentity(n)
    for col in xrange(n):
        diagonalRow = col
        assert(m[diagonalRow][col] != 0)
        k = 1.0/m[diagonalRow][col]
        m = multiplyRowOfSquareMatrix(m, diagonalRow, k)
        inverse = multiplyRowOfSquareMatrix(inverse, diagonalRow, k)
        sourceRow = diagonalRow
        for targetRow in xrange(n):
            if (sourceRow != targetRow):
                k = -m[targetRow][col]
                m = addMultipleOfRowOfSquareMatrix(m, sourceRow, k, targetRow)
                inverse = addMultipleOfRowOfSquareMatrix(inverse, sourceRow, k, targetRow)
    return inverse

# def make2dList(rows, cols):
#     a=[]
#     for row in xrange(rows): a += [[0]*cols]
#     return a

#coeff stores coefficients in decreasing order of powers
def computePolynomial(coeff,deg,x):
	val = 0
	for i in range(deg+1):
		val = (val + coeff[i])*x
	return val

def derivative(coeff,deg):
	ans = [0] * deg
	for i in range(deg):
		ans[i] = coeff[i]*(deg - i)
	return (ans,deg-1)

def multiplyPolys(coeff1,coeff2,deg1,deg2):
	deg = deg1 + deg2
	ans = [0] * (deg+1)
	for i in range(deg1+1):
		for j in range(deg2+1):
			ans[i+j] += coeff1[i] * coeff2[j]
	return (ans,deg)

def subPolys(coeff1,coeff2,deg1,deg2):
	#deg1>deg2 <- Assume
	ans = coeff1
	for i in range(deg2+1):
		ans[i+deg1-deg2] = coeff1[i+deg1-deg2] - coeff2[i]
	return (ans,deg1)

def multiplyPolyByNum(coeff,deg,x):
	for i in range(deg+1):
		coeff[i] *= x
	return (coeff,deg)

def newtonRaphson(coeff,deg,x0,err):
	x = x0
	while True:
		val = computePolynomial(coeff,deg,x)
		if abs(val)<=err:
			return x
		deriv = computePolynomial(derivative(coeff,deg)[0],deg-1,x)
		x = x - val/deriv

def solveDifferential(coeff,deg,x,y,x0):
	deriv = derivative(coeff,deg)
	product = multiplyPolys(deriv[0],coeff,deg-1,deg)
	temp = multiplyPolyByNum(deriv[0],deg-1,y)
	temp = subPolys(product[0],temp[0],product[1],temp[1])
	final = subPolys(temp[0],[-1,x],temp[1],1)
	return newtonRaphson(final[0],final[1],x0,0.0001)

def getDistance(coeff,deg,x,y):
	point = solveDifferential(coeff,deg,x,y,1954)
	dist = math.sqrt((x-point)**2 + (y-computePolynomial(coeff,deg,point))**2)
	return dist	

def read_sentences():
	#changes made to selected_indicators, please add "no_unit" to the file at the end of population and internet users relation separated by tab
	with open('sentences.tsv','r') as tsvin:
		tsvin = csv.reader(tsvin, delimiter='\t')
		for row in tsvin:
			key = row[0].strip()
			sent_dict[key] = []
			sent_dict[key].append(row[1])
			sent_dict[key].append(row[2])
			sent_dict[key].append(row[3])
		
	with open('selected_indicators','r') as tsvin:
		tsvin = csv.reader(tsvin, delimiter='\t')
		for row in tsvin:
			key = row[2].strip()
			relation_dict[key]=[]
			relation_dict[key].append(row[1])
			relation_dict[key].append(row[3])

def countryCodes():
	with open('countries_id_map.txt','r') as f:
		for line in f:
			x = line.strip().split('	')
			countryMap[x[0]] = x[1]

def read_knowledge_base():
	country_dict = {}
	with open('kb-facts-train_SI.tsv','r') as tsvin:
		tsvin = csv.reader(tsvin, delimiter='\t')
		count = 0
		for row in tsvin:
			count += 1
			key = row[0].strip() +"#"+ row[2].strip()
			try:
				country_dict[key].append(row[1])
			except :
				country_dict[key]=[]
				country_dict[key].append(row[1])

def keywords_relation():
	with open('keywords_relation.csv','r') as csvin:
		csvin = csv.reader(csvin)
		for row in csvin:
			key = row[1].strip()
			relation_keywords_dict[key] = []
			s = row[2]
			while (len(s) > 0):
				i = s.find(",")
				if(i == -1):
					s1 = s
					s = ""
				else:	
					s1 = s[0:i]
				s = s[i+1:len(s)]
				relation_keywords_dict[key].append(s1)

def find_which_relation():
	for x in sent_dict:
		l = sent_dict[x]
		sent = l[0]
		countlist = []
		frequency = 0
		relation = []
		for y in relation_keywords_dict:
			list1 = relation_keywords_dict[y]
			count = 0

			for i in list1:
				i = i .strip()
				score = int(i[len(i)-1])
				i = i[0:len(i)-1]
				sent1 = sent
				while (sent1.find(i) != -1):
					
					count = count + score
					sent1 = sent1[sent1.find(i)+1:len(sent1)]
				
			if(count > frequency):
				relation = []
				l = []
				l.append(y)
				l.append(count)
				frequency = count
				relation.append(l)
			elif(count == frequency):	
				l = []
				l.append(y)
				l.append(count)
				frequency = count	
				relation.append(l)

		if(len(relation) == 0 or relation[0][1] == 0):
			l = ["none"]
			relation = [l]

		s = ""
		for i in relation:
			s = s+","+i[0]
		sent_dict[x].append(s[1:len(s)])	

	#pls read this, relation is list of list,and its list element has the code of relation it is talking about at index 0
	# sent dictionary: sentenec, numbers,country,relation

class Country():
	def __init__(self,id):
		self.id = id
		self.name = countryMap[id]
		self.generateCoefficients()

	def generateDistribution(self):
		self.normalDist = {}
		for elm in targetRelations:
			self.normalDist[elm] = []
		with open('kb-facts-train_SI.tsv') as f:
			for line in f:
				x = line.strip().split('	')
				if x[0] == self.id:
					self.normalDist[x[2]].append(float(x[1]))
		self.meanList = {}
		self.stdList = {}
		for elm in self.normalDist:
			val = 0.0
			count = 0
			for x in self.normalDist[elm]:
				val += x
				count += 1
			try:
				self.meanList[elm] = val/count
			except ZeroDivisionError:
				self.meanList[elm] = 0
		for elm in self.normalDist: 
			val = 0.0
			count = 0
			for x in self.normalDist[elm]:
				val += (x - self.meanList[elm])**2
				count += 1
			try:
				self.stdList[elm] = math.sqrt(val/count)
			except ZeroDivisionError:
				self.stdList[elm] = 0
	
	def generateCoefficients(self):
		self.coeffDict = {}
		self.dataDict = {}
		for elm in targetRelations:
			self.coeffDict[elm] = []
			self.dataDict[elm] = []
		with open('kb-facts-train_SI.tsv') as f:
			for line in f:
				x = line.strip().split('	')
				if x[0] == self.id:
					self.dataDict[x[2]].append(float(x[1]))
		for elm in self.dataDict:
			if(len(self.dataDict[elm]) != 0):	
				self.dataDict[elm].sort()
				self.index=[]
				#self.index = [None]*(len(self.dataDict[elm]))
				for i in range(len(self.dataDict[elm])):
					l = []
					for j in range(2):
						l.append(i**j)
					self.index.append(l)
				self.tr = transpose(self.index)
				self.a = invertMatrix(multiplyMatrices(self.tr,self.index))     # Doing inverse(X' * X)
				self.dataDict[elm] = [self.dataDict[elm]]							
				self.b = multiplyMatrices(self.tr,transpose(self.dataDict[elm]))  # doing X' * Y
				self.coeffDict[elm] = multiplyMatrices(self.a,self.b)		# doing W = inverse(X' * X) * X'* Y

	def reverseColumnMatrix(self,a):
		b = list(a)
		l = len(a)
		for i in range(l):
			b[l-1-i] = a[i]
		return b

def populateClassObjects():
	for elm in countryMap:
		try:
			countryClassDict[countryMap[elm]] = Country(elm)
		except Exception as e:
			print e

def getIndex(l,x):
	for i in range(len(l)):
		if x>=l[i]: return i
	return len(l)

def assignConfidenceScores():
	for elm in sent_dict:
		try:
			numbers = sent_dict[elm][1].strip().split(',')
			numbers = [float(x.strip(',')) for x in numbers]
			countries = sent_dict[elm][2].split(',')
			countries = [x.strip(' ') for x in countries]
			relations = sent_dict[elm][3].split(',')
			relations = [x.strip(' ') for x in relations]
		except Exception as e:
			print "No number found. Cannot do anything."
		try:
			for t in countries:
				p = countryAliases[t]
				obj = countryClassDict[p].coeffDict
				for s in relations:
					coeffMatrix = obj[s]
					for r in numbers:
						index = getIndex(countryClassDict[p].dataDict[s],r)
						newMatrix = []
						for val in coeffMatrix:
							newMatrix.append(val[0])
						dist = getDistance(newMatrix,len(newMatrix)-1,index,r)
						score = math.exp(-(dist*0.000000001))
						#SentId [TAB] Country [TAB] Relation [TAB] Number [TAB] Score
						print(str(elm)+'	'+str(t)+'	'+str(s)+'	'+str(r)+'	'+str(score))
		except:
			pass

def doEverything():
	read_sentences()      # reads the sentences and stores them in sent_dict
	read_knowledge_base()	# creates the knowledge base for training later on
	keywords_relation() # stores all the key words associated with relations
	find_which_relation()	# sentences get their relations
	countryCodes()		 #reads the country codes
	generateCountryAliases()
	populateClassObjects()
	assignConfidenceScores()

doEverything()