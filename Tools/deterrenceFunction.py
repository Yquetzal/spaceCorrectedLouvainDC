import math
import matplotlib.pyplot as plt


#this file contains all useful functions for computing and ploting the deterrence function
def plotScatterFree(seriesToDisplay,maxX=-1,log=False,axesNames=("x","f(x)")):
	for serie in seriesToDisplay:
		plt.plot(serie[0],serie[1],label=serie[2])
	plt.legend(loc=2)
	#plt.xscale('log')
	if log:
		plt.yscale('log')

	if maxX>0:
		plt.xlim(plt.xlim()[0], maxX)
	ax = plt.axes()
	ax.xaxis.grid()
	ax.yaxis.grid()
	plt.ylabel(axesNames[1])
	plt.xlabel(axesNames[0])
	plt.show()


def fromDictionaryOutputOrderedKeysAndValuesByKey(aDict):
	sortedKeys = aDict.keys()
	sortedKeys = sorted(sortedKeys)
	sortedValues = []
	for key in sortedKeys:
		sortedValues.append(aDict[key])
	return(sortedKeys,sortedValues)


def convertWithPrecision(val,precision):
	#precision of rounding as a distance to the "." for instance :
	# 123.456 with precision 2 = 123.46
	# 123.456 with precision -2 = 100
	if precision >= 0:
		theDist = round(val, precision)
	else:
		theDist = int(val / math.pow(10, abs(precision))) * math.pow(10, abs(
			precision))  # if negative, round to closest dimension
	return theDist



#this class is used to compute a deterrence function. First call deterrenceFunctionEstimation and then getDeterrenceAtDistance
class deterrenceFunction():
	def __init__(self):
		self.distancesDic = {}
		self.roundDecimals=-1
		self.maximalDistance=-1

	def deterrenceFunctionEstimation(self,INs,OUTs,observedGraph, distances, roundDecimals, maximalDistance=100000, plot=False ):
		self.roundDecimals=roundDecimals
		self.maximalDistance=maximalDistance
		sumIn = sum(observedGraph.in_degree(weight="weight").values())


		byDistEstimated ={}
		byDistObserved ={}

		for e in observedGraph.edges(data=True):
			source = e[0]
			dest = e[1]


			theDist = distances(source,dest)
			#if theDist!=-1:
			theDist = convertWithPrecision(theDist,roundDecimals)

			if theDist<maximalDistance:
				byDistEstimated.setdefault(theDist,[])
				byDistObserved.setdefault(theDist,[])

				byDistEstimated[theDist].append(OUTs[source]*INs[dest]/sumIn)
				byDistObserved[theDist].append(e[2]["weight"])


		dicCoeffdistance = {}
		minVals = 3
		for d in byDistObserved:
			if len(byDistObserved[d])>minVals and sum(byDistEstimated[d])>0: #consider only if we have at least minVals values
				dicCoeffdistance[d]=sum(byDistObserved[d])/sum(byDistEstimated[d])


		if plot:
			(x, y) = fromDictionaryOutputOrderedKeysAndValuesByKey(dicCoeffdistance)
			plotScatterFree([(x[:maximalDistance],y[:maximalDistance],"deterrence function")])
		self.distancesDic = dicCoeffdistance

	def getDeterrenceAtDistance(self,dist):
		distNorm = convertWithPrecision(dist,self.roundDecimals)
		if distNorm>=self.maximalDistance or not distNorm in self.distancesDic:
			return min(self.distancesDic.values())/10
		return self.distancesDic[distNorm]
	#return distFunc
