import math
import matplotlib.pyplot as plt


#this file contains all useful functions for computing and ploting the deterrence function
def _plotScatterFree(seriesToDisplay, maxX=-1, log=False, axesNames=("x", "f(x)")):
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


def _fromDictionaryOutputOrderedKeysAndValuesByKey(aDict):
	sortedKeys = aDict.keys()
	sortedKeys = sorted(sortedKeys)
	sortedValues = []
	for key in sortedKeys:
		sortedValues.append(aDict[key])
	return(sortedKeys,sortedValues)


def _convertWithPrecision(val, precision):
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

	def deterrenceFunctionEstimation(self,INs,OUTs,observedGraph, distances, roundDecimals, minVals = 3, maximalDistance=100000, plot=False ):
		"""
		Compute the deterrence function
		:param INs: dictionary of in-degrees
		:param OUTs: dictionary of out-degrees
		:param observedGraph: a nx.Graph , the observed network
		:param distances: a function that return the distance betwee two provided nodes
			:param roundDecimal: the rounding used to compute bins of the deterrence function. for a distance d=123.456 :
		 if roundDecimal=2, binned value is 123.46.
		 if roundDecimal=-2, binned value is 100
		:param maximalDistance: ignore in most cases, parameter of the deterence function to set an upper bound on the considered distances
		:param minValsBin: parameter of the deterrence function, minimum number of observations in a bin to consider it. (avoid abherent values for rare distances)
		:param plot: if True, plot the deterrence function before the doubly constrained process and at the end of the process
		"""
		self.roundDecimals=roundDecimals
		self.maximalDistance=maximalDistance
		sumIn = sum(dict(observedGraph.in_degree(weight="weight")).values())


		byDistEstimated ={}
		byDistObserved ={}

		for e in observedGraph.edges(data=True):
			source = e[0]
			dest = e[1]


			theDist = distances(source,dest)
			#if theDist!=-1:
			theDist = _convertWithPrecision(theDist, roundDecimals)

			if theDist<maximalDistance:
				byDistEstimated.setdefault(theDist,[])
				byDistObserved.setdefault(theDist,[])

				byDistEstimated[theDist].append(OUTs[source]*INs[dest]/sumIn)
				byDistObserved[theDist].append(e[2]["weight"])


		dicCoeffdistance = {}
		for d in byDistObserved:
			if len(byDistObserved[d])>minVals and sum(byDistEstimated[d])>0: #consider only if we have at least minVals values
				dicCoeffdistance[d]=sum(byDistObserved[d])/sum(byDistEstimated[d])


		if plot:
			(x, y) = _fromDictionaryOutputOrderedKeysAndValuesByKey(dicCoeffdistance)
			_plotScatterFree([(x[:maximalDistance], y[:maximalDistance], "deterrence function")])
		self.distancesDic = dicCoeffdistance

	def getDeterrenceAtDistance(self,dist):
		"""
		for a provided distance, return the associated deterrence
		:param dist: a distance
		"""
		distNorm = _convertWithPrecision(dist, self.roundDecimals)
		if distNorm>=self.maximalDistance or not distNorm in self.distancesDic:
			#return min(self.distancesDic.values())/10
			return 0
		return self.distancesDic[distNorm]
	#return distFunc
