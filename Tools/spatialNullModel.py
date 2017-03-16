from .deterrenceFunction import *
from .utilSpatialNullM import *
from .GraphModels import *
import math
from scipy.stats.stats import pearsonr
import warnings


#This file contains the function to compute the gravity model and the degree constrained gravity model

def _checkSystematicBiasPosition(distances, graphModelref, graphModelComputed):
	excentricity = {}
	for node1 in graphModelref.getNodes():
		for node2 in graphModelref.getNodes():
			excentricity.setdefault(node1,0)
			excentricity[node1]+=distances(node1,node2)/len(graphModelref.getNodes())

	degreeBiasIn=[]
	degreeBiasOut=[]
	vectorExcentricity = []
	(InDegRef,OutDegRef) = graphModelref.getDegrees()
	(InDegComp,OutDegComp) = graphModelComputed.getDegrees()


	for node in graphModelref.getNodes():
		degreeBiasIn.append(round((InDegRef[node]-InDegComp[node])/InDegRef[node],ndigits=2))
		degreeBiasOut.append(round((OutDegRef[node]-OutDegComp[node])/InDegRef[node],ndigits=2))
		vectorExcentricity.append(excentricity[node])
	with warnings.catch_warnings():
		warnings.filterwarnings('error')
		try:
			corrIn = pearsonr(degreeBiasIn,vectorExcentricity)[0]
			corrOut = pearsonr(degreeBiasOut,vectorExcentricity)[0]
		except Warning: #if there is an empty array
			corrIn=0
			corrOut=0
	print("Correlations degrees/excentricity, IN: %s OUT:%s"%(corrIn,corrOut))



def _estimateEISsIN(EISsIN, EISsOUT, INs, OUTs, deterrencefunc, distances, normalized=False):
	#EIS : estimated Intrinsic Strenght
	newEISsIN = {}
	sumIn = sum(INs.values())

	#for each node
	for nodeDest in INs:
		#compute how many interaction it receives
		sumReceived =0.0
		for nodeSource in EISsOUT:
			sumReceived+=deterrencefunc(distances(nodeSource,nodeDest))*EISsOUT[nodeSource]*INs[nodeDest]/sumIn

		#modify its "intrinsic degree" to receive the right number according to reference
		newEISsIN[nodeDest] = INs[nodeDest] / sumReceived * INs[nodeDest]
		#print(nodeDest,sumReceived,INs[nodeDest],newEISsIN[nodeDest])


	return newEISsIN

def _estimateEISsOUT(EISsIN, EISsOUT, INs, OUTs, deterrencefunc, distances, normalized=False):
	# EIS : estimated Intrinsic Strenght
	newEISsOUT = {}
	sumOut = sum(OUTs.values())

	# for each node
	for nodeSource in OUTs:
		# compute how many interaction it receives
		sumSent = 0.0
		for nodeDest in EISsIN:
			sumSent += deterrencefunc(distances(nodeSource, nodeDest)) * EISsIN[nodeDest] * OUTs[nodeSource]/sumOut

		# modify its "intrinsic degree" to receive the right number according to reference
		newEISsOUT[nodeSource] = OUTs[nodeSource] / sumSent * OUTs[nodeSource]

	return newEISsOUT


def getSpatialNullModelExpertEtAl(originalNetwork,distances,roundDecimal,maximalDistance=100000,plot=False,printDebug=False):
	"""
	:param originalNetwork: the observed network for which we want the corresponding null model
	:param distances: a function that return the distance between two nodes of the provided graph
	:param roundDecimal: the rounding used to compute bins of the deterrence function. for a distance d=123.456 :
	 if roundDecimal=2, binned value is 123.46.
	 if roundDecimal=-2, binned value is 100
	:param maximalDistance: ignore in most cases, parameter of the deterence function to set an upper bound on the considered distances
	:param minValsBin: parameter of the deterrence function, minimum number of observations in a bin to consider it. (avoid abherent values for rare distances)
	:param plot: if True, plot the deterrence function before the doubly constrained process and at the end of the process
	:param printDebug: print at each step a trace of the current model: edit distance with original network, degree bias towards central nodes...
	slow down the process a lot, use only to understand or check that everything is going well.
	"""
	return getSpatialNullModel(originalNetwork,distances,roundDecimal,maximalDistance=maximalDistance,plot=plot,iterations=0,printDebug=printDebug)

def getSpatialNullModel(originalNetwork,distances,roundDecimal,maximalDistance=100000,minValsBin = 3,plot=False,iterations=5,printDebug=False):
	"""
	:param originalNetwork: the observed network for which we want the corresponding null model
	:param distances: a function that return the distance between two nodes of the provided graph
	:param roundDecimal: the rounding used to compute bins of the deterrence function. for a distance d=123.456 :
	 if roundDecimal=2, binned value is 123.46.
	 if roundDecimal=-2, binned value is 100
	:param maximalDistance: ignore in most cases, parameter of the deterence function to set an upper bound on the considered distances
	:param minValsBin: parameter of the deterrence function, minimum number of observations in a bin to consider it. (avoid abherent values for rare distances)
	:param plot: if True, plot the deterrence function before the doubly constrained process and at the end of the process
	:param iterations: number of iterations in the doubly constrained process
	:param printDebug: print at each step a trace of the current model: edit distance with original network, degree bias towards central nodes...
	slow down the process a lot, use only to understand or check that everything is going well.
	"""
	print("Computing the spatial null model with %s iterations"%(iterations))
	norm = False
	GraphModelOriginal = GraphModelAsnxGraph(originalNetwork)


	# compute in and out degrees
	INs = originalNetwork.in_degree(weight="weight")
	OUTs = originalNetwork.out_degree(weight="weight")



	# ---Print summary progress---------
	if printDebug:

		temporaryModel = ConfigurationModel(INs, OUTs)
		print("EDIT distance with a configuration model: %s" % temporaryModel.getEditDistWithGraphModel(GraphModelOriginal))

		(degreesConvergenceIN,degreesConvergenceOUT) = temporaryModel.getDegreesSimilarity(GraphModelOriginal)
		print("--convergence of degrees, IN: %s OUT: %s" % (degreesConvergenceIN, degreesConvergenceOUT))

		_checkSystematicBiasPosition(distances, GraphModelOriginal, temporaryModel)

	# --------------------------------------------


	print("computing the original deterrence function")
	deterrencefunc = deterrenceFunction()
	deterrencefunc.deterrenceFunctionEstimation(INs, OUTs, originalNetwork, distances,
													  roundDecimals=roundDecimal, maximalDistance=maximalDistance, plot=plot,minVals = minValsBin)

	# ---Print summary progress---------
	if printDebug or iterations==0:
		temporaryModel = GravityModel(INs,OUTs,deterrencefunc.getDeterrenceAtDistance,distances)
		print(
			"EDIT distance with Simple Gravity model: %s" % temporaryModel.getEditDistWithGraphModel(GraphModelOriginal))
		(degreesConvergenceIN, degreesConvergenceOUT) = temporaryModel.getDegreesSimilarity(GraphModelOriginal)
		print("--convergence of degrees, IN: %s OUT: %s" % (degreesConvergenceIN, degreesConvergenceOUT))
		_checkSystematicBiasPosition(distances, GraphModelOriginal, temporaryModel)

	# --------------------------------------------




	EISsIN=INs
	EISsOUT=OUTs

	for i in range(iterations):
		print("----------step %s" %(i))



		EISsIN = _estimateEISsIN(EISsIN, EISsOUT, INs, OUTs, deterrencefunc.getDeterrenceAtDistance, distances, normalized=norm)
		EISsOUT = _estimateEISsOUT(EISsIN, EISsOUT, INs, OUTs, deterrencefunc.getDeterrenceAtDistance, distances, normalized=norm)

		# ---Print summary progress---------
		if printDebug:

			temporaryModel = GravityModel(EISsIN, EISsOUT, deterrencefunc.getDeterrenceAtDistance, distances,
										  desiredInDegrees=INs, desiredOutDegrees=OUTs)
			print(
				"EDIT distance with DC gravity model: %s" % temporaryModel.getEditDistWithGraphModel(
					GraphModelOriginal))
			(degreesConvergenceIN, degreesConvergenceOUT) = temporaryModel.getDegreesSimilarity(GraphModelOriginal)
			print("--convergence of degrees, IN: %s OUT: %s" % (degreesConvergenceIN, degreesConvergenceOUT))
			_checkSystematicBiasPosition(distances, GraphModelOriginal, temporaryModel)

		# -------------------------------------
		deterrencefunc = deterrenceFunction()
		deterrencefunc.deterrenceFunctionEstimation(EISsIN, EISsOUT, originalNetwork, distances,
													roundDecimals=roundDecimal, maximalDistance=maximalDistance,
													plot=False,minVals = minValsBin)

	if plot:
		deterrencefunc = deterrenceFunction()
		deterrencefunc.deterrenceFunctionEstimation(EISsIN, EISsOUT, originalNetwork, distances,
													roundDecimals=roundDecimal, maximalDistance=maximalDistance,
													plot=plot,minVals = minValsBin)


	temporaryModel = GravityModel(EISsIN, EISsOUT, deterrencefunc.getDeterrenceAtDistance, distances,
								  desiredInDegrees=INs, desiredOutDegrees=OUTs)


	return(temporaryModel,deterrencefunc)
