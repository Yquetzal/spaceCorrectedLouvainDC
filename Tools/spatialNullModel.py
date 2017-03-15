from .deterrenceFunction import *
from .utilSpatialNullM import *
from .GraphModels import *
import math
from scipy.stats.stats import pearsonr
import warnings


#This file contains the function to compute the gravity model and the degree constrained gravity model

def checkSystematicBiasPosition(distances,graphModelref,graphModelComputed):
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



def estimateEISsIN(EISsIN,EISsOUT,INs,OUTs,deterrencefunc,distances,normalized=False):
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

def estimateEISsOUT(EISsIN,EISsOUT,INs,OUTs,deterrencefunc,distances,normalized=False):
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
	return getSpatialNullModel(originalNetwork,distances,roundDecimal,maximalDistance=maximalDistance,plot=plot,iterations=0,printDebug=printDebug)

def getSpatialNullModel(originalNetwork,distances,roundDecimal,maximalDistance=100000,plot=False,iterations=5,printDebug=False):
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

		checkSystematicBiasPosition(distances, GraphModelOriginal, temporaryModel)

	# --------------------------------------------


	print("computing the original deterrence function")
	deterrencefunc = deterrenceFunction()
	deterrencefunc.deterrenceFunctionEstimation(INs, OUTs, originalNetwork, distances,
													  roundDecimals=roundDecimal, maximalDistance=maximalDistance, plot=plot)

	# ---Print summary progress---------
	if printDebug or iterations==0:
		temporaryModel = GravityModel(INs,OUTs,deterrencefunc.getDeterrenceAtDistance,distances)
		print(
			"EDIT distance with Simple Gravity model: %s" % temporaryModel.getEditDistWithGraphModel(GraphModelOriginal))
		(degreesConvergenceIN, degreesConvergenceOUT) = temporaryModel.getDegreesSimilarity(GraphModelOriginal)
		print("--convergence of degrees, IN: %s OUT: %s" % (degreesConvergenceIN, degreesConvergenceOUT))
		checkSystematicBiasPosition(distances, GraphModelOriginal, temporaryModel)

	# --------------------------------------------




	EISsIN=INs
	EISsOUT=OUTs

	for i in range(iterations):
		print("----------step %s" %(i))



		EISsIN = estimateEISsIN(EISsIN, EISsOUT, INs, OUTs, deterrencefunc.getDeterrenceAtDistance,distances, normalized=norm)
		EISsOUT = estimateEISsOUT(EISsIN, EISsOUT, INs, OUTs, deterrencefunc.getDeterrenceAtDistance, distances,normalized=norm)

		# ---Print summary progress---------
		if printDebug:

			temporaryModel = GravityModel(EISsIN, EISsOUT, deterrencefunc.getDeterrenceAtDistance, distances,
										  desiredInDegrees=INs, desiredOutDegrees=OUTs)
			print(
				"EDIT distance with DC gravity model: %s" % temporaryModel.getEditDistWithGraphModel(
					GraphModelOriginal))
			(degreesConvergenceIN, degreesConvergenceOUT) = temporaryModel.getDegreesSimilarity(GraphModelOriginal)
			print("--convergence of degrees, IN: %s OUT: %s" % (degreesConvergenceIN, degreesConvergenceOUT))
			checkSystematicBiasPosition(distances, GraphModelOriginal, temporaryModel)

		# -------------------------------------
		deterrencefunc = deterrenceFunction()
		deterrencefunc.deterrenceFunctionEstimation(EISsIN, EISsOUT, originalNetwork, distances,
													roundDecimals=roundDecimal, maximalDistance=maximalDistance,
													plot=False)

	if plot:
		deterrencefunc = deterrenceFunction()
		deterrencefunc.deterrenceFunctionEstimation(EISsIN, EISsOUT, originalNetwork, distances,
													roundDecimals=roundDecimal, maximalDistance=maximalDistance,
													plot=plot)


	temporaryModel = GravityModel(EISsIN, EISsOUT, deterrencefunc.getDeterrenceAtDistance, distances,
								  desiredInDegrees=INs, desiredOutDegrees=OUTs)


	return(temporaryModel,deterrencefunc)
