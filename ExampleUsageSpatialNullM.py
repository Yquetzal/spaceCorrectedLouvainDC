from Tools.spatialNullModel_CBJ import *
from Tools.utilSpatialNullM import *

distanceFile = ""
flowFile = ""



# compute the spatial null model of a provided network with distances between nodes,
# and find spatial independant communities based on this NM
if __name__ == "__main__":

	plot=False
	#exemples of files to handle
	flowFile = "files/Global2011.txt"
	distanceFile = "files/distancesLyon.txt"

	#obtain the observed network as a networkx object
	observedNetwork = readVelovFileAsNetworkX(flowFile)

	#obtain the distance between nodes as a function. In little example like this one, precomputed, but could be computed on the fly
	distancesBetweenNodes = Distances()
	distancesBetweenNodes.getDistanceFunctionVelov(distanceFile)


	undirectedObservedGraph = createUndirectedGraphWithSumWeight(observedNetwork)


	###-------------- COMPUTE COMMUNITIES WITH THE PROPOSED DEGREE CORRECTED NULL MODEL
	#compute the null model of a spatial network. roundDecimal = -2 : bin every 100 unity (meters)
	(nullModel, distFunction) = getSpatialNullModel(observedNetwork,distancesBetweenNodes.getDistanceBetween,roundDecimal=-2,iterations=5,plot=plot)

	#run community detection using computed null model
	graphOfNullModel = createnxGraphFromGraphModel(nullModel)

	DCcommunitiesLevel1 = computeCommunityDetectionUsingRefNullModel(undirectedObservedGraph,graphOfNullModel,firstLevel=True)
	DCcommunitiesLastLevel = computeCommunityDetectionUsingRefNullModel(undirectedObservedGraph,graphOfNullModel,firstLevel=False)




	###-------------- COMPUTE COMMUNITIES WITH THE ORIGINAL DEGREE CORRECTED NULL MODEL by Expert et al.

	#run community detection using normal gravity model
	(nullModel, distFunction) = getSpatialNullModelExpertEtAl(observedNetwork,distancesBetweenNodes.getDistanceBetween,roundDecimal=-2,plot=plot)
	graphOfNullModel = createnxGraphFromGraphModel(nullModel)

	ExpertcommunitiesLevel1 = computeCommunityDetectionUsingRefNullModel(undirectedObservedGraph, graphOfNullModel,
																   firstLevel=True)
	ExpertcommunitiesLastLevel = computeCommunityDetectionUsingRefNullModel(undirectedObservedGraph, graphOfNullModel,
																	  firstLevel=False)



	#---------- COMPUTE COMMUNITIES WITH THE SIMPLE CONFIGURATION MODEL (ingore spatial effect)
	nullModel = ConfigurationModel(observedNetwork.in_degree(weight="weight"),observedNetwork.out_degree(weight="weight"))
	graphOfNullModel = createnxGraphFromGraphModel(nullModel)

	ConfigcommunitiesLevel1 = computeCommunityDetectionUsingRefNullModel(undirectedObservedGraph, graphOfNullModel,firstLevel=True)
	ConfigcommunitiesLastLevel = computeCommunityDetectionUsingRefNullModel(undirectedObservedGraph, graphOfNullModel,
																			firstLevel=False)


	printCommunitiesInVisualisationFormat([DCcommunitiesLevel1,DCcommunitiesLastLevel,ExpertcommunitiesLevel1,ExpertcommunitiesLastLevel,ConfigcommunitiesLevel1,ConfigcommunitiesLastLevel])


