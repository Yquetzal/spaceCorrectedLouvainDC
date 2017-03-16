import networkx as nx
from . import louvainModified as lvm

#this file contains a set of useful functions to handle graph transformations, community format transformations, reading files ...
def printCommunitiesInVisualisationFormat(communities):
	for n in communities[0]:
		toPrint=str(n)
		for c in communities:
			toPrint+="\t"+str(c[n])
		print(toPrint)


def createnxGraphFromGraphModel(graphModel):
	"""
	:param graphModel:
	:return:
	"""
	aGraph = nx.Graph()
	for source in graphModel.getNodes():
		for dest in graphModel.getNodes():
			aGraph.add_edge(source, dest, weight=0)

	for source in graphModel.getNodes():
		for dest in graphModel.getNodes():
			aGraph[source][dest]["weight"]+=graphModel.getExpectedEdges(source,dest)
	return aGraph

def createUndirectedGraphWithSumWeight(graphModel):
	aGraph = nx.Graph()
	for e in graphModel.edges():
		aGraph.add_edge(e[0],e[1],weight=0)

	for e in graphModel.edges(data=True):
		aGraph[e[0]][e[1]]["weight"]+=e[2]["weight"]
	return aGraph

def computeCommunityDetectionUsingRefNullModel(aGraph,aNullModel,firstLevel=True):

	print("----computing community detection----")
	if (firstLevel):
		communities = lvm.getPartitionAtSpecificLevel(aGraph, aNullModel)
	else:
		communities = lvm.best_partition(aGraph,aNullModel)

	#communities = {idToNode[k]:v for k,v in communities.iteritems()}

	return communities





def readVelovFileAsNetworkX(file):
	graph = nx.DiGraph()
	f = open(file)
	toAdd=[]
	for l in f:
		if l[0] != "#":
			elts = l.split("\t")
			(n1,n2) = elts[0].split("_")
			toAdd.append((n1,n2,{'weight':float(elts[1])}))
	graph.add_edges_from(toAdd) #add in batch
	nbtrips = graph.size(weight="weight")
	print("read graph with %s nodes, %s edges and %s trips" %(graph.number_of_nodes(),graph.number_of_edges(),nbtrips))
	return graph

class Distances():
	def __init__(self):
		self.allDistances = {}

	def getDistanceBetween(self,node1,node2):
		if not (node1,node2) in self.allDistances:
			raise Exception(" distance from %s to %s unknown" %(node1,node2))
		return self.allDistances[(node1,node2)]

	def getDistanceFunctionVelov(self,file):
		f = open(file)
		for l in f:
			if l[0] != "#":
				elts = l.split("\t")
				(n1, n2) = elts[0].split("_")
				self.allDistances[(n1,n2)]=float(elts[1])

