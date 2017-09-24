#this file contains an interface of a graph model and several implementations
#Convenient to handle frequent operations between typical graph models.

#A graph model can be defined by mathematic functions (e.g the configuration model, the gravity model ...
# for which the probability of an edge is defined according to the degree of nodes (and a deterrence functio for gravity model)

#A graph can also be defined edges by edges, in this case one can use a networkx object, or a dictionary of edges and their weights
#An observed graph is therefore a "particular case" of a graph model...



class GraphModel():
	def getSumEdges(self):
		raise Exception("Function not implemented")

	def getExpectedEdges(self, source, dest):
		raise Exception("Function not implemented")

	def getNodes(self):
		raise Exception("Function not implemented")

	def getEditDistWithGraphModel(self,referenceGraph,ignoreSelfLink=True):

		editDist = 0.
		for source in self.getNodes():
			for dest in self.getNodes():
				if not ignoreSelfLink or source!=dest:
					editDist+=abs(self.getExpectedEdges(source,dest)-referenceGraph.getExpectedEdges(source,dest))
		if referenceGraph.getSumEdges()!=self.getSumEdges():
			print("Be careful, edit distance between graph models with different nbEdges: %s and %s"%(self.getSumEdges(),referenceGraph.getSumEdges()))
		return editDist/2/referenceGraph.getSumEdges()

	def getDegrees(self):
		InDegreesSelf ={}
		OutDegreesSelf ={}
		for source in self.getNodes():
			for dest in self.getNodes():
				weightSelf = self.getExpectedEdges(source, dest)
				OutDegreesSelf.setdefault(source, 0)
				OutDegreesSelf[source] += weightSelf
				InDegreesSelf.setdefault(dest, 0)
				InDegreesSelf[dest] += weightSelf
		return(InDegreesSelf,OutDegreesSelf)

	def getDegreesSimilarity(self,referenceGraph):

		if self.getNodes()!=referenceGraph.getNodes():
			print(self.getNodes())
			print(referenceGraph.getNodes())
			raise Exception("Not the same nodes between self and compared graphs")

		(InDegreesSelf,OutDegreesSelf) = self.getDegrees()
		(InDegreesRef,OutDegreesRef) = referenceGraph.getDegrees()

		degreesSimilarityIN = 1 - sum([abs(InDegreesRef[n] - InDegreesSelf[n]) for n in self.getNodes()]) / 2 / referenceGraph.getSumEdges()
		degreesSimilarityOUT = 1 - sum([abs(OutDegreesRef[n] - OutDegreesSelf[n]) for n in self.getNodes()]) / 2 / referenceGraph.getSumEdges()
		return (degreesSimilarityIN,degreesSimilarityOUT)




class GraphModelAsnxGraph(GraphModel):
	def __init__(self,nxGraph):
		self.nxGraph=nxGraph

	def getNodes(self):
		return set(self.nxGraph.nodes())

	def getExpectedEdges(self, source, dest):
		#return self.nxGraph.get_edge_data(n1,n2,default=0)["weight"]
		try:
			return self.nxGraph[source][dest]["weight"]
		except KeyError:
			return 0

	def getSumEdges(self):
		return sum(dict(self.nxGraph.in_degree(weight="weight")).values())





class ConfigurationModel(GraphModel):
	def __init__(self,INdegrees,OUTdegrees):
		"""
		input are dictionaries
		"""
		self.inDegrees=INdegrees
		self.outDegrees=OUTdegrees
		self.totalDegrees = sum(dict(self.inDegrees).values())

	def getNodes(self):
		return set(self.inDegrees.keys())

	def getExpectedEdges(self, source, dest):
		return self.outDegrees[source]*self.inDegrees[dest]/self.totalDegrees

	def getSumEdges(self):
		return self.totalDegrees


class GravityModel(GraphModel):
	def __init__(self, INdegrees, OUTdegrees,deterrenceFunction,distances,desiredInDegrees=-1,desiredOutDegrees=-1):
		"""
		input are dictionaries
		"""
		self.inDegrees = INdegrees
		self.outDegrees = OUTdegrees
		self.desiredInDegrees=INdegrees
		self.desiredOutDegrees=OUTdegrees
		if desiredInDegrees!=-1:
			self.desiredInDegrees = desiredInDegrees
			self.desiredOutDegrees = desiredOutDegrees

		self.deterrence = deterrenceFunction
		self.distances=distances
		self.totalDesiredEdgesWeights = sum(self.desiredInDegrees.values())
		self.totalEdgesBeforeCorrection =0.
		for source in self.getNodes():
			for dest in self.getNodes():
				deterrenceAtDistance = deterrenceFunction(distances(source,dest))
				self.totalEdgesBeforeCorrection+=self.inDegrees[dest]*self.outDegrees[source]*deterrenceAtDistance

	def getNodes(self):
		return set(self.inDegrees.keys())

	def getExpectedEdges(self, source, dest):
		return self.outDegrees[source] * self.inDegrees[dest]* self.deterrence(self.distances(source,dest)) / self.totalEdgesBeforeCorrection*self.totalDesiredEdgesWeights

	def getSumEdges(self):
		return self.totalDesiredEdgesWeights


class GraphModelAsDictionary(GraphModel):
	"""
	This is a fast implementation for null models
	"""
	def __init__(self,expectedEdges,inDegrees,outDegrees):
		self.expectedEdges=expectedEdges
		self.inDegrees = inDegrees
		self.outDegrees=outDegrees

	def getNodes(self):
		return set(self.inDegrees.keys())

	def getExpectedEdges(self,source,dest):
		if not (source,dest) in self.expectedEdges:
			return 0
		return self.expectedEdges[(source,dest)]

	def getSumEdges(self):
		return sum(self.inDegrees.values())

	def getDegrees(self):
		return(self.inDegrees,self.outDegrees)