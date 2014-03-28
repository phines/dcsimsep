import sys
if sys.version_info < (2, 6):
	raise "must use Python 2.7 or greater"
	sys.exit(-1)
import math
import random
import os
from collections import defaultdict
import csv
import time
import json
import multiprocessing as mlt
import subprocess
import argparse
import networkx as nx
from networkx.readwrite import json_graph
#optional modules
try: 
	import powerlaw
	from pymongo import Connection 
	from pymongo.errors import ConnectionFailure 
except ImportError:
	print "Missing powerlaw and/or pymongo modules."

# python couplednetworks.py -1 -1 -1 config.json
parser = argparse.ArgumentParser(description='Coupled Network Model')
parser.add_argument('mpid', metavar='MATLAB_PID', type=int, nargs=1, help='process ID from MATLAB, '+
	'if running without MATLAB use -1 as a command line argument')
parser.add_argument('cfs_iter', metavar='CFS_iteration', type=int, nargs=1, help='iteration from MATLAB, '+
	'if running without MATLAB use -1 as a command line argument')
parser.add_argument('percent_removed', metavar='percent_removed', type=float, nargs=1, help='Percent of nodes removed '+
	'if running without MATLAB use -1 as a command line argument')
# parser.add_argument('config_name', metavar='config_name', type=str, nargs=1, help='Name and location of config.json file')
args = parser.parse_args()
mpid = args.mpid[0]
cfs_iter = args.cfs_iter[0]
percent_removed = args.percent_removed[0]
# config_name = args.config_name[0]

NUM_PROCS = mlt.cpu_count()
# config = json.load(open(config_name))
config = json.load(open('/Users/veneman/workspace/coupled-networks/source/config-cfs.json')) # read in the input configuration
k=config['k'] 	#average degree, <k>, must be a float, 4.0 is what Buldrev used
n=config['n']	#number of nodes, 50000 is what Buldyrev used, 2383 for Polish grid
#p is the fraction of nodes to keep, defined below
runs = config['runs']	# number of runs to make 
pValues = config['pValues']	# how many increments of p, looks like Buldyrev used 200 pValues in [0,1]
pMin = config['pMin']
pMax = config['pMax']
targeted = config['targeted'] # If this is true the highest degree nodes are removed first
outputRemovedNodesToDB = config['outputRemovedNodesToDB'] # For visualization, this writes a csv file with all the node removals at each p-value
outputRemovedNodes = config['outputRemovedNodes']
# Output files at each p-value or just a file at the end. This should probably be false... 
# unless you're concerned that you may need to restart a simulation mid run
logAllPValues = config['logAllPValues']
real = config['real'] # set to true when using CFS and not just topo networks
writeNetworksOut = config['writeNetworksOut'] # writes a json formatted copy of the networks to file for each run. Caution, will produce lots of files
showPlot = config['showPlot']
# if randomRewireProb is -1, base the comms network on a configuration model...
# otherwise the comms will be randomly rewired with this probability
# For scale free replication lambda/gamma = ~3 when randomRewireProb = 0.16...
# and ~2.7 when randomRewireProb = 0.31 and ~2.3 with rewireProb = 0.57
randomRewireProb = config['randomRewireProb'] 
# the fraction of the number of components in the giant component under which failure is declared
# for replication this should be 0 which gets changed to a threshold of 1 node
gCThreshold = config['gCThreshold']	
gridGCThreshold = config['gridGCThreshold'] # the fraction of the number of components in the giant component under which blackout is declared
shuffleNetworks = config['shuffleNetworks'] # randomly renumber the nodes on the networks, for replication this should be True
generateEachRun = config['generateEachRun'] # generate a new network for each run, True, or keep the same for all, False
hostOS = config['hostOS'] # for the compiled MATLAB this determines where to find cmp_dcsimsep, 0 for Mac, 1 for Linux
networkType = config['networkType']	# must be either 'SF', 'RR', 'ER', 'CFS-SW', and others TBD in createNetworks()
startWithComms = config['startWithComms']
verbose = config['verbose']
outagesFromFile = config['outagesFromFile']
minOutagePoint = config['minOutagePoint']
maxOutagePoint = config['maxOutagePoint']
numOutagesInFile = config['numOutagesInFile']
networkFromFile = config['networkFromFile']
findScalingExponent = config['findScalingExponent']

ep=k/(n-1)	#probability for edge creation
# print "Probability for edge creation: " + str(ep)
nodes = range(1,n,1)
pRange = pMax - pMin
home = os.path.expanduser("~")
threshold = gCThreshold * n
if gCThreshold == 0:
	threshold = 1

network_filename = '/tmp/comms_net_' + str(mpid) + '.edgelist'

if(mpid == -1):
	relpath = ".."
else:
	relpath = "/Users/veneman/workspace/coupled-networks"

def configType(shuffleNetworks, generateEachRun):
	gen = "0"
	shuffle = "0"
	if generateEachRun == True:
		gen = "1"
	if shuffleNetworks == True:
		shuffle = "1"
	return gen + shuffle


networkLabel = (networkType + '{nodes}nodes_{rw}rw_{thresh}thresh_{sruns}Avg'.format(nodes = n, sruns = runs, rw = randomRewireProb, thresh = threshold) + 
	'_Config_%s' % configType(shuffleNetworks, generateEachRun) + "_") # this gets put in the file names

if verbose == True: 
	print ("Network type: " + networkLabel + ", comms threshold: " + str(gCThreshold) +
	", grid threshold: " + str(gridGCThreshold))
	print "Number of processors " + str(NUM_PROCS)

# Only import matplotlib if you want to generate plots. It uses about 20MB RAM.
if showPlot == True:
	try:
		import matplotlib.pyplot as plt
	except ImportError:
		raise "matplotlib not found, can't show plots"
		sys.exit(-1)

# holds the p values used
x = [] 	
# holds whether the giant mutually connected component exists for each p value as: 
# [p,exists?] e.g. [0.01, 1] for p=0.01 and the GMCC existing
allY = []	

'''
Go through the networks and remove links per the process given in 
Buldyrev 2010, Catastrophic cascade of failures in interdependent networks.
"Each node in network A depends on one and only one node in network B, and 
vice versa. One node from network A is removed ('attack'). b, Stage 1: a 
dependent node in network B is also eliminated and network A breaks into 
three a1-clusters, namely a11, a12 and a13. c, Stage 2: B-links that link 
sets of B-nodes connected to separate a1-clusters are eliminated and network B
breaks into four b2-clusters, namely b21, b22, b23 and b24. d, Stage 3: 
A-links that link sets of A-nodes connected to separate b2-clusters are 
eliminated and network A breaks into four a3- clusters, namely a31, a32, a33 
and a34. These coincide with the clusters b21, b22, b23 and b24, and no 
further link elimination and network breaking occurs. Therefore, each 
connected b2-cluster/a3-cluster pair is a mutually connected cluster and the 
clusters b24 and a34, which are the largest among them, constitute the giant 
mutually connected component."
'''
def removeLinks(networkA, networkB, swapNetworks, iteration):
	# Get the connected components of each network.
	gMCCA = nx.connected_component_subgraphs(networkA)
	gMCCB = nx.connected_component_subgraphs(networkB)

	# Find which network has the fewest number of connected component subgraphs
	numSubnets = len(gMCCA) if len(gMCCA) < len(gMCCB) else len(gMCCB)

	identicalSubnets = 0
	swapAgain = swapNetworks
	nodesToRemove = []

	if numSubnets==0:
		swapAgain = 0
		return [swapAgain,nodesToRemove]

	# Go through each subnet
	for j in range(0,numSubnets,1):
		# check against threshold 
		# if len(gMCCA[0].nodes()) < threshold or len(gMCCB[0].nodes()) < threshold:
		# 	return [0,nodesToRemove]
		# Make the nodes sets so you can do .difference() on them	
		sA = set(gMCCA[j].nodes())	
		sB = set(gMCCB[j].nodes())
		# Find the elements in B but not in A
		outNodesB = sB.difference(sA)
		if outputRemovedNodes == True:
			if len(outNodesB) > 0:
				nodesToRemove.extend(list(outNodesB))
		if verbose == True: print str(len(nodesToRemove)) + " nodesToRemove in subnet " + str(j+1) + " of " + str(numSubnets) + "."
		if verbose == True: print "Number elements in B but not in A " + str(len(outNodesB))
		if len(outNodesB) == 0:
			if verbose == True: print "No differences in subnet " + str(j+1)
			identicalSubnets += 1
			# Check to see if the nodes in each subgraph are the same
			# If they are then no more link removal is needed and we can return
			if identicalSubnets == numSubnets:
				swapAgain = 0
				if verbose == True:
					print ("All the same" + ", iteration: " + str(iteration) + ", gMCCA size " + 
						str(len(gMCCA)) + " gMCCB size " + str(len(gMCCB)))
				return [swapAgain,nodesToRemove]
		for k in range(0,len(outNodesB),1):
			# For each node that's not in this MCC subnet of A and B networks
			nodeInSubgraph = list(outNodesB)[k]
			if verbose == True:
				print ("Node that's different in subgraph k:" + str(k) + " is: " + 
					str(nodeInSubgraph) + ", iteration: " + str(iteration))
			# Find their neighbors
			nodeInSubgraphNeighbors = nx.neighbors(gMCCB[j],nodeInSubgraph)
			# Add to removed nodes list (for visualization only)
			for m in range(0,len(nodeInSubgraphNeighbors),1):
				# Remove all the edges going to those neighbor nodes
				networkB.remove_edges_from([(nodeInSubgraph,nodeInSubgraphNeighbors[m])]) 
		# At the last subnet increment swapAgain so that the networks get swapped next time
		if j==numSubnets-1:
			swapAgain += 1
	if verbose == True:
		print "swapAgain " + str(swapAgain) + " netA size " + str(len(networkA.nodes())) + " netB size " +str(len(networkB.nodes()))
		print "gMCCA size " + str(len(gMCCA)) + " gMCCB size " + str(len(gMCCB))
	return [swapAgain,nodesToRemove]

def attackNetwork(run, networks):
	# Change the internal state of the random generator for each run
	newstate = random.randint(2,n*2)
	random.jumpahead(newstate)
	if real == False: print "Run number " + str(run+1) + " of " + str(runs) + ". Random state " + str(newstate)

	runstart = time.time()
	y = []
	generated = False

	# Create the networks 
	if generateEachRun == True:
		[networkA, networkB] = createNetworks(networkType)
	else: 
		networkA = networks[0]
		networkB = networks[1]
	# print "networks created, networkA: " + nx.info(networkA)
	# print "networks created, networkB: " + nx.info(networkB)
	
	dbh = -1 # initialize the database handle to something in case it's *not* used

	if writeNetworksOut == True:
		writeNetworks(run,networkA,networkB)
	if outputRemovedNodesToDB == True:
		try: # Connect to MongoDB
			c = Connection(host="localhost", port=27017) 
		except ConnectionFailure, e: 
			sys.stderr.write("Could not connect to MongoDB: %s" % e)
			# sys.exit(1)
		dbh = c["runs"]	# sets the db to "runs" 
		assert dbh.connection == c

	checkForFailure(networkA, networkB, dbh, run, y)

	d = defaultdict(list)
	for key, value in y:
		d[key].append(value)
	average = [float(sum(value)) / len(value) for key, value in d.items()]
	x = [key for key, value in d.items()]
	writeOutput(x, average, run)
	if verbose == True: print "Comms run time was " + str(time.time() - runstart) + " seconds"
	return y

def checkForFailure(networkA, networkB, dbh, run, y):
	if real == True:
		networkAcopy = networkA.copy()
		networkBcopy = networkB.copy()
		numNodesAttacked = int(math.floor(percent_removed * n))
		busessep = [] # bus separations read from MATLAB
		nodesAttacked = []
		if targeted == False and outagesFromFile == False:
			nodesAttacked = random.sample(networkA.nodes(),numNodesAttacked)
		elif targeted == True and outagesFromFile == False:
			nodesAttacked = targetedNodes(networkA,numNodesAttacked)
		elif targeted == False and outagesFromFile == True and startWithComms == True:
			# this is only used if startWithComms is true
			nodesAttacked = getNodesFromFile(run+1,1-percent_removed,networkA)	#TODO fix this
		elif targeted == False and outagesFromFile == True and startWithComms == False:
			nodesAttacked = []
		else:
			sys.stderr.write("Unknown configuration of inputs: targeted and outagesFromFile")
			sys.exit(-1) 
			

		blackout = False

		comm_status_filename = '/tmp/comm_status_' + str(mpid) + '.csv'
		grid_status_filename = '/tmp/grid_status_' + str(mpid) + '.csv'
		# grid_status_filename = '/tmp/grid_status_test.csv'  
		# print "\t !!!!!!! Using test grid status file !!!!!!!!!!"
		nodes = range(1,n+1,1) # for interfacing with MATLAB start at 1

		if startWithComms == True:
			if verbose == True: print "Starting with comms, percent_removed " + str(percent_removed) + ", numNodesAttacked " + str(numNodesAttacked)
			# remove nodes from the comms network
			networkAcopy.remove_nodes_from(nodesAttacked)
			#networkBcopy is the copy of the grid that's used in MATLAB, remove the bus separations on it for comparing with networkAcopy
			networkBcopy.remove_nodes_from(nodesAttacked)
		else:
			#read grid status and remove nodes accordingly
			try:
				with open(grid_status_filename, 'rb') as f:
					try:
						reader = csv.DictReader(f)
						for item in reader:
							try:
								# print str(item['status']) + " " + str(item['bus'])
								if int(item['status']) == 0:
									busessep.append(int(item['bus']))
							except ValueError:
								# print "No bus separations"
								pass
					except:
						print "CSV reader error. Check headers in file " + grid_status_filename + " and check for Unix line endings in that file."
			except:
				print "*************** -> Missing grid status file <- ****************"
			networkBcopy.remove_nodes_from(busessep)
			# Remove the links on networkA, comms, that go to the bus separations on networkB, power
			networkAcopy.remove_nodes_from(busessep)
			nodesAttacked = busessep
			numNodesAttacked = len(busessep)
			if verbose == True: print "Starting with grid, number of bus separations: " + str(len(busessep))
		if verbose == True:
			print ("Subgraphs in networkB: " + str(len(nx.connected_component_subgraphs(networkBcopy))) + 
				", subgraphs in networkA: " + str(len(nx.connected_component_subgraphs(networkAcopy))))
		
		nodes_out_pre = n-len(networkAcopy.nodes())
		if verbose == True: print "***Total nodes out PRE comms removal: " + str(nodes_out_pre)
		
		swapNetworks = 1
		
		# while nodes_out_pre != nodes_out_post:
		result = removeLinks(networkBcopy, networkAcopy, swapNetworks, cfs_iter)
		
		swapNetworks = result[0]
		if len(result[1]) == 0:
			if verbose == True: print "\t ######## No nodes removed #########"
			pass
		else:
			nodesAttacked.extend(result[1])
			nodesAttacked = set(nodesAttacked) # remove duplicates
			nodesAttacked = list(nodesAttacked) # convert back to a list
			# print ("Additional number of nodes out: " + str(len(result[1])) + 
			# 	", total nodes out: " + str(len(nodesAttacked)))
		networkAcopy.remove_nodes_from(nodesAttacked)
		nodes_out_post = n-len(networkAcopy.nodes())
		if verbose == True: print ">>>Total nodes out POST comms removal: " + str(nodes_out_post)


		with open(comm_status_filename, 'w') as f:
			try:
				writer = csv.writer(f)
				writer.writerow(['bus', 'status']) # Header
				status = 0
				for item in nodes:
					if item not in nodesAttacked:
						status = 1
					else:
						status = 0
					writer.writerow([item,status])
			except:
				print "CSV writer error"
		
		# with open(network_filename, 'w') as f:	#TODO, why is this being written?
		# 	try:
		# 		nx.write_edgelist(networkAcopy,f)
		# 		if verbose == True: print "Wrote to file at: " + network_filename
		# 	except Exception, e:
		# 		print "Network write to edgelist error"
		# 		raise
		del result # free up memory

	else:
		for i in range(0,pValues+1,1):
			# innerrunstart = time.time()
			#copy the networks so that all node removals are done on the same network layout
			#new network layouts will be generated for each run if generateEachRun is true
			networkAcopy = networkA.copy()
			networkBcopy = networkB.copy()

			p = i/float(pValues)
			p = p/(1/float(pRange))+pMin
			randomRemovalFraction=1-p	# Fraction of nodes to remove
			
			numNodesAttacked = int(math.floor(randomRemovalFraction * n))
			# print "p " + str(p) + ", numNodesAttacked " + str(numNodesAttacked)
			if targeted == False and outagesFromFile == False:
				nodesAttacked = random.sample(networkA.nodes(),numNodesAttacked)
			elif targeted == True and outagesFromFile == False:
				nodesAttacked = targetedNodes(networkA,numNodesAttacked)
			elif targeted == False and outagesFromFile == True:
				nodesAttacked = getNodesFromFile(run+1,p,networkAcopy)
			else:
				sys.stderr.write("Unknown configuration of inputs: targeted and outagesFromFile")
				sys.exit(1) 

			blackout = False
			#print "<><><><>Starting, iteration " + str(i) + ", " + str(len(nodesAttacked)) + " nodes removed<><><><>"
			# remove the attacked nodes from both networks
			networkAcopy.remove_nodes_from(nodesAttacked)
			networkBcopy.remove_nodes_from(nodesAttacked)
			if verbose == True:
				print ("Number of nodes attacked: " + str(numNodesAttacked) + 
				", fraction: " + str(randomRemovalFraction)) #+
				# ", nodes remaining in networkA: " + str(len(networkAcopy.nodes())))
			# Track whether the subnets are the same, if they're not we need to 
			# swap them and run removeLinks again.
			swapNetworks = 1
			if outputRemovedNodesToDB == True:
				writeRemovedNodes(run,p,swapNetworks,nodesAttacked,dbh,"A",networkAcopy)
				# print "Writing removed nodes: " + str(outputRemovedNodesToDB)
			while swapNetworks !=0:
				if swapNetworks%2 == 1:
					result = removeLinks(networkAcopy, networkBcopy, swapNetworks, i)
					swapNetworks = result[0]
					if swapNetworks != 0 and outputRemovedNodesToDB == True:
						newNodesAttacked = list(set(result[1]) - set(nodesAttacked)) # get only the new nodes removed
						writeRemovedNodes(run,p,swapNetworks,newNodesAttacked,dbh,"A",networkAcopy)
						nodesAttacked.extend(result[1])
						nodesAttacked = set(nodesAttacked) # remove duplicates
						nodesAttacked = list(nodesAttacked) # convert back to a list
				else:
					result = removeLinks(networkBcopy, networkAcopy, swapNetworks, i)
					swapNetworks = result[0]
					if swapNetworks != 0 and outputRemovedNodesToDB == True:
						newNodesAttacked = list(set(result[1]) - set(nodesAttacked)) # get only the new nodes removed
						writeRemovedNodes(run,p,swapNetworks,newNodesAttacked,dbh,"B",networkBcopy)
						nodesAttacked.extend(result[1])
						nodesAttacked = set(nodesAttacked) # remove duplicates
						nodesAttacked = list(nodesAttacked) # convert back to a list
			

			gMCCA = nx.connected_component_subgraphs(networkAcopy)
			
			# print "Number of CC sub-graphs: " + str(len(gMCCA))

			gMCCA0 = []
			if not gMCCA:	# check if gMCCA is empty
				gMCCASize = 0
			else:
				gMCCA0= gMCCA[0]	# gMCCA of the new network
				gMCCASize = len(gMCCA0)

			gMCCB = nx.connected_component_subgraphs(networkBcopy)
			if not gMCCB:	# check if gMCCB is empty
				gMCCBSize = 0
			else:
				gMCCB0= gMCCB[0]	# gMCCB of the new network
				gMCCBSize = len(gMCCB0)

			#print "Nodes in gMCCA: " + str(gMCCASize) + ", nodes in gMCCB: " + str(gMCCBSize) + ". Blackout? " + str(blackout)
			# does the giant MCC exist and was there no blackout?
			# if gMCCASize > 1 and blackout == False:
			if gMCCASize > threshold and blackout == False:
				y.append([p,1])
				# y=[p,1]
			else:
				y.append([p,0])
				# y=[p,0]
			# print "y is: " + str(y)
			# print "Inner run time was " + str(time.time() - innerrunstart) + " seconds"

def logResult(result):
	allY.extend(result)

def writeNetworks(run,networkA,networkB):
	data = json_graph.node_link_data(networkA)
	filename = 'networkA_' + networkLabel + '_' + str(n) + 'nodes_run' + str(run) + '.json'
	path = relpath + "/output/" + filename
	completeName = os.path.abspath(path)
	with open(completeName, 'w') as f:
		json.dump(data,f)
	data = json_graph.node_link_data(networkB)
	filename = 'networkB_' + networkLabel + '_' + str(n) + 'nodes_run' + str(run) + '.json'
	path = relpath + "/output/" + filename
	completeName = os.path.abspath(path)
	with open(completeName, 'w') as f:
		json.dump(data,f)

def writeRemovedNodes(run,p,swapNetworks,nodesAttacked,dbh,networkName,network):
	filename = 'nodesRemoved_' + networkLabel + '_' + str(n) + 'nodes_run' + str(run) +'.csv' # + '_p_%.2f'%p + '.csv'
	path = relpath + "/output/" + filename
	completeName = os.path.abspath(path)
	# print(nodesAttacked)
	writeToCSV = False
	if writeToCSV == True:
		try:
			with open(completeName):
				with open(completeName, 'a') as f:	#append to the file if it already exists
					writer = csv.writer(f)
					# create a list of tuples and append that to the file
					for row in [(swapNetworks, nodesAttacked[i], p) for i in range(len(nodesAttacked))]:
						writer.writerow(row)
					pass
		except IOError:	# create the file and write the header if it doesn't exist yet
			with open(completeName, 'w') as f:
				writer = csv.writer(f)
				writer.writerow(['Step','Node','p'])
				for row in [(swapNetworks, nodesAttacked[i], p) for i in range(len(nodesAttacked))]:
					writer.writerow(row)

	experiment = networkType + "_" + str(n) # This is the collection name

	nodes_remaining = list(set(nodes) - set(nodesAttacked))
	nodes_removed = {
		"run" : run,
		"network" : networkName,
		"step" : swapNetworks,
		"nodes" : nodes_remaining,
		"nodesAttacked" : nodesAttacked,
		"p" : p,
		"size" : len(nodesAttacked)
	}

	# If it's the first run put in the 
	if run == 0 and swapNetworks == 1:
		nodes_removed["initial_network"] = json_graph.node_link_data(network)
		nodes_removed["config"] = config

	dbh[experiment].insert(nodes_removed, safe=True) # sets the collection to experiment

def targetedNodes(networkA, numNodesAttacked):
	degree_sequence =nx.degree(networkA).values()
	node_by_deg=sorted(zip(degree_sequence,networkA),reverse=True)
	nodesAttacked = []

	for i in range(0,numNodesAttacked,1):
		nodesAttacked.append(node_by_deg[i][1])

	return nodesAttacked

def getNodesFromFile(run, p_point, networkA):
	nodesAttacked = []
	file_name = relpath + '/data/node-removals/bus_outage_' + str(run) + '.csv'
	p_column = int(numOutagesInFile - 10000*(p_point-minOutagePoint)/numOutagesInFile) # 50-10000*(0.995-0.75)/50 = 1
	# print "p_column " + str(p_column)
	if p_column > numOutagesInFile or p_column == 0:
		if verbose == True: print "No outages found for p_point " + str(p_point) + " at column " + str(p_column) + " creating outage instead"
		randomRemovalFraction=1-p_point	# Fraction of nodes to remove
		numNodesAttacked = int(math.floor(randomRemovalFraction * n))
		nodesAttacked = random.sample(networkA.nodes(),numNodesAttacked)
		if verbose == True: print "size nodesAttacked list " + str(len(nodesAttacked)) + " for run " + str(run) + ", p_point " + str(p_point) 
		return nodesAttacked
	try:
		with open(file_name, 'rb') as f:
			data = csv.DictReader(f)
			for next_row in data:
				item = int(next_row[str(p_column)])
				if item != 0:
					nodesAttacked.append(item)
			if verbose == True: print "size nodesAttacked list " + str(len(nodesAttacked)) + " for run " + str(run) + ", p_point " + str(p_point) 
	except Exception, e:
				print "Node outage read error: " + str(e)
				raise
	# print "nodesAttacked " + str(nodesAttacked) + " run " + str(run) + ", p_point " + str(p_point) 
	return nodesAttacked

def getNetworkFromFile():
	#network_ER_2383nodes_2.422deg_12345seed_0edges-skipped_1subgraphs.edgelist
	#network_SF_2383nodes_2.422deg_0.9rw_12345seed_0edges-skipped_1subgraphs.edgelist
	networkName = "network_" + networkType + "_" + str(n) + "nodes_2.422deg_"
	if networkType == "ER" or networkType == "RR":
		networkName = networkName + "12345seed_0edges-skipped_1subgraphs.edgelist"
	elif networkType == "SF":
		networkName = networkName + str(randomRewireProb) +  "rw_12345seed_0edges-skipped_1subgraphs.edgelist"
	elif networkType == "CFS-SW" and randomRewireProb == 0.5:
		networkName = networkName + str(randomRewireProb) +  "rw_12345seed_211edges-skipped_1subgraphs.edgelist"
	elif networkType == "CFS-SW" and randomRewireProb == 0.9:
		networkName = networkName + str(randomRewireProb) +  "rw_12345seed_1235edges-skipped_1subgraphs.edgelist"
	elif networkType == "CFS-SW" and randomRewireProb == -1.0:
		networkName = networkName + str(randomRewireProb) +  "rw_12345seed_0edges-skipped_32subgraphs.edgelist"
	elif networkType == "CFS-SW" and (randomRewireProb == 0.1 or randomRewireProb == 0.0):
		networkName = networkName + str(randomRewireProb) +  "rw_12345seed_0edges-skipped_1subgraphs.edgelist"
	else:
		print "Unknown network type for " + networkName + " and rw=" + str(randomRewireProb)
		sys.exit(-1)
		
		
	
	path = relpath + "/data/networks/generated/" + networkName
	try:
		network=nx.read_edgelist(path, nodetype=int)
	except Exception, e:
		print "Node outage read error: " + str(e)
		raise

	return network

def makeComms(networkB,copy):
	if randomRewireProb == -1:
		degree_sequence=sorted(nx.degree(networkB).values(),reverse=True)
		networkA=nx.configuration_model(degree_sequence)
	else:
		if copy == True:
			networkA=networkB.copy()
		else:
			networkA = networkB
		nodes = []
		targets = []
		for i,j in networkA.edges():
			nodes.append(i)
			targets.append(j)

		# rewire edges from each node, adapted from NetworkX W/S graph generator
		# http://networkx.github.io/documentation/latest/_modules/networkx/generators/random_graphs.html#watts_strogatz_graph
		# no self loops or multiple edges allowed
		for u,v in networkA.edges(): 
			if random.random() < randomRewireProb:
				w = random.choice(nodes)
				# Enforce no self-loops or multiple edges
				while w == u or networkA.has_edge(u, w): 
					w = random.choice(nodes)
				# print "node: " + str(u) + ", target: " + str(v)
				networkA.remove_edge(u,v)  
				networkA.add_edge(u,w)
	if copy == True:
		mapping=dict(zip(networkA.nodes(),range(1,2384)))#2384 for Polish, 4942 for western. relabel the nodes to start at 1 like networkA
		networkA=nx.relabel_nodes(networkA,mapping)
	return networkA

def createNetworks(networkType):
	if networkFromFile == True:
		networkA = getNetworkFromFile()
		path = relpath + "/data/power-grid/case2383wp.edgelist"
		networkB=nx.read_edgelist(path, nodetype=int)
		mapping=dict(zip(networkB.nodes(),range(1,2384))) #renumber the nodes to start at 1 for MATLAB
		networkB=nx.relabel_nodes(networkB,mapping)
	else:
		if networkType == 'ER':
			# Erdos-Renyi random graphs
			networkA=nx.gnp_random_graph(n,ep)
			networkB=nx.gnp_random_graph(n,ep)
		elif networkType == 'RR':
			# random regular graphs
			networkA=nx.random_regular_graph(4,n)
			networkB=nx.random_regular_graph(4,n)
		elif networkType == 'SF':
			# Scale free networks
			# m==2 gives <k>==4, for this lambda/gamma is always 3
			networkA=nx.barabasi_albert_graph(n,2)
			networkB=nx.barabasi_albert_graph(n,2)
			if randomRewireProb != -1:
				networkA=makeComms(networkA,False)
				networkB=makeComms(networkB,False)
		elif networkType == 'CFS-SW':
			'''networkB is a topological representation of the power network. networkA
			 is a generated	communication network created through either a configuration 
			 model of the power grid or by randomly rewiring the power grid topology.
			'''
			path = relpath + "/data/power-grid/case2383wp.edgelist"
			networkB=nx.read_edgelist(path, nodetype=int)
			mapping=dict(zip(networkB.nodes(),range(1,2384))) #renumber the nodes to start at 1 for MATLAB
			networkB=nx.relabel_nodes(networkB,mapping)
			'''make the comms network'''
			if cfs_iter == 1 or real == False:
				if verbose == True: print "^^^^^ Making comm network ^^^^^^"
				networkA = makeComms(networkB,True)
			else:
				try:
					if verbose == True: print "reading from comm network at: " + network_filename
					networkA = nx.read_edgelist(network_filename, nodetype=int)
				except Exception as e:
					print e
					print "Read edgelist error"
					raise
		else:
			print 'Invalid network type: ' + networkType
			return []

	# Order the nodes of the networks randomly
	if real == False and shuffleNetworks == True:
		randomlist=range(1,n+1) # node must names start at 1 for MATLAB, remapping to start at 1 done above
		random.shuffle(randomlist)
		mapping=dict(zip(networkA.nodes(),randomlist))
		networkA=nx.relabel_nodes(networkA,mapping)
		random.shuffle(randomlist)
		mapping=dict(zip(networkB.nodes(),randomlist))
		networkB=nx.relabel_nodes(networkB,mapping)
		del mapping
		del randomlist
	elif real == True and shuffleNetworks == True and verbose == True:
		print "\t ******* -> Not shuffling networks when using a real grid with CFS <- ********"
	if generateEachRun == False and findScalingExponent == True:
		degseq=sorted(nx.degree(networkA).values(),reverse=False)
		fit = powerlaw.Fit(degseq,xmin=2.0)
		if verbose == True: print "Scaling exponent networkA " + str(fit.power_law.alpha)
		degseq=sorted(nx.degree(networkB).values(),reverse=False)
		fit = powerlaw.Fit(degseq,xmin=2.0)
		if verbose == True: print "Scaling exponent networkB " + str(fit.power_law.alpha)

	return [networkA,networkB]

def main():
	manager = mlt.Manager()
	networks = []
	if generateEachRun == False and real == False:
		# allow the generated networks to be shared by all processes/threads
		networks = manager.list(createNetworks(networkType))
	else:
		networks = createNetworks(networkType)
	if real == False:
		runstart = time.time()
		pool = mlt.Pool()
		for i in range(0,runs,1):
			pool.apply_async(attackNetwork, args = (i, networks, ), callback = logResult)
		
		pool.close()
		pool.join()

		d = defaultdict(list)
		for key, value in allY:
			d[key].append(value)
		average = [float(sum(value)) / len(value) for key, value in d.items()]
		x = [key for key, value in d.items()]

		print "x is: " + str(x)
		print "average is: " + str(average)

		writeOutput(x,average,runs)

		print "Run time was " + str(time.time() - runstart) + " seconds"
		
		if showPlot == True:
			#plot the result
			plt.plot(x,average,'bo')	# 'rx' for red x 'g+' for green + marker
			# show the plot
			plt.show()	
	else:
		attackNetwork(0, networks)

def writeOutput(x, average, run):
	if logAllPValues == False and run != runs:
		return

	# transform the data to enable writing to csv
	output = []
	output.append(x)
	output.append(average)
	output = zip(*output)	
	numRun = ''
	if run == runs:
		numRun = 'run_final_'
	else:
		numRun = 'run_' + str(run) + '_' 

	# write the output to file
	filename = networkLabel + str(n) + 'nodes_' + numRun + 'rewire' + str(randomRewireProb) + '_GCThresh' + str(gCThreshold) + '_GridGCThresh' + str(gridGCThreshold) + '.csv'
	path = relpath + "/output/" + filename
	completeName = os.path.abspath(path)
	with open(completeName, 'wb') as f:
		writer = csv.writer(f)
		writer.writerow(['p','Pinf'])
		for row in output:
			writer.writerow(row)

if __name__ == '__main__':
	try:
		main()
	except Exception as e:
		print e
		raise
