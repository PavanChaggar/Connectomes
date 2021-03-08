# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#
#     Results parsing script for parcellations
#     constructed using CMTK (parcellation) 
#     and FSL (tractography)
#     
#  CMTK: https://hub.docker.com/r/sebastientourbier/multiscalebrainparcellator/
#        docker pull sebastientourbier/multiscalebrainparcellator:v1.1.1
#
#  FSL: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation
#       
#       -or-
#       sudo apt update
#       sudo apt-get install fsl
#       sudo apt-get upgrade
#
#
#  Authors:
#  ================================================  
#     Pavan Chaggar       chaggar@maths.ox.ac.uk
#                    ----
#     Travis B. Thompson  thompsont@maths.ox.ac.uk
#                    ----
#     Alain Goriely       
#
# ---------------------------------------------------------
import sys
from os import path
from os import walk
import numpy as np
import matplotlib.pyplot as plt

# A simple progress bar
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
	"""
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
	percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
	filledLength = int(length * iteration // total)
	bar = fill * filledLength + '-' * (length - filledLength)
	print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
	# Print New Line on Complete
	if iteration == total:
		print()


# ----------------------------------------------------------
# implements a node that can generate its own xml tag
# ----------------------------------------------------------

class Node:

	def __init__(self):
		self.number = -1
		self.volume = -1.0

		# d4 (region type i.e. cortical or subcortical)
		self.region = ""
		# d5 (freesurfer name)
		self.fsname = ""
		# d6 (node name)
		self.nname = ""
		# d7 (left / right hemisphere)
		self.hemisphere = ""

		self.x = -1.0
		self.y = -1.0
		self.z = -1.0

		# the final xml string will be generated from the contents of
		# the variables
		self.xmlstr = ""

	# === === === ===
	def setnodenumber(self, strNumber):
		self.number = int(strNumber)
	# === === === ===
	def setlabels(self, strLabel):
		# There are two types of labels
		# 1. those that start with ctx-
		# 2. those that start with Left- or Right-
		hyp = int(strLabel.find('-',0))
		prefix = (strLabel[0:hyp]).strip()


		if prefix == "ctx":
			# this is a cortical region (d4)
			self.region = "cortical"

			# the next important label is either `lh' or `rh' followed by
			# another hyphen.  i.e. ctx-rh-[the rest]
			hemi = (strLabel[hyp+1:hyp+3]).strip()
			rgn = (strLabel[hyp+4:]).strip()

			# the freesurfer name (d5) is the back half of the string
			self.fsname = rgn

			# the node name (d6) is hemi.rgn (i.e. separated by a dot
			# and not a hyphen
			self.nname = hemi + "." + rgn

			if hemi == "lh": # it is either lh or rh
				self.hemisphere = "left"
			else:
				self.hemisphere = "right"
		else:
			# this is a subcortical region
			self.region = "subcortical"

			if prefix == "Left":
				self.hemisphere = "left"
			else:
				self.hemisphere = "right"

			# Now, we essentially want to keep this whole string
			# but we need to replace any instance of _ with -
			# for example replacing Right-Accumbens_area with
			# Right-Accumbens-area
			# --> save the label off to both d5 and d6
			self.fsname = strLabel.replace('_', '-')
			self.nname = self.fsname
	# === === === ===
	def setvolume(self, strVol):
		self.volume = float(strVol)
	# === === === ===
	def setcoords(self,x,y,z):
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)

	def genxml(self):
		# node id
		self.xmlstr  = "    <node id=\"" + str(self.number) + "\">\n"
		# coordinate keys
		self.xmlstr += "      <data key =\"d0\">" + str(self.x) + "</data>\n"
		self.xmlstr += "      <data key =\"d1\">" + str(self.y) + "</data>\n"
		self.xmlstr += "      <data key =\"d2\">" + str(self.z) + "</data>\n"
		# node number (should be the same as node id)
		self.xmlstr += "      <data key =\"d3\">" + str(self.number) + "</data>\n"
		# node region (d4)
		self.xmlstr += "      <data key =\"d4\">" + str(self.region) + "</data>\n"
		# freesurfer name (d5)
		self.xmlstr += "      <data key =\"d5\">" + str(self.fsname) + "</data>\n"
		# node name (d6)
		self.xmlstr += "      <data key =\"d6\">" + str(self.nname) + "</data>\n"
		# left or right hemisphere (d7)
		self.xmlstr += "      <data key =\"d7\">" + str(self.hemisphere) + "</data>\n"
		#close the tag
		self.xmlstr += "    </node>\n"

	def getxml(self):
		self.genxml()
		return self.xmlstr
#----------------------------------------------------------

class Edge:

	def __init__(self):
		# node connection information
		self.source = -1
		self.target = -1

		# number of fibers (d9)
		self.nfibers = -1.0

		# Fiber length (d12)
		self.lfibers = -1.0

		# the node's xml string
		self.xmlstr = ""

		# -- These quantities are not currently available
		# -- in the Human c:onnectome data processesd via
		# -- FSL.  They are, however, in the braingraph.org
		# -- file format.  We fill in placeholder values for
		# -- posterity and consistency with existing software.

		# Fractional anisotropy mean value(d10)
		self.famean = -1.0

		# Fiber length std. deviation (d11)
		self.flstdev = -1.0

		# Fractional isotropy standard deviation (d13)
		self.fastdev = -1.0

	def setsourcetarget(self,source,target):
		s = int(source)
		t = int(target)

		if s == t or s<0 or t<0:
			print("Error in call to setsourcetarget(..): Invalid values for source and target")
			print("1. Self loops are not permitted (source = target)")
			print("2. source and target must be non-negative integers")
		else:
			self.source = s
			self.target = t

	def setfibernumbermean(self,fibn):
		fn = float(fibn)

		if fn <= 0:
			print("Error in call to setfibernumbermean(..)")
			print("1. fiber number cannot be negative")
			print("2. an edge with no fibers is not allowed (number of fibers = 0)")
		else:
			self.nfibers = fn

	def setfiberlengthmean(self,fibl):
		fl = float(fibl)

		if fl <= 0:
			print("Error in call to setfiberlengthmean(..)")
			print("1. fiber length cannot be negative")
			print("2. an edge with zer fiber length is not allowed (leangth of fibers = 0)")
		else:
			self.lfibers = fl

	def genxml(self, srctargetoffset):
		# node id: the srctargetoffset provides a shift for different numbering
		# schemes (i.e. zero based or one based, etc)
		self.xmlstr = "    <edge source=\"" + str(self.source + int(srctargetoffset)) + "\" target=\"" + str(self.target + int(srctargetoffset)) + "\">\n"
		# number of fibers (we use)
		self.xmlstr += "      <data key =\"d9\">" + str(self.nfibers) + "</data>\n"
		# Fractional anisotropy mean (for consistency, we do not use)
		self.xmlstr += "      <data key =\"d10\">" + str(self.famean) + "</data>\n"
		# Fiber length standard deviation (for consistency, we do not use)
		self.xmlstr += "      <data key =\"d11\">" + str(self.flstdev) + "</data>\n"
		# Fiber length mean (we use)
		self.xmlstr += "      <data key =\"d12\">" + str(self.lfibers) + "</data>\n"
		# Fractional isotropy standard deviation (for consistency, we do not use)
		self.xmlstr += "      <data key =\"d13\">" + str(self.fastdev) + "</data>\n"
		# close the tag
		self.xmlstr += "    </edge>\n"

	def getxml(self,srctargetoffset):
		self.genxml(srctargetoffset)
		return self.xmlstr

	# we can identify an edge by its source-->target
	# designation.  An important note here is that
	# 1. (a,b) and (b,a) are not the same
	# 2. We cannot index a dictionary by sets such as {a,b}
	def getidentifier(self):
		return (self.source, self.target)
#----------------------------------------------------------

class NodeList:
	def __init__(self):
		self.nodes = {}

	def addnode(self, number, node):
		self.nodes[number] = node

	def getnode(self,number):
		return self.nodes[number]

	def removenode(self, number):
		del self.nodes[number]

	def getnodenumbers(self):
		return self.nodes.keys()

	def getlength(self):
		return len(self.nodes)

class EdgeList:

	def __init__(self):
		self.edges = {}

	def addedge(self,edg):
		eid = edg.getidentifier()
		self.edges[eid] = edg

	# get an edge using the native identifier
	# for your edge type.
	def getedge(self, eid):
		return self.edges[eid]

	def removeedge(self,eid):
		del self.edges[eid]

	def getedgeidentifiers(self):
		return self.edges.keys()

	def removealledges(self):
		#kys = self.getedgeidentifiers()
		#for k in kys:
		#	self.removeedge(k)
		self.edges = {}

# ----------------------------------------------------------
# Implements a node container read from a TSV file
# ----------------------------------------------------------
class ParcellationNodes(NodeList):

	# pass in the fully qualified path to the tsv file to build the 
	#	initial nodelist from
	def __init__(self):
		super().__init__()
		self.tsv = False
		self.csv = False
		self.xmlgen = False
		self.minnodenumber = 10000000000000

	def csvcheck(self,csvfile):
		retv = False

		if path.isfile(csvfile) == True:
			retv = True
		else:
			print("File [" + csvfile + "] does not exist")
			print("You need to call readcoords(..) and pass in the parcellation csv coordinate filename")

		return retv

	# === === === ===
	def tsvcheck(self, tsvfile):
		retv = False

		if path.isfile(tsvfile) == True:
			retv = True
		else:
			print("File [" + tsvfile + "] does not exist")
			print("You will need to call readtsv(..) directly to initialize this object")
		return retv

	# === === === ===
	def readtsv(self,tsvfile):

		if self.tsvcheck(tsvfile) == True:
			f = open(tsvfile, "r")

			# skip the header
			f.readline()

			# begin reading
			line = f.readline()
			delimiter = ','
			while line:
				prvdelim = 0
				nxtdelim = 0
				fields = []
				anode = Node()


				# ----------------------------------
				# -- There are 4 fields per line --
				#         Get the first three
				for i in range(3):
					nxtdelim = line.find(delimiter, prvdelim)
					fields.append(line[prvdelim:nxtdelim])
					prvdelim = nxtdelim+1
				# now get the last one
				fields.append(line[prvdelim:])
				# ----------------------------------

				#-----------------------------------
				#   Process the individual fields
				for i in range(len(fields)):
					#strip off leading and trailing
					#  whitespace from the string
					fields[i] = fields[i].strip()

				# the first field is the node number
				anode.setnodenumber(fields[0])
				if self.minnodenumber > int(fields[0]):
					self.minnodenumber = int(fields[0])

				# the second field is used to set various labels for the node
				anode.setlabels(fields[1])
				# third field is cortical/subcortical (but we have already set this)
				# fourth field is the volume (in mm^3)
				anode.setvolume(fields[3])
				#-----------------------------------

				# add this node to the internal dictionary
				# (the key is the node number)
				self.nodes[int(fields[0])] = anode
				line = f.readline()
			f.close()
			self.tsv = True

	# === === === ===
	def readcoords(self,csvfile):
		# this function assumes that you are providing the coordinate csv file
		# for the parcellation.  This file can be created using the get_coordinates.py
		# source file on the appropriate raw parcellation output files

		if not self.tsv:
			print("Load nodal data from the .tsv file first by calling")
			print("the function readtst(\"/path/to/parcellation/tsv/file\"")
			return

		if self.csvcheck(tsvfile):
			f = open(csvfile, "r")
			# the first line of the CSV file corresponds to the node
			# with ID 1, the second line corresponds to the node with
			# ID 2, etc.
			inode = 1

			# there are no headers in this file and it is delimited by
			# spaces.
			delimiter = ' '
			line = f.readline()
			while line:
				prvdelim = 0
				nxtdelim = 0
				fields = []

				for i in range(2):
					nxtdelim = line.find(delimiter, prvdelim)
					fields.append(line[prvdelim:nxtdelim].strip())
					prvdelim = nxtdelim + 1
				# now get the last coordinate
				fields.append(line[prvdelim:].strip())

				(self.getnode(inode)).setcoords(float(fields[0]), float(fields[1]), float(fields[2]))

				inode += 1
				line = f.readline()

			f.close()
			self.csv = True

	def generatexmltags(self):
		if not (self.tsv and self.csv):
			print("You must first call readtsv(..) and then readcoords(..) before")
			print("calling generateXML")
		else:
			for k in self.nodes.keys():
				(self.nodes[k]).genxml()
			self.xmlgen = True

	def getlistxmlblock(self):
		xmlstr = ""

		if not (self.tsv and self.csv):
			print("You must first call readtsv(..) and then readcoords(..) before")
			print("you can call getlistxmlblock()")
		elif not self.xmlgen:
			self.generatexmltags()

		# we write the nodes out in order
		keyslist = list(self.nodes.keys())
		keyslist.sort()

		for k in keyslist:
			xmlstr += (self.nodes[k]).getxml()

		return xmlstr

	def getminnodeid(self):
		return self.minnodenumber
#----------------------------------------------------------
class ParcellationEdges(EdgeList):

	# this constructor expects a ParcellationNodes
	# object to be passed to it.  The edges are
	# assumed to be undirected
	def __init__(self, parcnodes):
		super().__init__()

		# get the number of nodes in the parcellation
		self.matsz = parcnodes.getlength()

		# the raw edges
		self.matnumbers = np.zeros((self.matsz, self.matsz))
		self.matlengths = np.zeros((self.matsz, self.matsz))

		# the flags for edge setting
		self.fdtnetworkmatrix  = False
		self.fdtnetworklengths = False

		# have the raw matrices benn symmetrized yet
		self.sym = False

	def reset(self):
		self.removealledges()
		self.matnumbers = np.zeros((self.matsz, self.matsz))
		self.matlengths = np.zeros((self.matsz, self.matsz))
		self.fdtnetworkmatrix = False
		self.fdtnetworklengths = False

	# an internal helper function used to check
	# the existence of a given file
	def checkexist(self, pathto):
		retv = False

		if path.isfile(csvfile) == True:
			retv = True
		else:
			print("File [" + pathto + "] does not exist")
		return retv

	# this function parses a line in the raw FSL tractography
	# patient data files.  E.g. either fdt_network_matrix or
	# fdt_network_matrix_lengths.  These files contain nROI
	# lines with each line having nROI entries; delimited
	# by double spaces (the default value for delim)
	def parsematline(self, linestr, nROI, delim='  '):
		prvdelim = 0
		nxtdelim = 0
		fields = []

		for i in range(nROI-1):
			nxtdelim = linestr.find(delim, prvdelim)
			fields.append(float(linestr[prvdelim:nxtdelim].strip()))
			prvdelim = nxtdelim + len(delim)
		# now get the last entry
		fields.append(float(linestr[prvdelim:].strip()))

		return fields


	# read the fiber numbers from a file. Pass in the fully qualified
	# path to the file (e.g. the patient fdt_network_matrix
	# file for a particular subject)
	def readfibernumbers(self,pathto,masked=False,maskObj=None):

		# construct a mask (which is all ones by default)
		mask = np.ones((self.matsz,self.matsz))

		if masked == True and type(maskObj) != None:
			mask = (maskObj.getMask().tolist())


		if not self.fdtnetworkmatrix and self.checkexist(pathto):
			# fill in the fiber number matrix (connectivity) from the raw file

			f = open(pathto, "r")
			line = f.readline()
			rowndx = 0

			while line:
				# process each line of the fiber number file
				# Each line of the fiber number file contains the number of connections
				# between the ROIs.  There are a number of entries, in each line, equal
				# to the number of ROIs; there are the same number of lines in each file
				# (square matrix of size nROI x nROI).  The file is double-space delimited
				entries = self.parsematline(line, self.matsz)

				for colndx in range(self.matsz):
					# update the internal matrix
					self.matnumbers[rowndx][colndx] = mask[rowndx][colndx]*entries[colndx]

				line = f.readline()
				rowndx += 1

			f.close()
			self.fdtnetworkmatrix = True
			self.sym = False # system matrices are no longer symmetric

			# uncomment to show a visual of the current numbers matrix
			# plt.imshow(self.matnumbers)
		else:
			print("Cannot load matrix numbers")
			print("Provide a valid file or call reset(..) before proceeding")

	# read the lengths from a file. Pass in the fully qualified
	# path to the file (e.g. the patient fdt_network_matrix_lengths
	# file for a particular subject)
	def readfiberlenghts(self,pathto,masked=False,maskObj=None):
		mask = np.ones((self.matsz,self.matsz))

		if masked == True and type(maskObj) != None:
			mask = (maskObj.getMask().tolist())

		if not self.fdtnetworklengths and self.checkexist(pathto):
			# fill in the fiber number lengths matrix from the raw file

			f = open(pathto, "r")
			line = f.readline()
			rowndx = 0

			while line:
				# process each line of the fiber number file
				# Each line of the fiber number file contains the number of connections
				# between the ROIs.  There are a number of entries, in each line, equal
				# to the number of ROIs; there are the same number of lines in each file
				# (square matrix of size nROI x nROI).  The file is double-space delimited
				entries = self.parsematline(line, self.matsz)

				for colndx in range(self.matsz):
					# update the internal matrix
					self.matlengths[rowndx][colndx] = mask[rowndx][colndx]*entries[colndx]

				line = f.readline()
				rowndx += 1

			f.close()
			self.fdtnetworklengths = True
			self.sym = False  # system matrices are no longer symmetric

			# uncomment to see a visual of the current length matrix
			# plt.imshow(self.matlengths)
		else:
			print("Cannot load matrix lengths")
			print("Provide a valid file or call reset(..) before proceeding")

	def symmetrize(self):
		if self.fdtnetworkmatrix and self.fdtnetworklengths and not self.sym:
			self.matnumbers = 0.5*(self.matnumbers + np.matrix.transpose(self.matnumbers))
			self.matlengths = 0.5*(self.matlengths + np.matrix.transpose(self.matlengths))
			self.sym = True

			# uncomment to see a visual of the symmetrized numbers and length matrices
			# plt.imshow(self.matnumbers)
			# plt.imshow(self.matlengths)
		else:
			print("Cannot symmetrize.  Please ensure that the network matrix")
			print("and matrix lengths have been read and that this object has")
			print("been constructed with sym=True (default)")

	# thesholds the internal matrices based on an input percentage
	# 0 < perc <= 1.  Entries where the number of connections are
	# at least (perc)% of the maximal connectivity matrix value will
	# be kept while everything else will be set to zero.
	# Example: calling threshold(0.1) will keep every entry which is at least 10%
	# of the maximum connectivity.  So if the maximum connectivity entry is 1050
	# and perc = 0.1 then every entry e with (connectivity score) e/1050 < 0.1
	# will be set to zero while entries with (connectivity score) e/1050 > 0.1
	# will be retained.
	# c.f. https://stackoverflow.com/questions/36719997/threshold-in-2d-numpy-array/36720130
	def thresholdbyconnectivity(self,perc):
		if self.fdtnetworkmatrix and self.fdtnetworklengths:
			if not self.sym:
				print("Matrices have not yet been symmetrized.  Calling symmetrize(..) and")
				print("then applying the requested threshold")
				self.symmetrize()

			maxn = self.matnumbers.max()

			# create the matrix masks.  These matrices are True
			# where the entry / maximum exceeds (perc) and False
			# where entry / maximum is below the (perc) threshold
			maskconn = (self.matnumbers / maxn >= perc)

			# threshold the internal matrices by the connectivity mask
			self.matnumbers = maskconn*self.matnumbers
			self.matlengths = maskconn*self.matlengths

			#uncomment to see a visual of the symmetrized and
			# thresholded number and length matrices
			# plt.imshow(self.matnumbers)
			# plt.imshow(self.matlengths)
		else:
			if not self.fdtnetworkmatrix:
				print("Please call readfibernumbers(..) first")
			if not self.fdtnetworklengths:
				print("Please call readfiberlengths(..) first")

	# this function generates an edge list from the internal matrices
	# if the matrices have not been symmetrized this function will symmetrize
	# them first. Note: this function only creates edge objects for edges with
	# non-zero connectivity.  Thus, any thresholding must be done before calling
	# the generateedges(..) function.
	def generateedges(self):

		if not self.fdtnetworkmatrix:
			print("Please call readfibernumbers(..) first")
			return

		if not self.fdtnetworklengths:
			print("Please call readfiberlengths(..) first")
			return

		if not self.sym:
			self.symmetrize()

		# clear out any existing edges before building a new list
		if len(self.edges) > 0:
			self.removealledges()

		# since this is an undirected graph, we need only walk
		# the upper triangular portion of the connectivity matrix
		# to add the edge objects.
		for i in range(self.matsz):
			for j in range(i, self.matsz):
				# we do not add self edges
				if i == j:
					continue
				conentry = self.matnumbers[i][j]
				lenentry = self.matlengths[i][j]
				if conentry > 0:
					if lenentry > 0:
						e = Edge()
						e.setsourcetarget(i, j)
						e.setfibernumbermean(conentry)
						e.setfiberlengthmean(lenentry)
						self.addedge(e)
					else:
						print("**MATRIX ERROR**: connectivity at ["+ str(i) + "][" + str(j) + "]")
						print("is non-zero but the corresponding length entry is zero")

	def getlistxmlblock(self,srctargetoffset=0):
		xmlstr = ""
		if len(self.edges) == 0:
			print("Please call generateedges(..) first")
			return xmlstr
		for ky in self.edges.keys():
			e = self.getedge(ky)
			xmlstr += e.getxml(srctargetoffset)

		return xmlstr
#----------------------------------------------------------
#----------------------------------------------------------


#----------------------------------------------------------
# This class creates a frequency mask that can be applied
# to a ParcellationEdges class.  This class assumes that
# you have a directory containing (at least one) subdirectories
# which correspond to individual subjects.  It is assumed that
# each directory contains a CSV file satisfying:
#    1. containing the patient connectivity matrix
#    2. each such file has the same name
# ---------------
# The `makeMask' function creates a True/False matrix that reflects
# whether or not an edge occurs in at least `0<percentage<=1' of the
# patient subject directories in `fullpath'
class groupMask():

	# Construct with
	#	1. The full path to the directory containing the subject subdirectories (path ending with /)
	#	2. The name of the static file `fname' that contains the connectivity matrix in each subject directory
	#	3. The number of vertices represented by the connectivity matrix
	def __init__(self,fullpath,fname,nV):
		self.reset(fullpath,fname,nV)

	# private function that provides an iterator for the different mask building options
	def __getpatientconnectivities(self,delimiter):
		for rootdir, subjectdirs, files in walk(self.path):
			totalsubj = len(subjectdirs)
			thissubj = 0
			for subj in subjectdirs:
				thissubj +=1
				subjectmatf = self.path + subj + "/" + self.fname
				csv = np.genfromtxt(subjectmatf, delimiter=delimiter)
				printProgressBar(thissubj, totalsubj, prefix='Creating group frequency mask:', suffix='Complete', length=50)


				# pass back to the calling generator
				# for further processing
				yield csv

	def makeFrequencyMask(self,freqPercentage,delimiter,medianFilter=False,averageFilter=False,filterPerc=1.0):
			if medianFilter == averageFilter:
				if medianFilter:
					print("Only one of medianFilter or averageFilter should be True")
					return
				else:
					print("Generating unfiltered frequency mask with occurence percentage {}%".format(freqPercentage*100))

			numpatients = 0
			for csvmat in self.__getpatientconnectivities(delimiter):

				# increment the number of patients
				numpatients += 1

				#------------------------------
				# walk the raw csv matrix and
				# increment the mask accordingly

				npcsv = np.array(csvmat)
				filter = 0

				if medianFilter:
					filter = np.median(npcsv)
				if averageFilter:
					filter = np.average(npcsv)

				for ndx, value in np.ndenumerate(npcsv):
					if value > filter*filterPerc:
						# record that the edge has passed
						# in the aggregation matrix
						self.tmat[ndx] += 1

			# now we walk the internal mask matrix and determine whether or not
			# to reset the entries to one or zero.  First, determine the number
			# of times that an entry must be aggregated in order to stay
			npat = int(freqPercentage*numpatients)

			# now determine the mask
			for ndx, value in np.ndenumerate(self.tmat):
				if value > npat:
					self.tmat[ndx] = 1.0
				else:
					self.tmat[ndx] = 0.0

			plt.imshow(self.tmat.tolist())

	def getMask(self):
		return self.tmat

	def reset(self,fullpath,fname,nV):
		self.path 	= fullpath
		self.fname 	= fname
		self.nV 	= nV
		self.tmat 	= np.zeros((self.nV, self.nV))


# Outputs a graphml file.  Requires a fully-constructed
# ParcellationNodes object (parcnodes) and ParcellationEdges
# object (parcedges) in addition to a fully-qualified path
# for output
class GraphMLFile:

	def __init__(self,parcnodes,parcedges,outpath):
		self.nodeList = parcnodes
		self.edgeList = parcedges
		self.outputfile = outpath

	def outputxml(self):
		fl = open(self.outputfile,"w")
		nodeidoffset = self.nodeList.getminnodeid()


		# print the header
		headerstr = "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" "
		headerstr += "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
		headerstr += "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns "
		headerstr += "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"

		keyblock  = "  <key attr.name=\"dn_position_x\" attr.type=\"double\" for=\"node\" id=\"d0\" />\n"
		keyblock += "  <key attr.name=\"dn_position_y\" attr.type=\"double\" for=\"node\" id=\"d1\" />\n"
		keyblock += "  <key attr.name=\"dn_position_z\" attr.type=\"double\" for=\"node\" id=\"d2\" />\n"
		keyblock += "  <key attr.name=\"dn_correspondence_id\" attr.type=\"string\" for=\"node\" id=\"d3\" />\n"
		keyblock += "  <key attr.name=\"dn_region\" attr.type=\"string\" for=\"node\" id=\"d4\" />\n"
		keyblock += "  <key attr.name=\"dn_fsname\" attr.type=\"string\" for=\"node\" id=\"d5\" />\n"
		keyblock += "  <key attr.name=\"dn_name\" attr.type=\"string\" for=\"node\" id=\"d6\" />\n"
		keyblock += "  <key attr.name=\"dn_hemisphere\" attr.type=\"string\" for=\"node\" id=\"d7\" />\n"
		keyblock += "  <key attr.name=\"number_of_fibers\" attr.type=\"double\" for=\"edge\" id=\"d9\" />\n"
		keyblock += "  <key attr.name=\"FA_mean\" attr.type=\"double\" for=\"edge\" id=\"d10\" />\n"
		keyblock += "  <key attr.name=\"fiber_length_std\" attr.type=\"double\" for=\"edge\" id=\"d11\" />\n"
		keyblock += "  <key attr.name=\"fiber_length_mean\" attr.type=\"double\" for=\"edge\" id=\"d12\" />\n"
		keyblock += "  <key attr.name=\"FA_std\" attr.type=\"double\" for=\"edge\" id=\"d13\" />\n"

		# get the graph tag
		keyblock += "  <graph edgedefault=\"undirected\">\n"

		nodestr = self.nodeList.getlistxmlblock()
		edgestr = self.edgeList.getlistxmlblock(nodeidoffset)

		footerstr = "  </graph>\n</graphml>"

		fl.write(headerstr)
		fl.write(keyblock)
		fl.write(nodestr)
		fl.write(edgestr)
		fl.write(footerstr)
		fl.close()
#----------------------------------------------------------


if __name__ == "__main__":

	# Set singlesubject to true to create connectome files for a single subject.
	# Set singlesubject to false to create connectome files for multiple subjects whose
	# files are all contained in directories located under a common root directory
	# Example: [path]/subjects/subj-1/  [path]/subjects/subj-2/ and so on
	singlesubject = False

	# the scale of the connectome(s) to be processed.  This should be an integer from 1 to 5
	scale = 1

	# You should set this value to the number of vertices (Regions of interests ROIs) for the
	# specific connectome parcellation (and scale) that you are working with.
	# The ROIs are as follows:
	# ---------
	# Scale 1:
	#....Standard - 83
	#....Hippocampal subfields - 121
	#....All subfields - 126
	# --------
	# Scale 2:
	#....Standard - 129
	#....Hippocampal subfields - 157
	#....All subfields - 172
	# ---------
	# Scale 3:
	#....Standard - 231
	#....Hippocampal subfields - 259
	#....All subfields - 274
	# ---------
	# Scale 4:
	#....Standard - 461
	#....Hippocampal subfields - 489
	#....All subfields - 504
	# ---------
	# Scale 5:
	#....Standard - 1017
	#....Hippocampal subfields - 1045
	#....All subfields - 1060
	# ---------
	nROI = 83

	# ==================== Parcellation Information ================================
	# Set the directories that determine the nodes.  These files are determined by the parcellation
	# strategy and are common for every tractography dataset constructed using this parcellation
	# -----
	# tsvfile: the .tsv file for the parcellation.  Example: [parcellationroot]/sub-01_label-L2018_desc-scale2_stats.tsv
	# -----
	# csvfile: the .csv file containing the coordinates for the parcellation.  This file can be created by running
	#	get_coordinates.py on the raw parcellation files (see get_coordinates.py for details)
	#	Example: [parcellationroot]/sub-01_label-L2018_desc-scale2_atlas_coordinates.csv

	#----
	# Standard connectome
	parcellationroot = "/home/zcxp/Documents/repos/oxford/Connectomes/standard_connectome/parcellation/parcellation-files/"

	# standard + hippocampal subfields connectome
	#parcellationroot = "/home/zcxp/Documents/repos/oxford/Connectomes/hippsubfields_connectome/parcellation/parcellation-files/"

	# standard + hippocampal + thalamic + brainstem subfields connectome
	#parcellationroot = "/home/zcxp/Documents/repos/oxford/Connectomes/full_connectome/parcellation/parcellation-files/"
	#----

	tsvfile = parcellationroot + f"sub-01_label-L2018_desc-scale{scale}_stats.tsv"
	csvfile = parcellationroot + f"sub-01_label-L2018_desc-scale{scale}_atlas_coordinates.csv"

	# =========================================================================

	# ==================== Subject Information ================================
	# Specify the top-level directory that contains your subject file(s)

	#----
	# Standard connectome
	subjectroot = f"/home/zcxp/Documents/repos/oxford/Connectomes/standard_connectome/scale{scale}/subjects/"

	# Standard connectome + hippocampal subfields connectome
	# subjectroot = f"/home/zcxp/Documents/repos/oxford/Connectomes/hippsubfields_connectome/scale{scale}/subjects/"

	# Standard + hippocampal + thalamic + brainstem subfields connectome
	#subjectroot = f"/home/zcxp/Documents/repos/oxford/Connectomes/full_connectome/scale{scale}/subjects/"
	#----


	# if singlesubject == True: set this string to the name of the directory
	# that contains the subject files you want to parse.  Example: subj-1/
	# if singlesubject == False, above, you can ignore the following as
	# it is determined automatically by the processing algorithm
	singlesubjectstr = "100206/"

	# Set the names of the connectivity files.  The standard FSL name for the
	# file that contains the number of connections between two ROIs in
	# a parcellation is `fdt_network_matrix' while the standard FSL name for
	# the file that contains the average fiber length is `fdt_network_matrix_lengths'
	fdtconnectivity = "fdt_network_matrix"
	fdtlengths = "fdt_network_matrix_lengths"
	# =========================================================================

	# ---- Debugging ----
	# Uncomment to do a test read on a single subject matrix to make sure dimensions agree
	#testm = np.genfromtxt(subjectroot + singlesubjectstr + fdtconnectivity,delimiter='  ')

	# ==================== Output-related Information ================================
	# Set the threshold you want to use.  This value should satisfy 0 < threshold <= 1.
	# The thresholding algorithm  does the following
	# 1. Determines the strongest connection (Cmax) in the connectivity matrix.
	# 2. For every edge with connection strength C, calculate C/Cmax
	# 3. If C/Cmax < threshold: set the edge connectivity strength to zero
	# This means that a lower values of (threshold) will set fewer edges to zero, while
	# a higher value of (threshold) will set more edges to zero

	# optimal (MND = 40) for standard connectome
	# threshold = 0.0079

	# optimal (MND = 40) for the hippocampal subfield connectome
	# threshold = 0.0062

	# optimal (MND=40) for the hippocampal + thalamic + brainstem subfield connectome
	threshold = 0.0060

	# Set the root directory for your desired output(s)
	#----
	# Standard connectome
	outputroot = f"/home/zcxp/Documents/repos/oxford/Connectomes/standard_connectome/scale{scale}/processed/"

	# Standard connectome + hippocampal subfields
	# outputroot = "/home/zcxp/Documents/repos/oxford/Connectomes/hippsubfields_connectome/scale1/processed/"

	# Standard + hippocampal + thalamic + brainstem subfields connectome
	#outputroot = "/home/zcxp/Documents/repos/oxford/Connectomes/full_connectome/scale1/processed/"
	#----


	# make a frequency mask based on a 50% frequency rate using average
	# filtering at the patient level (e.g. patient edges are retained,
	# for frequency matching, if their edge strength is 25% of average
	# or more).  We can use this to threshold the patient matrices
	# when reading in the individual files
	mask = groupMask(subjectroot, fdtconnectivity, nROI)
	#mask.makeFrequencyMask(0.5,'  ',averageFilter=True, filterPerc=0.25)
	mask.makeFrequencyMask(0.99, '  ', medianFilter=True, filterPerc=1.0)

	# specify a prefix for your output file.  The script will generate the appropriate
	# suffix and extensions (e.g. scale1.graphml or scale1-connectivity.csv etc)
	outputprefix = f"scale{scale}-standard-average-50pfreq"
	# =================================================================================

	if len(sys.argv) != 1:
		print("")
		print("***********************************************************************************************")
		print("This script parses the various parcellation output files produced by CMTK and FSL to construct")
		print("a connectome graph in the format of braingraph.org (used by BrainNet, PrYon, etc)")
		print("")
		print("Please modify the directory paths in the main function to parse your connectome files")
		print("***********************************************************************************************")
	else:
		#----------------------------------------------
		# Create connectome files for a single subject
		#----------------------------------------------
		if singlesubject:

			# Populate the parcellation nodes class.  We will use this
			# list of nodes to construct every patient connectome set
			# of edges
			parcNodes = ParcellationNodes()
			parcNodes.readtsv(tsvfile) 		# get the nodes and nodal information
			parcNodes.readcoords(csvfile)	# get the nodal coordinates (c.f. get_coordinates.py to create this file)
			parcNodes.generatexmltags()		# generate the XML tags for the individual nodes based on the file data

			# Create and populate the parcellation edges for the subject above
			parcEdges = ParcellationEdges(parcNodes)

			# We can apply a mask to the individual matrices, when reading them in, if desired.
			# For example, by using a groupMask() object.  The mask object must implement a function
			# called getMask() which returns a numpy 2D array (of 1s and 0s) to apply. A 1 means keep
			# and a 0 means discard
			parcEdges.readfibernumbers(subjectroot + singlesubjectstr + fdtconnectivity, masked=True, maskObj=mask)
			parcEdges.readfiberlenghts(subjectroot + singlesubjectstr + fdtlengths, masked=True, maskObj=mask)
			parcEdges.symmetrize()

			# Uncomment, below, to further threshold the matrices by percentage of maximal entry
			# i.e. all edges such that edge/maxentry < threshold will be discarded
			#parcEdges.thresholdbyconnectivity(threshold)
			parcEdges.generateedges()

			# the name(s) of the file(s) that will be generated for each subject
			graphmlfilename = outputprefix + ".graphml"

			# create the graphml file
			graphmlfile = GraphMLFile(parcNodes, parcEdges, outputroot + graphmlfilename)
			graphmlfile.outputxml()
		# ----------------------------------------------
		# Create connectome files for multiple subjects
		# ----------------------------------------------
		else:
			# if singlesubject is False, it is assumed that we are parsing multiple subjects.
			#
			# First, get the node list.  We only need to do this once
			parcNodes = ParcellationNodes()
			parcNodes.readtsv(tsvfile)     # get the nodes and nodal information
			parcNodes.readcoords(csvfile)  # get the nodal coordinates (c.f. get_coordinates.py to create this file)
			parcNodes.generatexmltags()    # generate the XML tags for the individual nodes based on the file data

			# our parcellation edges object. We will reset this object
			# whenever we switch subjects; we therefore only need one
			# instance of it.
			parcEdges = ParcellationEdges(parcNodes)

			# the name(s) of the file(s) that will be generated for each subject
			graphmlfilename = ""

			# get the directory list in the subject root
			for rootdir, subjectdirs, files in walk(subjectroot):
				totalsubj = len(subjectdirs)
				thissubj = 0
				for subj in subjectdirs:
					thissubj += 1

					printProgressBar(thissubj, totalsubj, prefix='Creating subject graphml connectomes:', suffix='Complete',length=50)
					subjectstr = subj + "/"
					graphmlfilename = outputprefix + "-" + subj + ".graphml"

					# the file generation block
					parcEdges.reset()

					# We can apply a mask to the individual matrices, when reading them in, if desired.
					# For example, by using a groupMask() object.  The mask object must implement a function
					# called getMask() which returns a numpy 2D array (of 1s and 0s) to apply. A 1 means keep
					# and a 0 means discard
					parcEdges.readfibernumbers(subjectroot + subjectstr + fdtconnectivity, masked=True, maskObj=mask)
					parcEdges.readfiberlenghts(subjectroot + subjectstr + fdtlengths, masked=True, maskObj=mask)
					parcEdges.symmetrize()

					# Uncomment, below, to further threshold the matrices by percentage of maximal entry
					# i.e. all edges such that edge/maxentry < threshold will be discarded
					#parcEdges.thresholdbyconnectivity(threshold)

					parcEdges.generateedges()

					# output the files for this subject
					graphmlfile = GraphMLFile(parcNodes, parcEdges, outputroot + graphmlfilename)
					graphmlfile.outputxml()