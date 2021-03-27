from os import path
import numpy as np

# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#
#   General python utilities for parsing FSL-formatted raw
#   output files and generating subsequent graphml
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


# ==========================================================
#              Basic node and Edge Objects
# ==========================================================

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

# -- -- -- -- -- --

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


# ==========================================================
#                   Node and Edge Lists
# ==========================================================

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

# -- -- -- -- -- --

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


# ==========================================================
#              Node and edge containers for
#                   use in generating
# ==========================================================

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

		if self.csvcheck(csvfile):
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

# -- -- -- -- -- --

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

		if path.isfile(pathto) == True:
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
