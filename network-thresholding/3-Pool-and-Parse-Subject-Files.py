import os
import shutil
from thresholding import parsing

# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#
#   Pool and parse subject files to graphml
#   -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
#   This file assumes that 2-Cosica-Thresholded-to-FSL.py
#   has been run to generate thresholded graphs for a set
#   of patient data and that these graphs have been
#   converted to the FSL format.
#
#   This script generates a comprehensive set of patient
#   graphml images - for each threshold method and for
#   each selected threshold value
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


#-----------------------------------------------------------
#                   Directory structure
#-----------------------------------------------------------

# -..-..-..-..-..-
# [Modify] > The relative directory (to this script) where you have stored
#          > the FSL connectome connectivity matrix data to be processed by this script
parcellationScale = 1
connectome = "standard_connectome"

#	[Do not modify]
#    > The directory of the current script
scriptdir = os.path.dirname(__file__)

# -..-..-..-..-..-
# [Modify]  > The relative directory (to this script) where you have stored
#           > the FSL connectome connectivity matrix data to be processed by this script
#           > ** This should be the same as `relative_finished_results_root' from the
#           > script 2-Cosica-Threholded-to-FSL
relative_thresholded_fslsubjectroot = "fslformatted/"


# -..-..-..-..-..-
# [No need to modify if running from within the repository structure]
#   > path to the raw fsl parcellation files for this connectome
relative_raw_parcellationroot = f"../{connectome}/parcellation/parcellation-files/"
#   > path to the raw FSL subject directory (un-thresholded)
relative_raw_fslsubjectroot = f"../{connectome}/scale{parcellationScale}/subjects/"

# -..-..-..-.-..-
# [Modify]  > The final (relative path) to the output directory for the finished (graphml) patient data file set
relative_final_subject_output_path = "completed/"

# -..-..-..-..-..-
#	[Do not modify]
#    > These are the file names of the connectivity and length matrices output by FSL.  This should
#    > not need to be modified (these are the names hard-coded in 1-FSL-to-Cosica.py)
fdtconnectivity = "fdt_network_matrix"
fdtlengths = "fdt_network_matrix_lengths"

parcellation_tsvfile = f"sub-01_label-L2018_desc-scale{parcellationScale}_stats.tsv"
parcellation_csvfile = f"sub-01_label-L2018_desc-scale{parcellationScale}_atlas_coordinates.csv"
# ---------------------------------------------------------

# --------------------------------------------------------
# Any paths input by the user will be formatted according to
# script expectations by encapsulating with this function.
def fixpath(pathstr):
    if pathstr[-1] != '/':
        pathstr = pathstr + '/'
    return pathstr

#--------------
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
#--------------



#--------
# Description:
# The directory structure generated from 2-Cosica-Thresholded-to-FSL.py
# is uniform.  Meaning that each subject directory has the same structure,
# and set of files, for every subject.  Therefore, we need only read the
# structure of a single directory in order to garner the full picture
#
# Input:
# The relative directory of the reformatted FSL files produced by
# 2-Cosica-Thresholded-to-FSL.
#
# Return values:
# Return a list of subjects, network sparsification methods, and
# the associated dictionary of thresholds found in the directory
# `reformattedFSLDir'
#--------
def getDirectoryStructure(reformattedFSLDir):
    fsldir = fixpath(reformattedFSLDir)

    #  > First, get the list of subject directories
    allsubjects = next(os.walk(fsldir))[1]

    #  > Now, get the list of thresholding methods used for each
    #    subject.  We can do this simply by looking at the first
    #    subject's results.
    firstsubj = fsldir + allsubjects[0]
    allmethods = next(os.walk(firstsubj))[1]

    #   > Now build a dictionary of the relevant thresholds for
    #   each method
    thresholds = {}
    for method in allmethods:
        methoddir = fixpath(firstsubj) + method
        thresholds[method] = next(os.walk(methoddir))[1]

    return allsubjects, allmethods, thresholds
# ---------------------------------------------------------



#----------------------------------------------------------
# Begin generation of thresholded patient graphml files
#----------------------------------------------------------
# the prefix for the graphml files we will generate
graphmlprefix = f"{connectome}-scale{parcellationScale}"

# Gather directory structure information
subjects, methods, thresholds = getDirectoryStructure(relative_thresholded_fslsubjectroot)

# Create and populate a single parcellation object
# this object interprets the connectivity matrices
# in terms of the parcellation used to generate
# them.
parcpath = fixpath(relative_raw_parcellationroot)
parcNodes = parsing.ParcellationNodes()
parcNodes.readtsv(parcpath + parcellation_tsvfile)
parcNodes.readcoords(parcpath + parcellation_csvfile)
parcNodes.generatexmltags()

# our parcellation edges object. We will reset this object
# whenever we switch subjects; we therefore only need one
# instance of it.
parcEdges = parsing.ParcellationEdges(parcNodes)


thissubj = 0
totalsubjs = len(subjects)
for subject in subjects:

    thissubj += 1
    print("Generating graphml files for subject: " + subject)
    printProgressBar(thissubj, totalsubjs, prefix='Writing graphml files:', suffix='Complete', length=100)

    # This is the path to the raw FSL length data for the subject
    fdt_subj_lengths = fixpath(relative_raw_fslsubjectroot) + subject + "/" + fdtlengths

    for method in methods:
        for threshold in thresholds[method]:
            # Read in the (reformmatted and thresholded)
            # parcellation edge object for this subject
            parcEdges.reset()

            # The path to our reformatted FSL matrix file generated
            # from for this specific method and at this specific threshold
            thispath = fixpath(relative_thresholded_fslsubjectroot) + subject + "/" + method + "/" + threshold
            thispath = fixpath(thispath)
            thisconnectivityfile = os.path.join(scriptdir, thispath + fdtconnectivity)

            # The fiber number connectivity files that we have generated
            # are already symmetric and have already been thresholded so
            # so we don't need to apply a further mask to these
            parcEdges.readfibernumbers(thisconnectivityfile, masked=False)

            # The fiber lengths matrices come from the raw data.  They are
            # not symmetric and have not been thresholded.  However, our
            # parcellation object only creates edges for entries with non-zero
            # connectivity.  Thus, we don't need to threshold the fiber lengths
            # explicitly.  However, we do need to symmetrize this matrix.
            parcEdges.readfiberlenghts(fdt_subj_lengths, masked=False)

            # Symmetrize the system matrices - this will result in the same
            # connectivity matrix (since it is already symmetric) but will
            # symmetrize the fiber lengths.
            parcEdges.symmetrize()

            # generate the edges (as noted above, we only create edges
            #   for the entries in the connectivity matrix, which has
            #   already been thresholded, that are non-zero)
            parcEdges.generateedges()

            #output the graphml file
            graphmlout = graphmlprefix + f"-{subject}" + ".graphml"

            # output the files for this subject
            graphmlfile = parsing.GraphMLFile(parcNodes, parcEdges, thispath + graphmlout)
            graphmlfile.outputxml()

#--------------------------------------------
# Begin gathering of results
#--------------------------------------------
relative_final_subject_output_path = fixpath(relative_final_subject_output_path)

if os.path.exists(relative_final_subject_output_path):
    shutil.rmtree(relative_final_subject_output_path)

os.mkdir(relative_final_subject_output_path)

#--
thisStatus = 0
totalStatus = len(methods)*len(thresholds)+1
for method in methods:
    methodpath = fixpath(relative_final_subject_output_path + method)
    os.mkdir(methodpath)

    for threshold in thresholds[method]:
        thisStatus += 1
        printProgressBar(thisStatus, totalStatus, prefix=f'Gathering final files in: {relative_final_subject_output_path} ', suffix='Complete', length=100)
        destinationDirectory = fixpath(methodpath + threshold)
        os.mkdir(destinationDirectory)

        # Now we are going to go through the subject list and aggregate
        # the individual graphml files for this path and threshold
        for subject in subjects:
            targetfile = fixpath(relative_thresholded_fslsubjectroot)
            targetfile = targetfile + subject + "/" + method + "/" + threshold + "/" + graphmlprefix + f"-{subject}" + ".graphml"
            shutil.copy(targetfile,destinationDirectory)


# Clear the intermediate directory
shutil.rmtree(relative_thresholded_fslsubjectroot)
