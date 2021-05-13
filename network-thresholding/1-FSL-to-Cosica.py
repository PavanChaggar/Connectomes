import os
import shutil
import pandas
import numpy as np
from scipy import stats


# --------------------------------------------------------
#
#  ***Oxford Mathematical Brain Modeling Group***
#
#  FSL raw connectivity matrices, symmetrization, and output
#  to Cosica format for further thresholding using 2-Cosica-Thresholded-to-FSL.py
#
#   Michel Cosica's thresholding library: https://www.michelecoscia.com/?page_id=287
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

# ---------------------------------------------------------
# This script reformats a standard FSL output connectivity matrix into the form required
# by Prof. Michele Coscia's network backboning toolbox (https://www.michelecoscia.com/?page_id=287)


# ----- Connectome and scale
parcellationScale = 5
connectome = "standard_connectome"

# ----- Directory structure -----
#    > The directory of the current script
scriptdir = os.path.dirname(__file__)
# [Modify] > The relative directory (to this script) where you have stored
#          > the FSL connectome connectivity matrix data to be processed by this script
relative_subjectroot = f"../{connectome}/scale{parcellationScale}/subjects/"

#          > Create the subject root directory
subjectroot = os.path.join(scriptdir, relative_subjectroot)

# [Modify] > The relative directory where you want outputs to be stored.  Appropriate subdirectories will
#            be created if they do not already exist
relative_output_root = "cosicaformatted/"

#    > These are the file names of the connectivity and length matrices output by FSL.  This should
#    > not need to be modified
fdtconnectivity = "fdt_network_matrix"
fdtlengths      = "fdt_network_matrix_lengths"
# ---------------------------------------------------------

# ----------------------------------------------------------------
#                         Helper Functions
#                         (Do Not Modify)
# ----------------------------------------------------------------

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

# -----
# Any paths input by the user will be formatted according to
# script expectations by encapsulating with this function.
def fixpath(pathstr):
    if pathstr[-1] != '/':
        pathstr = pathstr + '/'
    return pathstr
# -----

def mkoutdirs(relative_output_root):
    relative_output_root = fixpath(relative_output_root)

    if not(os.path.exists(relative_output_root)):
        os.mkdir(relative_output_root)



def readSubjectCSV(fslSubjMatPath):
    return np.genfromtxt(fslSubjMatPath, delimiter='  ')

def generateSymmetricCosicasTable(npsym):
    # build the dictionary
    tab = {'src':[], 'trg':[], 'nij':[]}
    for srcroi in range(len(npsym)):
        rowlen = len(npsym[srcroi])
        # generate an upper triangular table to save space
        # the matrix npsym is expected to be symmetric
        for trgroi in range(srcroi, rowlen):
            if srcroi == trgroi:
                # Do not permit self loops
                weight = 0.0
            else:
                weight = float(npsym[srcroi][trgroi])
            tab['src'].append(srcroi)
            tab['trg'].append(trgroi)
            tab['nij'].append(weight)

    return pandas.DataFrame(tab, columns=['src', 'trg', 'nij'])

def writeDataframeStats(df,dir):
    writeto = fixpath(dir)

    # compute the arithmetic mean, geometric mean, and ROI size of
    # the dataframe and store it in a statistics file so that we
    # can read it in later for further processing
    # the (expected) columns of the dataframe created are
    # 0: src ids
    # 1: trg ids
    # 2: weights

    # extract the rows - the maximum row number (plus one)
    # should be the number of ROIs in the graph
    rws = df.iloc[:,0]
    nROI = np.max(rws) + 1

    # extract the weights.  We want the arithmetic and
    # geometric means of the non-zero entries
    wgts = df.iloc[:, 2]
    nzwgts = wgts.loc[wgts != 0]

    amean = np.mean(nzwgts)
    gmean = stats.gmean(nzwgts)

    # create and write the stats dataframe
    dfstats = pandas.DataFrame(np.array([[nROI, amean, gmean]]), columns=['nROI', 'amean', 'gmean'])
    dfstats.to_csv(writeto + "fslstats-cosica", sep='\t', index=False)

# -----


# ----------------------------------------------------------------
#               Begin reformatting procedure
# ----------------------------------------------------------------

if os.path.exists(relative_output_root):
    shutil.rmtree(relative_output_root)

mkoutdirs(relative_output_root)

for rootdir, subjectdirs, files in os.walk(subjectroot):
    totalsubj = len(subjectdirs)
    thissubj = 0
    for subj in subjectdirs:
        thissubj += 1
        printProgressBar(thissubj, totalsubj, prefix='Reformatting FSL subject files to Coscia format:', suffix='Complete', length=100)

        fslSubj = os.path.join(subjectroot,subj)
        fslSubjMatPath = fslSubj + "/" + fdtconnectivity

        # Read the raw csv into a numpy subject array
        subjArr = readSubjectCSV(fslSubjMatPath)

        # symmetrize the FSL data (which is not a priori symmetric)
        sym = 0.5 * (subjArr + np.matrix.transpose(subjArr))

        # generate the Pandas table that forms the basis of
        # the output format needed by Coscia's thresholding routines
        dataframe = generateSymmetricCosicasTable(sym)

        # Make the subject output directory
        outsubj = relative_output_root + subj
        if not (os.path.exists(outsubj)):
            os.mkdir(outsubj)

        # Write the reformatted file with tab delimiters
        # (as expected by Cosima's file format)
        targetFile = outsubj + "/" + fdtconnectivity + "-cosica"
        dataframe.to_csv(targetFile, sep='\t', index=False)

        writeDataframeStats(dataframe,outsubj)
