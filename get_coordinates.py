import nibabel as nib
from nilearn.plotting import find_parcellation_cut_coords
import numpy as np

# SET PATH TO PARCELLATIO NIFTI IMAGE 
root_path='/home/fs0/exet5033/Connectomes/standardconnectome/parcellation/parcellation-files/'

for i in range(1,6):
    parcellation_path = root_path + 'sub-01_label-L2018_desc-scale%d_atlas.nii.gz' %i

    # gets x, y, z coordinates from input nifti image
    coordinates = find_parcellation_cut_coords(parcellation_path)

    output_path = root_path + 'sub-01_label-L2018_desc-scale%d_atlas_coordinates.csv' %i
    np.savetxt(output_path, coordinates)
