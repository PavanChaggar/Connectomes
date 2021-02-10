import nibabel as nib
from nilearn.plotting import find_parcellation_cut_coords
import numpy as np

# SET PATH TO PARCELLATIO NIFTI IMAGE 
input_path='<insert here>'

# SET OUPTUT PATH 
output_path='<insert here>'

coordinates = find_parcellation_cut_coords(parcellation)

np.savetxt(output_path, coordinates)
