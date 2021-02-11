from fsl.data.image import Image
from fsl.utils.image import resample
import numpy as np
import os

# Split atlas into individual rois

dirname  = '/home/fs0/exet5033/Connectomes/hippsubfields_connectome'
filename = os.path.join(dirname,'hippoc-subfields/parcellation-files/sub-01_label-L2018_desc-scale1_atlas.nii.gz')
atlas    = Image(filename)
ref      = Image('/opt/fmrib/fsl/data/standard/MNI152_T1_2mm.nii.gz')


if not os.path.exists(os.path.join(dirname,'rois')):
    os.makedirs(os.path.join(dirname,'rois'))


for i in np.unique(atlas.data):    
    if i>0:
        #print(i)
        x = np.zeros(atlas.shape)
        x[atlas.data==i] = 1
        img = Image(x,xform=atlas.voxToWorldMat)
        data,xform = resample.resampleToReference(img,ref,order=0)
        img = Image(data,xform=xform)
        img.save(os.path.join(dirname,'rois',f'roi{int(i):03d}.nii.gz'))
        
