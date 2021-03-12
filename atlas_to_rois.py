from fsl.data.image import Image
from fsl.utils.image import resample
import numpy as np
import os

# Split atlas into individual rois

dirname  = '/home/fs0/exet5033/Connectomes/standard_connectome'
scale = 5
filename = os.path.join(dirname,'parcellation/parcellation-files/sub-01_label-L2018_desc-scale%d_atlas.nii.gz' %scale)
out_dir = os.path.join(dirname, 'scale%d/rois' %scale)
atlas    = Image(filename)
ref      = Image('/opt/fmrib/fsl/data/standard/MNI152_T1_2mm.nii.gz')


if not os.path.exists(out_dir):
    os.makedirs(out_dir)


for i in np.unique(atlas.data):    
    if i>0:
        #print(i)
        x = np.zeros(atlas.shape)
        x[atlas.data==i] = 1
        img = Image(x,xform=atlas.voxToWorldMat)
        data,xform = resample.resampleToReference(img,ref,order=0)
        img = Image(data,xform=xform)
        img.save(os.path.join(out_dir,f'roi{int(i):04d}.nii.gz'))
        
