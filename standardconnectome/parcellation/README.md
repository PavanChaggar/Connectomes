This directory contains the input and outputs from the multiscale brain 
parcellator. 

`cmtk-commands` contains the command line scripts needed to generate the
parcellations. 

`parcellation-file` contains the main ouput from the multiscale brain parcellator. 
Within it, one will find the following for each parcellation scale: 
* .graphml (contains hemisphere, name, number information for parcels)
* .nii.gz (contains nifti image representation of the parcellation)
* FreeSurferColourLUT.txt (contains region number, label and RGBA colour spec)
* coordinates.csv (x, y, z coordinates in starndard space)
* stats.tsv (containts index, label, type and mm3 volume)


