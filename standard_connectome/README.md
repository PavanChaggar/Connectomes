# Standard Connectome

This directory contains the connectomes generated using the standard Lausanne atlas 
and HCP data.  

At present, connectomes at scales 1 and 5 have been processed and can be found in their respective directories. Within these directories, individual subject dat can be found in the subdirectory `subjects`. 

Inside each of the subject directories, there exist: 
* fdt_network_matrix
* fdt_network_matrix_lengths
* fdt_network_matrix.nii.gz
* fdt_network_matrix_lengths.nii.gz
* NumSeeds_of_ROIs 
* probtrackx.log
* waytotal 

Of particular importance are fdt_network_matrix and fdt_network_matrix_lengths
which each contain an adjacency matricies weighted by number of connections and 
connection length, respectively. Their corresponding nifti images are provided 
in compressed form. 
