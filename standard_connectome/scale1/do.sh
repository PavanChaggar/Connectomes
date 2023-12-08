#!/bin/bash

hcp_folder=/vols/Scratch/HCP # dir containing all the subjects and bedpostx 
out_folder=/home/fs0/exet5033/Connectomes/standard_connectome/scale1 # dir to save output
filename=/vols/Scratch/HCP/Diffusion/Q1200/all_subjects # .csv file pointing to subject ID numbers, e.g. 100307

#rm $out_folder/commands/command.txt
#mkdir $out_folder/commands


for ((i = 170 ; i < 171 ; i++)); do
	subj=$(sed -n $i'p' $filename)

    #100307  119025  137431  155635  175136  195445  214019  341834  486759  599671  744553  878877

    do_split=False  


    # step 0 : turn parcellation to individual rois
    # Need to run python script to generate nifti files for each roi.

    if [ "$do_split" == True  ];then
        python3 atlas_to_rois.py
    fi

    # create text file with list of seed masks
    echo $out_folder/rois/roi???.nii.gz > $out_folder/roilist
    seeds=$out_folder/roilist


    # step 1 : run probtrackx
    # main input files needed to run probtrackx
    bpxdir=$hcp_folder/Diffusion/Q1200/$subj/T1w/Diffusion.bedpostX
    xfm=$hcp_folder/Structural/Q1200/$subj/MNINonLinear/xfms/standard2acpc_dc.nii.gz
    invxfm=$hcp_folder/Structural/Q1200/$subj/MNINonLinear/xfms/acpc_dc2standard.nii.gz

    # compile options
    o=" --opd --loopcheck --forcedir" 
    o=" $o -P 10000 --sampvox=1 --randfib=1 --network"
    o=" $o --xfm=$xfm --invxfm=$invxfm"
    o=" $o --samples=$bpxdir/merged --mask=$bpxdir/nodif_brain_mask"
    o=" $o --seed=$seeds --ompl -V 1 --dir=$out_folder/subjects/$subj"

    # run

    # On GPU
    echo $FSLDIR/bin/probtrackx2_gpu $o >> $out_folder/commands/command.txt

    # On CPU
    #fsl_sub -q veryshort.q $FSLDIR/bin/probtrackx2 $o 

done

fsl_sub -q cuda.q -t $out_folder/commands/command.txt
