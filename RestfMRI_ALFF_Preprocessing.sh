#!/bin/bash

# Despike
# Nuisance regression
# Global normalization
# Detrend
# Register to MNI
# Spatial smoothing


#[ $# -lt 1 ] && start_loc=1 || start_loc=$1
#[ $# -eq 2 ] && end_loc=$2 || end_loc=0

start_loc=0
end_loc=26

afni_dir=~/Desktop/Toolbox/afni
wkpath=/Volumes/iEEG/RestOlfAnaly/SubjectData
script_dir=~/Desktop/Utilities/Projects/RestfMRI
subjlist=/Volumes/iEEG/RestOlfAnaly/nasal_oral_subject_list.txt

subjlist=~/Desktop/Data/Results/RestfMRI_VBM/eLife_subj_list.txt


# spatial smoothing paramters
sigma=3

# band-pass filtering frequencies, in Hz
lp=0.1
hp=0.008

export PATH=$PATH:${afni_dir}
stdbrain=${FSLDIR}/data/standard/MNI152_T1_2mm_brain

subjs=($(cat ${subjlist} 2>/dev/null))

nb_subj=${#subjs[@]}
[[ ${end_loc} -le 0 || ${end_loc} -gt ${nb_subj} ]] && end_loc=${nb_subj}
start_loc=$(( ${start_loc} - 1))
[ ${start_loc} -le 0 ] && start_loc=0
[ ${start_loc} -gt ${end_loc} ] && start_loc=${end_loc}

for ((i=${start_loc}; i<${end_loc}; i++)); do
    subj=${subjs[$i]}
    for task in Nasal10Min Mouth10Min; do
        featdir=${wkpath}/${subj}/${task}/${task}.feat        
        [ -d ${featdir} ] || continue

        fsmaksdir=${wkpath}/${subj}/anat/NuiMask

        # Motion corrected data (FSL Feat)
        mc_vol=${featdir}/filtered_func_data

        reg=${featdir}/reg/example_func2standard.mat
        warp=${featdir}/reg/example_func2standard_warp.nii.gz

        # Despike
        dspk_vol=${mc_vol}_dspk
        [ -f ${dspk_vol}.nii.gz ] || ${afni_dir}/3dDespike -prefix ${dspk_vol}.nii.gz ${mc_vol}.nii.gz

        mc_file=${featdir}/mc/prefiltered_func_data_mcf.par
        mc_file24=${featdir}/mc/prefiltered_func_data_mcf.par_24motion.txt
        tsdir=${featdir}/TimeSeries/Nuisance_alff
        # white matter, csf, whole brain, 24 motion
        mc_surfix=WmCSFWb24Motion
        nuicorr_design=${tsdir}/nuisance_regressor_${mc_surfix}.txt

        nuicorr_vol=${dspk_vol}_${mc_surfix}
        if [[ -f ${nuicorr_vol}.nii.gz ]]; then
            pass
        else
            ${script_dir}/Extract_Brain_Nuisance.sh ${featdir} ${dspk_vol} ${fsmaksdir} Nuisance_alff
            paste -d ' ' ${mc_file24} ${tsdir}/wholebrain_wm_csf.txt > ${nuicorr_design}
            fsl_glm -i ${dspk_vol} -d ${nuicorr_design} -o ${nuicorr_vol}_beta --out_res=${nuicorr_vol}_res --demean
            fslmaths ${dspk_vol} -Tmean -add ${nuicorr_vol}_res ${nuicorr_vol} -odt float
            #[ -f ${nuicorr_vol}_res.nii.gz ] && rm ${nuicorr_vol}_res.nii.gz
        fi


        vol4norm=${dspk_vol}_${mc_surfix}

        # Global intensity normalization
        printf "Global normalization\n"
        [ -f ${vol4norm}_globalnorm.nii.gz ] || fslmaths ${vol4norm} -ing 10000 ${vol4norm}_globalnorm -odt float

        # detrend
        detrend_vol=${vol4norm}_globalnorm_rlt
        printf "Detrending\n"
        [ -f ${detrend_vol}.nii.gz ] || 3dTcat -rlt+ -prefix ${detrend_vol}.nii.gz ${vol4norm}_globalnorm.nii.gz

        printf "Registering to MNI\n"
        applywarp -i ${detrend_vol} -r ${stdbrain} -o ${detrend_vol}2std_nonlinear -w ${warp}
        #flirt -in ${detrend_vol} -ref ${stdbrain} -applyxfm -init ${reg} -out ${detrend_vol}2std

        printf "Spatial smoothing\n"
        fslmaths ${detrend_vol}2std_nonlinear -kernel gauss ${sigma} -fmean ${detrend_vol}2std_nonlinear_s${sigma} -odt float

       # printf "0.25 Hz low pass\n"
       # [ -f ${detrend_vol}_filtLP025.nii.gz ] || 3dFourier -lowpass 0.25 -prefix ${detrend_vol}_filtLP025.nii.gz ${detrend_vol}.nii.gz

       # printf "Registering low-pass filtered image to MNI\n"
       # #flirt -in ${detrend_vol}_filtLP025 -ref ${stdbrain} -applyxfm -init ${reg} -out ${detrend_vol}_filtLP0252std
       # applywarp -i ${detrend_vol}_filtLP025 -r ${stdbrain} -o ${detrend_vol}_filtLP0252std_nonlinear -w ${warp}

       # printf "Spatial smoothing low-pass filtered image\n"
       # fslmaths ${detrend_vol}_filtLP0252std_nonlinear -kernel gauss ${sigma} -fmean ${detrend_vol}_filtLP0252std_nonlinear_s${sigma} -odt float

       # [ -f ${detrend_vol}_filtLP025.nii.gz ] && rm ${detrend_vol}_filtLP025.nii.gz
       # [ -f ${detrend_vol}_filtLP0252std_nonlinear.nii.gz ] && rm ${detrend_vol}_filtLP0252std_nonlinear.nii.gz

        # band-pass filtering and smooth for functional connectivity analysis
        filt_vol=${detrend_vol}_filtBP00801
        [ -f ${filt_vol}.nii.gz ] && rm ${filt_vol}.nii.gz
        3dFourier -lowpass ${lp} -highpass ${hp} -prefix ${filt_vol}.nii.gz ${detrend_vol}.nii.gz

        # registration to standard brain
        applywarp -i ${filt_vol} -r ${stdbrain} -o ${filt_vol}2std_nonlinear -w ${warp}
        fslmaths ${filt_vol}2std_nonlinear -kernel gauss ${sigma} -fmean ${filt_vol}2std_nonlinear_s${sigma} -odt float

        [ -f ${filt_vol}.nii.gz ] && rm ${filt_vol}.nii.gz

    done # Condition loop
done # subject loop
