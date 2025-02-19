#!/usr/bin/env bash
#
# Resample and combine PET images by target.
#
# Download raw PET images from https://github.com/netneurolab/hansen_receptors
#
# GZ


src_petdir=~/Desktop/Toolbox/hansen_receptors-main/data/PET_nifti_images

# saving directory
pet_dir=~/Desktop/test/REACT/PET_Atlas/Downsampled


## Changes below does not seem to be necessary
std_brain=${FSLDIR}/data/standard/MNI152_T1_2mm_brain

# Combine by receptor
d=$(dirname ${pet_dir})
comb_petdir=${d}/CombineByReceptor

printf "Downsampling raw PET images ...\n"
[ -d ${pet_dir} ] || mkdir -p ${pet_dir}
pet_nam=($( find ${src_petdir} -type f -name "[0-9|a-z|A-Z]*" ))
ids=()
for ((i=0; i<${#pet_nam[@]}; i++)); do
    vol_nam=$( basename ${pet_nam[$i]})
    nam=$( echo ${vol_nam} | cut -d '_' -f1)
    nam=$( echo ${nam} | cut -d '-' -f1)
    ids+=(${nam})
    # Resample to 2 mm
    vol_nam=$(remove_ext ${vol_nam})
    atlas_file=${pet_dir}/${vol_nam}_atlas_2mm.nii.gz
    [ -f "${atlas_file}" ] || flirt -in ${pet_nam[$i]} -ref ${std_brain} -applyxfm -usesqform -out ${atlas_file}
done

IFS=" " read -r -a ids <<< "$(echo "${ids[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')"

# Check if it's correct
printf "Unique IDs (%s): " ${#ids[@]}
echo "${ids[@]}"
# 5HT1a 5HT1b 5HT2a 5HT4 5HT6 5HTT A4B2 CB1 D1 D2 DAT FDOPA GABAa H3 M1 MU NAT NMDA VAChT mGluR5

printf "Merging PET images by receptor ...\n"
[ -d ${comb_petdir} ] || mkdir -p ${comb_petdir}
for ((i=0; i<${#ids[@]}; i++)); do
    save_file=${comb_petdir}/Combined_${ids[${i}]}_2mm
    fslmerge -t ${save_file} $(ls ${pet_dir}/${ids[$i]}*_atlas_2mm.nii.gz)
done
