
#!/bin/bash
# Code for pre-processing HRV-ER fMRI data
# on Vanderbilt cluster (ACCRE)

# "main" directory with raw data
maindir_raw=/data1/neurdylab/datasets/HRV-ER/HRV-ER_raw

# "main" directory that will contain our processed data
maindir_proc=/data1/neurdylab/datasets/HRV-ER/HRV-ER_proc

# paths on platypus
MNI_T1_2mm_template="/data1/neurdylab/MNI152_T1_2mm_brain.nii.gz"
scripts_path="/data1/neurdylab/scripts/vu_meica_pipeline"
afni_init="singularity exec --bind /data1:/data1 ${scripts_path}/afni_cmake_build_AFNI_21.1.03.sif" 

# module load for FSL
module load GCC/8.2.0  OpenMPI/3.1.4
module load FSL/6.0.1-Python-3.7.2

# for ants:
pyscripts_path="/data1/neurdylab/eegfmri_vu_pipeline/scripts/ants_python"

# start loop over subjects

#for sub_id in `cat /data1/neurdylab/datasets/HRV-ER/HRV-ER_raw/subs.txt` # RS: make a text file containing all participants for HRV-ER
for sub_id in `cat /data1/neurdylab/datasets/HRV-ER/derivates/rerun.txt` #RS: rerun on all files that different echo1 name

do
echo $sub_id

# dirs where subject's processed fmri files will be written:
proc_func_dir=${maindir_proc}/${sub_id}/ses-pre/func #RS: changed to ses-pre instead of ses-BAS1
proc_anat_dir=${maindir_proc}/${sub_id}/ses-pre/anat
mkdir -p $proc_func_dir $proc_anat_dir

# subject's raw functional & anatomic dirs:
raw_func_dir=${maindir_raw}/${sub_id}/ses-pre/func
raw_anat_dir=${maindir_raw}/${sub_id}/ses-pre/anat
anat_img=${raw_anat_dir}/${sub_id}_ses-pre_acq-wholebrain_T1w.nii.gz #RS: either this or_ses-pre_acq-wholebrain_T1w.nii.gz 

func_prefix=${sub_id}_ses-pre_task-rest_bold #RS: might need this for echo2 and echo3 right? 

log_path="${proc_func_dir}/logs"
mkdir -p $log_path
log_file="${log_path}/${sub_id}.log"
printf "Starting processing for ${sub_id}\n\n" > ${log_file}
echo "Log file created at ${log_file}"


# anatomic (T1) image processing
# ------------------------------

if ! [ -e ${proc_anat_dir}/T1_bet.nii.gz ]; then

        echo "Bet and unifize T1 anatomical, done once per subject..."

        # Brain Extraction (bet)
        cmd="bet ${anat_img} ${proc_anat_dir}/T1_bet.nii.gz -f 0.5 -o -m" #RS: check 2 or 3 subjects to see if these work 
        printf "Running ${cmd}\n" >> $log_file
        $cmd
fi

T1_bet_img=${proc_anat_dir}/T1_bet.nii.gz

if ! [ -e ${proc_anat_dir}/T1_unifize.nii.gz ]; then

        # unifize
        ${afni_init} 3dcopy ${T1_bet_img} ${proc_anat_dir}/T1_bet
        ${afni_init} 3dUnifize -prefix ${proc_anat_dir}/T1_unifize ${proc_anat_dir}/T1_bet+orig
        ${afni_init} 3dcopy ${proc_anat_dir}/T1_unifize+orig ${proc_anat_dir}/T1_unifize.nii.gz

        rm -f ${proc_anat_dir}/T1_bet+orig*
        rm -f ${proc_anat_dir}/T1_unifize+orig*

fi

# Motion correction old (re-alignment)
#----------------------------------

# if ! [ -e ${proc_func_dir}/${func_prefix}_mo.nii.gz ]; then
#         echo "mcflirt done once per subject"

#         cmd="mcflirt -plots -in ${raw_func_dir}/${func_prefix}.nii.gz \
#         -out ${proc_func_dir}/${func_prefix}_mo.nii.gz"
#         printf "Running ${cmd}\n" >> $log_file
#         echo $cmd
#         $cmd

# fi

# Motion correction new (RS taken from pre-meica)
#-----------------------------------
echo "Applying motion coregistration..."
printf "Applying motion coregistration\n----------------------------------------------\n" >> $log_file

curr_in="${raw_func_dir}/${sub_id}_ses-pre_task-rest_echo-2_bold"
curr_out="${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-2_bold_mo"

# First run motion correction on echo2
if ! [ -f ${proc_func_dir}/*.volreg_mats.aff12.1D ]; then
        cmd="${afni_init} 3dvolreg -verbose -1Dmatrix_save ${curr_out}.volreg_mats -1Dfile ${curr_out}.volreg_par -base ${curr_in}.nii[0] -prefix ${curr_out}.nii -input ${curr_in}.nii"

        printf "Running ${cmd}\n" >> $log_file
        $cmd
fi 

# Use motion paramters from echo2 on echo1 (note: echo1 naming inconsistent across participants)
if ! [-e ${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-1_bold_mo]; then
        mot_co_check_path=$(ls ${proc_func_dir}/*.volreg_mats.aff12.1D) # assumes that there is only one file with that extension

        curr_in="${raw_func_dir}/${sub_id}_ses-pre_task-rest_bold"
        curr_out="${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-1_bold_mo"
        cmd="${afni_init} 3dcopy ${curr_in}.nii ${proc_func_dir}/volreg_tmp"
        printf "Running ${cmd}\n" >> $log_file
        $cmd

        tmp_in="${proc_func_dir}/volreg_tmp+orig.BRIK"
        tmp_out="${proc_func_dir}/volreg_tmp_mo"

        cmd="${afni_init} 3dAllineate -1Dmatrix_apply ${mot_co_check_path} -master ${tmp_in} -prefix ${tmp_out} -source ${tmp_in}"
        printf "Running ${cmd}\n" >> $log_file
        $cmd

        cmd="${afni_init} 3dcopy ${tmp_out}+orig.BRIK ${curr_out}.nii"
        printf "Running ${cmd}\n" >> $log_file
        $cmd

        rm -f ${proc_func_dir}/volreg_tmp*
fi

# Use motion parameters from echo2 on echo3
if ! [-e ${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-3_bold_mo]; then
        curr_in="${raw_func_dir}/${sub_id}_ses-pre_task-rest_echo-3_bold"
        curr_out="${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-3_bold_mo"
        cmd="${afni_init} 3dcopy ${curr_in}.nii ${proc_func_dir}/volreg_tmp"
        printf "Running ${cmd}\n" >> $log_file
        $cmd

        tmp_in="${proc_func_dir}/volreg_tmp+orig.BRIK"
        tmp_out="${proc_func_dir}/volreg_tmp_mo"

        cmd="${afni_init} 3dAllineate -1Dmatrix_apply ${mot_co_check_path} -master ${tmp_in} -prefix ${tmp_out} -source ${tmp_in}"
        printf "Running ${cmd}\n" >> $log_file
        $cmd

        cmd="${afni_init} 3dcopy ${tmp_out}+orig.BRIK ${curr_out}.nii"
        printf "Running ${cmd}\n" >> $log_file
        $cmd

        rm -f ${proc_func_dir}/volreg_tmp*
fi

# RS: Slice timing correction on all echos (note: echo1 name inconsistent across participants)
# ----------------------------------------------
tpattern="alt+z2" # RS for VU Scans, might have to do somethign different for HRV-ER 
echo "Correcting slice-timing..."
printf "Correcting slice-timing\n----------------------------------------------\n" >> $log_file

echo1_in="${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-1_bold_mo"
echo1_out="${echo1_in}_st"
cmd="${afni_init} 3dTshift -tpattern ${tpattern} -prefix ${echo1_out}.nii ${echo1_in}.nii"

printf "Running ${cmd}\n" >> $log_file
$cmd

echo2_in="${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-2_bold_mo"
echo2_out="${echo2_in}_st"
cmd="${afni_init} 3dTshift -tpattern ${tpattern} -prefix ${echo2_out}.nii ${echo2_in}.nii"

printf "Running ${cmd}\n" >> $log_file
$cmd

echo3_in="${proc_func_dir}/${sub_id}_ses-pre_task-rest_echo-3_bold_mo"
echo3_out="${echo3_in}_st"
cmd="${afni_init} 3dTshift -tpattern ${tpattern} -prefix ${echo3_out}.nii ${echo3_in}.nii"

printf "Running ${cmd}\n" >> $log_file
$cmd

printf ".----------------------------------------------\n" >> $log_file

# RS: Obtaining auto masks for Meica
# -------------------------------
out_dir="${proc_func_dir}/meica"
mkdir -p ${out_dir}
nframes=175
# Extract reference EPI middle volume - once per subject scan
if ! [ -f "${proc_func_dir}/oneVol-${sub_id}.nii" ]; then
    echo "Extracting one EPI reference vol, done once per subject scan..."

    curr_in="${echo2_out}.nii"
    voln=$((nframes/2))
    out_fn="${proc_func_dir}/oneVol-${sub_id}"

    ${afni_init} 3dcalc -a ${curr_in}[${voln}] -expr a -prefix ${out_fn}.nii
fi

# Automask middle volume of mot-corr data - once per subject scan
if ! [ -f "${out_dir}/autoMask-${sub_id}.nii" ] && [ -f "${proc_func_dir}/oneVol-${sub_id}.nii" ]; then
    echo "Middle volume automask, done once per subject scan..."

    in_fn="${proc_func_dir}/oneVol-${sub_id}"

    ${afni_init} 3dAutomask -dilate 2 -prefix ${out_dir}/autoMask-${sub_id}.nii ${in_fn}.nii
fi

# RS: ME-ICA using tedana on HRV-ER TE: 18, 35.59, 53.18
echo "Running tedana on $sub_id"
current_path=${out_dir}
source /data1/neurdylab/software/py3env/bin/activate
tedana -d "${echo1_out}.nii" "${echo2_out}.nii" "${echo3_out}.nii" -e 18 35.59 53.18 --mask ${current_path}/autoMask-${sub_id}.nii --out-dir ${current_path} 
deactivate

# register one BOLD volume to same subject's T1

if ! [ -e ${proc_func_dir}/oneVol_EPI2T1.nii.gz ]; then

 echo "epi_reg done once per subject"

        cmd="epi_reg --epi=${proc_func_dir}/oneVol-${sub_id}.nii.gz \
        --t1=${anat_img} \
        --t1brain=${T1_bet_img} \
        --out=${proc_func_dir}/oneVol_EPI2T1.nii.gz"
 echo $cmd
	printf "Running ${cmd}\n" >> $log_file
        $cmd

fi


# ANTS VERSION
# ------------------------------
if ! [ -e ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2T1.nii.gz ]; then

echo "ants reg..."
mkdir -p $proc_func_dir/ants_out

# check playtpus tedana version 
applywarp --ref=${proc_anat_dir}/T1_unifize.nii.gz --in=${out_dir}/dn_ts_OC.nii.gz --out=${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2T1.nii.gz \
--postmat=${proc_func_dir}/oneVol_EPI2T1.mat

fi

# module change for ANTs
module load GCC/6.4.0-2.28 OpenMPI/2.1.1
module load ANTs/2.3.0-Python-2.7.14


if ! [ -e ${proc_anat_dir}/ants_reg_tforms/tform1_0GenericAffine.mat ]; then

    echo "running ANTS!!"
    mkdir -p ${proc_anat_dir}/ants_reg_tforms

    # T1 to MNI
${pyscripts_path}/ants_doReg_accre.py \
                 --in_nii ${proc_anat_dir}/T1_unifize.nii.gz \
                 --ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
                 --out_dir ${proc_anat_dir}/ants_reg_tforms/


# EPI to MNI one Vol
${pyscripts_path}/ants_applyReg_accre.py \
                 --in_nii  ${proc_func_dir}/oneVol_EPI2T1.nii.gz \
                 --out_nii ${proc_func_dir}/ants_out/oneVol_EPI2MNI.nii.gz  \
                 --ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
                 --tform_dir ${proc_anat_dir}/ants_reg_tforms/
fi

# EPI to MNI for whole fMRI 4D
${pyscripts_path}/ants_applyReg_accre.py \
		--in_nii ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2T1.nii.gz \
		--out_nii ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI.nii.gz \
		--ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
		--tform_dir ${proc_anat_dir}/ants_reg_tforms/

# original
# ${pyscripts_path}/ants_applyReg_accre.py \
#                 --in_nii ${current_path}/${sub_id}-${scan_num}_EPI2T1.nii.gz \
#                 --out_nii ${current_path}/${sub_id}-${scan_num}_EPI2MNI.nii.gz \
#                 --ref_nii /data1/neurdylab/MNI152_T1_2mm_brain.nii.gz \
#                 --tform_dir ${base_dir}/reg_tforms/


# module change back to FSL
module load GCC/8.2.0  OpenMPI/3.1.4
module load FSL/6.0.1-Python-3.7.2

# ------------------------------

#Post-registration smoothing and nuisance regression
# --------------------------

scripts_path="/data1/neurdylab/scripts/vu_meica_pipeline"
afni_init="singularity exec --bind /data1:/data1 ${scripts_path}/afni_cmake_build_AFNI_21.1.03.sif"

#Spatial Blurring
if ! [ -e ${proc_func_dir}/${func_prefix}_mo_EPI2MNI_sm.nii.gz ]; then

        echo "spatial blurring for subject"
        cmd="${afni_init} 3dmerge -1blur_fwhm 3.0 -doall -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm.nii.gz" "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI.nii.gz""
        echo $cmd
        printf "Running ${cmd}\n" >> $log_file
        $cmd

fi

#Calculate Mean

if ! [ -e ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_mean.nii.gz ]; then

	echo "calculating mean for subject"

	cmd="${afni_init} 3dTstat -mean -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_mean.nii.gz" "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm.nii.gz""
	echo $cmd
	printf "Running ${cmd}\n" >> $log_file
	$cmd

fi

#Nuisance Regression

if ! [ -e ${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_nr_demean.nii.gz ]; then
	echo "completing nuisance regression for subject"
	cmd="${afni_init} 3dDetrend -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_nr_demean.nii.gz" -polort 4 "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm.nii.gz""
	echo $cmd
	printf "Running ${cmd}\n" >> $log_file
	$cmd

fi


#Add Back Mean

echo "adding back the mean"
${afni_init} 3dcalc -a "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_mean.nii.gz" \
             -b "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_nr_demean.nii.gz" \
             -expr 'a+b' -prefix "${proc_func_dir}/ants_out/${func_prefix}_mo_EPI2MNI_sm_nr.nii.gz"


# QA for ants:
mkdir -p ${proc_func_dir}/ants_out/QA/imgs

slicer "${proc_func_dir}/ants_out/oneVol_EPI2MNI.nii.gz" "/data1/neurdylab/MNI152_T1_2mm_brain.nii.gz"  -s 2 \
        -x 0.35 "${proc_func_dir}/ants_out/QA/imgs/sla.png" -x 0.45 "${proc_func_dir}/ants_out/QA/imgs/slb.png" -x 0.55 "${proc_func_dir}/ants_out/QA/imgs/slc.png" -x 0.65 "${proc_func_dir}/ants_out/QA/imgs/sld.png" \
        -y 0.35 "${proc_func_dir}/ants_out/QA/imgs/sle.png" -y 0.45 "${proc_func_dir}/ants_out/QA/imgs/slf.png" -y 0.55 "${proc_func_dir}/ants_out/QA/imgs/slg.png" -y 0.65 "${proc_func_dir}/ants_out/QA/imgs/slh.png" \
        -z 0.35 "${proc_func_dir}/ants_out/QA/imgs/sli.png" -z 0.45 "${proc_func_dir}/ants_out/QA/imgs/slj.png" -z 0.55 "${proc_func_dir}/ants_out/QA/imgs/slk.png" -z 0.65 "${proc_func_dir}/ants_out/QA/imgs/sll.png"

pngappend "${proc_func_dir}/ants_out/QA/imgs/sla.png" + "${proc_func_dir}/ants_out/QA/imgs/slb.png" + "${proc_func_dir}/ants_out/QA/imgs/slc.png" + "${proc_func_dir}/ants_out/QA/imgs/sld.png" + "${proc_func_dir}/ants_out/QA/imgs/sle.png" \
	+ "${proc_func_dir}/ants_out/QA/imgs/slf.png" + "${proc_func_dir}/ants_out/QA/imgs/slg.png" + "${proc_func_dir}/ants_out/QA/imgs/slh.png" + "${proc_func_dir}/ants_out/QA/imgs/sli.png" + "${proc_func_dir}/ants_out/QA/imgs/slj.png" \
	+ "${proc_func_dir}/ants_out/QA/imgs/slk.png" + "${proc_func_dir}/ants_out/QA/imgs/sll.png" "${proc_func_dir}/ants_out/QA/${func_prefix}_EPI2MNI_ants.png"


# create write permissions
chmod g+w -R $proc_func_dir $proc_anat_dir

done
