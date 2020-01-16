#!/bin/bash

source ~/libtbx_env.sh


RUN=222 #$1
TRIAL=0
RG=3 #$3 
BASE_PATH=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3
EXP=mfxls4916
PHOTON_ENERGY=0 #$5

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"

mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out/
mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/stdout/

cp ${BASE_PATH}/input/jungfrau_spotfinding_1.phil ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil
cp ${BASE_PATH}/input/new_jungfrau_mask_panel13.pickle ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/mask.pickle

cp  ${BASE_PATH}/submit_jungfrau_spotfinding.sh  ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/submit_jungfrau_spotfinding.sh


locator_str="experiment=${EXP}\nrun=${RUN}\ndetector_address=MfxEndstation.0:Jungfrau.0\ncalib_dir=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/input/calib"
echo -e $locator_str > ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/data.loc


t_start=`date +%s`


echo ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${PHOTON_ENERGY} 
#exit
bash ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/submit_jungfrau_spotfinding.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${PHOTON_ENERGY}  
#srun -n 16 -c 4 ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/submit_jungfrau_spotfinding.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${PHOTON_ENERGY}  

t_end=`date +%s`
sleep 5

