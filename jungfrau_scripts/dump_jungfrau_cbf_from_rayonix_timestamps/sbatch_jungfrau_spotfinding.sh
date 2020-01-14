#!/bin/bash

source /reg/neh/home/asmit/ExaFEL/dials/newconda/etc/profile.d/conda.sh
conda activate myEnv
source /reg/neh/home/asmit/ExaFEL/dials/build/setpaths.sh
QUEUE=$1

for i in `seq $2 $3`
do
  echo 'Processing LCLS run ' $i

  RUN=$i #$1
  TRIAL=1
  RG=3 #$3 
  BASE_PATH=/reg/d/psdm/mfx/mfxls4916/scratch/asmit/LS49_2019
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


  locator_str="experiment=${EXP}\nrun=${RUN}\ndetector_address=MfxEndstation.0:Jungfrau.0\ncalib_dir=/reg/data/ana04/mfx/mfxls4916/scratch/asmit/LS49_2019/input/calib"
  echo -e $locator_str > ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/data.loc


  t_start=`date +%s`
  bsub -n 12 -o ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/stdout/log.out -q ${QUEUE} -J ${RUN_DIR} ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/submit_jungfrau_spotfinding.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${PHOTON_ENERGY}
  t_end=`date +%s`
  sleep 5

done
