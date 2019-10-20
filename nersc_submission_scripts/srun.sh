#! /bin/sh
export SIT_PSDM_DATA=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3
export SIT_DATA=/global/cscratch1/sd/asmit/LS10/psdm/psdm/data
#export OMP_PROC_BIND=true
#export OMP_NUM_THREADS=4
#export OMP_PLACES=threads

source /build/setpaths.sh

RUN=$1
TRIAL=$2
RG=$3
BASE_PATH=$4
EXP=$5
GEOM=$6


RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"


START_XTC=$(date +"%s")


cctbx_args="input.experiment=${EXP} input.run_num=${RUN} output.logging_dir=/tmp output.output_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out mp.method=mpi input.reference_geometry=${GEOM} output.composite_output=False format.cbf.invalid_pixel_mask=/tmp/mask.pickle /tmp/params_1.phil"

cctbx_args="${cctbx_args} indexing.real_space_grid_smart_search.use_openmp=True"

libtbx.python /tmp/xtc_process.py ${cctbx_args}


#libtbx.python ~/dials/modules/LS49_project/SAD/xtc_process_sad.py input.trial=${TRIAL} output.output_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out output.logging_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/stdout input.experiment=${EXP} input.run_num=${RUN} input.cfg=None input.xtc_dir=None input.use_ffb=False ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil input.rungroup=${RG} mp.method=mpi input.reference_geometry=${GEOM} output.composite_output=False

END_XTC=$(date +"%s")
ELAPSED=$((END_XTC-START_XTC))
echo IOTA_XTC_SingleRank_TimeElapsed ${ELAPSED} ${START_XTC} ${END_XTC} ${RUN}
