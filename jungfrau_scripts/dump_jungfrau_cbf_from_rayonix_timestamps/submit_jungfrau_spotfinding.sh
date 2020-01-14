RUN=$1
TRIAL=$2
RG=$3
BASE_PATH=$4
EXP=$5
PHOTON_ENERGY=$6

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"
OUT_DIR=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out

##cctbx_args="$cctbx_args iota.timeout_cutoff_sec=2400"
cctbx_args="spotfinder.lookup.mask=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/mask.pickle output.output_dir=${OUT_DIR} output.logging_dir=${OUT_DIR}"
echo $cctbx_args
t_start=`date +%s`
mpirun libtbx.python /reg/neh/home5/asmit/ExaFEL/dials/modules/exafel_project/ADSE13_25/spot_finding/find_spots.py ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/data.loc ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil ${cctbx_args}
