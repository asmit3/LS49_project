#!/bin/bash -l 

#SBATCH -C knl
#SBATCH -N 384
#SBATCH -q debug
#SBATCH --image=docker:asmit3/ls49_gcc48:latest
#SBATCH -J L1419
#SBATCH -t 0:15:00 
#SBATCH -A m2859 
#SBATCH --module=mpich-cle6
#SBATCH --mail-user=abhowmick@lbl.gov

NODES=384
NUM_RANKS=$((NODES*68))
echo $NUM_RANKS

RUN=999 # Does not matter, runs specified manually downstairs
TRIAL=0
RG=2
BASE_PATH=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3
EXP=mfxls4916
GEOM=${BASE_PATH}/flat_refined_geometry/t011_temp2_141.9.json


for i in `seq 175 255`
do
  RUN_DIR="r$(printf "%04d" ${i})"
  TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
  RG_3d="rg_$(printf "%03d" ${RG})"

  mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out/
  mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/stdout/
done

sbcast -p ${BASE_PATH}/phil_files/${RG_3d}.phil /tmp/params_1.phil
sbcast -p ${BASE_PATH}/input/mask_r4.pickle /tmp/mask.pickle
sbcast -p /global/u2/a/asmit/dials/modules/LS49_project/SAD/xtc_process_iota_ransac_sad.py /tmp/xtc_process.py
sbcast -p /global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/input/FEE_r143_v1.pickle /tmp/FEE_r143_v1.pickle

#Optimal, getting 4x savings on KNL with omp
export OMP_PROC_BIND=true
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export PMI_MMAP_SYNC_WAIT_TIME=600
module unload craype-hugepages2M

t_start=`date +%s`
# this is where the fun begins


srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 175 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 176 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 177 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 178 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 179 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 183 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 13 -n 884 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 184 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 7 -n 476 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 185 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 10 -n 680 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 186 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 3 -n 204 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 187 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 3 -n 204 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 188 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 3 -n 204 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 189 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 190 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 3 -n 204 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 191 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 2 -n 136 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 192 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 193 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 194 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 195 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 199 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 7 -n 476 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 200 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 201 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 7 -n 476 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 202 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 9 -n 612 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 203 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 3 -n 204 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 204 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 7 -n 476 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 208 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 2 -n 136 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 209 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 210 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 211 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 212 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 213 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 214 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 215 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 216 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 217 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 3 -n 204 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 218 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 219 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 7 -n 476 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 220 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 221 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 1 -n 68 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 222 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 12 -n 816 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 223 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 12 -n 816 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 224 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 9 -n 612 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 225 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 9 -n 612 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 226 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 9 -n 612 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 227 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 11 -n 748 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 228 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 8 -n 544 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 229 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 9 -n 612 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 230 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 10 -n 680 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 231 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 2 -n 136 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 232 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 233 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 234 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 235 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 236 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 7 -n 476 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 237 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 3 -n 204 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 238 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 240 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 8 -n 544 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 241 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 242 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 243 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 244 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 8 -n 544 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 245 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 246 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 19 -n 1292 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 247 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 250 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 251 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 4 -n 272 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 252 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 5 -n 340 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 253 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 8 -n 544 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 254 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
srun -N 6 -n 408 -c 4 --cpu_bind=cores shifter ${BASE_PATH}/srun.sh 255 ${TRIAL} ${RG} ${BASE_PATH} ${EXP} ${GEOM} &
wait

t_end=`date +%s`
echo IOTA_XTC_JobCompleted_TimeElapsed $((t_end-t_start)) $t_start $t_end ${RUN}
