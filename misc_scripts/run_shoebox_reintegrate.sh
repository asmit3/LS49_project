export SIT_PSDM_DATA=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3
export SIT_DATA=/global/cscratch1/sd/asmit/LS10/psdm/psdm/data

for run in `seq 175 255`
do
  echo "Processing run # " ${run}
  #rm -r ${PWD}/r0${run}/rayonix_expt 
  libtbx.python ~/dials/modules/LS49_project/SAD/reintegrate_with_shoebox_info.py output.output_dir=${PWD}/r0${run}/rayonix_expt reintegrate.refined_expt_fnames=${PWD}/timestamps_for_LS49_regression_${run}.dat reintegrate.run_number=${run} ${PWD}/phil_files/reintegrate_000_rg002.phil &
done

#mp.method=mpi
