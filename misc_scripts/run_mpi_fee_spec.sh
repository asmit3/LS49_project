export SIT_PSDM_DATA=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3
export SIT_DATA=/global/cscratch1/sd/asmit/LS10/psdm/psdm/data

for run in `seq 176 255`
do
  echo "Processing run # " ${run}
  #rm -r ${PWD}/r0${run}/rayonix_expt 
  srun -n 136 -c 4 libtbx.python ~/dials/modules/LS49_project/SAD/mpi_fee_spec_dumper.py output.output_dir=${PWD}/r0${run}/rayonix_expt fee_spec.refined_expt_fnames=${PWD}/timestamps_for_LS49_regression_${run}.dat fee_spec.run_number=${run} ${PWD}/phil_files/reintegrate_000_rg002.phil mp.method=mpi
done

#mp.method=mpi
