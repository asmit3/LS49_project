trial_rg=009_rg002

for run in `seq 227 255`
do
  echo "Processing run # " ${run}
  
  srun -n 56 libtbx.python /global/common/software/m2859/asmit/dials/modules/LS49_project/SAD/reintegrate_with_specified_mosaicity.py output.output_dir=${PWD}/r0${run}/${trial_rg}/out_reint reintegrate.initial_directory=${PWD}/r0${run}/${trial_rg}/out reintegrate.set_domain_size_ang_value=1093.0 reintegrate.set_mosaic_half_deg_value=0.0001 mp.method=mpi ${PWD}/phil_files/reintegrate_000_rg002.phil
done
