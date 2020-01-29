# Run on psana folder in the LS49_project/jungfrau_scripts/jungfrau_spotfinding folder
# Run on all the files in folder
mpirun -n 4 libtbx.python /reg/neh/home5/asmit/ExaFEL/dials/modules/exafel_project/ADSE13_25/spot_finding/find_spots.py /reg/neh/home5/asmit/ExaFEL/dials/modules/LS49_regression/diffBragg_work/iota_r0222_cori/jungfrau_cbf jungfrau_spotfinding_1.phil


# Run on a specific image for debugging
#img=jungfrauhit_20180501143547150.cbf # Debugging 2 spots showing up super close to each other in spotfinder
#img=jungfrauhit_20180501143650626.cbf
#libtbx.python /reg/neh/home5/asmit/ExaFEL/dials/modules/exafel_project/ADSE13_25/spot_finding/find_spots.py /reg/neh/home5/asmit/ExaFEL/dials/modules/LS49_regression/diffBragg_work/iota_r0222_cori/jungfrau_cbf/${img} jungfrau_spotfinding_1.phil

