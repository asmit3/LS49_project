# Run this script after logging into a data transfer node on nersc
# ssh asmit@dtn01.nersc.gov
#./bbcp -T "ssh -x -a -oFallBackToRsh=no %I -l %U %H bbcp" asmit@psexport.slac.stanford.edu:/reg/d/psdm/mfx/mfxls4916/xtc "/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/MFX/mfxls4916/xtc/"

for run in `seq 247 247`
do 
  for stream in `seq 0 3`
  do
    fname=e1182-r0${run}-s0${stream}-c02.xtc
    echo "Copying File " ${fname}
    bbcp -z -T "ssh -x -a -oFallBackToRsh=no %I -l %U %H bbcp" asmit@psexport.slac.stanford.edu:/reg/d/psdm/mfx/mfxls4916/xtc/${fname} "/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/MFX/mfxls4916/xtc/"
  done
  # s80 stream
  fname=e1182-r0${run}-s80-c00.xtc
  echo "Copying File " ${fname}
  #bbcp -z -T "ssh -x -a -oFallBackToRsh=no %I -l %U %H bbcp" asmit@psexport.slac.stanford.edu:/reg/d/psdm/mfx/mfxls4916/xtc/${fname} "/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/MFX/mfxls4916/xtc/"
  
done
