from __future__ import absolute_import, division, print_function
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.spot_prediction import StillsReflectionPredictor
from dials.array_family import flex
import copy
from dxtbx.model.experiment_list import ExperimentListDumper
from libtbx.easy_pickle import dump
from scitbx.matrix import col
import dxtbx
import os
import sys

def displace_jungfrau(old_expt, deltaX, deltaY, deltaZ):
  """ Returns a new detector model with specificed displacement """

  expt = copy.deepcopy(old_expt)

  pg0=expt[0].detector.hierarchy()
  x0,y0,z0=pg0.get_origin()
  x=x0+deltaX
  y=y0+deltaY
  z=z0+deltaZ
  pg0.set_frame(pg0.get_fast_axis(), pg0.get_slow_axis(), (x,y,z))
  '''
  oriD0=col((x,y,z))
  for quad_num in range(2):
    pg1=pg0[quad_num]
    x0,y0,z0=pg1.get_origin()
    x=x0+deltaX
    y=y0+deltaY
    z=z0+deltaZ
    oriQx=col((x,y,z))-oriD0
    pg1.set_local_frame(pg1.get_fast_axis(), pg1.get_slow_axis(), oriQx.elems)
    for asic_num in range(8):
      p=pg1[asic_num]
      x0,y0,z0 = p.get_origin()
      x=x0+deltaX
      y=y0+deltaY
      z=z0+deltaZ
      ori = col((x,y,z))-oriQx
      p.set_local_frame(p.get_fast_axis(), p.get_slow_axis(), ori.elems)
  from dxtbx.model.experiment_list import ExperimentListDumper
  dumper=ExperimentListDumper(expt)
  old_dumper=ExperimentListDumper(old_expt)
  dumper.as_json('new.json')
  old_dumper.as_json('old.json')
  exit()
  '''
  return expt

def displace_jungfrau_v2(old_expt, deltaX, deltaY, deltaZ):
  expt=copy.deepcopy(old_expt)
  for panel in expt[0].detector:
    x0,y0,z0=panel.get_origin()
    x=x0+deltaX
    y=y0+deltaY
    z=z0+deltaY
    panel.set_local_frame(panel.get_fast_axis(), panel.get_slow_axis(), (x,y,z))
  from dxtbx.model.experiment_list import ExperimentListDumper
  dumper=ExperimentListDumper(expt)
  old_dumper=ExperimentListDumper(old_expt)
  dumper.as_json('new.json')
  old_dumper.as_json('old.json')
  exit()
  return expt


# Specify timestamp to process
#ts='20180501102146130'
ts=sys.argv[1]
#ts='20180501143547150'

# Load rayonix integrated experiment
rayonix_fname=os.path.join('rayonix_expt', 'idx-%s_integrated_experiments.json'%ts)
rayonix_expt=ExperimentListFactory.from_json_file(rayonix_fname)
jungfrau_det=ExperimentListFactory.from_json_file('rotated_plus_90.json')

# Load jungfrau image and imageset
jungfrau_fname=os.path.join('jungfrau_cbf', 'jungfrauhit_%s.cbf'%ts)
jungfrau_img=dxtbx.load(jungfrau_fname)
jungfrau_iset=jungfrau_img.get_imageset([jungfrau_fname])

# Replace the detector model in the rayonix expt and the imageset as well with the jungfrau specifics
# Thus we have in the experiment
# 1. Crystal model = from rayonix indexing
# 2. Beam model = from rayonix processing
# 3. Detector model = from jungfrau rotated model
# 4. Imageset = from jungfrau actual cbf

rayonix_expt[0].detector=jungfrau_det[0].detector
rayonix_expt[0].imageset=jungfrau_iset

dx = .0 #mm displacement
dy = .0 #mm displacement
dz = .0 #mm displacement
N_scan = 1 

for i in range(N_scan):
  for j in range(N_scan):
    for k in range(N_scan):

      deltaX=(i-N_scan/2.0)*dx
      deltaY=(j-N_scan/2.0)*dy
      deltaZ=(k-N_scan/2.0)*dz
      new_expt = displace_jungfrau(rayonix_expt, deltaX, deltaY, deltaZ)
      predictor=StillsReflectionPredictor(new_expt[0])
      #from IPython import embed; embed(); exit()
      ubx=predictor.for_ub(new_expt[0].crystal.get_A())
      ubx['id']=flex.int(len(ubx),0)
      should_dump=False
      if len(ubx) > 0:
        print ('NUmber of preds and deltas = ', deltaX, deltaY, deltaZ, len(ubx), ' ', ts)
        should_dump=True
      #print ('NUMBER OF PREDS=', len(ubx))
      if should_dump:
        #fname='my_prediction_%d_%d_%d'%(i,j,k)
        #fname='my_pred_%s'%ts
        fname='my_pred_%s_%s_%s_%s'%(i,j,k,ts)
        dumper=ExperimentListDumper(new_expt)
        dumper.as_json(fname+'.json')
        dump(fname+'.pickle',ubx)

print ('======================================================================================')
