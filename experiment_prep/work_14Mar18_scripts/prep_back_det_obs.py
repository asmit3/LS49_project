from dxtbx.datablock import DataBlockFactory, DataBlockDumper
from dxtbx.model.detector import DetectorFactory
from dials.array_family import flex
import numpy as np
import sys

assert len(sys.argv) == 4, 'Need to supply dx,dy,dz'
dx = float(sys.argv[1])
dy = float(sys.argv[2])
dz = float(sys.argv[3])

db = DataBlockFactory.from_json_file('../../idx-step5_MPIbatch_000000_jungfrau_datablock.json', check_format=False)
datab=db[0]
d = datab.unique_detectors()[0]
b = datab.unique_beams()[0]
p = d[0]
p.get_beam_centre_px(b.get_s0())
error_in_dim = (dx,dy,dz) # in mm
#error_in_dim = (50., 50., 0.) # in mm
jungfrau_z_offset = error_in_dim[2] # in mm 

#from IPython import embed; embed(); exit()
# Reference: https://www.psi.ch/detectors/jungfrau
true_jungfrau = DetectorFactory.simple(            # beam center
                  'PAD', 331, (0,0), '+x', '-y',
                                  # pixel size image size
                  (0.075, 0.075), (4*256,4*256),
                                  # trusted range
                   (-1, 2e6))
multiplier = 4.35
true_panel = true_jungfrau[0]
pix_size = true_panel.get_pixel_size()[0]
hole_size = 7.0833333
true_panel.set_frame(true_panel.get_fast_axis(), true_panel.get_slow_axis(),
                (multiplier*((-pix_size*true_panel.get_image_size()[0])-hole_size), pix_size*true_panel.get_image_size()[1]/2, -331))
print true_jungfrau

image_size_fast = true_panel.get_image_size()[0]
image_size_slow = true_panel.get_image_size()[1]

import copy
wrong_detector = copy.deepcopy(true_jungfrau)
wrong_panel = wrong_detector[0]
wrong_panel.set_frame(true_panel.get_fast_axis(), true_panel.get_slow_axis(), (true_panel.get_origin()[0]-error_in_dim[0], true_panel.get_origin()[1]-error_in_dim[1], -331.0+error_in_dim[2] ))
print wrong_panel

from libtbx import easy_pickle
for i in xrange(1000):
  data = easy_pickle.load('../../idx-step5_MPIbatch_%06d_jungfrau_strong.pickle'%i)
  data = data.select((data['intensity.sum.value']/flex.sqrt(data['intensity.sum.variance']))>100.0)
  if len(data) == 0:
    continue
  obs = flex.vec3_double()
  for j in xrange(len(data)):
    c = p.get_pixel_lab_coord(data['xyzobs.px.value'][j][0:2])
    c_true = true_panel.get_ray_intersection_px(c)
    obs.append((c_true[0], c_true[1], 0.0))
  data['xyzobs.px.value'] = obs 
  easy_pickle.dump('wrong_panel_spots_%06d.pickle'%i,data)
datab.extract_imagesets()[0].set_detector(wrong_detector)
dump = DataBlockDumper(db)
dump.as_file('wrong_detector_datablock.json')
