from dxtbx.datablock import DataBlockFactory
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.spot_prediction import StillsReflectionPredictor
from dxtbx.model.detector import DetectorFactory
from dials.array_family import flex
import numpy as np

db = DataBlockFactory.from_json_file('./wrong_detector_datablock.json', check_format=False)
datab=db[0]
d = datab.unique_detectors()[0]
b = datab.unique_beams()[0]
p = d[0]
p.get_beam_centre_px(b.get_s0())
padding = 133 # in pix
#p.set_frame(p.get_fast_axis(), p.get_slow_axis(), (-370.71644483097697, 38.77388127287489, -338.00842424245764))
p.set_frame(p.get_fast_axis(), p.get_slow_axis(), (-423.45834800551455, -18.11826751359136, -323.1693454250518))
#p.set_frame(p.get_fast_axis(), p.get_slow_axis(), (-364.892, 38.4, -331.0 ))


image_size_fast, image_size_slow = p.get_image_size()
big_jungfrau = DetectorFactory.simple(
      sensor = 'UNKNOWN',
      distance = -p.get_origin()[2],
      beam_centre = (0., 0.),
      fast_direction = '+x',
      slow_direction = '-y',
      pixel_size = (p.get_pixel_size()[0], p.get_pixel_size()[0]),
      image_size = (image_size_fast+2*padding, image_size_slow+2*padding ),
      trusted_range = (-1, 2e6),
      mask = [])

big_panel = big_jungfrau[0]
big_jungfrau[0].set_frame(big_jungfrau[0].get_fast_axis(), big_jungfrau[0].get_slow_axis(), (p.get_pixel_lab_coord((-padding, -padding))[0], p.get_pixel_lab_coord((-padding,-padding))[1],p.get_origin()[2]))

print p, big_panel
xobs = []
yobs = []
xcal = []
ycal = []

from libtbx import easy_pickle
for i in xrange(1000):
  try:
    experiments = ExperimentListFactory.from_json_file('../../../rayonix/idx-step5_MPIbatch_%06d_rayonix_refined_experiments.json'%i, check_format=False)
  except Exception:
    continue
  experiments[0].detector = big_jungfrau
  predictor = StillsReflectionPredictor(experiments[0] )
  ubx = predictor.for_ub(experiments[0].crystal.get_A())
  try:
    data = easy_pickle.load('./wrong_panel_spots_%06d.pickle'%i)
  except Exception:
    continue
  data = data.select((data['intensity.sum.value']/flex.sqrt(data['intensity.sum.variance']))>100.0)
  if len(ubx) == 0 or len(data) == 0:
    continue
  for j in xrange(len(ubx)):
    c = experiments[0].detector[0].get_pixel_lab_coord(ubx['xyzcal.px'][j][0:2])
    xcal.append(c[0])
    ycal.append(c[1])

  for j in xrange(len(data)):
    c = p.get_pixel_lab_coord(data['xyzobs.px.value'][j][0:2])
    xobs.append(c[0])
    yobs.append(c[1])

  print 'PRED=%d, OBS = %d'%(len(ubx), len(data)), list(data['xyzobs.px.value'])

import matplotlib.pyplot as plt
plt.scatter(xobs,yobs)
plt.scatter(xcal,ycal)
plt.show()



