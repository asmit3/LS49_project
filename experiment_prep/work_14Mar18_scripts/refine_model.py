from dxtbx.datablock import DataBlockFactory
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.spot_prediction import StillsReflectionPredictor
from dxtbx.model.detector import DetectorFactory
from dials.array_family import flex
import numpy as np
from scitbx.matrix import col
from libtbx import easy_pickle
from scitbx.math import five_number_summary
import matplotlib.pyplot as plt
from scitbx.simplex import simplex_opt
import copy

db = DataBlockFactory.from_json_file('./wrong_detector_datablock.json', check_format=False)
datab=db[0]
d = datab.unique_detectors()[0]
b = datab.unique_beams()[0]
p = d[0]
p.get_beam_centre_px(b.get_s0())
padding = 133 # in pix
UPPER_DIST = 30.36 # obtained from tukey's rule of thumb



class SimplexMinimizer(object):
  def __init__(self, model, seed=None,plot=False):
    self.model = model
    self.plot = plot
    assert seed is not None, 'seed should not be None'
    self.n = 2 # Number of params in shift (x,y,z)
    self.x = flex.double(list(model[0].get_origin()[0:self.n]))
    self.starting_simplex = []
    mt = flex.mersenne_twister(seed)
    random_scale = 5.0
    self.nimages = 1000
    self.experiments = []
    self.data = []
    for i in xrange(self.nimages):
      try:
        self.experiments.append(ExperimentListFactory.from_json_file('../../../rayonix/idx-step5_MPIbatch_%06d_rayonix_refined_experiments.json'%i, check_format=False))
      except Exception:
        self.experiments.append(None)
      try:
        data = easy_pickle.load('./wrong_panel_spots_%06d.pickle'%i)
      except Exception:
        self.data.append(None)
        continue
      self.data.append(data.select((data['intensity.sum.value']/flex.sqrt(data['intensity.sum.variance']))>100.0))
    for i in xrange(self.n+1):
      self.starting_simplex.append(random_scale*(((mt.random_double(self.n))/2.0)-1.0)+self.x)

    self.optimizer = simplex_opt(dimension=self.n, matrix=self.starting_simplex, evaluator=self, tolerance=1e-3)
    self.x = self.optimizer.get_solution()

  def target(self, vector):
    model = copy.deepcopy(self.model)
    panel = model[0]
    true_det_dist = -331.0
    panel.set_frame(panel.get_fast_axis(), panel.get_slow_axis(), tuple(list(vector)+[true_det_dist]))

    image_size_fast, image_size_slow = panel.get_image_size()
    # Deliberately constructing a bigger detector to get predictions that may otherwise
    # be missed due to wrong detector configuration
    big_jungfrau = DetectorFactory.simple(
          sensor = 'UNKNOWN',
          distance = -panel.get_origin()[2],
          beam_centre = (0., 0.),
          fast_direction = '+x',
          slow_direction = '-y',
          pixel_size = (panel.get_pixel_size()[0], panel.get_pixel_size()[0]),
          image_size = (image_size_fast+2*padding, image_size_slow+2*padding ),
          trusted_range = (-1, 2e6),
          mask = [])
    big_panel = big_jungfrau[0]
    big_jungfrau[0].set_frame(big_jungfrau[0].get_fast_axis(), big_jungfrau[0].get_slow_axis(), (panel.get_pixel_lab_coord((-padding, -padding))[0], panel.get_pixel_lab_coord((-padding,-padding))[1],panel.get_origin()[2]))

    xobs = flex.double()
    yobs = flex.double()
    xcal = flex.double()
    ycal = flex.double()
    all_dist = flex.double()
    for i in xrange(self.nimages):
      if self.experiments[i] is None or self.data[i] is None: continue
      self.experiments[i][0].detector = copy.deepcopy(big_jungfrau)
      predictor = StillsReflectionPredictor(self.experiments[i][0] )
      ubx = predictor.for_ub(self.experiments[i][0].crystal.get_A())
      if len(ubx) == 0 or len(self.data[i]) == 0:
        continue
      assert len(self.data[i]) == 1, 'If more than one, need ANN'
      dist = []
      cobs = panel.get_pixel_lab_coord(self.data[i]['xyzobs.px.value'][0][0:2]) # was p.get_pixel
    
      for j in xrange(len(ubx)):
        c = self.experiments[i][0].detector[0].get_pixel_lab_coord(ubx['xyzcal.px'][j][0:2])
        dist.append((col(cobs)-col(c)).length())
      dist_min = min(dist)
      minidx = dist.index(dist_min)
      ccalmin = self.experiments[i][0].detector[0].get_pixel_lab_coord(ubx['xyzcal.px'][minidx][0:2])
      xobs.append(cobs[0])
      yobs.append(cobs[1])
      xcal.append(ccalmin[0])
      ycal.append(ccalmin[1])
      all_dist.append(min(dist))
      #print 'Minimum distance = ',min(dist)
      #print 'PRED=%d, OBS = %d'%(len(ubx), len(data)), list(data['xyzobs.px.value'])
    min_dist, q1_dist, med_dist, q3_dist, max_dist = five_number_summary(all_dist)
    sel = all_dist < ((q3_dist-q1_dist)*1.5/2) + med_dist
    all_dist = all_dist.select(sel)
    xobs = xobs.select(sel)
    yobs = yobs.select(sel)
    xcal = xcal.select(sel)
    ycal = ycal.select(sel)
    if self.plot:
      plt.figure()
      plt.scatter(xobs,yobs)
      plt.scatter(xcal,ycal)
      plt.figure()
      plt.hist(all_dist, bins=10)
      plt.show()
    f = flex.sum(all_dist**2)
    #print 'TARGET STATS = ',f,list(vector),five_number_summary(all_dist)  
    print 'AVERAGE_ERROR_IN_SPOT_DIST = ',np.sqrt(f)/len(all_dist)
    return f

x = SimplexMinimizer(d,seed=42).x
print list(x)

