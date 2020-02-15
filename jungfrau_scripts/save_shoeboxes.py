from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.spot_prediction import StillsReflectionPredictor
import libtbx.load_env
import os, time

message = ''' Script for saving shoeboxes on Jungfrau given a rayonix experiment and jungfrau cbf. Hard codes some
              paths. Please be careful '''


class myShoeboxWriter(object):
  def __init__ (self, reflections, experiments, executor=None):
    self.reflections=reflections
    self.experiments=experiments
    self.executor=executor


def save_shoeboxes(ts='20180501132541396', outdir='.'):

  print ('Working on timestamp = %s'%ts)
  # First specify some paths
  diffBragg_dir = os.path.join('/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_grid_search_4_or_more_regression')
  rayonix_expt=ExperimentListFactory.from_json_file(os.path.join(diffBragg_dir,'rayonix_images_4_or_more_spots_r183_255/idx-'+ts+'_integrated.expt'), check_format=False)
  jungfrau_det=ExperimentListFactory.from_json_file('jungfrau_geom_run_183_255_v1.json', check_format=False)
  jungfrau_cbf_expt=ExperimentListFactory.from_json_file('jungfrau_4_or_more_spots/idx-jungfrauhit_%s_imported.expt'%ts)
  
  # Make sure jungfrau cbf has right detector model. Also be careful about beam 
  jungfrau_cbf_expt[0].detector = jungfrau_det[0].detector
  jungfrau_cbf_expt[0].crystal = rayonix_expt[0].crystal
  jungfrau_cbf_expt[0].beam = rayonix_expt[0].beam
  
  # Set domain size to a low value because this will cover more spots
  jungfrau_cbf_expt[0].crystal.set_domain_size_ang(93.0)
  jungfrau_cbf_expt[0].crystal.set_half_mosaicity_deg(0.0001)
  
  predictor=StillsReflectionPredictor(jungfrau_cbf_expt[0])
  ubx=predictor.for_ub(jungfrau_cbf_expt[0].crystal.get_A())
  from libtbx.easy_pickle import dump
  from dials.array_family import flex
  ubx['id']=flex.int(len(ubx),0)
  print ('NUMBER OF PREDS=', len(ubx))
  
  if False:
    fname='my_prediction2'
    jungfrau_cbf_expt.as_file(fname+'.expt')
    dump(fname+'.pickle',ubx)
  ################
  # Hard code some sigma values here. 
  sigma_b=0.15
  sigma_m=0.0
  n_sigma=1.3
  
  from dials.algorithms.profile_model.factory import ProfileModelFactory
  from dials.algorithms.integration.integrator import IntegratorFactory
  from dials.command_line.stills_process import phil_scope as stills_process_phil_scope

  # Set params here
  params=stills_process_phil_scope.extract()
  params.profile.gaussian_rs.parameters.sigma_b=sigma_b
  params.profile.gaussian_rs.parameters.sigma_m=sigma_m
  params.profile.gaussian_rs.parameters.n_sigma=n_sigma
  params.profile.gaussian_rs.centroid_definition='com'
  params.integration.debug.output=True
  params.integration.debug.separate_files=False
  
  experiments = ProfileModelFactory.create(params, jungfrau_cbf_expt, ubx)
  bbox=ubx.compute_bbox(experiments)
  ubx['bbox']=bbox

  # Assign shoebox below but does not care about the image pixels yet
  ubx['shoebox']=flex.shoebox(ubx['panel'], ubx['bbox'], allocate=False, flatten=False)
  ubx.compute_partiality(experiments)
  imageset=experiments[0].imageset
  from dials_algorithms_integration_integrator_ext import ShoeboxProcessor
  from dials.algorithms.integration.integrator import IntegratorExecutor
  from dials.model.data import make_image

  executor = IntegratorExecutor(experiments, None)
  executor.initialize(0, 1, ubx)
  frame0=0
  frame1=1
  save_shoebox=True
  processor=ShoeboxProcessor(ubx, len(imageset.get_detector()), frame0, frame1, save_shoebox)
  image=imageset.get_corrected_data(0)
  
  myshoeboxwriter=myShoeboxWriter(ubx, experiments, executor)
  
  # Loop through the imageset, extract pixels and process reflections
  for i in range(len(imageset)):
      image = imageset.get_corrected_data(i)
      if imageset.is_marked_for_rejection(i):
          mask = tuple(flex.bool(im.accessor(), False) for im in image)
      else:
          mask = imageset.get_mask(i)
      processor.next(make_image(image, mask), myshoeboxwriter.executor)
      del image
      del mask
  assert processor.finished(), "Data processor is not finished"
  
  experiments=myshoeboxwriter.experiments
  reflections=myshoeboxwriter.reflections
  reflections.as_file(os.path.join(outdir, 'jungfrau_shoeboxes_%s.refl'%ts))
  experiments.as_file(os.path.join(outdir, 'jungfrau_shoeboxes_%s.expt'%ts))


def get_timestamps(strong_directory=None):
  assert strong_directory is not None, 'Need to provide a directory of strong.refl to get timestamps'
  all_ts=[]
  for f in os.listdir(strong_directory):
    if 'strong.refl' not in f: continue
    ts=f[16:16+17]
    all_ts.append(ts)
  return all_ts

if __name__ == '__main__':
  print ('Start saving shoeboxes')
  timestamps = get_timestamps('jungfrau_4_or_more_spots')
  for ts in timestamps:
    save_shoeboxes(ts=ts, outdir='out_jungfrau_shoeboxes')
  print ('Done saving shoeboxes')
