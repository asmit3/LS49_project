from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.spot_prediction import StillsReflectionPredictor
import libtbx.load_env
import os, time


class myShoeboxWriter(object):
  def __init__ (self, reflections, experiments, executor=None):
    self.reflections=reflections
    self.experiments=experiments
    self.executor=executor


#ts='20180501143632300' # Used last in year 2019
ts='20180501132541396'
LS49_regression= libtbx.env.find_in_repositories(
        relative_path="LS49_regression",
        test=os.path.isdir)
diffBragg_dir = os.path.join('/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_grid_search_4_or_more_regression')
rayonix_expt=ExperimentListFactory.from_json_file(os.path.join(diffBragg_dir,'rayonix_images_4_or_more_spots_r183_255/idx-'+ts+'_integrated.expt'), check_format=False)
jungfrau_det=ExperimentListFactory.from_json_file('jungfrau_geom_run_183_255_v1.json', check_format=False)
#jungfrau_cbf_expt = ExperimentListFactory.from_filenames([os.path.join(diffBragg_dir, 'jungfrau_cbf/jungfrauhit_'+ts+'.cbf')])
jungfrau_cbf_expt=ExperimentListFactory.from_json_file('jungfrau_4_or_more_spots/idx-jungfrauhit_%s_imported.expt'%ts)
#rayonix_expt[0].detector=jungfrau_det[0].detector

jungfrau_cbf_expt[0].detector = jungfrau_det[0].detector
jungfrau_cbf_expt[0].crystal = rayonix_expt[0].crystal
jungfrau_cbf_expt[0].beam = rayonix_expt[0].beam

# Set mosaicity value to 1093.0 as estimated as median value from work done on psana
jungfrau_cbf_expt[0].crystal.set_domain_size_ang(293.0)
jungfrau_cbf_expt[0].crystal.set_half_mosaicity_deg(0.0001)


#from IPython import embed; embed(); exit()

predictor=StillsReflectionPredictor(jungfrau_cbf_expt[0])
ubx=predictor.for_ub(jungfrau_cbf_expt[0].crystal.get_A())
from libtbx.easy_pickle import dump
from dials.array_family import flex
ubx['id']=flex.int(len(ubx),0)
print ('NUMBER OF PREDS=', len(ubx))

#from dxtbx.model.experiment_list import ExperimentListDumper
#dumper=ExperimentListDumper(jungfrau_det)
fname='my_prediction2'
#dumper.as_json(fname+'.json')
jungfrau_cbf_expt.as_file(fname+'.expt')
dump(fname+'.pickle',ubx)
################
sigma_b=0.15
sigma_m=0.0
n_sigma=1.3

from dials.algorithms.profile_model.factory import ProfileModelFactory
from dials.command_line.stills_process import phil_scope as stills_process_phil_scope
params=stills_process_phil_scope.extract()
params.profile.gaussian_rs.parameters.sigma_b=sigma_b
params.profile.gaussian_rs.parameters.sigma_m=sigma_m
params.profile.gaussian_rs.parameters.sigma_m=n_sigma
#params.profile.gaussian_rs.centroid_definition='com'
#from IPython import embed; embed();exit()
experiments = ProfileModelFactory.create(params, jungfrau_cbf_expt, ubx)
bbox=ubx.compute_bbox(experiments)
ubx['bbox']=bbox
from IPython import embed; embed(); exit()
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
read_time = 0.0
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
reflections.as_file('try4.refl')
experiments.as_file('try4.expt')

#from IPython import embed; embed(); exit()
