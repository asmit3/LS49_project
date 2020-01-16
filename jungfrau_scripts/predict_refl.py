from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.spot_prediction import StillsReflectionPredictor
import libtbx.load_env
import os

#ts='20180501143632300' # Used last in year 2019
ts='20180501143546717'
LS49_regression= libtbx.env.find_in_repositories(
        relative_path="LS49_regression",
        test=os.path.isdir)
diffBragg_dir = os.path.join(LS49_regression, 'diffBragg_work', 'iota_LS49_regression_r0222')
rayonix_expt=ExperimentListFactory.from_json_file(os.path.join(diffBragg_dir,'rayonix_expt/idx-'+ts+'_integrated.expt'), check_format=False)
jungfrau_det=ExperimentListFactory.from_json_file('rotated_plus_90.json', check_format=False)
jungfrau_cbf_expt = ExperimentListFactory.from_filenames([os.path.join(diffBragg_dir, 'jungfrau_cbf/jungfrauhit_'+ts+'.cbf')])
#rayonix_expt[0].detector=jungfrau_det[0].detector

jungfrau_cbf_expt[0].detector = jungfrau_det[0].detector
jungfrau_cbf_expt[0].crystal = rayonix_expt[0].crystal

# Set mosaicity value to 1093.0 as estimated as median value from work done on psana
jungfrau_cbf_expt[0].crystal.set_domain_size_ang(493.0)
jungfrau_cbf_expt[0].crystal.set_half_mosaicity_deg(0.01)


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
