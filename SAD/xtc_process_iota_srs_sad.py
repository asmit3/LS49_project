from __future__ import division
from six.moves import range
#
PSANA2_VERSION = 0
try:
  import psana
  PSANA2_VERSION = psana.__version__
except ImportError:
  pass # for running at home without psdm build
except AttributeError:
  pass

import time
from LS49_project.FEE_spec.fee_spectra import FEE_spec
from xfel.cftbx.detector import cspad_cbf_tbx
from xfel.cxi.cspad_ana import cspad_tbx, rayonix_tbx
import pycbf, os, sys, copy, socket
import libtbx.load_env
from libtbx.utils import Sorry, Usage
from dials.util.options import OptionParser
from libtbx.phil import parse
from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx.array_family import flex
import numpy as np
from libtbx import easy_pickle
from LS49_project.FEE_spec.fee_spectra import FEE_spec

# phil string imports
iota_phil_str = '''
  iota {
    timeout_cutoff_sec = None
      .type = float
      .help = Timeout cutoff to end iota loop (seconds)
    method = off *random_sub_sampling
      .type = choice
      .help = Type of IOTA processing to be done. \
              off : No IOTA processing is done. \
              random-sub-sampling : randomly sub-sample observed bragg spots and index. Can be done multiple times. See options for random-sub-sampling if this is used.
    random_sub_sampling {
      ntrials = 50
        .type = int
        .help = Number of random sub-samples to be selected
      fraction_sub_sample = 0.8
        .type = float
        .help = fraction of sample to be sub-sampled. Should be between 0 and 1
      consensus_function = *unit_cell
        .type = choice
        .help = choose the type of consensus function to be employed for random_sub_sampling. More details \
                in the functions themselves
      show_plot = False
        .type = bool
        .help = Flag to indicate whether plots for clustering are to be displayed. Useful for debugging
      no_outlier_rejection_and_candidates_refinement=False
        .type = bool
        .help = Flag to indicate if candidate basis vectors should be refined and whether \
                outlier rejectionis needed
  }
}

'''
from LS49_project.SAD.xtc_process_sad import LS49_phil_str
from xfel.command_line.xtc_process import xtc_phil_str, extra_dials_phil_str
from dials.command_line.stills_process import dials_phil_str, program_defaults_phil_str
from xfel.ui import db_phil_str
from xfel.command_line.xfel_process import radial_average_phil_str
phil_scope = parse(iota_phil_str + LS49_phil_str + xtc_phil_str + dials_phil_str + extra_dials_phil_str + db_phil_str + radial_average_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))
# Function imports
from xfel.command_line.xtc_process import filter, run_psana2
# Class imports
from xfel.command_line.xfel_process import Script as DialsProcessScript
from xfel.ui.db.frame_logging import DialsProcessorWithLogging
from xfel.command_line.xtc_process import EventOffsetSerializer
from LS49_project.SAD.xtc_process_sad import InMemScript_LS49_SAD


class IOTA_TimeoutError(Exception):
  """Raises a timeout error if IOTA indexing takes too long """
  pass


class InMemScript_iota(InMemScript_LS49_SAD):

  def index(self, experiments, observed):
    ''' IOTA-SRS indexing. Goes through ntrials and indexes subsamples'''
    # index and refine
    if self.params.iota.method == 'random_sub_sampling':
      from scitbx.array_family import flex
      experiments_list = []
      #No outlier rejection or refinement should be done for the candidate basis vectors
      self.known_crystal_models = None
      outlier_rejection_flag=self.params.indexing.stills.candidate_outlier_rejection
      refine_all_candidates_flag=self.params.indexing.stills.refine_all_candidates
      if self.params.iota.random_sub_sampling.no_outlier_rejection_and_candidates_refinement:
        self.params.indexing.stills.candidate_outlier_rejection=False
        self.params.indexing.stills.refine_all_candidates=False

      # Adding timeout option for IOTA
      initial_time = time.time()
      try:
        experiments_tmp, indexed_tmp = self.index_with_iota(experiments, observed_sample)
        return experiments_tmp, indexed_tmp
      else:
        print ('gotta do iota sampling')
      for trial in range(self.params.iota.random_sub_sampling.ntrials):
        curr_time = time.time()
        if self.params.iota.timeout_cutoff_sec is not None:
          if curr_time - initial_time > self.params.iota.timeout_cutoff_sec:
            raise IOTA_TimeoutError('IOTA_TIMEOUT ',curr_time-initial_time)
        flex.set_random_seed(trial+1001)
        observed_sample = observed.select(flex.random_selection(len(observed), int(len(observed)*self.params.iota.random_sub_sampling.fraction_sub_sample)))
        try:
          print('IOTA:SUM_INTENSITY_VALUE=%d',sum(observed_sample['intensity.sum.value']),' ', trial)
          experiments_tmp, indexed_tmp = self.index_with_iota(experiments, observed_sample)
          #XXX
          return experiments_tmp, indexed_tmp
          experiments_list.append(experiments_tmp)
        except Exception as e:
          print('Indexing failed for some reason', str(e))
      if self.params.iota.random_sub_sampling.consensus_function == 'unit_cell':
        #from IPython import embed; embed(); exit()
        from exafel_project.ADSE13_25.clustering.old_consensus_functions import get_uc_consensus as get_consensus
          #known_crystal_models = get_consensus(experiments_list, show_plot=self.params.iota.random_sub_sampling.show_plot, return_only_first_indexed_model = False)
        if len(experiments_list) > 0:
          known_crystal_models, clustered_experiments_list = get_consensus(experiments_list, show_plot=False, return_only_first_indexed_model=True, finalize_method=None, clustering_params=None)
          self.known_crystal_models = known_crystal_models
        print ('IOTA: Reindexing with best chosen crystal model')
          # Set back whatever PHIL parameter was supplied by user for outlier rejection and refinement
        self.params.indexing.stills.candidate_outlier_rejection=outlier_rejection_flag
        self.params.indexing.stills.refine_all_candidates=refine_all_candidates_flag
    #
        experiments, indexed = self.index_with_known_orientation(experiments, observed)
        return experiments,indexed
    return
    #else:
    #  experiments, indexed = self.index_with_iota(experiments, observed)
    #  return experiments,indexed
#
  def index_with_iota(self, experiments, reflections):
    ''' Copy of index method from DialsFrameLogging class in xfel/ui/frame'''
    experiments, indexed = super(DialsProcessorWithLogging, self).index(experiments, reflections)
    return experiments, indexed

  def index_with_known_orientation(self, experiments, reflections):
    ''' Copy of the index function from stills_process to force IOTA to use stills_indexer during known_orientation '''
    from dials.algorithms.indexing.stills_indexer import stills_indexer
    from time import time
    import copy
    st = time()

    imagesets = experiments.imagesets()
    params = copy.deepcopy(self.params)
    # don't do scan-varying refinement during indexing
    params.refinement.parameterisation.scan_varying = False

    if hasattr(self, 'known_crystal_models'):
      known_crystal_models = self.known_crystal_models
    else:
      known_crystal_models = None
    if params.indexing.stills.method_list is None:
      idxr = stills_indexer.from_parameters(
        reflections, imagesets, known_crystal_models=known_crystal_models,
        params=params)
      idxr.index()
    else:
      indexing_error = None
      for method in params.indexing.stills.method_list:
        params.indexing.method = method
        try:
          idxr = stills_indexer.from_parameters(
            reflections, imagesets,
            params=params)
          idxr.index()
        except Exception as e:
          logger.info("Couldn't index using method %s"%method)
          if indexing_error is None:
            if e is None:
              e = Exception("Couldn't index using method %s"%method)
            indexing_error = e
        else:
          indexing_error = None
          break
      if indexing_error is not None:
        raise indexing_error

    indexed = idxr.refined_reflections
    experiments = idxr.refined_experiments

    if known_crystal_models is not None:
      from dials.array_family import flex
      filtered = flex.reflection_table()
      for idx in set(indexed['miller_index']):
        sel = indexed['miller_index'] == idx
        if sel.count(True) == 1:
          filtered.extend(indexed.select(sel))
      print("Filtered duplicate reflections, %d out of %d remaining"%(len(filtered),len(indexed)))
      indexed = filtered

    return experiments, indexed


if __name__ == "__main__":
  # Fix to make init method work for InMemScript_iota
  import xfel.command_line.xtc_process
  xfel.command_line.xtc_process.phil_scope=phil_scope
  from dials.util import halraiser
  try:
    script = InMemScript_iota()
    script.run()
  except Exception as e:
    halraiser(e)
