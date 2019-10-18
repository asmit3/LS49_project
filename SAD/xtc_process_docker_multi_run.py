from __future__ import absolute_import, division, print_function
from six.moves import range
from six.moves import zip
# -*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.xtc_process
#
PSANA2_VERSION = 0
try:
  import psana
  PSANA2_VERSION = psana.__version__
except ImportError:
  pass # for running at home without psdm build
except AttributeError:
  pass

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

from LS49_project.SAD.xtc_process_sad import LS49_phil_str
from xfel.command_line.xtc_process import xtc_phil_str, extra_dials_phil_str
# Replace the int run_num with str. Fail program if not present
original_str='''run_num = None\n      .type = int\n      .help = Run number or run range to process\n'''
new_str='''run_num = None\n      .type = str\n      .help = Run number or run range to process2\n'''
if original_str not in xtc_phil_str: raise Sorry("Cannot replace xtc_phil_str: Please check code carefully")
xtc_phil_str=xtc_phil_str.replace(original_str, new_str)
from dials.command_line.stills_process import dials_phil_str, program_defaults_phil_str
from xfel.ui import db_phil_str
from xfel.command_line.xfel_process import radial_average_phil_str
from exafel_project.ADSE13_25.command_line.stills_process_filter_spots import iota_phil_str
phil_scope = parse(LS49_phil_str + iota_phil_str + xtc_phil_str + dials_phil_str + extra_dials_phil_str + db_phil_str + radial_average_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))
# Function imports
from xfel.command_line.xtc_process import filter, run_psana2
# Class imports
#from xfel.command_line.xfel_process import Script as DialsProcessScript
#from xfel.ui.db.frame_logging import DialsProcessorWithLogging
from xfel.command_line.xtc_process import EventOffsetSerializer
from exafel_project.ADSE13_25.command_line.stills_process_filter_spots import Script_iota, Processor_iota
from LS49_project.SAD.xtc_process_sad_multiple_runs import InMemScript_LS49_SAD_multiple_runs

class InMemScript_LS49_RANSAC_multiple_runs(InMemScript_LS49_SAD_multiple_runs, Processor_iota):
  '''
      Mixin class trying to use IOTA code 
      Stuff to change in this script and use InMemScript instead
      1. init func 
      2. finalize
      3. debug_write
      4. debug_start '''
  def __init__(self):
    InMemScript_LS49_SAD_multiple_runs.__init__(self)

  def finalize(self):
    InMemScript_LS49_SAD_multiple_runs.finalize(self)

  def debug_write(self, string, state=None):
    InMemScript_LS49_SAD_multiple_runs.debug_write(self, string, state)

  def debug_start(self, ts):
    InMemScript_LS49_SAD_multiple_runs.debug_start(self, ts)

  def run(self):
    """ Process all images assigned to this thread """

    try:
      params, options = self.parser.parse_args(
        show_diff_phil=True, quick_parse=True)
    except Exception as e:
      if "Unknown command line parameter definition" in str(e) or \
          "The following definitions were not recognised" in str(e):
        deprecated_params = ['mask_nonbonded_pixels','gain_mask_value','algorithm','custom_parameterization']
        deprecated_strs = ['%s','%s','common_mode.%s','common_mode.%s']
        for i in range(len(deprecated_params)):
          if deprecated_params[i] in str(e):
            print("format.cbf.%s"%(deprecated_strs[i]%deprecated_params[i]), "has changed to format.cbf.cspad.%s"%(deprecated_strs[i]%deprecated_params[i]))
      raise

    # Check inputs
    if params.input.experiment is None or \
       params.input.run_num is None or \
       (params.input.address is None and params.format.file_format != 'pickle'):
      raise Usage(self.usage)

    if params.format.file_format == "cbf":
      if params.format.cbf.detz_offset is None:
        raise Usage(self.usage)
    elif params.format.file_format == "pickle":
      if params.input.cfg is None:
        raise Usage(self.usage)
    else:
      raise Usage(self.usage)

    if not os.path.exists(params.output.output_dir):
      raise Sorry("Output path not found:" + params.output.output_dir)

    if params.format.file_format == "cbf":
      if params.output.tmp_output_dir == "(NONE)":
        tmp_dir = params.output.tmp_output_dir
      else:
        #Environment variable redirect for CBFLib temporary CBF_TMP_XYZ file output
        if params.output.tmp_output_dir is None:
          tmp_dir = os.path.join(params.output.output_dir, '.tmp')
        else:
          tmp_dir = os.path.join(params.output.tmp_output_dir, '.tmp')
        if not os.path.exists(tmp_dir):
          try:
            os.makedirs(tmp_dir)
          except Exception as e:
            # Can fail if running multiprocessed, which is ok if the tmp folder was created
            if not os.path.exists(tmp_dir):
              halraiser(e)
      os.environ['CBF_TMP_DIR'] = tmp_dir

    for abs_params in params.integration.absorption_correction:
      if abs_params.apply and abs_params.algorithm == "fuller_kapton":
        if not (params.integration.debug.output and not params.integration.debug.separate_files):
          raise Sorry('Shoeboxes must be saved to integration intermediates to apply an absorption correction. '\
            +'Set integration.debug.output=True and integration.debug.separate_files=False to save shoeboxes.')

    self.params = params
    self.load_reference_geometry()

    if params.output.composite_output:
      from dxtbx.model.experiment_list import ExperimentList
      from dials.array_family import flex
      #self.all_strong_reflections = flex.reflection_table() # no composite strong pickles yet
      self.all_indexed_experiments = ExperimentList()
      self.all_indexed_reflections = flex.reflection_table()
      self.all_integrated_experiments = ExperimentList()
      self.all_integrated_reflections = flex.reflection_table()
    else:
      # The convention is to put %s in the phil parameter to add a time stamp to
      # each output datafile. Save the initial templates here.
      self.strong_filename_template                 = params.output.strong_filename
      self.indexed_filename_template                = params.output.indexed_filename
      self.refined_experiments_filename_template    = params.output.refined_experiments_filename
      self.integrated_filename_template             = params.output.integrated_filename
      self.integrated_experiments_filename_template = params.output.integrated_experiments_filename
      self.reindexedstrong_filename_template        = params.output.reindexedstrong_filename

    # Don't allow the strong reflections to be written unless there are enough to
    # process
    params.output.strong_filename = None

    # Save the paramters
    self.params_cache = copy.deepcopy(params)
    self.options = options

    if params.mp.method == "mpi":
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
      size = comm.Get_size() # size: number of processes running in this job
    elif params.mp.method == "sge" and \
        'SGE_TASK_ID'    in os.environ and \
        'SGE_TASK_FIRST' in os.environ and \
        'SGE_TASK_LAST'  in os.environ:
      if 'SGE_STEP_SIZE' in os.environ:
        assert int(os.environ['SGE_STEP_SIZE']) == 1
      if os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined' or os.environ['SGE_TASK_ID'] == 'undefined':
        rank = 0
        size = 1
      else:
        rank = int(os.environ['SGE_TASK_ID']) - int(os.environ['SGE_TASK_FIRST'])
        size = int(os.environ['SGE_TASK_LAST']) - int(os.environ['SGE_TASK_FIRST']) + 1
    else:
      rank = 0
      size = 1
    self.composite_tag = "%04d"%rank

    # Configure the logging
    if params.output.logging_dir is None:
      info_path = ''
      debug_path = ''
    elif params.output.logging_dir == 'DevNull':
      print("Redirecting stdout, stderr and other DIALS logfiles to /dev/null")
      sys.stdout = open(os.devnull,'w', buffering=0)
      sys.stderr = open(os.devnull,'w',buffering=0)
      info_path = os.devnull
      debug_path = os.devnull
    else:
      log_path = os.path.join(params.output.logging_dir, "log_rank%04d.out"%rank)
      error_path = os.path.join(params.output.logging_dir, "error_rank%04d.out"%rank)
      print("Redirecting stdout to %s"%log_path)
      print("Redirecting stderr to %s"%error_path)
      sys.stdout = open(log_path,'a', buffering=0)
      sys.stderr = open(error_path,'a',buffering=0)
      print("Should be redirected now")

      info_path = os.path.join(params.output.logging_dir, "info_rank%04d.out"%rank)
      debug_path = os.path.join(params.output.logging_dir, "debug_rank%04d.out"%rank)

    from dials.util import log
    log.config(params.verbosity, info=info_path, debug=debug_path)

    debug_dir = os.path.join(params.output.output_dir, "debug")
    if not os.path.exists(debug_dir):
      try:
        os.makedirs(debug_dir)
      except OSError as e:
        pass # due to multiprocessing, makedirs can sometimes fail
    assert os.path.exists(debug_dir)

    if params.debug.skip_processed_events or params.debug.skip_unprocessed_events or params.debug.skip_bad_events:
      print("Reading debug files...")
      self.known_events = {}
      for filename in os.listdir(debug_dir):
        # format: hostname,timestamp_event,timestamp_now,status,detail
        for line in open(os.path.join(debug_dir, filename)):
          vals = line.strip().split(',')
          if len(vals) != 5:
            continue
          _, ts, _, status, detail = vals
          if status in ["done", "stop", "fail"]:
            self.known_events[ts] = status
          else:
            self.known_events[ts] = "unknown"

    self.debug_file_path = os.path.join(debug_dir, "debug_%d.txt"%rank)
    write_newline = os.path.exists(self.debug_file_path)
    if write_newline: # needed if the there was a crash
      self.debug_write("")

    if params.mp.method != 'mpi' or params.mp.mpi.method == 'client_server':
      if rank == 0:
        self.mpi_log_file_path = os.path.join(debug_dir, "mpilog.out")
        write_newline = os.path.exists(self.mpi_log_file_path)
        if write_newline: # needed if the there was a crash
          self.mpi_log_write("\n")

    # FIXME MONA: psana 2 has pedestals and geometry hardcoded for cxid9114.
    # We can remove after return code when all interfaces are ready.
    if PSANA2_VERSION:
        print(("PSANA2_VERSION", PSANA2_VERSION))
        run_psana2(self, params, comm)
        return

    # set up psana
    if params.input.cfg is not None:
      psana.setConfigFile(params.input.cfg)
    # all cores in stripe mode and the master in client-server mode read smd
    if params.dispatch.datasource is None:
      datasource = "exp=%s:run=%s:%s"%(params.input.experiment,params.input.run_num,'smd')
      if params.input.xtc_dir is not None:
        if params.input.use_ffb:
          raise Sorry("Cannot specify the xtc_dir and use SLAC's ffb system")
        datasource += ":dir=%s"%params.input.xtc_dir
      elif params.input.use_ffb:
      # as ffb is only at SLAC, ok to hardcode /reg/d here
        datasource += ":dir=/reg/d/ffb/%s/%s/xtc"%(params.input.experiment[0:3],params.input.experiment)
      if params.input.stream is not None and len(params.input.stream) > 0:
        datasource += ":stream=%s"%(",".join(["%d"%stream for stream in params.input.stream]))
      if params.input.calib_dir is not None:
        psana.setOption('psana.calib-dir',params.input.calib_dir)
      if params.mp.method == "mpi" and params.mp.mpi.method == 'client_server' and size > 2:
        dataset_name_client = datasource.replace(":smd",":rax")
      # for client-server, master reads smd - clients read rax
        if rank == 0:
          ds = psana.DataSource(datasource)
        else:
          ds = psana.DataSource(dataset_name_client)

      else:
      # for stripe, all cores read smd
        ds = psana.DataSource(datasource)
    else:
      datasource = params.dispatch.datasource
      ds = psana.DataSource(datasource)

    if params.format.file_format == "cbf":
      self.psana_det = psana.Detector(params.input.address, ds.env())

    # set this to sys.maxint to analyze all events
    if params.dispatch.max_events is None:
      max_events = sys.maxsize
    else:
      max_events = params.dispatch.max_events


    ##########################################
    # Get the total number of runs to process by parsing the run_num phil parameter
    # Note that the format hopefully will be as given here
    # https://confluence.slac.stanford.edu/display/PSDM/Manual#Manual-Datasetspecification
    run_count = 0
    all_runs = []
    xrun_num = params.input.run_num.strip().split(',')
    for x in xrun_num:
      z = x.split('-')
      if len(z) == 1:
        all_runs.append(int(z[0])) # can't extend a list with an integer !
      elif len(z) == 2:
        all_runs.extend(range(int(z[0]), int(z[1])+1))
      else:
        print ('cant extend run(s)', str(x))
    total_runs = len(list(set(all_runs)))
    print ('TOTAL_RUNS',total_runs)
    #########################################
    # Set up db connection if being used
    if self.params.experiment_tag is not None:
      self.db_app = dxtbx_xfel_db_application(params)

    for run in ds.runs():
      if params.format.file_format == "cbf":
        if params.format.cbf.mode == "cspad":
          # load a header only cspad cbf from the slac metrology
          try:
            self.base_dxtbx = cspad_cbf_tbx.env_dxtbx_from_slac_metrology(run, params.input.address)
          except Exception as e:
            raise Sorry("Couldn't load calibration file for run %d, %s"%(run.run(), str(e)))
        elif params.format.cbf.mode == "rayonix":
          # load a header only rayonix cbf from the input parameters
          detector_size = rayonix_tbx.get_rayonix_detector_dimensions(ds.env())
          self.base_dxtbx = rayonix_tbx.get_dxtbx_from_params(params.format.cbf.rayonix, detector_size)

        if self.base_dxtbx is None:
          raise Sorry("Couldn't load calibration file for run %d"%run.run())

        if params.format.file_format == 'cbf':
          if params.format.cbf.cspad.common_mode.algorithm == "custom":
            self.common_mode = params.format.cbf.cspad.common_mode.custom_parameterization
            assert self.common_mode is not None
          else:
            self.common_mode = params.format.cbf.cspad.common_mode.algorithm # could be None or default

        if params.format.cbf.invalid_pixel_mask is not None:
          self.dials_mask = easy_pickle.load(params.format.cbf.invalid_pixel_mask)
          if params.format.cbf.mode == "cspad":
            assert len(self.dials_mask) == 64
            if self.params.format.cbf.cspad.mask_nonbonded_pixels:
              psana_mask = self.psana_det.mask(run.run(),calib=False,status=False,edges=False,central=False,unbond=True,unbondnbrs=True)
              dials_mask = self.psana_mask_to_dials_mask(psana_mask)
              self.dials_mask = [self.dials_mask[i] & dials_mask[i] for i in range(len(dials_mask))]
        else:
          if params.format.cbf.mode == "cspad":
            psana_mask = self.psana_det.mask(run.run(),calib=True,status=True,edges=True,central=True,unbond=True,unbondnbrs=True)
            self.dials_mask = self.psana_mask_to_dials_mask(psana_mask)
          else:
            self.dials_mask = None
      if self.params.spotfinder.lookup.mask is not None:
        self.spotfinder_mask = easy_pickle.load(self.params.spotfinder.lookup.mask)
      else:
        self.spotfinder_mask = None
      if self.params.integration.lookup.mask is not None:
        self.integration_mask = easy_pickle.load(self.params.integration.lookup.mask)
      else:
        self.integration_mask = None

      # prepare fractions of process_percent, if given
      process_fractions = None
      if params.dispatch.process_percent:
        import fractions
        percent = params.dispatch.process_percent / 100
        process_fractions = fractions.Fraction(percent).limit_denominator(100)
      # list of all events
      # only cycle through times in client_server mode
      if params.mp.mpi.method != 'client_server' and total_runs > 1:
        raise Sorry("This script currently is supported only for client-server MPI.")
      if params.mp.method == "mpi" and params.mp.mpi.method == 'client_server' and size > 2:
        run_count +=1
        print ('RUN_COUNT_COMBO = ',run_count,' ', run.run())
        # process fractions only works in idx-striping mode
        if params.dispatch.process_percent:
          raise Sorry("Process percent only works in striping mode.")
        print("Using MPI client server in rax mode")

        # use a client/server approach to be sure every process is busy as much as possible
        # only do this if there are more than 2 processes, as one process will be a server
        try:
          if rank == 0:
            # server process
            self.mpi_log_write("MPI START\n")
            for nevt, evt in enumerate(run.events()):
              if nevt == max_events: break
              self.mpi_log_write("Getting next available process\n")
              offset = evt.get(psana.EventOffset)
              rankreq = comm.recv(source=MPI.ANY_SOURCE)
              t = evt.get(psana.EventId).time()
              ts = cspad_tbx.evt_timestamp((t[0],t[1]/1e6))
              self.mpi_log_write("Process %s is ready, sending ts %s\n"%(rankreq, ts))
              comm.send(EventOffsetSerializer(offset),dest=rankreq)
            # send a stop command to each process
            if run_count == total_runs:
              self.mpi_log_write("MPI DONE, sending stops\n")
              for rankreq in range(size-1):
                self.mpi_log_write("Getting next available process\n")
                rankreq = comm.recv(source=MPI.ANY_SOURCE)
                self.mpi_log_write("Sending stop to %d\n"%rankreq)
                comm.send('endrun',dest=rankreq)
              self.mpi_log_write("All stops sent.")
            # rank 0 has finished iterating through this run: tell other folks to move on who are waiting
            else:
              self.mpi_log_write("Rank_0 done with run %d :: Sending transition signal\n"%run.run())
              for rankreq in range(size-1):
                self.mpi_log_write("Getting next available process\n")
                rankreq = comm.recv(source=MPI.ANY_SOURCE)
                self.mpi_log_write("Sending transition signal to %d\n"%rankreq)
                comm.send('transitionrun',dest=rankreq)
              self.mpi_log_write("All transitions sent.\n")
          else:
            # client process
            while True:
              # inform the server this process is ready for an event
              print("Rank %d getting next task"%rank)
              comm.send(rank,dest=0)
              print("Rank %d waiting for response"%rank)
              offset = comm.recv(source=0)
              if offset == 'endrun':
                print("Rank %d recieved endrun"%rank)
                break
              if offset == 'transitionrun':
                print("Rank %d recieved transitionrun"%rank)
                break
              evt = ds.jump(offset.filenames, offset.offsets, offset.lastBeginCalibCycleDgram)
              print("Rank %d beginning processing"%rank)
              try:
                self.process_event(run, evt)
              except Exception as e:
                print("Rank %d unhandled exception processing event"%rank, str(e))
              print("Rank %d event processed"%rank)
        except Exception as e:
          print("Error caught in main loop")
          print(str(e))
        print("Synchronizing rank %d"%rank)
        #comm.Barrier()
        print("Rank %d done with main loop"%rank)
      else:
        import resource
        # chop the list into pieces, depending on rank.  This assigns each process
        # events such that the get every Nth event where N is the number of processes
        print("Striping events")

        nevent = mem = first = last = 0
        if process_fractions:
          def process_this_event(nevent):
            # nevent modulo the denominator gives us which fraction we're in
            n_mod_denom = nevent % process_fractions.denominator
            # compare the 0-indexed modulo against the 1-indexed numerator (intentionally not <=)
            n_accept = n_mod_denom < process_fractions.numerator
            return n_accept
        for nevent, evt in enumerate(run.events()):
          if nevent%size != rank: continue
          if nevent >= max_events: break
          if process_fractions and not process_this_event(nevent): continue

          self.process_event(run, evt)

          mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
          if nevent < 50:
            #print "Mem test rank %03d"%rank, i, mem
            continue
          #print "Mem test rank %03d"%rank, 'Cycle %6d total %7dkB increase %4dkB' % (i, mem, mem - last)
          if not first:
            first = mem
          last = mem
        print('Total memory leaked in %d cycles: %dkB' % (nevent+1-50, mem - first))

    comm.Barrier()
    print("Rank %d finalizing"%rank)
    try:
      self.finalize()
    except Exception as e:
      print("Rank %d, exception caught in finalize"%rank)
      print(str(e))

    if params.format.file_format == "cbf" and params.output.tmp_output_dir == "(NONE)":
      try:
        os.rmdir(tmp_dir)
      except Exception as e:
        pass

    if params.joint_reintegration.enable:
      if params.output.composite_output:
        raise NotImplementedError("Joint reintegration not implemented for composite output yet")
      assert self.params.dispatch.dump_indexed, "Cannot do joint reintegration unless indexed files were dumped"
      if rank == 0:
        reint_dir = os.path.join(params.output.output_dir, "reint")
        if not os.path.exists(reint_dir):
          os.makedirs(reint_dir)
        images = []
        experiment_jsons = []
        indexed_tables = []
        for filename in os.listdir(params.output.output_dir):
          if not filename.endswith("_indexed.refl"):
            continue
          experiment_jsons.append(os.path.join(params.output.output_dir, filename.split("_indexed.refl")[0] + "_refined.expt"))
          indexed_tables.append(os.path.join(params.output.output_dir, filename))
          if params.format.file_format == "cbf":
            images.append(os.path.join(params.output.output_dir, filename.split("_indexed.refl")[0] + ".cbf"))
          elif params.format.file_format == "pickle":
            images.append(os.path.join(params.output.output_dir, filename.split("_indexed.refl")[0] + ".pickle"))

        if len(images) < params.joint_reintegration.minimum_results:
          pass # print and return

        # TODO: maximum_results_per_chunk = 500
        combo_input = os.path.join(reint_dir, "input.phil")
        f = open(combo_input, 'w')
        for json, indexed in zip(experiment_jsons, indexed_tables):
          f.write("input {\n")
          f.write("  experiments = %s\n"%json)
          f.write("  reflections = %s\n"%indexed)
          f.write("}\n")
        f.close()

        combined_experiments_file = os.path.join(reint_dir, "combined.expt")
        combined_reflections_file = os.path.join(reint_dir, "combined.refl")
        command = "dials.combine_experiments reference_from_experiment.average_detector=True %s output.reflections=%s output.experiments=%s"% \
          (combo_input, combined_reflections_file, combined_experiments_file)
        print(command)
        from libtbx import easy_run
        easy_run.fully_buffered(command).raise_if_errors().show_stdout()

        from dxtbx.model.experiment_list import ExperimentListFactory

        combined_experiments = ExperimentListFactory.from_json_file(combined_experiments_file, check_format=False)
        combined_reflections = flex.reflection_table.from_file(combined_reflections_file)

        from dials.algorithms.refinement import RefinerFactory

        refiner = RefinerFactory.from_parameters_data_experiments(
          params.joint_reintegration, combined_reflections, combined_experiments)

        refiner.run()
        experiments = refiner.get_experiments()
        reflections = combined_reflections.select(refiner.selection_used_for_refinement())

        from dxtbx.model.experiment_list import ExperimentListDumper
        from dxtbx.model import ExperimentList
        dump = ExperimentListDumper(experiments)
        dump.as_json(os.path.join(reint_dir, "refined.expt"))
        reflections.as_pickle(os.path.join(reint_dir, "refined.refl"))

        for expt_id, (expt, img_file) in enumerate(zip(experiments, images)):
          try:
            refls = reflections.select(reflections['id'] == expt_id)
            refls['id'] = flex.int(len(refls), 0)
            base_name = os.path.splitext(os.path.basename(img_file))[0]
            self.params.output.integrated_filename = os.path.join(reint_dir, base_name + "_integrated.refl")

            expts = ExperimentList([expt])
            self.integrate(expts, refls)
            dump = ExperimentListDumper(expts)
            dump.as_json(os.path.join(reint_dir, base_name + "_refined.expt"))
          except Exception as e:
            print("Couldn't reintegrate", img_file, str(e))
    print("Rank %d signing off"%rank)

if __name__ == "__main__":
  # Fix to make init method work for InMemScript_iota
  import xfel.command_line.xtc_process
  xfel.command_line.xtc_process.phil_scope=phil_scope
  
  from dials.util import halraiser
  try:
    script = InMemScript_LS49_RANSAC_multiple_runs()
    script.run()
  except Exception as e:
    halraiser(e)
