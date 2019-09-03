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

LS49_phil_str='''
  LS49 {
    photon_energy_from_FEE = True
      .type = bool
      .help = read in photon energy from FEE spectrometer instead of ebeam.
    fee_calibration_file = None
      .type = str
      .help = fee calibration file that will be used to read off slope/intercept values \
              to map fee pixels to ev axis

}
'''
from xfel.command_line.xtc_process import xtc_phil_str, extra_dials_phil_str
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
from xfel.command_line.xtc_process import EventOffsetSerializer, InMemScript
from exafel_project.ADSE13_25.command_line.stills_process_filter_spots import Script_iota, Processor_iota

class InMemScript_LS49_RANSAC(InMemScript, Processor_iota):
  '''
      Mixin class trying to use IOTA code 
      Stuff to change in this script and use InMemScript instead
      a. init func 
      b. finalize
      c. debug_write
      d. debug_start '''
  def __init__(self):
    InMemScript.__init__(self)

  def finalize(self):
    InMemScript.finalize(self)

  def debug_write(self, string, state=None):
    InMemScript.debug_write(self, string, state)

  def debug_start(self, ts):
    InMemScript.debug_start(self, ts)

  def process_event(self, run, evt):
    """
    Process a single event from a run
    @param run psana run object
    @param timestamp psana timestamp object
    """
    if PSANA2_VERSION:
      sec  = evt.seconds
      nsec = evt.nanoseconds
    else:
      time = evt.get(psana.EventId).time()
      fid = evt.get(psana.EventId).fiducials()
      sec  = time[0]
      nsec = time[1]

    ts = cspad_tbx.evt_timestamp((sec,nsec/1e6))
    if ts is None:
      print("No timestamp, skipping shot")
      return

    if len(self.params_cache.debug.event_timestamp) > 0 and ts not in self.params_cache.debug.event_timestamp:
      return
    self.run = run

    if self.params_cache.debug.skip_processed_events or self.params_cache.debug.skip_unprocessed_events or self.params_cache.debug.skip_bad_events:
      if ts in self.known_events:
        if self.known_events[ts] not in ["stop", "done", "fail"]:
          if self.params_cache.debug.skip_bad_events:
            print("Skipping event %s: possibly caused an unknown exception previously"%ts)
            return
        elif self.params_cache.debug.skip_processed_events:
          print("Skipping event %s: processed successfully previously"%ts)
          return
      else:
        if self.params_cache.debug.skip_unprocessed_events:
          print("Skipping event %s: not processed previously"%ts)
          return

    self.debug_start(ts)

    # FIXME MONA: below will be replaced with filter() callback
    if not PSANA2_VERSION:
      if evt.get("skip_event") or "skip_event" in [key.key() for key in evt.keys()]:
        print("Skipping event",ts)
        self.debug_write("psana_skip", "skip")
        return

    print("Accepted", ts)
    self.params = copy.deepcopy(self.params_cache)

    # the data needs to have already been processed and put into the event by psana
    if self.params.format.file_format == 'cbf':
      if self.params.format.cbf.mode == "cspad":
        # get numpy array, 32x185x388
        data = cspad_cbf_tbx.get_psana_corrected_data(self.psana_det, evt, use_default=False, dark=True,
                                                    common_mode=self.common_mode,
                                                    apply_gain_mask=self.params.format.cbf.cspad.gain_mask_value is not None,
                                                    gain_mask_value=self.params.format.cbf.cspad.gain_mask_value,
                                                    per_pixel_gain=self.params.format.cbf.cspad.per_pixel_gain,
                                                    additional_gain_factor=self.params.format.cbf.cspad.additional_gain_factor)
      elif self.params.format.cbf.mode == "rayonix":
        data = rayonix_tbx.get_data_from_psana_event(evt, self.params.input.address)
      if data is None:
        print("No data")
        self.debug_write("no_data", "skip")
        return
      if self.params.format.cbf.override_distance is None:
        if self.params.format.cbf.mode == "cspad":
          distance = cspad_tbx.env_distance(self.params.input.address, run.env(), self.params.format.cbf.detz_offset)
        elif self.params.format.cbf.mode == "rayonix":
          distance = self.params.format.cbf.detz_offset
        if distance is None:
          print("No distance, skipping shot")
          self.debug_write("no_distance", "skip")
          return
      else:
        distance = self.params.format.cbf.override_distance


      if self.params.format.cbf.override_energy is None and not self.params.LS49.photon_energy_from_FEE:
        if PSANA2_VERSION:
          wavelength = 12398.4187/self.psana_det.photonEnergy(evt)
        else:
          wavelength = cspad_tbx.evt_wavelength(evt)
        if wavelength is None:
          print ("No wavelength, skipping shot")
          self.debug_write("no_wavelength", "skip")
          return
      else:
        if self.params.LS49.photon_energy_from_FEE:
          from libtbx.easy_pickle import load
          assert self.params.LS49.fee_calibration_file is not None, 'Please supply FEE calibration pickle'
          FEE=load(os.path.join(self.params.LS49.fee_calibration_file))
          #FEE=load('/reg/d/psdm/mfx/mfxls4916/scratch/asmit/LS49_2019/input/FEE_r143_v0.pickle')
          photon_energy=FEE.get_photon_energy(run=self.params.input.run_num, evt=evt)
          if photon_energy is None:
            print ("No wavelength from FEE, skipping shot")
            self.debug_write("no_wavelength_FEE", "skip")
            return
          print ('LS49_PR', photon_energy)
          wavelength=12398.4187/photon_energy
        else:
          wavelength = 12398.4187/self.params.format.cbf.override_energy

    if self.params.format.file_format == 'pickle':
      image_dict = evt.get(self.params.format.pickle.out_key)
      data = image_dict['DATA']

    self.timestamp = timestamp = t = ts
    s = t[0:4] + t[5:7] + t[8:10] + t[11:13] + t[14:16] + t[17:19] + t[20:23]
    print ("Processing shot", s)

    def build_dxtbx_image():
      if self.params.format.file_format == 'cbf':
        # stitch together the header, data and metadata into the final dxtbx format object
        if self.params.format.cbf.mode == "cspad":
          dxtbx_img = cspad_cbf_tbx.format_object_from_data(self.base_dxtbx, data, distance, wavelength, timestamp, self.params.input.address, round_to_int=False)
        elif self.params.format.cbf.mode == "rayonix":
          dxtbx_img = rayonix_tbx.format_object_from_data(self.base_dxtbx, data, distance, wavelength, timestamp, self.params.input.address)

        if self.params.input.reference_geometry is not None:
          from dxtbx.model import Detector
          # copy.deep_copy(self.reference_detctor) seems unsafe based on tests. Use from_dict(to_dict()) instead.
          dxtbx_img._detector_instance = Detector.from_dict(self.reference_detector.to_dict())
          if self.params.format.cbf.mode == "cspad":
            dxtbx_img.sync_detector_to_cbf() #FIXME need a rayonix version of this??

      elif self.params.format.file_format == 'pickle':
        from dxtbx.format.FormatPYunspecifiedStill import FormatPYunspecifiedStillInMemory
        dxtbx_img = FormatPYunspecifiedStillInMemory(image_dict)
      return dxtbx_img

    dxtbx_img = build_dxtbx_image()
    for correction in self.params.format.per_pixel_absorption_correction:
      if correction.apply:
        if correction.algorithm == "fuller_kapton":
          from dials.algorithms.integration.kapton_correction import all_pixel_image_data_kapton_correction
          data = all_pixel_image_data_kapton_correction(image_data=dxtbx_img, params=correction.fuller_kapton)()
          dxtbx_img = build_dxtbx_image() # repeat as necessary to update the image pixel data and rebuild the image

    self.tag = s # used when writing integration pickle

    if self.params.dispatch.dump_all:
      self.save_image(dxtbx_img, self.params, os.path.join(self.params.output.output_dir, "shot-" + s))

    self.cache_ranges(dxtbx_img, self.params.input.override_spotfinding_trusted_min, self.params.input.override_spotfinding_trusted_max)

    from dxtbx.imageset import ImageSet, ImageSetData, MemReader
    imgset = ImageSet(ImageSetData(MemReader([dxtbx_img]), None))
    imgset.set_beam(dxtbx_img.get_beam())
    imgset.set_detector(dxtbx_img.get_detector())

    if self.params.dispatch.estimate_gain_only:
      from dials.command_line.estimate_gain import estimate_gain
      estimate_gain(imgset)
      return


    # FIXME MONA: radial avg. is currently disabled
    if not PSANA2_VERSION:
      # Two values from a radial average can be stored by mod_radial_average. If present, retrieve them here
      key_low = 'cctbx.xfel.radial_average.two_theta_low'
      key_high = 'cctbx.xfel.radial_average.two_theta_high'
      tt_low = evt.get(key_low)
      tt_high = evt.get(key_high)

    if self.params.radial_average.enable:
      if tt_low is not None or tt_high is not None:
        print("Warning, mod_radial_average is being used while also using xtc_process radial averaging. mod_radial_averaging results will not be logged to the database.")

    from dxtbx.model.experiment_list import ExperimentListFactory
    experiments = ExperimentListFactory.from_imageset_and_crystal(imgset, None)

    try:
      self.pre_process(experiments)
    except Exception as e:
      self.debug_write("preprocess_exception", "fail")
      return

    if not self.params.dispatch.find_spots:
      self.debug_write("data_loaded", "done")
      return

    # before calling DIALS for processing, set output paths according to the templates
    if not self.params.output.composite_output:
      if self.indexed_filename_template is not None and "%s" in self.indexed_filename_template:
        self.params.output.indexed_filename = os.path.join(self.params.output.output_dir, self.indexed_filename_template%("idx-" + s))
      if "%s" in self.refined_experiments_filename_template:
        self.params.output.refined_experiments_filename = os.path.join(self.params.output.output_dir, self.refined_experiments_filename_template%("idx-" + s))
      if "%s" in self.integrated_filename_template:
        self.params.output.integrated_filename = os.path.join(self.params.output.output_dir, self.integrated_filename_template%("idx-" + s))
      if "%s" in self.integrated_experiments_filename_template:
        self.params.output.integrated_experiments_filename = os.path.join(self.params.output.output_dir, self.integrated_experiments_filename_template%("idx-" + s))
      if "%s" in self.reindexedstrong_filename_template:
        self.params.output.reindexedstrong_filename = os.path.join(self.params.output.output_dir, self.reindexedstrong_filename_template%("idx-" + s))

    if self.params.input.known_orientations_folder is not None:
      expected_orientation_path = os.path.join(self.params.input.known_orientations_folder, os.path.basename(self.params.output.refined_experiments_filename))
      if os.path.exists(expected_orientation_path):
        print("Known orientation found")
        from dxtbx.model.experiment_list import ExperimentListFactory
        self.known_crystal_models = ExperimentListFactory.from_json_file(expected_orientation_path, check_format=False).crystals()
      else:
        print("Image not previously indexed, skipping.")
        self.debug_write("not_previously_indexed", "stop")
        return
    # Load a dials mask from the trusted range and psana mask
    from dials.util.masking import MaskGenerator
    generator = MaskGenerator(self.params.border_mask)
    mask = generator.generate(imgset)
    if self.params.format.file_format == "cbf" and self.dials_mask is not None:
      mask = tuple([a&b for a, b in zip(mask,self.dials_mask)])
    if self.spotfinder_mask is None:
      self.params.spotfinder.lookup.mask = mask
    else:
      self.params.spotfinder.lookup.mask = tuple([a&b for a, b in zip(mask,self.spotfinder_mask)])

    self.debug_write("spotfind_start")
    try:
      observed = self.find_spots(experiments)
    except Exception as e:
      import traceback; traceback.print_exc()
      print(str(e), "event", timestamp)
      self.debug_write("spotfinding_exception", "fail")
      return

    print("Found %d bright spots"%len(observed))

    if self.params.dispatch.hit_finder.enable and len(observed) < self.params.dispatch.hit_finder.minimum_number_of_reflections:
      print("Not enough spots to index")
      self.debug_write("not_enough_spots_%d"%len(observed), "stop")
      return
    if self.params.dispatch.hit_finder.maximum_number_of_reflections is not None:
      if self.params.dispatch.hit_finder.enable and len(observed) > self.params.dispatch.hit_finder.maximum_number_of_reflections:
        print("Too many spots to index - Possibly junk")
        self.debug_write("too_many_spots_%d"%len(observed), "stop")
        return

    self.restore_ranges(dxtbx_img)
    # save cbf file
    if self.params.dispatch.dump_strong:
      self.save_image(dxtbx_img, self.params, os.path.join(self.params.output.output_dir, "hit-" + s))

      # save strong reflections.  self.find_spots() would have done this, but we only
      # want to save data if it is enough to try and index it
      if self.strong_filename_template:
        if "%s" in self.strong_filename_template:
          strong_filename = self.strong_filename_template%("hit-" + s)
        else:
          strong_filename = self.strong_filename_template
        strong_filename = os.path.join(self.params.output.output_dir, strong_filename)

        from dials.util.command_line import Command
        Command.start('Saving {0} reflections to {1}'.format(
            len(observed), os.path.basename(strong_filename)))
        observed.as_pickle(strong_filename)
        Command.end('Saved {0} observed to {1}'.format(
            len(observed), os.path.basename(strong_filename)))

    if not self.params.dispatch.index:
      self.debug_write("strong_shot_%d"%len(observed), "done")
      return

    # index and refine
    self.debug_write("index_start")
    try:
      experiments, indexed = self.index(experiments, observed)
    except Exception as e:
      import traceback; traceback.print_exc()
      print(str(e), "event", timestamp)
      self.debug_write("indexing_failed_%d"%len(observed), "stop")
      return

    if self.params.dispatch.dump_indexed:
      img_path = self.save_image(dxtbx_img, self.params, os.path.join(self.params.output.output_dir, "idx-" + s))
      imgset = ExperimentListFactory.from_filenames([img_path]).imagesets()[0]
      assert len(experiments.detectors()) == 1;   imgset.set_detector(experiments[0].detector)
      assert len(experiments.beams()) == 1;       imgset.set_beam(experiments[0].beam)
      assert len(experiments.scans()) <= 1;       imgset.set_scan(experiments[0].scan)
      assert len(experiments.goniometers()) <= 1; imgset.set_goniometer(experiments[0].goniometer)
      for expt_id, expt in enumerate(experiments):
        expt.imageset = imgset

    self.debug_write("refine_start")

    try:
      experiments, indexed = self.refine(experiments, indexed)
    except Exception as e:
      import traceback; traceback.print_exc()
      print(str(e), "event", timestamp)
      self.debug_write("refine_failed_%d"%len(indexed), "fail")
      return

    if self.params.dispatch.reindex_strong:
      self.debug_write("reindex_start")
      try:
        self.reindex_strong(experiments, observed)
      except Exception as e:
        import traceback; traceback.print_exc()
        print(str(e), "event", timestamp)
        self.debug_write("reindexstrong_failed_%d"%len(indexed), "fail")
        return

    if not self.params.dispatch.integrate:
      self.debug_write("index_ok_%d"%len(indexed), "done")
      return

    # integrate
    self.debug_write("integrate_start")
    self.cache_ranges(dxtbx_img, self.params.input.override_integration_trusted_min, self.params.input.override_integration_trusted_max)

    if self.cached_ranges is not None:
      # Load a dials mask from the trusted range and psana mask
      imgset = ImageSet(ImageSetData(MemReader([dxtbx_img]), None))
      imgset.set_beam(dxtbx_img.get_beam())
      imgset.set_detector(dxtbx_img.get_detector())
      from dials.util.masking import MaskGenerator
      generator = MaskGenerator(self.params.border_mask)
      mask = generator.generate(imgset)
      if self.params.format.file_format == "cbf" and self.dials_mask is not None:
        mask = tuple([a&b for a, b in zip(mask,self.dials_mask)])
    if self.integration_mask is None:
      self.params.integration.lookup.mask = mask
    else:
      self.params.integration.lookup.mask = tuple([a&b for a, b in zip(mask,self.integration_mask)])

    try:
      integrated = self.integrate(experiments, indexed)
    except Exception as e:
      import traceback; traceback.print_exc()
      print(str(e), "event", timestamp)
      self.debug_write("integrate_failed_%d"%len(indexed), "fail")
      return
    self.restore_ranges(dxtbx_img)

    self.debug_write("integrate_ok_%d"%len(integrated), "done")


if __name__ == "__main__":
  # Fix to make init method work for InMemScript_iota
  import xfel.command_line.xtc_process
  xfel.command_line.xtc_process.phil_scope=phil_scope
  
  from dials.util import halraiser
  try:
    script = InMemScript_LS49_RANSAC()
    script.run()
  except Exception as e:
    halraiser(e)
