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
from LS49_project.SAD.xtc_process_sad import InMemScript_LS49_SAD

class InMemScript_LS49_RANSAC(InMemScript_LS49_SAD, Processor_iota):
  '''
      Mixin class trying to use IOTA code 
      Stuff to change in this script and use InMemScript instead
      1. init func 
      2. finalize
      3. debug_write
      4. debug_start '''
  def __init__(self):
    InMemScript.__init__(self)

  def finalize(self):
    InMemScript.finalize(self)

  def debug_write(self, string, state=None):
    InMemScript.debug_write(self, string, state)

  def debug_start(self, ts):
    InMemScript.debug_start(self, ts)

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
