#
# MPI version of load_ls49
#
from __future__ import absolute_import, print_function, division 
from libtbx.utils import Abort, Sorry
from libtbx.phil import parse
import sys, os

from LS49_project.diffBragg_work.load_ls49 import strong_spot_mask, process_ls49_image_real

help_message = ''' MPI version of load_ls49 script to distribute diffBragg jobs over several cores/nodes using MPI'''
usage = "usage: libtbx.python script_name params.phil"


LS49_diffBragg_phil_str=''' 
  LS49_diffBragg {
    timestamps_file = None
      .type = str
      .help = Path to filename that contains all timestamps to be processed
    rayonix_expt_path = None
      .type = str
      .help = Path to folder that contains all rayonix expt/refl files
    output_dir = None
      .type = str
      .help = Path to output directory 
    plot = False
      .type = bool
      .help = Flag to indicate if plotting will be used.
    Deff = 1000.0
      .type = float
      .help = Mosaic domain size initial estimate, will be refined eventually
}
'''
phil_scope = parse(LS49_diffBragg_phil_str)


def params_from_phil(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse(file_name=arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  return params

class Script(object):
  def __init__(self):
    """Initialise the script."""
    from dials.util.options import OptionParser
    import libtbx.load_env
    self.params=params_from_phil(sys.argv[1:])
    if self.params.LS49_diffBragg.output_dir is not None:
      self.output_dir = self.params.LS49_diffBragg.output_dir
    else:
      self.output_dir = './'

  def run(self):
    """Execute the script."""
    from dials.util import log
    from time import time
    from libtbx import easy_mp
    import copy

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Set up logger for each rank
    log_path = os.path.join(self.output_dir, "log_rank%04d.out"%rank)
    error_path = os.path.join(self.output_dir, "error_rank%04d.out"%rank)
    print("Redirecting stdout to %s"%log_path)
    print("Redirecting stderr to %s"%error_path)
    sys.stdout = open(log_path,'a', buffering=0)
    sys.stderr = open(error_path,'a',buffering=0)
    print("Should be redirected now")

    # Let rank 0 read in the timestamps to process and distribute the work
    if rank == 0:
      all_timestamps = []
      with open(self.params.LS49_diffBragg.timestamps_file, 'r') as fin:
        for line in fin:
          if line !='/n':
            all_timestamps.append(line.strip())
      print ('Total number of timestamps = %d'%len(all_timestamps))
      iterable = all_timestamps
    comm.Barrier()
    if rank == 0:
      # server process
      for item in iterable:
        print("Getting next available process")
        rankreq = comm.recv(source=MPI.ANY_SOURCE)
        print("Process %s is ready, sending %s\n" % (rankreq, item))
        comm.send(item, dest=rankreq)
        # send a stop command to each process
      print("MPI DONE, sending stops\n")
      for rankreq in range(size - 1):
        rankreq = comm.recv(source=MPI.ANY_SOURCE)
        print("Sending stop to %d\n" % rankreq)
        comm.send("endrun", dest=rankreq)
      print("All stops sent.")
    else:
      # client process
      while True:
      # inform the server this process is ready for an event
        print("Rank %d getting next task" % rank)
        comm.send(rank, dest=0)
        print("Rank %d waiting for response" % rank)
        item = comm.recv(source=0)
        if item == "endrun":
          print("Rank %d received endrun" % rank)
          break
        print("Rank %d beginning processing" % rank)
        try:
          self.do_work(rank, [item])
        except Exception as e:
          print("Rank %d unhandled exception processing event" % rank,str(e))
        print("Rank %d event processed" % rank)

  def do_work(self, rank, item_list=None):
    from simtbx.diffBragg.refiners import RefineAll
    from simtbx.diffBragg.sim_data import SimData
    import numpy as np
    from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
    from simtbx.diffBragg import nanoBragg_crystal, nanoBragg_beam
    from copy import deepcopy
    from dxtbx.model.experiment_list import Experiment, ExperimentList, ExperimentListFactory

    ts = item_list[0] 
    data = process_ls49_image_real(tstamp=ts,Nstrongest=10, resmin=2.0, resmax=13.5, ls49_data_dir=self.params.LS49_diffBragg.rayonix_expt_path)
 

    C = data["dxcrystal"]
    D = data["dxdetector"]
    B = data["dxbeam"]
    exp=data['experiment']
    indexed_reflections = deepcopy(data['indexed_reflections'])
    dump_exp = Experiment(imageset=data['cbf_imageset'], 
                         beam=B,
                         detector=D,
                         goniometer=exp.goniometer,
                         scan=exp.scan,
                         crystal=C)
    dump_explist = ExperimentList([dump_exp])
    dump_explist.as_file(os.path.join(self.output_dir, 'before_refinement_%s.expt'%ts))
    indexed_reflections.as_file(os.path.join(self.output_dir, 'before_refinement_%s.refl'%ts))

    # Some global variables here for LS49
    mos_spread_deg=0.01
    n_mos_domains=1
    a,b,c,_,_,_ = C.get_unit_cell().parameters()
    Deff = self.params.LS49_diffBragg.Deff
    Ncells_abc = np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
    nbcryst = nanoBragg_crystal.nanoBragg_crystal()
    nbcryst.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
    nbcryst.mos_spread_deg = mos_spread_deg
    nbcryst.n_mos_domains = n_mos_domains
    nbcryst.thick_mm = 0.005
    nbcryst.miller_array = data["sfall"]
    nbcryst.dxtbx_crystal = C

    nbbeam = nanoBragg_beam.nanoBragg_beam()
    nbbeam.unit_s0 = B.get_unit_s0()
    nbbeam.spectrum = data["spectrum"]

    SIM = SimData()
    SIM.crystal = nbcryst
    SIM.detector = D
    SIM.beam = nbbeam

    SIM.instantiate_diffBragg(adc_offset=0,
                              oversample=0,
                              interpolate=0,
                              verbose=0)
    SIM.D.show_params()

    ucell = C.get_unit_cell().parameters()
    UcellMan = MonoclinicManager(a=ucell[0], b=ucell[1], c=ucell[2], beta=ucell[4]*np.pi / 180.)

    RUC = RefineAll(
        spot_rois=data["bboxes_x1x2y1y2"],
        spot_resolution=data['resolution'],
        abc_init=data["tilt_abc"],
        img=data["data_img"],
        SimData_instance=SIM,
        plot_images=False,
        ucell_manager=UcellMan,
        init_gain=1,
        init_scale=1.0)
    RUC.trad_conv = True
    RUC.trad_conv_eps = 1e-5
    RUC.max_calls = 250
    RUC.refine_background_planes = False
    RUC.refine_Umatrix = False
    RUC.refine_Bmatrix = False
    RUC.verbose=True
    RUC.plot_stride=1
    RUC.refine_background_planes = False
    RUC.refine_ncells=False #
    RUC.refine_crystal_scale = False # scale factor
    RUC.refine_gain_fac = False
    RUC.refine_detdist=False
    RUC.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
    RUC.use_curvatures=False #
    RUC.calc_curvatures=True # if set to False, never uses curvatures
    RUC.refine_with_restraints=True
    #RUC.run()
    #if RUC.hit_break_to_use_curvatures:
    #  RUC.num_positive_curvatures=0
    #  RUC.use_curvatures=True
    #  RUC.run(setup=False) # Now it will use curvatures

    # First refine scale factors and ncells, then refine everything else
    RUC.refine_ncells=True
    RUC.refine_crystal_scale=True
    RUC.run()
    refined_ncells = RUC.x[7]
    refined_scale = RUC.x[10]

    # Now set the other stuff to refine
    print (' ================== now refining all the variables ======================= ')
    Ncells_abc = refined_ncells #np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
    nbcryst = nanoBragg_crystal.nanoBragg_crystal()
    nbcryst.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
    nbcryst.mos_spread_deg = mos_spread_deg
    nbcryst.n_mos_domains = n_mos_domains
    nbcryst.thick_mm = 0.005
    nbcryst.miller_array = data["sfall"]
    nbcryst.dxtbx_crystal = C

    nbbeam = nanoBragg_beam.nanoBragg_beam()
    nbbeam.unit_s0 = B.get_unit_s0()
    nbbeam.spectrum = data["spectrum"]

    SIM2 = SimData()
    SIM2.crystal = nbcryst
    SIM2.detector = D
    SIM2.beam = nbbeam

    SIM2.instantiate_diffBragg(adc_offset=0,
                              oversample=0,
                              interpolate=0,
                              verbose=0)
    SIM2.D.show_params()
    RUC2 = RefineAll(
        spot_rois=data["bboxes_x1x2y1y2"],
        spot_resolution=data['resolution'],
        abc_init=data["tilt_abc"],
        img=data["data_img"],
        SimData_instance=SIM2,
        plot_images=False,
        ucell_manager=UcellMan,
        init_gain=1,
        init_scale=refined_scale)
    RUC2.trad_conv = True
    RUC2.trad_conv_eps = 1e-5
    RUC2.max_calls = 250
    RUC2.refine_background_planes = False
    RUC2.refine_Umatrix = True
    RUC2.refine_Bmatrix = True
    RUC2.verbose=True
    RUC2.plot_stride=1
    RUC2.refine_background_planes = False
    RUC2.refine_ncells=True #
    RUC2.refine_crystal_scale = True # scale factor
    RUC2.refine_gain_fac = False
    RUC2.refine_detdist=False
    RUC2.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
    RUC2.use_curvatures=False #
    RUC2.calc_curvatures=True # if set to False, never uses curvatures
    RUC2.run()

    print("Done.")
    print("Refined scale =%f", RUC.x[-1])

    best=RUC2.best_image
    C2 = deepcopy(C)
    ang, ax = RUC2.get_correction_misset(as_axis_angle_deg=True)
    C2.rotate_around_origin(ax, ang)
    C2.set_B(RUC2.get_refined_Bmatrix())
      

    # refined unit cell parameters
    ucell_ref = C2.get_unit_cell().parameters()

    print("ground truth unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell)
    print("refined unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell_ref)
    print("")

    C2.show()
    dump_exp = Experiment(imageset=data['cbf_imageset'], 
                         beam=B,
                         detector=D,
                         goniometer=exp.goniometer,
                         scan=exp.scan,
                         crystal=C2)
    dump_explist = ExperimentList([dump_exp])
    dump_explist.as_file(os.path.join(self.output_dir, 'after_refinement_%s.expt'%ts))
    # Dump refl file as well based on prediction from refined model
    if True:
      from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
      ref_predictor = ExperimentsPredictorFactory.from_experiments(
                      dump_explist,
                      force_stills=True,
                      spherical_relp=False)
      ref_predictor(indexed_reflections)
      indexed_reflections.as_file(os.path.join(self.output_dir, 'after_refinement_%s.refl'%ts))

if __name__ == "__main__":
  try:
    script = Script()
    script.run()
  except Exception as e:
    print (str(e))
