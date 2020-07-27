#
# MPI version of load_ls49
#
from __future__ import absolute_import, print_function, division 
from libtbx.utils import Abort, Sorry
from libtbx.phil import parse
import sys, os

from LS49_project.diffBragg_work.mcmc_JF1M import strong_spot_mask, process_ls49_image_real, run_all_refine_ls49_JF1M

from dials.command_line.stills_process import control_phil_str, dials_phil_str, program_defaults_phil_str

help_message = ''' MPI version of load_ls49 script to distribute diffBragg jobs over several cores/nodes using MPI'''
usage = "usage: libtbx.python script_name params.phil"


LS49_diffBragg_phil_str=''' 
  LS49_diffBragg {
    timestamps_file = None
      .type = str
      .help = Path to filename that contains all timestamps to be processed
    ls49_data_dir = None
      .type = str
      .help = Path to folder that contains all jungfrau expt/refl files
    output_dir = None
      .type = str
      .help = Path to output directory 
    plot = False
      .type = bool
      .help = Flag to indicate if plotting will be used.
    Deff = 1000.0
      .type = float
      .help = Mosaic domain size initial estimate, will be refined eventually
    max_seed = 0
      .type = int
      .help = maximum seed value to do mcmc choosing fp/fdp displacement using random number with seed between 0 and max_seed  
    short_circuit_dir = None
      .type = str
      .help = Directory for RUC_info_timestamp.pickle file
    swap_spectra_timestamp = False
      .type = bool
      .help = If true, swap spectra of timestamp with another timestamp as read in from a pre-specified file
}
'''
#phil_scope = parse(LS49_diffBragg_phil_str)
phil_scope = parse(LS49_diffBragg_phil_str+control_phil_str + dials_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))

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

      # Define iterable here
      iterable = []
      seeds = range(1, self.params.LS49_diffBragg.max_seed+1) 
      for seed in seeds:
        for ts in all_timestamps:
          iterable.append((ts, seed))
      #iterable = all_timestamps
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
    import time
    ts, seed = item_list[0][0], item_list[0][1]
    print ('Send TYT', ts,seed)
    print ('Inside do_work for rank %d'%rank)
    t_start = time.time()
    final_likelihood = run_all_refine_ls49_JF1M(ts=ts, ls49_data_dir=self.params.LS49_diffBragg.ls49_data_dir, show_plotted_images=False, outdir=self.params.LS49_diffBragg.output_dir, params=self.params, seed=seed, short_circuit_dir=self.params.LS49_diffBragg.short_circuit_dir, swap_spectra_timestamp=self.params.LS49_diffBragg.swap_spectra_timestamp)
    t_end = time.time()
    delta_time = t_end - t_start
    print ('DiffBragg_LS49_timing %s  = %d'%(ts,delta_time))
    print ('Diff_Scattering_LogLikelihood = %d, %s, %.3f'%(seed, ts, final_likelihood))

if __name__ == "__main__":
  try:
    script = Script()
    script.run()
  except Exception as e:
    print (str(e))
