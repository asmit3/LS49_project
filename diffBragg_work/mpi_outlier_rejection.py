#
# MPI version of load_ls49
#
from __future__ import absolute_import, print_function, division 
from libtbx.utils import Abort, Sorry
from libtbx.phil import parse
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
import sys, os

from LS49_project.diffBragg_work.outlier_rejection import outlier_rejection_ls49

help_message = ''' MPI version of load_ls49 script to distribute outlier rejection with diffBragg jobs over several cores/nodes using MPI'''
usage = "usage: libtbx.python script_name params.phil"


outlier_diffBragg_phil_str=''' 
  outlier_diffBragg {
    timestamps_file = None
      .type = str
      .help = Path to filename that contains all timestamps to be processed
    rayonix_expt_path = None
      .type = str
      .help = Path to folder that contains all rayonix expt/refl files
    output_dir = None
      .type = str
      .help = Path to output directory 
    dump_outlier_output_files = False
      .type = bool
      .help = Determine if outlier rejection filtered integration pickles/json be dumped
}
'''
phil_scope = parse(outlier_diffBragg_phil_str)


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
    if self.params.outlier_diffBragg.output_dir is not None:
      self.output_dir = self.params.outlier_diffBragg.output_dir
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
      with open(self.params.outlier_diffBragg.timestamps_file, 'r') as fin:
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
    import time

    ts = item_list[0] 
    ls49_data_dir=self.params.outlier_diffBragg.rayonix_expt_path
    print ('Inside do_work for rank %d'%rank)
    t_start = time.time()
    experiments_path=os.path.join(ls49_data_dir,'idx-%s_integrated.expt'%ts)
    reflections_path=os.path.join(ls49_data_dir,'idx-%s_integrated.refl'%ts)
    experiments=ExperimentListFactory.from_json_file(experiments_path, check_format=False)[0:1]
    reflections=flex.reflection_table.from_file(reflections_path)
    reflections=reflections.select(reflections['id']==0)
    if ls49_data_dir is None:
      LS49_regression = libtbx.env.find_in_repositories(
        relative_path="LS49_regression",
        test=os.path.isdir)
      if LS49_regression is None:
        raise Sorry('LS49_regression folder needs to be present or else specify ls49_data_dir')
      ls49_data_dir = os.path.join(LS49_regression, 'diffBragg_work', 'iota_r0222_cori', 'rayonix_expt')
    outlier_rejection_ls49(experiments, reflections,ls49_data_dir=ls49_data_dir, ts=ts, show_plotted_images=False, outdir=self.params.outlier_diffBragg.output_dir, dump_output_files=self.params.outlier_diffBragg.dump_outlier_output_files)
    t_end = time.time()
    delta_time = t_end - t_start
    print ('Outlier_DiffBragg_timing %s  = %d'%(ts,delta_time))

if __name__ == "__main__":
  try:
    script = Script()
    script.run()
  except Exception as e:
    print (str(e))
