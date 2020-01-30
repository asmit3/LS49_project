from __future__ import print_function, division, absolute_import
import sys, os, copy
from libtbx.phil import parse
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory, ExperimentList
from dials.algorithms.spot_prediction import StillsReflectionPredictor
from scitbx.matrix import col

phil_scope = parse(''' 
  grid_search_jungfrau {
    spotfinding_directory = None
      .type = str
      .help = Provide directory with results from spotfinding i.e strong.refl and imported.expt
    rayonix_expt_directory = None
      .type = str
      .help = Provide directory with results from rayonix indexing. Will need integrated.expt
    jungfrau_geometry_file = None
      .type = str
      .help = path to jungfrau geometry file
    minimum_number_of_spots = 3
      .type = int
      .help = minimum number of spots found on jungfrau to pursue further processing
    mosaicity_Deff = 1093.0
      .type = float
      .help = Set crystal Deff value if provided
    mosaicity_eta_half = 0.0001
      .type = float
      .help = Set crystal eta half value if provided
    grid_search_bins_x = 10
      .type = int
      .help = number of increments in x-direction on either side of mean value i.e +/- grid_search_bins_x*dx_mm 
    grid_search_bins_y = 10
      .type = int
      .help = number of increments in y-direction on either side of mean value i.e +/- grid_search_bins_y*dy_mm 
    grid_search_bins_z = 10
      .type = int
      .help = number of increments in z-direction on either side of mean value i.e +/- grid_search_bins_z*dz_mm 
    dx_mm = 0.177
      .type = float
      .help = increment length in grid search x-direction
    dy_mm = 0.177
      .type = float
      .help = increment length in grid search y-direction
    dz_mm = 0.177
      .type = float
      .help = increment length in grid search z-direction
    use_mpi = False
      .type = bool
      .help = Whether to use MPI or not (striping mode)
    outdir = .
      .type = str
      .help = specify directory to dump output
    dump_intermediate_files = False
      .type = bool
      .help = flag to indicate whether you want to dump intermediate grid search files aka \
              refl/expt files from each prediction of each perturbed detector distance
    critical_pixel_dist = 10
      .type = float
      .help = critical pixel distance used to evaluate whether prediction is close enough to a strong spot \
}
''')

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

class GridSearch_Jungfrau(object):
  ''' Class to do a grid search of the jungfrau detector model
       Reads in a bunch of imported.expt and strong.refl files
       Filters if enough spots are there to do the grid search
       Then predicts spots based on rayonix crystal model for each point on grid'''

  def __init__(self):
    # Read phil stuff
    self.params = params_from_phil(sys.argv[1:])
    self.trial_result = [] # a list to store trial result evaluation

  def setup_grid(self):
    dx_mm = self.params.grid_search_jungfrau.dx_mm 
    dy_mm = self.params.grid_search_jungfrau.dy_mm 
    dz_mm = self.params.grid_search_jungfrau.dz_mm 
    self.dx = [0.0]
    self.dy = [0.0]
    self.dz = [0.0]
    for ix in range(1, self.params.grid_search_jungfrau.grid_search_bins_x+1):
      self.dx.append(ix*dx_mm)
      self.dx.append(-ix*dx_mm)
    for iy in range(1, self.params.grid_search_jungfrau.grid_search_bins_y+1):
      self.dy.append(iy*dy_mm)
      self.dy.append(-iy*dy_mm)
    for iz in range(1, self.params.grid_search_jungfrau.grid_search_bins_z+1):
      self.dz.append(iz*dz_mm)
      self.dz.append(-iz*dz_mm)

  def do_work(self, rank, iterable):
    for item in iterable:
      ts, strong_refl_f, imported_expt_f, integrated_rayonix_expt_f = item 

      strong_refl = flex.reflection_table.from_file(strong_refl_f)
      if len(strong_refl) <self.params.grid_search_jungfrau.minimum_number_of_spots: continue 
      imported_expt = ExperimentListFactory.from_json_file(imported_expt_f)
      integrated_rayonix_expt = ExperimentListFactory.from_json_file(integrated_rayonix_expt_f, check_format=False)

      imported_expt[0].detector = self.jungfrau_geom[0].detector
      imported_expt[0].crystal = integrated_rayonix_expt[0].crystal
      imported_expt[0].crystal.set_domain_size_ang(self.params.grid_search_jungfrau.mosaicity_Deff)
      imported_expt[0].crystal.set_half_mosaicity_deg(self.params.grid_search_jungfrau.mosaicity_eta_half)

      # This is the fun part aka grid search
      original_detector = copy.deepcopy(imported_expt[0].detector)
      d0 = original_detector.hierarchy()
      f0 = col(d0.get_local_fast_axis())
      s0 = col(d0.get_local_slow_axis())
      r0 = col(d0.get_local_origin())
      for dx in self.dx:
        print ('Stage dx')
        for dy in self.dy:
          for dz in self.dz:
            imported_expt[0].detector=original_detector
            displacement = col((dx,dy,dz))
            det0 = imported_expt[0].detector.hierarchy()
            df = col((f0[0]*displacement[0], f0[1]*displacement[1], f0[2]*displacement[2]))
            ds = col((s0[0]*displacement[0], s0[1]*displacement[1], s0[2]*displacement[2]))
            dorigin = col((0,0,dz)) 
            r_new = r0+df+ds+dorigin
            det0.set_local_frame(det0.get_fast_axis(), det0.get_slow_axis(), r_new)
            
            predictor=StillsReflectionPredictor(imported_expt[0])
            ubx=predictor.for_ub(imported_expt[0].crystal.get_A())
            ubx['id']=flex.int(len(ubx),0)
            #print ('NUMBER OF PREDS=', len(ubx), ts, dx,dy,dz)
#            if len(ubx) == len(strong_refl):
#              print ('NUMBER OF PREDS=', len(ubx), ts, dx,dy,dz)
            self.evaluate_trial(strong_refl, ubx, imported_expt[0].detector)
            if self.params.grid_search_jungfrau.dump_intermediate_files:
              imported_expt.as_file(os.path.join(self.params.grid_search_jungfrau.outdir, 'prediction_%s_%.2f_%.2f_%.2f.expt'%(ts,dx,dy,dz)))
              ubx.as_file(os.path.join(self.params.grid_search_jungfrau.outdir, 'prediction_%s_%.2f_%.2f_%.2f.refl'%(ts,dx,dy,dz)))
      from IPython import embed; embed(); exit()

  def evaluate_grid_search(self):
    pass

  def evaluate_trial(self, strong_refl, predicted_refl, detector):
    pixel_size = detector[0].get_pixel_size()[0]
    count = 0
    for ii, strong in enumerate(strong_refl.rows()):
      for jj, predicted in enumerate(predicted_refl.rows()):
        pstrong = detector[strong['panel']]
        ppredicted = detector[predicted['panel']]
        xyobs_px=strong['xyzobs.px.value'][0:2]
        xypred_px = predicted['xyzcal.px'][0:2]
        r_strong = col(pstrong.get_pixel_lab_coord(xyobs_px))
        r_predicted = col(ppredicted.get_pixel_lab_coord(xypred_px))
        dpx = (r_strong - r_predicted).length()/pixel_size
        if dpx < self.params.grid_search_jungfrau.critical_pixel_dist:
          count +=1
    self.trial_result.append(count)
    

  def run(self):
    # Figure out what data to process
    # First read in the spotfinding directory and create a list of timestamps and folders to read in
    # Set up the iterable variable which will be distributed to all ranks in case of MPI
    iterable = []
    self.jungfrau_geom = ExperimentListFactory.from_json_file(self.params.grid_search_jungfrau.jungfrau_geometry_file, check_format=False)
    for f in os.listdir(self.params.grid_search_jungfrau.spotfinding_directory):
      if 'strong.refl' in f:
        strong_refl_f = os.path.join(self.params.grid_search_jungfrau.spotfinding_directory, f)
        ts = f[16:16+17]  
        
        imported_expt_f = os.path.join(self.params.grid_search_jungfrau.spotfinding_directory, 'idx-jungfrauhit_%s_imported.expt'%ts)
        # Load up the rayonix integrated model
        integrated_rayonix_expt_f = os.path.join(self.params.grid_search_jungfrau.rayonix_expt_directory,'idx-%s_integrated.expt'%ts)
        #print (ts, strong_refl, imported_expt, integrated_rayonix_expt)
        # Make sure these files exist before proceeding
        if not os.path.exists(imported_expt_f) or not os.path.exists(integrated_rayonix_expt_f): continue
        iterable.append((ts, strong_refl_f, imported_expt_f, integrated_rayonix_expt_f))

    # Set up the grid 
    self.setup_grid()
     
    # Process the data
    if self.params.grid_search_jungfrau.use_mpi:
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
      size = comm.Get_size()  # size: number of processes running in this job

        # Configure the logging
      if params.output.logging_dir is None:
        logfile = None
      else:
        log_path = os.path.join(
            params.output.logging_dir, "log_rank%04d.out" % rank)
        error_path = os.path.join(
            params.output.logging_dir, "error_rank%04d.out" % rank)
        print("Redirecting stdout to %s" % log_path)
        print("Redirecting stderr to %s" % error_path)
        sys.stdout = open(log_path, "a", buffering=0)
        sys.stderr = open(error_path, "a", buffering=0)
        print("Should be redirected now")

        if True: # Just using the striping mode for MPI 
          subset = [item for i, item in enumerate(iterable) if (i + rank) % size == 0]
          self.do_work(rank, subset)
    else:
      self.do_work(0, iterable)

if __name__ == '__main__':
  try:
    script = GridSearch_Jungfrau()
    script.run()
  except Exception as e:
    print (str(e))
