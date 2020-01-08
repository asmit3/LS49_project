from __future__ import print_function, division, absolute_import
import glob
import os
from dials.array_family import flex
from scitbx.matrix import col
from libtbx.easy_pickle import load
import math
from scitbx.math import five_number_summary
from dxtbx.model.experiment_list import ExperimentListFactory
from cctbx import uctbx
import sys

message = ''' This is a script used to create the files in iota_LS49_regression_r0222. Following are the files you will be needing in that folder.
              Not all the files are copied by this script. In fact this script just writes out a text file with full paths of the following files
              with *. The remaining have to be manually copied over/created
              (1) All indexed refl files*
              (2) All refined expt files*
              (3) All integrated refl files*
              (4) All integrated expt files*
              (5) All CBF files*
              (6) mask_r4.pickle file
              (7) fee pickle file for all the frames. Needs to be done using LS49_project/FEE_spec/fee_spectra.py
              (8) MTZ file after initial merging.''' 



num_pixels=2.0
num_spots_cutoff=-1
all_idx_ts = {}
photon_energy_min = -1
photon_energy_max = 9999.0
wave_conv_factor=12398.4187

reference_cell=uctbx.unit_cell((63.6,28.8,35.6,90,106.5,90))
runs = [int(sys.argv[1])]#range(175, 255+1)
trial_rg = ['000_rg002', '001_rg002', '002_rg002', '003_rg002', '004_rg002', '005_rg002', '006_rg002', '007_rg002', '008_rg002', '009_rg002']
for run in runs:
  print ('run=',run)
  for trial in trial_rg:
    print ('trial = ',trial)
    folder = os.path.join('r%04d'%run, trial, 'out_reint')
    out_folder = os.path.join('r%04d'%run, trial, 'out')
    if not os.path.exists(folder): continue
    for f in os.listdir(folder):
      if 'indexed.refl' not in f: continue
      ts = f[4:21]
      tmp_folder=folder
      fpickle=os.path.join(folder, f)
      integrated_file = os.path.join(folder, f.replace('indexed', 'integrated')) 
      integrated_file_out = os.path.join(out_folder, f.replace('indexed', 'integrated')) 
      # Make sure unit cells are OK. No point considering it if they will be thrown out by cxi.merge
      integrated_expt_f = os.path.join(folder, f.replace('indexed.refl', 'integrated.expt'))
      integrated_expt = ExperimentListFactory.from_json_file(integrated_expt_f, check_format=False)[0] # considering only 0 case right now
      photon_energy = wave_conv_factor/integrated_expt.beam.get_wavelength()
      if photon_energy < photon_energy_min: continue
      if photon_energy > photon_energy_max: continue
      unit_cell=integrated_expt.crystal.get_unit_cell()
      if not unit_cell.is_similar_to(other=reference_cell,relative_length_tolerance=0.1, absolute_angle_tolerance=2.0 ): continue
      if not os.path.exists(integrated_file): continue
      #
      reflections = load(fpickle)
      refl_now = reflections.select(reflections['id'] == 0) # careful specific to 0
      if len(refl_now) < num_spots_cutoff: continue
      integrated_refl=flex.reflection_table.from_file(integrated_file)
      num_of_int_spots = len(integrated_refl.select(integrated_refl['id']==0)) # careful specific to 0
      # Load the out folder integrated pickle stuff to see which one to consider
      if os.path.exists(integrated_file_out):
        integrated_refl_out=flex.reflection_table.from_file(integrated_file_out)
        num_of_int_spots_out = len(integrated_refl_out.select(integrated_refl_out['id']==0)) # careful specific to 0
        if num_of_int_spots_out > num_of_int_spots:
          #print ('FILE DECISION',num_of_int_spots_out, num_of_int_spots, ts)
          tmp_folder=out_folder
      count =0
      strong = [] # dummy array 
      dR = refl_now['xyzcal.px'] - refl_now['xyzobs.px.value']
      rmsd = dR.norm()/math.sqrt(len(dR))
      path_name = os.path.abspath(tmp_folder)
      if rmsd > num_pixels: continue
      information=(path_name, rmsd, len(dR), num_of_int_spots, ts)
      if ts not in all_idx_ts:
        all_idx_ts[ts]=[]
      all_idx_ts[ts].append(information)

all_rmsd=flex.double()
total_rmsd=0.0
indexed_spots_count = []
integrated_spots_count = []
dump_stuff=[]
print ('Done collecting info: Now select top rankers')
fout=open('timestamps_for_LS49_regression_%s.dat'%run, 'w')
for ts in all_idx_ts:
  #print (ts)
  ts_info = all_idx_ts[ts]
  num_of_spots = flex.int()
  rmsd_of_spots = flex.double()
  min_rmsd=999999.9
  ranking = []
  for info in all_idx_ts[ts]:
    if ts=='20180501102049119': print (info)
    num_of_spots.append(info[2])
    rmsd_of_spots.append(info[1])
    #if info[1] < min_rmsd:
    #  min_rmsd=info[1]


  perm_rmsd = flex.sort_permutation(rmsd_of_spots, reverse=False)
  perm_spots = flex.sort_permutation(num_of_spots, reverse=True)
  for i in range(len(num_of_spots)):
    idx_spots = list(perm_spots).index(i)
    idx_rmsd = list(perm_rmsd).index(i)
    ranking.append(idx_spots+idx_rmsd)

  top_ranker = ranking.index(min(ranking))
  stuff = ts_info[top_ranker]
  #print (stuff)
  #print('\n')
  min_rmsd = rmsd_of_spots[top_ranker]
  #print ('top ranker = ', stuff)
  dump_stuff.append(stuff)
  # careful only doing 0 id here
  int_file_0 = os.path.join(stuff[0], 'int-0-%s.pickle'%stuff[-1])
  integrated_expt=os.path.join(stuff[0], 'idx-%s_integrated.expt'%stuff[-1])
  integrated_refl=os.path.join(stuff[0], 'idx-%s_integrated.refl'%stuff[-1])
  refined_expt=os.path.join(stuff[0], 'idx-%s_refined.expt'%stuff[-1])
  indexed_refl=os.path.join(stuff[0], 'idx-%s_indexed.refl'%stuff[-1])
  if 'out_reint' in stuff[0]:
    new_stuff=stuff[0].replace('out_reint', 'out')
    cbf_file=os.path.join(new_stuff, 'idx-%s.cbf'%stuff[-1])
  else:
    cbf_file=os.path.join(stuff[0], 'idx-%s.cbf'%stuff[-1])
  
  #int_file_1 = os.path.join(stuff[0], 'int-1-%s.pickle'%stuff[-1])
  if os.path.exists(int_file_0):
    fout.write(integrated_expt+'\n') 
    fout.write(integrated_refl+'\n') 
    fout.write(refined_expt+'\n') 
    fout.write(indexed_refl+'\n') 
    fout.write(cbf_file+'\n')
    print (stuff[-1])
    #fout.write(int_file_0+'\n')
  #if os.path.exists(int_file_1):
  #  fout.write(int_file_1+'\n')
  
  
  all_rmsd.append(float(min_rmsd))
  indexed_spots_count.append(stuff[2])
  integrated_spots_count.append(stuff[3])
  #print ('--------------------------------------')
fout.close()
print ('NUmber of images indexed  =', len(all_idx_ts))
print ('Five number summary of best case scenario RMSDs (px)= ', five_number_summary(all_rmsd))
from libtbx.easy_pickle import dump
#dump('choose_refl_stats_%s.pickle'%run,dump_stuff)
#print ('Mean RMSD over all spots in dataset indexed (um) = ', math.sqrt(total_rmsd/total_count))
