from __future__ import print_function, division, absolute_import
import glob
import os
from dials.array_family import flex
from scitbx.matrix import col
from libtbx.easy_pickle import load
import math
from scitbx.math import five_number_summary

num_pixels=1.0
all_idx_ts = {}
for folder in glob.glob('r0*/00*_rg002/out_reint'):
  #from IPython import embed; embed(); exit()
  for f in os.listdir(folder):
    if 'indexed.refl' not in f: continue
    ts = f[4:21]
    fpickle=os.path.join(folder, f)
    reflections = load(fpickle)
    refl_now = reflections.select(reflections['id'] >= 0)
    count =0
    strong = [] # dummy array 
    #for refl in refl_now.rows():
    #    dR.append((col(refl['xyzcal.mm']) - col(refl['xyzobs.mm.value'])).length())
    #    count +=1
    #print ('RMSD(um) = ', filename,1000.0*math.sqrt(dR.dot(dR)/len(dR)), len(dR), len(strong))
    dR = refl_now['xyzcal.mm'] - refl_now['xyzobs.mm.value']
    #rmsd=int(1000.0*math.sqrt(dR.dot(dR)/len(dR)))
    #from IPython import embed; embed()
    rmsd = 1000.0*dR.norm()/len(dR)
    path_name = os.path.abspath(folder)
    if rmsd > 178.0*num_pixels: continue
    information=(path_name, rmsd, len(dR), ts)
    if ts not in all_idx_ts:
      all_idx_ts[ts]=[]
    all_idx_ts[ts].append(information)

all_rmsd=flex.double()
total_rmsd=0.0
total_count=0
print ('Done collecting info: Now select top rankers')
fout=open('timestamps_to_merge.dat', 'w')
for ts in all_idx_ts:
  #print (ts)
  ts_info = all_idx_ts[ts]
  num_of_spots = flex.int()
  rmsd_of_spots = flex.double()
  min_rmsd=999999.9
  ranking = []
  for info in all_idx_ts[ts]:
    #print (info)
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
  print (stuff)
  print('\n')
  min_rmsd = rmsd_of_spots[top_ranker]
  #print ('top ranker = ', stuff)
  int_file_0 = os.path.join(stuff[0], 'int-0-%s.pickle'%stuff[-1])
  int_file_1 = os.path.join(stuff[0], 'int-1-%s.pickle'%stuff[-1])
  if os.path.exists(int_file_0):
    fout.write(int_file_0+'\n')
  if os.path.exists(int_file_1):
    fout.write(int_file_1+'\n')
  
  
  all_rmsd.append(float(min_rmsd))
  #total_rmsd+=min_rmsd*min_rmsd*count
  total_count+=stuff[2]
  #print ('--------------------------------------')
print ('NUmber of images indexed  =', len(all_idx_ts))
print ('Five number summary of best case scenario RMSDs (um)= ', five_number_summary(all_rmsd))
#print ('Mean RMSD over all spots in dataset indexed (um) = ', math.sqrt(total_rmsd/total_count))
