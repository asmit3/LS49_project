import os, glob

message = ''' Figure out which timestamps were left unprocessed at the end of the KNL jobs for IOTA. Will be used to process again using the CBFs''' 

def psana_to_cctbx_ts_converter(psana_ts):
  cctbx_ts = psana_ts[0:4]+psana_ts[5:7]+psana_ts[8:10] +psana_ts[11:13]+psana_ts[14:16]+psana_ts[17:19]+psana_ts[20:23]
  return cctbx_ts

trial_rg = ['001_rg002']

runs = range(175, 255+1)
crucial_text = 'Sending stop to'
unfinished_business = [] 
cbf_to_process = [] # we care about this at the end


for trial_rg_folder in trial_rg:
  for run in runs:
    ranks_finished = []   
    ranks_unfinished = []   
    run_folder = 'r%04d'%run
    debug_folder=os.path.join(run_folder, trial_rg_folder, 'out', 'debug')
    # First find number of ranks by figuring out number of debug files
    nranks = len(glob.glob1(debug_folder, 'debug_*.txt'))
    #print ('NRANKS = %d'%nranks)
    # Let us first go through the mpilog.out and scrape ranks that returned stuff ?
    mpilog_f=os.path.join(debug_folder, 'mpilog.out')
    if not os.path.exists(mpilog_f): continue
    with open(mpilog_f, 'r') as fin:
      for line in fin:
        if line !='\n' and crucial_text in line:
          ranks_finished.append(int(line.split(crucial_text)[-1]))
    for rank in range(1,nranks+1):
      if rank not in ranks_finished:
        ranks_unfinished.append(rank)
    # Now go through the last line of those ranks and figure out the time stamp
    #print (ranks_unfinished)
    # Dumb way of going through entire file and seeking last line
    for rank in ranks_unfinished:
      debug_f = os.path.join(debug_folder, 'debug_%d.txt'%rank)
      with open(debug_f, 'r') as fin:
        for line in fin:
          if line !='\n':
            last_line = line
      unfinished_business.append(last_line.strip())
for business in unfinished_business:
  psana_ts = business.split(',')[1]
  cctbx_ts = psana_to_cctbx_ts_converter(psana_ts)
  cbf_to_process.append('hit-%s.cbf'%cctbx_ts)
  #print (cctbx_ts, psana_ts)
print (cbf_to_process)
from IPython import embed; embed(); exit()
