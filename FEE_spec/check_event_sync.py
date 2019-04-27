from __future__ import absolute_import, division, print_function
from psana import *
from Detector.EpicsDetector import EpicsDetector
from Detector.AreaDetector import AreaDetector

message=""" script to check whether fee_spec, fee_gas, ebeam, rayonix and jungfrau1M events are in sync for LS49
            It was noticed that fee_spec data is present for many events but no fee_gas/rayonix etc. Subsequently
            it was deduced that those were background signals in the fee_spec (magnitude ~ 1e5). Compare that with 
            events which had signal (i.e ebeam/fee_gas/rayonix) which had magnitude on fee_spec ~ 1e7
            
            **This function only uses psana functions so should run outside CCTBX**"""

# Specify all psana related interfaces here
ds = DataSource('exp=mfxls4916:run=223:smd')
rayonix = Detector('MfxEndstation.0:Rayonix.0')
env = ds.env()
fee_spec = AreaDetector('FEE_Spec', env)
fee_gas_det = Detector('FEEGasDetEnergy')
jungfrau1M = Detector('MfxEndstation.0:Jungfrau.0')
ebeam = Detector('EBeam')

# variables to keep track of
evt_count_with_wavelength=0
fee_spec_evt_count=0
fee_gas_evt_count=0
evt_with_rayonix =0
evt_with_jungfrau1M =0
total_events=0

for nevent,evt in enumerate(ds.events()):
  if nevent > 50: break
  total_events +=1
  ev_ebeam=ebeam.get(evt)
  print ('NEVENT = ',nevent)
  fee_gas = fee_gas_det.get(evt)
  calib_array = fee_spec.calib(evt)
  rayonix_img = rayonix.calib(evt)
  jungfrau1M_img = jungfrau1M.calib(evt)
  if ev_ebeam is not None:
    print ('event with ebeam',nevent)
    evt_count_with_wavelength +=1
  if calib_array is not None:
    #if nevent > 40 and fee_gas is not None:
    #  from IPython import embed; embed(); exit()
    fee_spec_evt_count +=1
    print ('fee spec event',nevent)
  if fee_gas is not None:
    fee_gas_evt_count +=1
    print ('FEE gas event',nevent)
  if rayonix_img is not None:
    print ('Rayonix event',nevent)
    evt_with_rayonix +=1
  if jungfrau1M_img is not None:
    evt_with_jungfrau1M +=1
    print ('jungfrau1M event',nevent)
  print ('-----------------------------------------------')

print ('Total events = ',total_events )
print ('Event count with wavelengths=',evt_count_with_wavelength)
print ('fee_spec events = ', fee_spec_evt_count)
print ('fee gas events = ', fee_gas_evt_count)
print ('Number of rayonix events = ', evt_with_rayonix)
print ('Number of Jungfrau1M events = ', evt_with_jungfrau1M)
