from __future__ import print_function, division, absolute_import

from psana import *
from Detector.EpicsDetector import EpicsDetector
from Detector.AreaDetector import AreaDetector
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import itertools
import math
from scipy.optimize import curve_fit
import h5py
from scipy.interpolate import CubicSpline

message="""Script to process, calibrate and correct FEE spectra data for LCLS experiment LS49
           Work done by Franklin Fuller (SLAC) and Asmit Bhowmick (LBNL)"""

#############################################################################################
#											    #
#                       Define all helper functions at the top				    #	
#											    #
#############################################################################################

def imshow_robust(i,min_p = 1, max_p = 99):
  """often times events do not have any data in it. This function is to ensure those events are properly handled"""
  vmin, vmax = np.percentile(np.ravel(i),[min_p, max_p])
  figure()
  imshow(i, vmin=vmin, vmax=vmax)

def robust_mean(x):
  return x[np.invert(np.isnan(x))].mean()

def gaus(x,a,x0,sigma):
  return a*np.exp(-(x-x0)**2/(2*sigma**2))

def fwhm(X,Y):
  half_max = Y.max() / 2.
  #find when function crosses line half_max (when sign of diff flips)
  #take the 'derivative' of signum(half_max - Y[])
  d = Y - half_max
  indexes = np.where(d > 0)[0]
  return abs(X[indexes[-1]] - X[indexes[0]])

def robust_mean_asmit(x):
  if len(x)==0: return None
  if np.isnan(x).any(): return None
  return x[np.invert(np.isnan(x))].mean()

def all_gaus(x, all_params):
  y = np.zeros(np.array(x).shape)
  for params in all_params:
    y += gaus(x,*params)
  return y

def find_nearest_idx(nparray, value):
  """courtesy Stack Overflow
  https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array"""
  idx = (np.abs(nparray - value)).argmin()
  array_value=nparray[idx]
  return array_value,idx

def flag_func(xp, a,comparison='less'):
  """Function compares each element of list x with value a and returns a boolean array
       of 1/0 depending on whether comparison is satisfied"""
  from scitbx.array_family import flex
  if comparison=='less':
    return flex.double([int(xp-a<0.0) for xp in x])
  elif comparison=='greater':
    return flex.double([int(xp-a>0.0) for xp in x])
  else:
    return flex.double([0]*len(xp))

def exgaussian(x, a, mu, sigma, tau):
  import scipy.special as sp
  erfc_input=[(mu+(sigma*sigma/tau)-xp)/(math.sqrt(2)*sigma) for xp in x]
  erfc=sp.erfc(erfc_input)
  arg_2=[(2*mu+(sigma*sigma/tau)-2*xp)/(2*tau) for xp in x]
  arg_3=np.exp(arg_2)
  return a*arg_3*erfc/(2*tau)

def plateau_func_1(x,A, mu, sigma):
  """ See https://stats.stackexchange.com/questions/203629/is-there-a-plateau-shaped-distribution"""
  import scipy.special as sp
  erf_input_1=[(xp+mu)/(math.sqrt(2)*sigma) for xp in x]
  erf_input_2=[(xp-mu)/(math.sqrt(2)*sigma) for xp in x]
  erf_1=sp.erf(erf_input_1)
  erf_2=sp.erf(erf_input_2)
  return A*(erf_1-erf_2)/(4*mu)

def plateau_func_2(x, A, mu, sigma, beta):
  """This is the generalized normal subbotin distribution
  beta=2.0 is the normal distribution
  See https://stats.stackexchange.com/questions/203629/is-there-a-plateau-shaped-distribution"""
  from scitbx.math import gamma_complete
  from scitbx.array_family import flex
  xp=flex.double(x)
  return 1.0+A*beta*flex.exp(-(flex.pow(flex.abs(xp-mu)/sigma,beta)))/(2*sigma*gamma_complete(1./beta))

def plateau_func_3(x, A, mu, sigma, w):
  from scitbx.array_family import flex
  xp=flex.double(x)
  h=1./(1+w/(math.sqrt(2*math.pi)*sigma))
  pre_exp=A*h/(math.sqrt(2*math.pi)*sigma)
  term_1=pre_exp*flex.exp(-(xp-mu+w/2)*(xp-mu+w/2)/(2*sigma*sigma))*flag_func(xp,mu-w/2,comparison='less')
  term_2=pre_exp*flag_func(xp,mu-w/2,comparison='less')*flag_func(xp,mu+w/2,comparison='greater')
  term_3=pre_exp*flex.exp(-(xp-mu-w/2)*(xp-mu-w/2)/(2*sigma*sigma))*flag_func(xp,mu+w/2,comparison='greater')
  return term_1+term_2+term_3


def flag_func2(x,lower, upper):
  """Returns a boolean of 1 between upper and lower limit. 0 outside limit"""
  n=len(x)
  y=[]
  for entry in range(n):
    if entry > lower and entry < upper:
      y.append(1.0)
    else:
      y.append(0.0)
  return y

def clamp_function(x, lower, upper):
  """Clamps a function to be 0 outside upper and lower limit"""
  y=[]
  n=len(x)
  for entry in range(n):
    if entry > lower and entry < upper:
      y.append(1.0/x[entry])
    else:
      y.append(0.0)
  return y

#############################################################################################
#											    #
#                       End defining helper functions; start real work;                     #
#			Define FEE_spec class and do all the work here	   		    #	
#											    #
#############################################################################################

class FEE_spec(object):
  def __init__ (self, run=223, Fe_foil_run=96, dithered_run=105):
    """ ---------- Explanation of variables ---------- 
        run = run number to do XRD processing, say run 223 in LS49
        Fe_foil_run = run number where Fe foil data was collected (96 for shift 2, 141 for shift 3) 
        dithered_run = vernier dithering of ebeam for beam profile correction. Data collected in run 105
    """
    self.run=run
    self.Fe_foil_run=Fe_foil_run
    self.dithered_run=dithered_run

  def vignetting_correction_h5(self, plot=False):
    """Vignetting correction using h5 files generated by Franklin
     There was no vignetting correction done in shift 2. Taking the data from runs 53-65 in shift 1 and hoping
     the correction factors are the same in shift 3
     Exgaussian fits looking better than gaussian"""

    runs = np.arange(53,65)
    pscam = []
    gauss=False
    exgauss=True
    if plot:
      plt.figure()
    for run in runs:
      with h5py.File('/reg/d/psdm/mfx/mfxls4916/hdf5/smalldata/ana-run{:03d}.h5'.format(run)) as fid:
        y = np.array(fid['FEE-SPEC0/hproj']).mean(axis=0)
        y0 = robust_mean(np.array(fid['gas_detector/f_11_ENRC']))
        x = range(len(y))
        #print ('FEE_GAS=',y0,np.array(fid['gas_detector/f_11_ENRC']))
        if gauss:
          popt,pcov = curve_fit(gaus,x,y/y0,p0=[y.max(),x[np.argmax(y)],fwhm(x,y)/2.355])
        if exgauss:
          popt,pcov = curve_fit(exgaussian,x,y/y0,p0=[y.max(),x[np.argmax(y)],fwhm(x,y)/2.355, fwhm(x,y)/2.355])
          popt[0]=max(exgaussian(x,*popt))
        pscam.append(popt)
        if plot:
          plt.plot(x,y/y0)
          if gauss:
            plt.plot(x,gaus(x,*popt))
          if exgauss:
            plt.plot(x,exgaussian(x,*popt), label='run= '+str(run))
    
    # Take the vignetting plot from above and fit a smooth line
    # Sweep the camera arm under constant central energy conditions
    axis = np.array([p[1] for p in pscam])
    if gauss:
        vals = np.array([p[0] for p in pscam])
    if exgauss:
        vals = np.array([p[0] for p in pscam])
        #vals=np.array([np.max(exgaussian(x,*p)) for p in pscam])
    s = np.argsort(axis)
    if plot:
      plt.figure()
      plt.plot(axis[s],vals[s],'-o',c='m')
      plt.xlabel('Energy Dispersion Axis')

    x=axis[s]
    y=vals[s]
    # Add satellite data points for consistency
    x=np.insert(x,0,0)
    y=np.insert(y,0,600000)
    x=np.append(x,2050)
    y=np.append(y,600000.0)
    self.vignetting_factor_cs=CubicSpline(x,y)
    xs=np.arange(min(x),max(x)-10,1.0)
    self.vignetting_factor=list(self.vignetting_factor_cs(xs))
    if plot:
      plt.plot(list(xs),vignetting_factor)
      plt.show()

  def beam_profile_correction(self, 
                              max_energy=7115.0, 
                              min_energy=7043.0, 
                              num_bins=18, 
                              apply_vignetting_corrections=True, 
                              plot=False):
    """ Code to do beam profile correction. Applies the vignetting correction if data is available and the
        flag is turned on. """
    ds = DataSource('exp=mfxls4916:run=%d:smd'%self.dithered_run)
    env = ds.env()
    fee_spec = AreaDetector('FEE_Spec', env)
    fee_gas_det = Detector('FEEGasDetEnergy')
    ebeam = Detector('EBeam')
    count=0
    max_event = 35000
    min_event = 0
    # Set up the arrays for binning
    delta_E=(max_energy-min_energy)/num_bins
    all_img_avg={}
    all_fee_gas={}
    central_energies=[]
    stored_central_energies=[]
    # Setup arrays in the dictionaries 
    for i_bin in range(1,num_bins+1):
      central_energy=min_energy+(i_bin-0.5)*delta_E
      # FEE spectrometer size hardcoded from shift 2
      all_img_avg[central_energy]=np.zeros((392, 2050),dtype=np.float32)
      all_fee_gas[central_energy]=[]
      central_energies.append(central_energy)
    
    for nevent,evt in enumerate(ds.events()):
      e = ebeam.get(evt)
      if nevent > max_event: break
      if nevent < min_event: continue
      if ebeam is None or e is None: continue
      # Filter further based on max and min energies otherwise outliers 
      # will mess things up 
      ev_ebeam = e.ebeamPhotonEnergy()
      if ev_ebeam > max_energy or ev_ebeam < min_energy: continue
      calib_array = fee_spec.calib(evt)
      if calib_array is not None:
        try:
          # Determine which bin ev_ebeam falls under
          fee_gas_data=fee_gas_det.get(evt).f_11_ENRC()
          # based on plot of sum(fee_intensity) vs fee_gas_intensity
          # there seems to be real signal only for fee_gas_intensity above 0.2 (linear regime)
          # Otherwise it looks like junk data
          if fee_gas_data < 0.2: continue
          central_energy,idx=find_nearest_idx(np.array(central_energies),ev_ebeam)
          all_fee_gas[central_energy].append(fee_gas_data)
          all_img_avg[central_energy]+=calib_array
        except Exception as e:
          print ('Error in calib_array',str(e))

    if True:
      psprof=[]
      all_fee_I_max=[]
      all_fee_I_px=[]
      if plot:
        plt.figure()
      # This is only for identifying the extremes of the distribution
      if False:
        #central_energies=[7112.5, 7057.5]
        central_energies=[7113.0, 7057.0]

      for i_energy, central_energy in enumerate(central_energies):
        img_avg=all_img_avg[central_energy]
        fee_gas=all_fee_gas[central_energy]
        fee_gas_correction=robust_mean_asmit(np.array(fee_gas))
        if fee_gas_correction is None: continue
        y=(img_avg[:,:-10].mean(axis=0) - img_avg[25:50,:-10].mean(axis=0))/fee_gas_correction
        x=range(len(y))
        if apply_vignetting_corrections and self.vignetting_factor is not None:
            y=y/self.vignetting_factor

        popt,pcov = curve_fit(gaus,x[:-10],y[:-10],p0=[y[:-10].max(),x[:-10][np.argmax(y[:-10])],fwhm(x[:-10],y[:-10])/2.355])
        psprof.append(popt)
        all_fee_I_max.append(np.max(y))
        all_fee_I_px.append(x[np.argmax(y)])
        #plot(x[300:1800],y[300:1800])
        if plot:
          plt.plot(x,y)
          plt.xlabel('Energy dispersion axis')
          plt.ylabel('FEE spectrometer intensity')
          plt.plot(x,gaus(x,*popt))
        #if central_energy in [7112.5,7057.5]:
        if central_energy in [7113.0, 7057.0]:
          print ('3 sigma boundaries ',popt[1]+3*popt[2],popt[1]-3*popt[2])
      # Gaussian fit to the sum of all 6 gaussians. Will help us obtain the 60 eV window on the spectrometer
      # Visually inspect this plot to get the reading for 60 eV window 
      if False:
          z = all_gaus(x,psprof)
          plt.plot(x, z, color='orange',label='sum of each of the gaussians at each point')

      if plot:
        plt.figure()
      axis = np.array([p[1] for p in psprof])
      vals = np.array([p[0] for p in psprof])
      s = np.argsort(axis)
      if plot:
        plt.plot(axis[s],vals[s],'o-',label='Max of each eV curve fit',c='y')
        plt.xlabel('Energy Dispersion Axis')
      x=axis[s]
      y=vals[s]#/gaus(x,*vignetting_params)
      # Add satellite data points for consistency
      y_min=1.0*min(y)
      x=np.insert(x,0,0)
      y=np.insert(y,0,0.5*y_min)
      x=np.insert(x,1,100)
      y=np.insert(y,1,0.55*y_min)
      #
      x=np.append(x,1648)
      y=np.append(y,0.54*y_min)
      x=np.append(x,1748)
      y=np.append(y,0.53*y_min)
      x=np.append(x,1848)
      y=np.append(y,0.52*y_min)
      x=np.append(x,1948)
      y=np.append(y,0.51*y_min)
      x=np.append(x,2048)
      y=np.append(y,0.5*y_min)
      self.beam_profile_cs=CubicSpline(x,y)
      #plot(x,beam_profile_cs(x),'*-',c='blue')
      # !!! Beam profile correction params saved here !!!
      if False:
       plt.figure()
       x=range(2048)
       xs=np.arange(min(x),max(x)-10,1.0)
       ys=self.beam_profile_cs(xs)
       plot(xs[500:1500],ys[500:1500],c='b')
    
    else:
      if plot:
        plt.figure()
        fee_gas_mean=robust_mean_asmit(np.array(all_fee_gas[central_energies[0]]))
        plt.imshow(all_img_avg[central_energies[0]][:,:-10]/fee_gas_mean)
    if plot:
      plt.show()
    
  def calibrate_fee_spec(self, plot=False):
    """Calibration
     Fe-edge position is known (x0,y0)=(1124,7112) (1141 from shift 2; 1124 from shift 3)
     Use delta_y = y2-y1 = 7115-7055 = 60eV
     Use delta_x = x2-x1= 139-1923 # Reading it off the gaussian fit (all_gaus; look at the extremes)
     Now fit a straight line using the equation
     (y-y0)/(x-x0)=(y2-y1)/(x2-x1) """

    #FIXME
    # !!! There are some hardcoded numbers in here right now - needs improvement !!!

    x0=self.Fe_edge_pixel # Accounting for 10 pixels ignored.  
    y0=7112 # Fe-edge 
    slope=60.0/(139-1923) # Ignore 10 pixels because they cancel out
    
    # defining slope using vignetting corrected beam profile data: use nbins=18
    # 1. using the mean of 7113 and 7057 instead of 3sigma
    #slope=(7113-7057)/(566.-1402.)
    # 2. use the 3 sigma values of 139 and 2055
    #slope=(72)/(139.0-2055.0)
    
    ##slope=60.0/(97-1963)
    intercept=y0-slope*x0
    print (slope, intercept)
    # Try Fe-edge value for consistency check
    x=1141
    y=slope*x+intercept
    print (x,y)
    boundary_pixels = []
    # Find pixel for 7122 eV
    y=7122
    x=(y-intercept)/slope
    print (x,y)
    boundary_pixels.append(int(x))
    # Find bounding pixels for 7122-30 & 7122+30 eV
    y=7092
    x=(y-intercept)/slope
    print (x,y)
    if x > 2000: x=2000
    boundary_pixels.append(int(x))
    y=7152
    x=(y-intercept)/slope
    print (x,y)
    if x <100: x=100
    boundary_pixels.append(int(x))
    print (boundary_pixels)
     
    # Store the slope and intercept
    self.slope=slope
    self.intercept=intercept
    
    if plot:
      x=[550,650,750,826,850,950,1050,1141,1250,1350,1450,1550,1650]
      y=[xx*slope+intercept for xx in x]
      plt.figure()
      plt.plot(x,y,'*')
      plt.show()

  def get_photon_energy(self, run):
    pass 

  def get_Fe_edge(self, plot=False):
    """ Return the pixel on the FEE spectrometer that is the Fe edge i.e 7112 eV. Using 0-indexing for pixels 
        Use the concept of it being the inflection point to determine it from the Fe_foil_run 
        Hard coding the spectrometer size that is in the stream from shift 2.
        spectrometer sizes
        shift 2 = 392, 2050
        shift 3 = 392, 2048
        full spectrometer recording area from the website= 2048, 2048 """

    ds = DataSource('exp=mfxls4916:run=%d:smd'%self.Fe_foil_run)
    env = ds.env()
    fee_spec = AreaDetector('FEE_Spec', env)
    img_avg = np.zeros((392, 2050),dtype=np.float32)
    count = 0
    max_events=4000 # Just putting a cutoff because sometimes towards the end of a run, ebeam does weird things
    for nevent,evt in enumerate(ds.events()):
      calib_array = fee_spec.calib(evt)
      #print ('Calib array shape = ',calib_array.shape)
      if calib_array is not None:
        img_avg += calib_array
        count += 1
        if count > max_events: break 

    # Get pixel position of Fe edge
    # Take the mid point as the inflection point
    # According to NKS, the curve is inverted and energy increases from right to left (i.e 2050--> 0 pixels)
    # Further update --> inflection point is the midpoint of the slower rising curve, not the faster one
    
    Fe_edge_data = img_avg[:,:-10].mean(axis=0) - img_avg[50:100,:-10].mean(axis=0)
    Fe_gradient=list(np.gradient(Fe_edge_data))
    from scitbx.array_family import flex
    Fe_gradient_sorted_idx=flex.sort_permutation(flex.double(Fe_gradient), reverse=True)
    # print out the top 5 indexes; the first one is an outlier at 1660 for run 96
    #print ('Top 5 indexes of Fe-gradient=', list(Fe_gradient_sorted_idx[:5]))
    # !!!!!!!! Take the second one if run 96 !!!!!!!!!!!!!!!!
    self.Fe_edge_pixel = Fe_gradient_sorted_idx[1]
    print ('Fe edge value intensity and pixel value with 0-indexing=', Fe_edge_data[self.Fe_edge_pixel], self.Fe_edge_pixel)

    # This is the plot of the spectrometer data
    # Background subtraction is done here
    # Run this for run 96 and see sharp iron edge
    if plot:
      plt.figure()
      fee_data=img_avg[:,:-10].mean(axis=0) - img_avg[50:100,:-10].mean(axis=0)
      plt.plot(fee_data[:])
      plt.xlabel('pixels on spectrometer')
      plt.ylabel('Intensity recorded on spectrometer')
      #Visualize the raw data on the spectrometer
      plt.figure()
      fee_data=img_avg[:,:-10]-img_avg[25:50,:-10].mean(axis=0)
      plt.imshow(fee_data)
      plt.show()

if __name__ == '__main__':
  FEE = FEE_spec()
  FEE.get_Fe_edge(plot=False)
  #FEE.vignetting_correction_h5(plot=False)
  #FEE.beam_profile_correction(plot=False)
  FEE.calibrate_fee_spec(plot=True)
