from __future__ import print_function, division, absolute_import

from psana import *
from Detector.EpicsDetector import EpicsDetector
from Detector.AreaDetector import AreaDetector
import numpy as np
#import matplotlib.pyplot as plt
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
  """ Exgaussian functional; See https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution"""
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
      # Save the gaussian mean positions. Might be used to calibrate fee_spectra
      self.gauss_window_beam_params=psprof
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


  def beam_profile_horizontal(self,
                              max_energy=7115.0,
                              min_energy=7043.0,
                              num_bins=18,
                              apply_vignetting_corrections=True,
                              plot=False):
    """Try to get the beam profile correction by getting profile from
       the horizontal direction. Try to be smart about running the binning calculations"""



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
#
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
        print ('CENTRAL', central_energy)
        img_avg=all_img_avg[central_energy]
        fee_gas=all_fee_gas[central_energy]
        fee_gas_correction=robust_mean_asmit(np.array(fee_gas))
        if fee_gas_correction is None: continue
        y=(img_avg[:,:-10].mean(axis=1) - img_avg[:,25:50].mean(axis=1))/fee_gas_correction
        x=range(len(y))
        v=np.argmax(img_avg[:,:-10].max(axis=0) - img_avg[25:50,:-10].mean(axis=0))
        # Not sure how to apply vignetting corrections to horizontal data
        if False: #apply_vignetting_corrections and self.vignetting_factor is not None:
            y=y/self.vignetting_factor
        popt,pcov = curve_fit(gaus,x[:-10],y[:-10],p0=[y[:-10].max(),x[:-10][np.argmax(y[:-10])],fwhm(x[:-10],y[:-10])/2.355])
        popt[1]=v*1.0
        psprof.append(popt)
        disp_axis=range(2038)
        all_fee_I_max.append(np.max(y))
        all_fee_I_px.append(x[np.argmax(y)])
        if plot:
          plt.xlabel('Energy dispersion axis-Horizontal')
          plt.ylabel('FEE spectrometer intensity')
          plt.plot(disp_axis, gaus(disp_axis, *popt))
        #if central_energy in [7112.5,7057.5]:
        if central_energy in [7113.0, 7057.0]:
          print ('3 sigma boundaries ',popt[1]+3*popt[2],popt[1]-3*popt[2])
      if plot:
        plt.show()

    
  def calibrate_fee_spec(self, plot=False):
    """Calibration
     Fe-edge position is known (x0,y0)=(1124,7112) (1141 from shift 2; 1124 from shift 3)
     Use delta_y = y2-y1 = 7115-7055 = 60eV
     Use delta_x = x2-x1= 139-1923 # Reading it off the gaussian fits (look at the extremes)
     Now fit a straight line using the equation
     (y-y0)/(x-x0)=(y2-y1)/(x2-x1) """

    #FIXME
    # !!! There are some hardcoded numbers in here right now - needs improvement !!!

    x0=self.Fe_edge_pixel # Accounting for 10 pixels ignored.  
    y0=7112 # Fe-edge 
    slope=60.0/(139-1923) # Ignore 10 pixels because they cancel out
    # attempt to take a mean/median slope instead of taking it between 2 extreme datapoints
    # Note that the diffeence is 4 eV windows.
    all_slopes=[]
    window_means=[p[1] for p in self.gauss_window_beam_params]
    for ii in range(len(window_means)-1):
      tmp_slope=-4.0/(window_means[ii]-window_means[ii+1])
      all_slopes.append(tmp_slope)
    
    # defining slope using vignetting corrected beam profile data: use nbins=18
    # 1. using the mean of 7113 and 7057 instead of 3sigma
    #slope=(7113-7057)/(566.-1402.)
    # 2. use the 3 sigma values of 139 and 2055
    #slope=(72)/(139.0-2055.0)
    # 3. take average slope between 7115 and 7043 means over 4ev windows
    #slope=-0.07295042
    #slope=np.mean(all_slopes)
    # 4. take median slope between 7115 and 7043 means over 4eV windows
    # take median
    slope=np.median(all_slopes)
    
    ##slope=60.0/(97-1963)
    intercept=y0-slope*x0
    print ('Slope=%.3f, Intercept=%.3f'%(slope, intercept))
    # Try Fe-edge value for consistency check
    x=1141
    y=slope*x+intercept
    print (x,y)
    self.boundary_pixels = []
    # Find pixel for 7122 eV
    y=7122
    x=(y-intercept)/slope
    print (x,y)
    self.boundary_pixels.append(int(x))
    # Find bounding pixels for 7122-30 & 7122+30 eV
    y=7092
    x=(y-intercept)/slope
    print (x,y)
    if x > 2000: x=2000
    self.boundary_pixels.append(int(x))
    y=7152
    x=(y-intercept)/slope
    print (x,y)
    if x <100: x=100
    self.boundary_pixels.append(int(x))
    print (self.boundary_pixels)
     
    # Store the slope and intercept
    # 
    self.slope=slope
    self.intercept=intercept
    
    if plot:
      x=[550,650,750,826,850,950,1050,1141,1250,1350,1450,1550,1650]
      y=[xx*slope+intercept for xx in x]
      plt.figure()
      plt.plot(x,y,'*')
      plt.show()

  def get_center_of_ev_spec(self,fee_spec, mode='weighted_mean'):
    """ Function to calculate a single representative photon energy from the FEE spectrometer reading.
        Various choices available. LCLS for example uses gaussian 1st moment"""
    numer=0.0
    denom=0.0
    # somehow the spectrometer size changed again during the myoglobin runs
    assert self.slope is not None, 'ev-axis calibration slope is None. please check code '
    assert self.intercept is not None, 'ev-axis calibtation intercept is None. please check code'
    ev_axis=[x*self.slope+self.intercept for x in range(2038)]

    if mode=='gauss_first_moment':
      y=fee_spec
      x=range(len(fee_spec))
      popt,pcov = curve_fit(gaus,x,y,p0=[y.max(),x[np.argmax(y)],fwhm(x,y)/2.355])
      x_at_first_moment=popt[1]
      print ('GAUSS_1ST_MOMENT')
      return x_at_first_moment*self.slope+self.intercept

    # Return photon energy at which FEE reading is maximum.
    if mode=='max_fee':
      return ev_axis[np.argmax(fee_spec)]

    if mode=='fixed_value':
      return 7122.0 # Return a fixed value; this is mostly for debugging

    # Return weighted mean of photon energy. Weighted by FEE intensity data
    # Reference: Franson et.al Biochemistry (2018), 57, 4629
    # Eq (1) 
    else:
      for ii, fee_data in enumerate(fee_spec):
        numer += float(ev_axis[ii])*fee_data
        denom += fee_data
      return numer/denom

  def plot_correlation_ev_ebeam(self):
    """  Function to see the correlation between eV reading from FEE spectrometer vs 
         ebeam reading from psana. They should be highly correlated (i.e > 0.9)"""

  # Plot correlation between ebeam data and the eV data. Should be positively correlated

    from xfel.cxi.cspad_ana import cspad_tbx
    ebeam_energy=[]
    fee_energy=[]
    xfel_ts = []

    ds = DataSource('exp=mfxls4916:run=%d:smd'%self.run)
    env = ds.env()
    fee_spec = AreaDetector('FEE_Spec', env)
    ebeam = Detector('EBeam')
    plt.figure()
    xs=range(2038)
    scale_factor = self.beam_profile_cs(xs) 
    ev_axis=[x*self.slope+self.intercept for x in range(2038)]
    for nevent,evt in enumerate(ds.events()):
      calib_array = fee_spec.calib(evt)
      e=ebeam.get(evt)
      if calib_array is not None and e is not None:
        data=calib_array[:,:-10].sum(axis=0) - calib_array[25:50, :-10].sum(axis=0)
        scaled_data=(data-min(data))/(scale_factor)
        ebeam_energy.append(e.ebeamPhotonEnergy())
        fee_energy.append(self.get_center_of_ev_spec(scaled_data))
        ts=cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt))
        xfel_ts.append(ts)
        #if nevent > 2000: break
    
    plt.plot(ebeam_energy,fee_energy,'*')
    print ('correlation coefficient', np.corrcoef(ebeam_energy, fee_energy))
    plt.xlabel('EBEAM energy')
    plt.ylabel('FEE weighted central energy')
    plt.show()

  def save_spectra_of_select_events(self, timestamps, fname='fee_spectra_info.pickle', save_event=True):
    ''' Save spectra of select timestamps in the form of a pickle file
        @params: timestamps is a list containing cctbx style timestamp eg. ['20180501143559313'] ''' 

    from xfel.cxi.cspad_ana import cspad_tbx
    from scitbx.array_family import flex
    fee_energy={}
    fee_spectra = {}

    ds = DataSource('exp=mfxls4916:run=%d:smd'%self.run)
    env = ds.env()
    fee_spec = AreaDetector('FEE_Spec', env)
    xs=range(2038)
    scale_factor = self.beam_profile_cs(xs)
    ev_axis=[x*self.slope+self.intercept for x in range(2038)]
    for nevent,evt in enumerate(ds.events()):
      ts=cspad_tbx.evt_timestamp(cspad_tbx.evt_time(evt))
      xfel_ts = ts[0:4] + ts[5:7] + ts[8:10] + ts[11:13] + ts[14:16] + ts[17:19] + ts[20:23]
      if xfel_ts not in timestamps: continue
      calib_array = fee_spec.calib(evt)
      print ('PLEASE SAVE THIS', xfel_ts)
      if calib_array is not None:
        data=calib_array[:,:-10].sum(axis=0) - calib_array[25:50, :-10].sum(axis=0)
        scaled_data=(data-min(data))/(scale_factor)
        fee_data_for_this_shot = flex.vec2_double()
        for x,y in zip(ev_axis, scaled_data):
          fee_data_for_this_shot.append((x,y))
        fee_spectra[xfel_ts] = fee_data_for_this_shot
        fee_energy[xfel_ts] = self.get_center_of_ev_spec(scaled_data)
    if save_event:
      from libtbx.easy_pickle import dump
      dump(fname, fee_spectra) 


  def plot_fee_spectra_of_event(self, event_timestamp=None):
    """ Plot fee spectra of a certain timestamp. Also print out the ebeam value. Useful for sanity checks"""
    
    from xfel.cxi.cspad_ana import cspad_tbx
    assert event_timestamp is not None, 'event_timestamp cannot be None. Please supply event timestamp you are interested in'
    ds = DataSource('exp=mfxls4916:run=%d:smd'%self.run)
    env = ds.env()
    fee_spec = AreaDetector('FEE_Spec', env)
    ebeam = Detector('EBeam')
    xs=range(2038)
    scale_factor = self.beam_profile_cs(xs)
    ev_axis=[x*self.slope+self.intercept for x in range(2038)]
    for nevent,evt in enumerate(ds.events()):
      tx = evt.get(EventId).time()
      sec=tx[0]
      nsec=tx[1]
      ts = cspad_tbx.evt_timestamp((sec,nsec/1e6))
      if ts!=event_timestamp: continue
      calib_array = fee_spec.calib(evt)
      e=ebeam.get(evt)
      if calib_array is not None and e is not None:
        data=calib_array[:,:-10].sum(axis=0) - calib_array[25:50, :-10].sum(axis=0)
        scaled_data=(data-min(data))/(scale_factor)
        ebeam_energy=e.ebeamPhotonEnergy()
        fee_energy=self.get_center_of_ev_spec(scaled_data)
        print ('Mean energy of event = %.3f (fee), %.3f (ebeam)'%(fee_energy, ebeam_energy))
        import numpy as np
        #import matplotlib.pyplot as plt
        y=np.arange(min(scaled_data), max(scaled_data), 1000)
        plt.figure()
        plt.plot(ev_axis, scaled_data)
        plt.plot([fee_energy]*len(y), y, '-r')
        plt.show()

  def get_photon_energy(self, run=199, evt=None,mode='weighted_mean'):
    """ Provided an event and a run number, returns the photon energy of that event 
        The photon energy is a weighted mean of the reading from the FEE spectrometer
        yielding a single number"""
    ds = DataSource('exp=mfxls4916:run=%d:smd'%run)
    env = ds.env()
    fee_spec = AreaDetector('FEE_Spec', env)
    #import pdb; pdb.set_trace()

    rayonix = Detector('MfxEndstation.0:Rayonix.0')
    fee_gas_det = Detector('FEEGasDetEnergy')
    #jungfrau1M = Detector('MfxEndstation.0:Jungfrau.0')
    ebeam = Detector('EBeam')

    calib_array = fee_spec.calib(evt)
    ev_ebeam_reading=ebeam.get(evt)
    fee_gas_reading=fee_gas_det.get(evt)
    if calib_array is None or ev_ebeam_reading is None or fee_gas_reading is None:
      return None

    xs=range(2038)
    scale_factor=self.beam_profile_cs(xs)
    ev_axis=[x*self.slope+self.intercept for x in range(2038)]
    if calib_array is not None:
      data=calib_array[:,:-10].sum(axis=0) - calib_array[25:50, :-10].sum(axis=0)
      scaled_data=(data-min(data))/(scale_factor)
      photon_energy=self.get_center_of_ev_spec(scaled_data, mode=mode)
      return photon_energy
    return None 
      
  def analyze_photon_energy(self, plot=False):
    """ This function can be used to analyze and plot scaled and unscaled FEE spec
         reading. Useful to see the effect if vignetting/beam profile corrections/calibration"""
    ds = DataSource('exp=mfxls4916:run=%d:smd'%self.run)
    env = ds.env()
    ebeam = Detector('EBeam')
    fee_spec = AreaDetector('FEE_Spec', env)
    ebeam_energy=[]
    fee_energy=[]
    if plot:
      fig1=plt.figure()
      fig2=plt.figure()
      ax1=fig1.gca()
      ax2=fig2.gca()
    ev_axis=[x*self.slope+self.intercept for x in range(2038)]
    #print max(ev_axis)
    #print scaled_data.shape
    for nevent,evt in enumerate(ds.events()):
      if nevent < 20: continue
      calib_array = fee_spec.calib(evt)
      e=ebeam.get(evt)
      #if e is not None:
        #print ('EBEAM photon energy = ',e.ebeamPhotonEnergy())
      if calib_array is not None:
        data=calib_array[:,:-10].sum(axis=0) - calib_array[25:50, :-10].sum(axis=0)
        xs=range(2038)
        scale_factor=self.beam_profile_cs(xs)
        if plot:
          ax1.plot(ev_axis,data-min(data))
          #ax1.plot(ev_axis,np.log10(1+data))
          #ylim((0,10.0))
        # scale data properly by beam profile correction factor
        scaled_data=(data-min(data))/(scale_factor)
        #scaled_data=data/scale_factor
        if plot:
          ax2.plot(ev_axis, scaled_data)
          #ax2.plot(ev_axis, np.log10(1.0+scaled_data))
        if e is not None:
          ebeam_energy.append(e.ebeamPhotonEnergy())
          fee_energy.append(self.get_center_of_ev_spec(scaled_data))
          #print ('EBEAM photon energy=',e.ebeamPhotonEnergy(), self.get_center_of_ev_spec(scaled_data))
        if plot:
          #xlim((500,1500))
          #print max(scaled_data)
          ax1.plot([ev_axis[self.boundary_pixels[0]]]*2,[0,max(data)],'-',linewidth=5, color='yellow', alpha=0.3)
          ax1.plot([ev_axis[self.boundary_pixels[1]]]*2,[0,max(data)],'-',linewidth=5, color='grey', alpha=0.1)
          ax1.plot([ev_axis[self.boundary_pixels[2]]]*2,[0,max(data)],'-',linewidth=5, color='grey', alpha=0.1)
          ax2.plot([ev_axis[self.boundary_pixels[0]]]*2,[0,max(scaled_data)],'-',linewidth=5, color='yellow', alpha=0.3)
          ax2.plot([ev_axis[self.boundary_pixels[1]]]*2,[0,max(scaled_data)],'-',linewidth=5, color='grey', alpha=0.1)
          ax2.plot([ev_axis[self.boundary_pixels[2]]]*2,[0,max(scaled_data)],'-',linewidth=5, color='grey', alpha=0.1)
          ax1.set_xlabel('eV')
          ax1.set_ylabel('unscaled intensity')
          ax2.set_xlabel('eV')
          ax2.set_ylabel('scaled intensity')
          #ax1.set_ylim((0,8.0))
          #print ('Calib array shape = ',calib_array.shape)

        #if nevent > 50: break
    print ('Printing summary statistics')
    print ('Mean energies = %.2f (fee), %.2f (ebeam)'%(np.mean(fee_energy), np.mean(ebeam_energy)))
    if plot:
      plt.show()

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
    if self.Fe_foil_run==96:
      img_avg = np.zeros((392, 2050),dtype=np.float32)
    elif self.Fe_foil_run==143:
      img_avg = np.zeros((392, 2048),dtype=np.float32)
      
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
    if self.Fe_foil_run==96:
      self.Fe_edge_pixel = Fe_gradient_sorted_idx[1]
    elif self.Fe_foil_run==143:
      self.Fe_edge_pixel = Fe_gradient_sorted_idx[0]
      
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


def cctbx_to_psana_ts_converter(cctbx_ts):
  ''' Given a cctbx style timestamp, convert it to a psana style timestamp
      Eg. 20180501143559313 gets converted to 2018-05-01T14:35Z59.313''' 
  return cctbx_ts[0:4]+'-'+cctbx_ts[4:6]+'-'+cctbx_ts[6:8]+'T'+cctbx_ts[8:10]+':'+\
         cctbx_ts[10:12]+'Z'+cctbx_ts[12:14]+'.'+cctbx_ts[14:17]


if __name__ == '__main__':
  short_circuit=True
  dump_pickle=False
  if short_circuit and dump_pickle :
    FEE = FEE_spec(run=219,Fe_foil_run=143)
    FEE.get_Fe_edge(plot=False)
    FEE.vignetting_correction_h5(plot=False)
    FEE.beam_profile_correction(plot=False)
    FEE.calibrate_fee_spec(plot=False)
    from libtbx.easy_pickle import dump
    dump('FEE.pickle', FEE)
    exit()
  if not dump_pickle:
    from libtbx.easy_pickle import load
    # On psana
    #FEE=load('/reg/d/psdm/mfx/mfxls4916/scratch/asmit/LS49_SAD_v3/input/FEE_r143_v1.pickle')
    # On NERSC
    FEE=load('/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/input/FEE_r143_v1.pickle')
    print ('Should be plotting stuff')
    FEE.run=222
    #FEE.analyze_photon_energy(plot=False)
    #FEE.beam_profile_horizontal(plot=False)
    #FEE.plot_correlation_ev_ebeam()
    #timestamps_of_interest = ['20180501143559313']
    #timestamps_of_interest = ['20180501143533988',
    #                          '20180501143546717',
    #                          '20180501143546817',
    #                          '20180501143547150',
    #                          '20180501143548650',
    #                          '20180501143549416',
    #                          '20180501143549949',
    #                          '20180501143551715',
    #                          '20180501143555114',
    #                          '20180501143559313',
    #                          '20180501143602713',
    #                          '20180501143606545',
    #                          '20180501143620206',
    #                          '20180501143625171',
    #                          '20180501143628702',
    #                          '20180501143628902',
    #                          '20180501143631168',
    #                          '20180501143632300',
    #                          '20180501143640763',
    #                          '20180501143643462',
    #                          '20180501143643662',
    #                          '20180501143652325',
    #                          '20180501143701853']
    # Timestamps of events indexed by IOTA-RANSAC for LS49
    timestamps_of_interest = [
                              '20180501143547850',
                              '20180501143545684',
                              '20180501143641996',
                              '20180501143545184',
                              '20180501143628303',
                              '20180501143602713',
                              '20180501143546817',
                              '20180501143651959',
                              '20180501143631101',
                              '20180501143546717',
                              '20180501143620206',
                              '20180501143546017',
                              '20180501143552248',
                              '20180501143651292',
                              '20180501143626270',
                              '20180501143652325',
                              '20180501143625171',
                              '20180501143549949',
                              '20180501143648827',
                              '20180501143620239',
                              '20180501143546983',
                              '20180501143606545',
                              '20180501143549216',
                              '20180501143533988',
                              '20180501143643462',
                              '20180501143548650',
                              '20180501143628902',
                              '20180501143643662',
                              '20180501143546883',
                              '20180501143535254',
                              '20180501143547150',
                              '20180501143549416',
                              '20180501143701853',
                              '20180501143603279',
                              '20180501143657722',
                              '20180501143632300',
                              '20180501143651092',
                              '20180501143640763',
                              '20180501143628702',
                              '20180501143631568',
                              '20180501143551715',
                              '20180501143533955',
                              '20180501143631201',
                              '20180501143547483',
                              '20180501143549183',
                              '20180501143648193',
                              '20180501143559313',
                              '20180501143630268',
                              '20180501143641130',
                              '20180501143555114',
                              '20180501143645328',
                              '20180501143650626',
                              '20180501143631168',
                              '20180501143553215']
    #for ts in timestamps_of_interest:
    #  FEE.plot_fee_spectra_of_event(cctbx_to_psana_ts_converter(ts))
    FEE.save_spectra_of_select_events(timestamps_of_interest, fname='fee_data_r0222_iota.pickle')
    # Timestamp used in run 2019 as an example '2018-05-01T14:25Z24.942' 
    #FEE.plot_fee_spectra_of_event(event_timestamp='2018-05-01T14:25Z24.942')
