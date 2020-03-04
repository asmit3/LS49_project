#http://homepages.inf.ed.ac.uk/rbf/HIPR2/fourier.htm
#http://journals.iucr.org/j/issues/2007/06/00/kk5014/kk5014.pdf
import sys, os
import numpy as np
import dxtbx
from scitbx.array_family import flex
from numpy.fft import *
import math
from scitbx.smoothing import convolve
#from scipy import signal
#from scipy.ndimage.filters import convolve
from dials.algorithms.image.filter import convolve
import libtbx
from dxtbx.format.cbf_writer import FullCBFWriter
import pycbf
seed = 22 #int(sys.argv[1])
np.random.seed(seed)


class PSF_deconvolve(object):
  ''' A class to deconvolve the PSF on a rayonix image
      Options present to use either power law PSF or gaussian PSF as proposed in Holton 2012 (Journal of Synchotron Rad)
      Also options present to use fft or Lucy-Richardson method to deconvolve PSF. FFT method is not perfect and needs an offset 
      when writing out the final image file'''

  def __init__(self):
    pass

  def convolve_fft(self, Input, psf, epsilon):
      InputFFT = fftn(Input)
      psfFFT = fftn(psf)+epsilon
      convolved = ifftn(InputFFT*psfFFT)
      convolved = np.abs(convolved)
      return convolved
  
  def deconvolve_fft(self, Input, psf, epsilon):
      InputFFT = fftn(Input)
      psfFFT = fftn(psf)+epsilon
      deconvolved = ifftn(InputFFT/psfFFT)
      deconvolved = np.abs(deconvolved)
      return deconvolved
  
  def makeGaussPSF(self, fwhm_pixel, sizeX, sizeY): 
      ''' Do not use this for detectors as integral form is more appropriate'''
      x,y = np.mgrid[-sizeY/2.:sizeY/2., -sizeX/2.:sizeX/2.]
      psf = np.exp(-np.log(16)*(x**2+y**2)/(fwhm_pixel*fwhm_pixel))
      psf = psf/psf.sum()
      psf = psf.tolist()
      psf = flex.double(psf)
      return psf
  
  def makeMoffatPSF(self, fwhm_pixel, sizeX, sizeY):
      ''' Do not use this for detectors as integral form is more appropriate'''
      g = fwhm_pixel*0.65238
      psf = np.zeros((sizeX, sizeY))
      sx = int(sizeX/2)
      sy = int(sizeY/2)
      for y in range(sy, -(sy+1), -1):
          for x in range(-sx, sx+1):
              psf[x+sx,-y+sy] = g/(2*np.pi*np.power((g*g+(1.0*x*x)+(1.0*y*y)),1.5))
      psf = psf/psf.sum()
      psf = psf.tolist()
      psf = flex.double(psf)
      return psf
  
  def fiber2D_integ(self, x,y,g):
      return math.atan((x*y)/(g*math.sqrt(g*g + x*x + y*y)))/(2.0*math.pi);
  
  def makeMoffat_integPSF(self, fwhm_pixel, sizeX, sizeY):
      ''' Integral form of contribution of moffat PSF to signal recorded in a square pixel'''
      g = fwhm_pixel*0.65238
      psf = np.zeros((sizeX, sizeY))
      sx = int(sizeX/2)
      sy = int(sizeY/2)
      for y in range(sy, -(sy-1), -1):
      #for y in range(sy, -(sy+1), -1):
        for x in range(-(sx), sx):
          #for x in range(-sx, sx+1):
          # Holton 2012 paper says x,y should be pixel center; this does not seem right ?
          #psf[x+sx, -y+sy] = self.fiber2D_integ(x+1./2,y+1./2,g)-self.fiber2D_integ(x+1./2,y-1./2,g)-self.fiber2D_integ(x-1./2,y+1./2,g)+self.fiber2D_integ(x-1./2,y-1./2,g)
          # Trying to get pixel center instead
          psf[x+sx, -y+sy] = self.fiber2D_integ(x+1,y+1,g)-self.fiber2D_integ(x+1,y,g)-self.fiber2D_integ(x,y+1,g)+self.fiber2D_integ(x,y,g)
      psf = psf/psf.sum()
      psf = psf.tolist()
      psf = flex.double(psf)
      return psf
  
  #====================GAUSSIAN STUFF=======================
  
  def gauss2D_integ(self, x_in,y_in,fwhm):
      x=x_in/fwhm
      y=y_in/fwhm
      return 0.125*(math.erf(2.0*x*np.sqrt(np.log(2.0)))*math.erf(y*np.sqrt(np.log(16.0)))*np.sqrt(np.log(16.0)/np.log(2.0))); 
  
  def makeGauss_integPSF(self, fwhm_pixel, sizeX, sizeY):
      ''' Integral form of contribution of gaussian PSF to signal recorded in a square pixel'''
      g = fwhm_pixel*1.0
      psf = np.zeros((sizeX, sizeY))
      sx = int(sizeX/2)
      sy = int(sizeY/2)
      for y in range(sy, -(sy-1), -1):
      #for y in range(sy, -(sy+1), -1):
        for x in range(-(sx), sx):
          #for x in range(-sx, sx+1):
          #psf[x+sx, -y+sy] = self.gauss2D_integ(x+1,y+1,g)-self.gauss2D_integ(x+1,y-1,g)-self.gauss2D_integ(x,y+1,g)+self.gauss2D_integ(x,y,g)
          psf[x+sx, -y+sy] = self.gauss2D_integ(x+1./2,y+1./2,g)-self.gauss2D_integ(x+1./2,y-1./2,g)-self.gauss2D_integ(x-1./2,y+1./2,g)+self.gauss2D_integ(x-1./2,y-1./2,g)
      psf = psf/psf.sum()
      psf = psf.tolist()
      psf = flex.double(psf)
      return psf
  
  #===========================================
  
  def exportArrayAsImage(self, img, filename, header=None):
      img = np.abs(img)
      img = np.where(img > 65536, 65536, img) 
      img = img.astype(np.uint16)
      fout = open(filename,'wb')
      if header is None:
          fout.write('{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;\nTYPE=unsigned_short;\n')
          fout.write('SIZE1=960;\nSIZE2=960;\nPIXEL_SIZE=0.177;\nDISTANCE=141.7;\nWAVELENGTH=1;\n')
          fout.write('BEAM_CENTER_X=-17.88;\nBEAM_CENTER_Y=66.40;\nORGX=511.5;\nORGY=512.5;\nCLOSE_DISTANCE=300;\nPHI=0;\n')
          fout.write('OSC_START=0;\nOSC_RANGE=0;\nTWOTHETA=0;\nDETECTOR_SN=000;\nBEAMLINE=fake;\n}\f')
      else:
          fout.write(header)
      row,col = img.shape 
      while fout.tell() < 512:
          fout.write(' ')
      img = (np.reshape(img,row*col))
      fout.write(img)
  
  def gaussnoise(self,image,mean,sigma):
      row,col = image.shape
      gauss = np.random.normal(mean, sigma, (row,col))
      gauss = gauss.reshape(row,col)
      return gauss
  
  def poissonnoise(self, image):
      row,col = image.shape
      poidev = np.random.poisson(image)
      return poidev
  
  def divergence(self, dF_dx, dF_dy):
      """ compute the divergence of n-D scalar field `F` """
      print 'd2Fdx2 info ',((np.gradient(dF_dx))[0].ndim)
      d2F_dx2 = (np.gradient(dF_dx))[0]
      d2F_dy2 = (np.gradient(dF_dy))[1]
      return d2F_dx2+d2F_dy2
  #    return np.sum(d2F_dx2,d2F_dy2)
  #    return reduce(np.add,np.gradient(F))
  
  
  def l1_norm(self, dF_dx, dF_dy):
     return np.linalg.norm(dF_dx,1)+np.linalg.norm(dF_dy,1)
  
  
  def tv_func(self, img, lambd):
      img_np = img.as_numpy_array()
      do_dx, do_dy = np.gradient(img_np)
      l1_grad_img = self.l1_norm(do_dx, do_dy)
      div_img = self.divergence(do_dx/l1_grad_img,do_dy/l1_grad_img)
      tvf_fact =  (1-lambd*div_img)
      tvf_fact = tvf_fact.tolist()
      tvf_fact = flex.double(tvf_fact)
      return tvf_fact
  
  
  
  
  def rl_standard(self, raw_image, psf, niter):
      """ Standerd lucy-richardson convolution
      arXiv 2002 Lauer
      https://stackoverflow.com/questions/9854312/how-does-richardson-lucy-algorithm-work-code-example
      """
     
  #    psf /= psf.sum() # PSF should always be normalized !
      xpsf, ypsf = psf.all()
      psf_hat = flex.double(range(xpsf*ypsf))
      psf_hat.reshape(flex.grid(xpsf, ypsf))
      for n in range(xpsf):
          for m in range(ypsf):
              for i in range(n,xpsf):
                  for j in range(m,ypsf):
                      psf_hat[n,m] = psf[i-n,j-m]
  # Take an initial guess of the image
      #lucy = 0.5*np.ones( raw_image.shape )
      lucy = 1.0*raw_image
      for i in xrange( niter ):
          print 'iteration # = ',i
          estimate = convolve(lucy, psf) + 1e-300 # Denominator (u(t)*p)
          #estimate[ np.isnan(estimate) ] = 0
          relative_blur = raw_image/estimate # (d/(u(t)*p)
          correction = convolve(relative_blur, psf_hat) #+ 1e-300 #([d/(u(t)*p)]*p_hat)
          #correction[ np.isnan(correction) ] = 0
          lucy *= correction   # u(t+1) = u(t)*[{d/(u(t)*p)}*p_hat]
      return lucy
  
  def rl_regularized_1(self, raw_image, psf, niter):
      """ Regularized lucy-richardson convolution
          http://cens.ioc.ee/~pearu/theses/martin_msc_2009.pdf
      """
  
  #    psf /= psf.sum() # PSF should always be normalized !
      xpsf, ypsf = psf.all()
      psf_hat = flex.double(range(xpsf*ypsf))
      psf_hat.reshape(flex.grid(xpsf, ypsf))
      for n in range(xpsf):
          for m in range(ypsf):
              for i in range(n,xpsf):
                  for j in range(m,ypsf):
                      psf_hat[n,m] = psf[i-n,j-m]
  # Take an initial guess of the image
  #    lucy = 0.5*np.ones( raw_image.shape )
      lucy = 1.0*raw_image
      lambd = 0.2
      for i in xrange( niter ):
          print 'iteration # = ',i
          estimate = convolve(lucy, psf) + 1e-15 # Denominator (u(t)*p)
  #        estimate[ np.isnan(estimate) ] = 0
          relative_blur = raw_image/estimate # (d/(u(t)*p)
          correction = convolve(relative_blur, psf_hat) + 1e-15 #([d/(u(t)*p)]*p_hat)
  #        correction[ np.isnan(correction) ] = 0
          tv_factor = self.tv_func(lucy, lambd)
          lucy /= tv_factor
          lucy *= (correction)   # u(t+1) = u(t)*[{d/(u(t)*p)}*p_hat]
  #        lucynp *= (correction/(1-lambd*(divO/1.)))   # u(t+1) = u(t)*[{d/(u(t)*p)}*p_hat]
  
      return lucy
  
  
  def find_spots(self, image,bg):
      # Construct a sub-image of dimension 3*3 (can be changed if needed)
      row,col = image.shape
      spot_location = []
      spot_intensity = []
      alt_image = np.lib.pad(image, ((1,1),(1,1)), 'constant', constant_values=0)
      subimage = np.zeros((3,3))
      for i in range(1,row) :
          for j in range(1,col):
              subimage = alt_image[i-1:i+2,j-1:j+2]
              if np.argmax(subimage) == 4 and subimage[1,1] > bg:
                  spot_intensity.append(np.sum(subimage))
                  spot_location.append([i-1,j-1])
      print spot_intensity
  
  def find_spots_1x1(self, image,bg):
      # Construct a sub-image of dimension 3*3 (can be changed if needed)
      row,col = image.shape
      spot_location = []
      spot_intensity = []
      alt_image = np.lib.pad(image, ((1,1),(1,1)), 'constant', constant_values=0) # 2x2
      subimage = np.zeros((3,3)) # 5x5
      for i in range(1,row) :    # 1 or 2
          for j in range(1,col): # 1 or 2
              subimage = alt_image[i-1:i+2,j-1:j+2]
              if np.argmax(subimage) == 4 and subimage[1,1] > bg: # 4, 1x1 
                  spot_intensity.append(np.max(subimage))
                  spot_location.append([i-1,j-1])
      return spot_intensity,spot_location
  
   
  def find_spots_3x3(self, image,bg):
      # Construct a sub-image of dimension 3*3 (can be changed if needed)
      row,col = image.shape
      spot_location = []
      spot_intensity = []
      alt_image = np.lib.pad(image, ((1,1),(1,1)), 'constant', constant_values=0) # 2x2
      subimage = np.zeros((3,3)) # 5x5
      for i in range(1,row) :    # 1 or 2
          for j in range(1,col): # 1 or 2
              subimage = alt_image[i-1:i+2,j-1:j+2]
              if np.argmax(subimage) == 4 and subimage[1,1] > bg: # 4, 1x1 
                  spot_intensity.append(np.sum(subimage))
                  spot_location.append([i-1,j-1])
      return spot_intensity,spot_location
  
  
  def find_spots_5x5(self, image,bg):
      # Construct a sub-image of dimension 3*3 (can be changed if needed)
      row,col = image.shape
      spot_location = []
      spot_intensity = []
      alt_image = np.lib.pad(image, ((2,2),(2,2)), 'constant', constant_values=0) # 2x2
      subimage = np.zeros((5,5)) # 5x5
      for i in range(2,row) :
          for j in range(2,col):
  #            subimage = alt_image[i-1:i+2,j-1:j+2]
              subimage = alt_image[i-2:i+3,j-2:j+3]
              if np.argmax(subimage) == 12 and subimage[2,2] > bg: # 4, 1x1 
                  spot_intensity.append(np.sum(subimage))
                  spot_location.append([i-1,j-1])
      return spot_intensity,spot_location
  
  def find_spots(self, image,bg):
      return find_spots_1x1(image,bg)
  


############################---MAIN---#########################
if __name__=='__main__':

  deconvolver=PSF_deconvolve()

  # Main parameters go in here
  epsilon = 0.0 # Slight fudge factor used in FFT method
  fwhm = 76.8 # fwhm of PSF in microns
  lpix = 177.8 # pixel size in microns for LS49, 4x4 binning
  fwhm_pixel = fwhm/lpix # Convert fwhm to pixels
  xpsf = 25 # PSF box shape x-dim
  ypsf = 25 # PSF box shape y-dim
  niter=50 # Number of iterations for Lucy Richardson algorithm
  fft_flag = False # Whether to use FFT deconvolver or Lucy Richardson
  
  # Load up image and raw data here
  ts='20180501143545184' # Good example image with some PSF effects clearly visible


  ls49_data_dir=None
  if ls49_data_dir is None:
    LS49_regression = libtbx.env.find_in_repositories(
      relative_path="LS49_regression",
      test=os.path.isdir)
    if LS49_regression is None:
      raise Sorry('LS49_regression folder needs to be present or else specify ls49_data_dir')
    ls49_data_dir = os.path.join(LS49_regression, 'diffBragg_work', 'iota_r0222_cori', 'rayonix_expt')


  image_filename=os.path.join(ls49_data_dir, 'idx-%s.cbf'%ts)
  img = dxtbx.load(image_filename)
  imgx = img.get_raw_data().as_double() # flex array
  writer=FullCBFWriter(image_filename) 
  cbf_handle=writer.get_cbf_handle()
  Input = imgx
  if (fft_flag):
  #Visualize the PSF
      sizeX=960
      sizeY=960
      #psf = makeGaussPSF(fwhm_pixel,sizeX, sizeY)
      psf = deconvolver.makeGauss_integPSF(fwhm_pixel,sizeX, sizeY)
      #psf = makeMoffat_integPSF(fwhm_pixel, sizeX, sizeY)
  #    exportArrayAsImage(psf, 'PSF.img')
      InputDeconv = np.abs(deconvolver.deconvolve_fft(Input.as_numpy_array(), psf.as_numpy_array(), epsilon))
  #    exportArrayAsImage(fftshift(InputNoisyConv), 'NoisyConvolved.img')
      deconvolver.exportArrayAsImage(InputDeconv, 'Deconvolved_real.img')

  ############################
  #  Lucy-Richardson stuff   #
  ############################
  else:
      psf = deconvolver.makeMoffat_integPSF(fwhm_pixel,xpsf,ypsf)
      #psf = deconvolver.makeGauss_integPSF(fwhm_pixel,xpsf,ypsf)
      print 'deconvolving'
      InputDeconv = ((deconvolver.rl_standard(Input, psf, niter)))
      InputDeconv = InputDeconv.as_numpy_array()
      header=None
      deconvolver.exportArrayAsImage(InputDeconv,'deconvolved_psf_00000.img',header)
      #from IPython import embed; embed(); exit()
      #writer.add_data_to_cbf(cbf_handle, data=flex.double(InputDeconv))
      #cbf_handle.write_widefile(
      #      'deconv.cbf', pycbf.CBF, pycbf.MIME_HEADERS | pycbf.MSG_DIGEST | pycbf.PAD_4K, 0 )
      #writer.write_cbf('deconvolved_psf.cbf', index=0)
