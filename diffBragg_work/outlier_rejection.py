from load_ls49 import strong_spot_mask

def process_ls49_image_real(tstamp='20180501143555114', #tstamp='20180501143559313',
                            Nstrongest = 30,
                            resmax=12.0, resmin=3.0,
                            #mtz_file='5cmv_Iobs.mtz',
                            #mtz_file='anom_ls49_oxy_2.3_t3_gentle_pr_s0_mark0.mtz',
                            mtz_file='anom_ls49_oxy_2.3_unit_pr_lorentz_primeref_m008_s0_mark0.mtz',
                            ls49_data_dir=None):
    import os, pickle, numpy as np
    from scipy.interpolate import interp1d
    import dxtbx
    from dxtbx.model.experiment_list import ExperimentListFactory
    from simtbx.diffBragg import utils
    from dials.array_family import flex
    from iotbx import mtz
    import libtbx.load_env
    from dials.util import Sorry

    if ls49_data_dir is None:
      LS49_regression = libtbx.env.find_in_repositories(
        relative_path="LS49_regression",
        test=os.path.isdir)
      if LS49_regression is None:
        raise Sorry('LS49_regression folder needs to be present or else specify ls49_data_dir')
      ls49_data_dir = os.path.join(LS49_regression, 'diffBragg_work', 'iota_r0222_cori', 'rayonix_expt')

    GAIN = 0.75
    loader = dxtbx.load(os.path.join(ls49_data_dir,'idx-%s.cbf'%tstamp))
    cbf_imageset = loader.get_imageset([os.path.join(ls49_data_dir,'idx-%s.cbf'%tstamp)])
    
    img = loader.get_raw_data().as_numpy_array() / GAIN
    exp_list = ExperimentListFactory.from_json_file(os.path.join(ls49_data_dir,'idx-%s_refined.expt'%tstamp), check_format=False)
    exp = exp_list[0]
    C = exp.crystal
    B = exp.beam
    D = exp.detector
    #refls = flex.reflection_table.from_file('idx-%s_indexed.refl' % tstamp)
    refls = flex.reflection_table.from_file(os.path.join(ls49_data_dir,'idx-%s_integrated.refl'%tstamp))
    Nbefore = len(refls)
    refls = refls.select(flex.bool([resmin < d < resmax for d in refls['d']]))
    print("Kept %d out of %d refls in the res range %2.2f to %2.2f"
          % (len(refls), Nbefore, resmin, resmax))

    bboxes = [list(refls['shoebox'][i].bbox)[:4] for i in range(len(refls)) ]
    resolutions = [D[0].get_resolution_at_pixel(B.get_s0(), refls[i]['xyzobs.px.value'][0:2]) for i in range(len(refls))]
    bboxes = np.array(bboxes)
    bboxes[bboxes > 960] = 959
    bboxes[bboxes < 0] = 0
    mill_idx = [ list(refls['miller_index'][i]) for i in range(len(refls)) ]
    R2 = flex.reflection_table.from_file(os.path.join(ls49_data_dir, 'idx-%s_indexed.refl'%tstamp))
    strong_mask = strong_spot_mask(refl_tbl=R2, img_size=img.shape)
    is_bg_pixel = np.logical_not(strong_mask)
    is_BAD_pixel = np.logical_not(pickle.load(open(os.path.join(ls49_data_dir,'../','mask_r4.pickle'), 'r'))[0].as_numpy_array())
    is_bg_pixel[is_BAD_pixel] = False
    num_spots = len(refls)
    tilt_abc = np.zeros((num_spots, 3))
    tilts = []
    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes):
        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]
        tilt, bgmask, coeff, _ = utils.tilting_plane(shoebox_img, mask=shoebox_mask, zscore=2)
        tilts.append(tilt)
        tilt_abc[i_spot] = (coeff[1], coeff[2], coeff[0])

    fee_file='idx-fee_data_%s.pickle'%tstamp
    chann_lambda, channI = np.array(pickle.load(open(os.path.join(ls49_data_dir,fee_file), 'r'))[tstamp]).T
    I = interp1d(chann_lambda, channI)
    max_energy  = chann_lambda[np.argmax(channI)]
    min_energy_interpol = max_energy - 3
    max_energy_interpol = max_energy + 3
    print ('INTERPOLATION ENERGIES = ', min_energy_interpol, max_energy_interpol)
    interp_energies = np.arange(min_energy_interpol, max_energy_interpol, 0.5)
    interp_fluxes = I(interp_energies)
    interp_fluxes /= interp_fluxes.sum()
    interp_fluxes *= 1000000000000.0
    spectrum = zip(12398.419739640716 / interp_energies, interp_fluxes)


    #M = mtz.object('ls49_oxy_2.5_s0_mark0.mtz')
    #sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs')]

    M = mtz.object(os.path.join(ls49_data_dir,'../',mtz_file))
    #sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs')]
    sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'IMEAN')]
    sfall = sfall.as_amplitude_array()
    return {'dxcrystal': C, 'dxdetector': D, 'dxbeam': B, 'mill_idx': mill_idx, 'data_img': img, 'bboxes_x1x2y1y2': bboxes, 
       'tilt_abc': tilt_abc, 'spectrum': spectrum, 'sfall': sfall,
            'mask': is_BAD_pixel, 'experiment':exp, 'indexed_reflections': R2, 'resolution': resolutions, 'cbf_imageset':cbf_imageset}


def outlier_rejection_ls49(ts=None, ls49_data_dir=None, show_plotted_images=False, params=None):
    from simtbx.diffBragg.refiners import RefineAll
    from simtbx.diffBragg.sim_data import SimData
    import numpy as np
    from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
    from simtbx.diffBragg import nanoBragg_crystal, nanoBragg_beam
    from copy import deepcopy
    from dxtbx.model.experiment_list import Experiment, ExperimentList, ExperimentListFactory
    from dials.array_family import flex
    import pylab as plt

    data = process_ls49_image_real(tstamp=ts,Nstrongest=10, resmin=2.3, resmax=13.5, ls49_data_dir=ls49_data_dir)
    refine_with_psf=True
    plot_images=True # This is a lie. Mostly need this to store model_Lambda for statistics etc

    C = data["dxcrystal"]
    D = data["dxdetector"]
    B = data["dxbeam"]
    exp=data['experiment']

    # Some global variables here for LS49
    mos_spread_deg=0.001
    n_mos_domains=1
    a,b,c,_,_,_ = C.get_unit_cell().parameters()
    Deff = 1000
    Ncells_abc_0 = np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
    n_spots = len(data['tilt_abc'])

    # Set up nbcryst and nbbeam
    nbcryst = nanoBragg_crystal.nanoBragg_crystal()
    nbcryst.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
    nbcryst.mos_spread_deg = mos_spread_deg
    nbcryst.n_mos_domains = n_mos_domains
    nbcryst.thick_mm = 0.005
    nbcryst.miller_array = data["sfall"]
    nbcryst.dxtbx_crystal = C
    nbbeam = nanoBragg_beam.nanoBragg_beam()
    nbbeam.unit_s0 = B.get_unit_s0()
    nbbeam.spectrum = data["spectrum"]
    # Set up simdata
    SIM = SimData()
    SIM.crystal = nbcryst
    SIM.detector = D
    SIM.beam = nbbeam

    SIM.instantiate_diffBragg(adc_offset=0,
                              oversample=0,
                              interpolate=0,
                              verbose=0)
    SIM.D.show_params()

    #SIM.D.spot_scale = 1e6 # to guide your eye
    import pylab as plt
    plt.clf()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
    ax1.imshow([[0, 1, 1], [0, 1, 2]])
    ax2.imshow([[0, 1, 1], [0, 1, 2]])
    #ax3.imshow([[0, 1, 1], [0, 1, 2]])
    scale=1.0
    for i_spot in range(n_spots):
      bbox=data['bboxes_x1x2y1y2'][i_spot]
      a,b,c=data['tilt_abc'][i_spot]
      x1,x2,y1,y2=bbox
      yr, xr = np.indices((y2-y1+1, x2-x1+1))
      tilt_plane=xr*a+yr*b+c

      SIM.D.region_of_interest=((bbox[0], bbox[1]), (bbox[2], bbox[3]))
      SIM.D.add_diffBragg_spots()
      img=SIM.D.raw_pixels.as_numpy_array()
      img=img[y1:y2+1, x1:x2+1]
      #img=tilt_plane+scale*img
      data_img=data['data_img'][y1:y2+1, x1:x2+1]-tilt_plane
      data_img[data_img<0]=0
      m = img[img>1.e-9].mean()
      s = img[img>1.e-9].std()
      vmin = 1
      vmax = m+10*s
      if vmax < vmin: vmax=vmin+10
      m0 = data_img[data_img>1.e-9].mean()
      s0 = data_img[data_img>1.e-9].std()
      #print ('vmax0 %d  %d'%(i_spot, m0+5*s0)) 
      #print ('vmax %d  %d'%(i_spot, m+5*s)) 
      vmin0 = 1
      vmax0 =  m0+10*s0
      normalized_data_img=(data_img-np.mean(data_img))/np.std(data_img)
      normalized_img=(img-np.mean(img))/np.std(img)
    
      summed_I_data=np.sum(data_img)
      summed_I_model=np.sum(img)
      print ('Summed Intensities: Data: %d     Model=%d'%(summed_I_data, summed_I_model))
      diff_normalized_img=normalized_img-normalized_data_img
       
      
      ax1.images[0].set_data(img)
      ax1.images[0].set_clim(vmin, vmax)
      ax2.images[0].set_data(data_img)
      ax2.images[0].set_clim(vmin0, vmax0)
      #ax3.images[0].set_data(diff_normalized_img)
      plt.suptitle("Spot Number %d  ||   I_Model = %d   || I_Data= %d || Res= %.2f"%(i_spot, summed_I_model, summed_I_data, data['resolution'][i_spot]))
      fig.canvas.draw()
      plt.pause(2.5)
   
if __name__ == "__main__":
    # Initial r0222 regression
    timestamps_of_interest = ['20180501143533988', # bad
                              '20180501143546717', # bad --> blows up
                              '20180501143546817', # looks ok
                              '20180501143547150', # does not work
                              '20180501143548650', # blows up
                              '20180501143549416', # does not seem to work
                              '20180501143549949', # Look ok !
                              '20180501143551715', # does not work
                              '20180501143555114', # seems to get 50% of the spots
                              '20180501143559313', # does not work
                              '20180501143602713', # does not work
                              '20180501143606545', # diverges
                              '20180501143620206', # diverges it seems although does not crash 
                              '20180501143625171', # diverges, not sure why ?
                              '20180501143628702', # fails with assertion error in curvatures. Did not look good till then
                              '20180501143628902', # Looks good actually
                              '20180501143631168', # Failes with assertion error in curvatures like above. did not look good
                              '20180501143632300', # Looks good 
                              '20180501143640763', # Looks good
                              '20180501143643462', # meehhh
                              '20180501143643662', # Looks OK
                              '20180501143652325', # curvature assertion error, looked good till then
                              '20180501143701853'] # Does not work

    ts = timestamps_of_interest[8]
    #ls49_data_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/all_files/rayonix_expt'
    #ls49_data_dir='/Users/abhowmick/Desktop/software/dials/modules/LS49_regression/diffBragg_work/jungfrau_grid_search_4_or_more_regression/rayonix_images_4_or_more_spots_r183_255'
    #ts='20180501114703722' # Image used in blog to compare on jungfrau
    #ts='20180501120317142'
    ls49_data_dir=None
    outlier_rejection_ls49(ts=ts, ls49_data_dir=ls49_data_dir)


