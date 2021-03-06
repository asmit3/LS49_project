#
# File courtesy Derek Mendez
# Copying it from cctbx_project/simtbx/diffBragg to make more local edits here instead of in cctbx
#
#


message = """ Refinement strategy 1 that seems to give good results
              Much better than just the 2 stage refinement previously
              1. ncells+s_i
              2. bg+U+B+s_i+S
              3. ncells+s_i
              4. U+B+(bg)"""



def strong_spot_mask(refl_tbl, img_size):
    """note only works for strong spot reflection tables
    img_size is slow-scan, fast-scan"""
    import numpy as np
    from dials.algorithms.shoebox import MaskCode
    Nrefl = len(refl_tbl)
    masks = [ refl_tbl[i]['shoebox'].mask.as_numpy_array() for i in range(Nrefl)
            ]
    code = MaskCode.Foreground.real
    x1, x2, y1, y2, z1, z2 = zip(*[ refl_tbl[i]['shoebox'].bbox for i in range(Nrefl)
                                  ])
    spot_mask = np.zeros(img_size, bool)
    for i1, i2, j1, j2, M in zip(x1, x2, y1, y2, masks):
        slcX = slice(i1, i2, 1)
        slcY = slice(j1, j2, 1)
        spot_mask[(slcY, slcX)] = M & code == code

    return spot_mask


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
    snr = refls['intensity.sum.value'] / flex.sqrt(refls['intensity.sum.variance'])
    order = np.argsort(snr)[::-1]
    min_snr=3.0
    refls=refls.select(snr>min_snr)
    #refls = refls.select(snr > snr[order[Nstrongest]])
    snr2 = refls['intensity.sum.value'] / flex.sqrt(refls['intensity.sum.variance'])

    # Now select only the highest resolution spots upto Nstrongest
    sorted_indices_by_resolution = flex.sort_permutation(refls['d'], reverse=False)
    refls = refls.select(sorted_indices_by_resolution[0:Nstrongest])

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
    min_energy_interpol = max_energy - 15
    max_energy_interpol = max_energy + 15
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


if __name__ == "__main__":
    from argparse import ArgumentParser
    from simtbx.diffBragg.refiners import RefineAll
    from simtbx.diffBragg.sim_data import SimData
    import numpy as np
    from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
    from simtbx.diffBragg import nanoBragg_crystal, nanoBragg_beam
    from copy import deepcopy
    from dxtbx.model.experiment_list import Experiment, ExperimentList, ExperimentListFactory
    from dials.array_family import flex

    parser = ArgumentParser()
    parser.add_argument("--plot", action='store_true')
    parser.add_argument("--scaleonly", action='store_true')
    args = parser.parse_args()


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

    ts = timestamps_of_interest[-2]
    #ls49_data_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/all_files/rayonix_expt'
    #ls49_data_dir='/Users/abhowmick/Desktop/software/dials/modules/LS49_regression/diffBragg_work/jungfrau_grid_search_4_or_more_regression/rayonix_images_4_or_more_spots_r183_255'
    #ts='20180501114703722'
    #ts='20180501120317142'
    ls49_data_dir=None
    data = process_ls49_image_real(tstamp=ts,Nstrongest=10, resmin=2.3, resmax=13.5, ls49_data_dir=ls49_data_dir)

    refine_with_psf=True

    C = data["dxcrystal"]
    D = data["dxdetector"]
    B = data["dxbeam"]
    exp=data['experiment']
    indexed_reflections = deepcopy(data['indexed_reflections'])
    dump_exp = Experiment(imageset=data['cbf_imageset'], 
                         beam=B,
                         detector=D,
                         goniometer=exp.goniometer,
                         scan=exp.scan,
                         crystal=C)
    dump_explist = ExperimentList([dump_exp])
    dump_explist.as_file('before_refinement_%s.expt'%ts)
    indexed_reflections.as_file('before_refinement_%s.refl'%ts)

    # Some global variables here for LS49
    mos_spread_deg=0.001
    n_mos_domains=1
    

    a,b,c,_,_,_ = C.get_unit_cell().parameters()
    Deff = 1000
    Ncells_abc = np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
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

    if False:
      D_jung = ExperimentListFactory.from_json_file('../jungfrau_scripts/fake_jungfrau_from_rayonix.expt', check_format=False)[0].detector
      #D_jung = ExperimentListFactory.from_json_file('../jungfrau_scripts/rotated_plus_90.json', check_format=False)[0].detector

      nbcryst_final = nanoBragg_crystal.nanoBragg_crystal()
      nbcryst_final.Ncells_abc = 5,5,5 #Ncells_abc, Ncells_abc, Ncells_abc
      nbcryst_final.mos_spread_deg = mos_spread_deg
      nbcryst_final.n_mos_domains = n_mos_domains
      nbcryst_final.thick_mm = 0.005
      nbcryst_final.miller_array = data["sfall"]
      nbcryst_final.dxtbx_crystal = C

      SIM_jung = SimData()
      SIM_jung.crystal = nbcryst
      SIM_jung.detector = D #D_jung
      SIM_jung.beam = nbbeam

      SIM_jung.instantiate_diffBragg(adc_offset=0,
                                oversample=0,
                                interpolate=0,
                                verbose=0)
      SIM_jung.D.show_params()
      SIM_jung.D.spot_scale = 1e6 # to guide your eye
      SIM_jung.panel_id = 0
      SIM_jung.D.add_diffBragg_spots()
      img=SIM_jung.D.raw_pixels.as_numpy_array()
      import matplotlib.pyplot as plt
      m = img[img>1.e-9].mean()
      s = img[img>1.e-9].std()
      vmin = 1
      vmax = 1000 #m+5*s
      plt.imshow(img, vmax=vmax, vmin=vmin)
      plt.show()
      plt.imshow(img, vmax=vmax, vmin=vmin)
      plt.show()
      exit()


    SIM = SimData()
    # maybe make Ncells very small e.g. 5x5x5
    #nbcryst.Ncells_abc (5,5,5)
    SIM.crystal = nbcryst
    SIM.detector = D
    SIM.beam = nbbeam

    SIM.instantiate_diffBragg(adc_offset=0,
                              oversample=0,
                              interpolate=0,
                              verbose=0)
    SIM.D.show_params()
    #SIM.D.spot_scale = 1e6 # to guide your eye
    #SIM.panel_id = 0
    #SIM.D.add_diffBragg_spots()
    #img=SIM.D.raw_pixels.as_numopy_array()

    ucell = C.get_unit_cell().parameters()
    UcellMan = MonoclinicManager(a=ucell[0], b=ucell[1], c=ucell[2], beta=ucell[4]*np.pi / 180.)
    import pylab as plt
    if args.plot:
        I = data["data_img"].copy()
        M = ~data["mask"]
        I*=M
        m = I[I > 0].mean()
        s = I[I > 0].std()
        vmin = m-s
        vmax=m+2.5*s

        plt.imshow(I, vmin=vmin, vmax=vmax)
        #for x1, x2, y1, y2 in data["bboxes_x1x2y1y2"]:
        #    patch = plt.Rectangle(
        #        width=x2-x1,
        #        height=y2-y1,
        #        xy=(x1, y1),
        #        ec='r', fc='none')
        #    plt.gca().add_patch(patch)
        #plt.show()
    n_spots = len(data['tilt_abc'])
    init_local_spotscale = flex.double([1.0]*n_spots)
    RUC = RefineAll(
        spot_rois=data["bboxes_x1x2y1y2"],
        spot_resolution=data['resolution'],
        abc_init=data["tilt_abc"],
        img=data["data_img"],
        SimData_instance=SIM,
        plot_images=args.plot,
        ucell_manager=UcellMan,
        init_gain=1,
        init_scale=1.0,
        init_local_spotscale = init_local_spotscale)
    RUC.trad_conv = True
    RUC.trad_conv_eps = 1e-5
    RUC.max_calls = 250
    RUC.refine_background_planes = False
    RUC.refine_Umatrix = False
    RUC.refine_Bmatrix = False
    RUC.verbose=True
    RUC.plot_stride=1
    RUC.refine_ncells=False #
    RUC.refine_crystal_scale = False # scale factor
    RUC.refine_gain_fac = False
    RUC.refine_detdist=False
    RUC.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
    RUC.use_curvatures=False #
    RUC.calc_curvatures=True # if set to False, never uses curvatures
    RUC.refine_with_restraints=False
    RUC.refine_with_psf=refine_with_psf
    #RUC.run()
    #if RUC.hit_break_to_use_curvatures:
    #  RUC.num_positive_curvatures=0
    #  RUC.use_curvatures=True
    #  RUC.run(setup=False) # Now it will use curvatures

    # First refine scale factors and ncells, then refine everything else
    RUC.refine_ncells=True
    RUC.refine_crystal_scale=False
    RUC.refine_local_spotscale=True

    RUC.run()
    refined_ncells = RUC.x[-4]
    refined_scale = RUC.x[-1]
    refined_gain = RUC.x[-2]
    refined_local_spotscale = RUC.x[RUC.n_bg:RUC.n_bg+RUC.n_spots]
    refined_tilt_abc=[]
    for i in range(RUC.n_spots):
      refined_tilt_abc.append([RUC.x[i], RUC.x[RUC.n_spots+i], RUC.x[2*RUC.n_spots+i]])
    refined_tilt_abc=np.array(refined_tilt_abc)

    # Refinement analysis: use a linear correlation and linear regression model ?
    if True:
      all_x = flex.double()
      all_y = flex.double()
      for i_spot in range(RUC.n_spots):
        x = RUC.store_model_Lambda[i_spot]
        y = RUC.store_Imeas[i_spot]
        all_x.extend(flex.double(x.flatten()))
        all_y.extend(flex.double(y.flatten()))
      z_regression=flex.linear_regression(all_x, all_y)
      z_correlation=flex.linear_correlation(all_x, all_y)
      print ('--------- Linear Regression Summary----------')
      z_regression.show_summary()
      print ('----------------------------------------------')
      print ('--------- Linear Correlation Summary----------')
      z_correlation.show_summary()


    if False:
      for i_spot in range(RUC.n_spots):
        fig, axs = plt.subplots(1,2)
        axs[0].imshow([[0, 1, 1], [0, 1, 2]])
        axs[1].imshow([[0, 1, 1], [0, 1, 2]])
        x = RUC.store_model_Lambda[i_spot]
        y = RUC.store_Imeas[i_spot]
        vmin = RUC.store_vmin[i_spot]
        vmax = RUC.store_vmax[i_spot]
        axs[0].images[0].set_data(x)
        axs[1].images[0].set_data(y)
        axs[0].images[0].set_clim(vmin, vmax)
        axs[1].images[0].set_clim(vmin, vmax)
        plt.suptitle("Spot number = %d"%i_spot)
      plt.show()

    #from IPython import embed; embed(); exit() 
    # Now set the other stuff to refine
    print (' ================== Stage 2 : Now refining all the variables ======================= ')
    Ncells_abc = refined_ncells #np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
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

    SIM2 = SimData()
    SIM2.crystal = nbcryst
    SIM2.detector = D
    SIM2.beam = nbbeam

    SIM2.instantiate_diffBragg(adc_offset=0,
                              oversample=0,
                              interpolate=0,
                              verbose=0)
    SIM2.D.show_params()
    RUC2 = RefineAll(
        spot_rois=data["bboxes_x1x2y1y2"],
        spot_resolution=data['resolution'],
        abc_init=refined_tilt_abc,
        img=data["data_img"],
        SimData_instance=SIM2,
        plot_images=args.plot,
        ucell_manager=UcellMan,
        init_gain=refined_gain,
        init_scale=refined_scale,
        init_local_spotscale=refined_local_spotscale)
    RUC2.trad_conv = True
    RUC2.trad_conv_eps = 1e-5
    RUC2.max_calls = 250
    RUC2.refine_background_planes = True
    RUC2.refine_Umatrix = True
    RUC2.refine_Bmatrix = True
    RUC2.verbose=True
    RUC2.plot_stride=1
    RUC2.refine_ncells=False #
    RUC2.refine_crystal_scale = True # scale factor
    RUC2.refine_gain_fac = False
    RUC2.refine_detdist=False
    RUC2.refine_local_spotscale=True
    RUC2.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
    RUC2.use_curvatures=False # using curvatures asmit
    RUC2.calc_curvatures=True # if set to False, never uses curvatures
    RUC2.refine_with_psf=refine_with_psf 
    RUC2.run()

    print("Done.")
    print("Refined scale =%f", RUC.x[-1])
    # Refinement analysis: use a linear correlation and linear regression model ?
    if True:
      all_x = flex.double()
      all_y = flex.double()
      for i_spot in range(RUC2.n_spots):
        x = RUC2.store_model_Lambda[i_spot]
        y = RUC2.store_Imeas[i_spot]
        all_x.extend(flex.double(x.flatten()))
        all_y.extend(flex.double(y.flatten()))
      z_regression=flex.linear_regression(all_x, all_y)
      z_correlation=flex.linear_correlation(all_x, all_y)
      print ('--------- Linear Regression Summary----------')
      z_regression.show_summary()
      print ('----------------------------------------------')
      print ('--------- Linear Correlation Summary----------')
      z_correlation.show_summary()
    show_pixel_values=False
    if False:
      for i_spot in range(RUC2.n_spots):
        fig, axs = plt.subplots(4,2)
        axs[0][0].imshow([[0, 1, 1], [0, 1, 2]])
        axs[1][0].imshow([[0, 1, 1], [0, 1, 2]])
        axs[0][1].imshow([[0, 1, 1], [0, 1, 2]])
        axs[1][1].imshow([[0, 1, 1], [0, 1, 2]])
        axs[2][0].imshow([[0, 1, 1], [0, 1, 2]])
        axs[2][1].imshow([[0, 1, 1], [0, 1, 2]])
        axs[3][0].imshow([[0, 1, 1], [0, 1, 2]]) # before
        axs[3][1].imshow([[0, 1, 1], [0, 1, 2]]) # after

        x = RUC2.store_model_Lambda[i_spot]
        y = RUC2.store_Imeas[i_spot]
        x1 = RUC.store_model_Lambda[i_spot]
        y1 = RUC.store_Imeas[i_spot]
        x0 = RUC.store_init_model_Lambda[i_spot]
        y0 = RUC.store_init_Imeas[i_spot]
        deltaI0 = np.abs(x0-y0)
        deltaI = np.abs(x-y)
      
        vmin = RUC2.store_vmin[i_spot]
        vmax = RUC2.store_vmax[i_spot]
        vmin1 = RUC.store_vmin[i_spot]
        vmax1 = RUC.store_vmax[i_spot]
        vmin0 = RUC.store_init_vmin[i_spot]
        vmax0 = RUC.store_init_vmax[i_spot]
        axs[0][0].images[0].set_data(x0)
        axs[1][0].images[0].set_data(x1)
        axs[2][0].images[0].set_data(x)
        axs[3][0].images[0].set_data(deltaI0)

        axs[0][1].images[0].set_data(y0)
        axs[1][1].images[0].set_data(y1)
        axs[2][1].images[0].set_data(y)
        axs[3][1].images[0].set_data(deltaI)

        axs[0][0].images[0].set_clim(vmin0, vmax0)
        axs[1][0].images[0].set_clim(vmin1, vmax1)
        axs[2][0].images[0].set_clim(vmin, vmax)
        # Stuff just for deltaI
        m = deltaI0[deltaI0 > 1e-9].mean()
        s = deltaI0[deltaI0 > 1e-9].std()
        dvmax = m+5*s
        dvmin = m-s
        axs[3][0].images[0].set_clim(dvmin, dvmax)
        m = deltaI[deltaI > 1e-9].mean()
        s = deltaI[deltaI > 1e-9].std()
        dvmax = m+5*s
        dvmin = m-s
        axs[3][1].images[0].set_clim(dvmin, dvmax)

        axs[0][1].images[0].set_clim(vmin0, vmax0)
        axs[1][1].images[0].set_clim(vmin1, vmax1)
        axs[2][1].images[0].set_clim(vmin, vmax)

        axs[0][0].set_title('Before Refinement: Calc')
        axs[1][0].set_title('Stage 1 Refinement: Calc')
        axs[2][0].set_title('Final Refinement: Calc')
        axs[3][0].set_title('Initial difference image')

        axs[0][1].set_title('Observation')
        axs[1][1].set_title('Observation')
        axs[2][1].set_title('Observation')
        axs[3][1].set_title('Final difference image')

        plt.suptitle("Spot number = %d at %.2f resolution"%(i_spot,RUC2.spot_resolution[i_spot]))
        if show_pixel_values:
          ds=3
          df=3
          smax=np.argmax(x0)//x0.shape[0]
          fmax=np.argmax(x0)%x0.shape[0]
          if smax-ds < 0 or smax +ds >= x0.shape[0]:
            ds=0
          if fmax-df <0 or fmax+df >= x0.shape[1]:
            df=0
          for s in range(smax-ds, smax+ds):
            for f in range(fmax-df, fmax+df):
              text=axs[0][0].text(f,s, int(x0[s,f]), ha="center", va="center", color="w")
      plt.show()

    best=RUC2.best_image
    C2 = deepcopy(C)
    ang, ax = RUC2.get_correction_misset(as_axis_angle_deg=True)
    C2.rotate_around_origin(ax, ang)
    C2.set_B(RUC2.get_refined_Bmatrix())
      

    # refined unit cell parameters
    ucell_ref = C2.get_unit_cell().parameters()

    print("ground truth unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell)
    print("refined unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell_ref)
    print("")

    C2.show()
    dump_exp = Experiment(imageset=data['cbf_imageset'], 
                         beam=B,
                         detector=D,
                         goniometer=exp.goniometer,
                         scan=exp.scan,
                         crystal=C2)
    dump_explist = ExperimentList([dump_exp])
    dump_explist.as_file('after_refinement_%s.expt'%ts)
    # Dump refl file as well based on prediction from refined model
    if True:
      from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
      ref_predictor = ExperimentsPredictorFactory.from_experiments(
                      dump_explist,
                      force_stills=True,
                      spherical_relp=False)
      ref_predictor(indexed_reflections)
      indexed_reflections.as_file('after_refinement_%s.refl'%ts)


########################################################################################
    if True:
      print (' ================== Stage 3 refinement ======================= ')
      refined_ncells = RUC2.x[-4]
      refined_scale = RUC2.x[-1]
      refined_gain = RUC2.x[-2]
      refined_local_spotscale = RUC2.x[RUC2.n_bg:RUC2.n_bg+RUC2.n_spots]
      refined_tilt_abc=[]
      for i in range(RUC2.n_spots):
        refined_tilt_abc.append([RUC2.x[i], RUC2.x[RUC2.n_spots+i], RUC2.x[2*RUC2.n_spots+i]])
      refined_tilt_abc=np.array(refined_tilt_abc)

      Ncells_abc = refined_ncells #np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
      nbcryst = nanoBragg_crystal.nanoBragg_crystal()
      nbcryst.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
      nbcryst.mos_spread_deg = mos_spread_deg
      nbcryst.n_mos_domains = n_mos_domains
      nbcryst.thick_mm = 0.005
      nbcryst.miller_array = data["sfall"]
      nbcryst.dxtbx_crystal = C2

      nbbeam = nanoBragg_beam.nanoBragg_beam()
      nbbeam.unit_s0 = B.get_unit_s0()
      nbbeam.spectrum = data["spectrum"]

      SIM3 = SimData()
      SIM3.crystal = nbcryst
      SIM3.detector = D
      SIM3.beam = nbbeam

      SIM3.instantiate_diffBragg(adc_offset=0,
                                oversample=0,
                                interpolate=0,
                                verbose=0)
      SIM3.D.show_params()
      ucell3 = C2.get_unit_cell().parameters()
      UcellMan3 = MonoclinicManager(a=ucell3[0], b=ucell3[1], c=ucell3[2], beta=ucell3[4]*np.pi / 180.)
      RUC3 = RefineAll(
          spot_rois=data["bboxes_x1x2y1y2"],
          spot_resolution=data['resolution'],
          abc_init=refined_tilt_abc,
          img=data["data_img"],
          SimData_instance=SIM3,
          plot_images=args.plot,
          ucell_manager=UcellMan3,
          init_gain=refined_gain,
          init_scale=refined_scale,
          init_local_spotscale=refined_local_spotscale)
      RUC3.trad_conv = True
      RUC3.trad_conv_eps = 1e-5
      RUC3.max_calls = 250
      RUC3.refine_background_planes = False
      RUC3.refine_Umatrix = False
      RUC3.refine_Bmatrix = False
      RUC3.verbose=True
      RUC3.plot_stride=1
      RUC3.refine_ncells=True #
      RUC3.refine_crystal_scale = False # scale factor
      RUC3.refine_gain_fac = False
      RUC3.refine_detdist=False
      RUC3.refine_local_spotscale=True
      RUC3.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
      RUC3.use_curvatures=False # using curvatures asmit
      RUC3.calc_curvatures=True # if set to False, never uses curvatures
      RUC3.refine_with_psf=refine_with_psf 
      RUC3.run()

      print("Done.")
      print("Refined scale =%f", RUC3.x[-1])
      # Refinement analysis: use a linear correlation and linear regression model ?
      if True:
        all_x = flex.double()
        all_y = flex.double()
        for i_spot in range(RUC3.n_spots):
          x = RUC3.store_model_Lambda[i_spot]
          y = RUC3.store_Imeas[i_spot]
          all_x.extend(flex.double(x.flatten()))
          all_y.extend(flex.double(y.flatten()))
        z_regression=flex.linear_regression(all_x, all_y)
        z_correlation=flex.linear_correlation(all_x, all_y)
        print ('--------- Linear Regression Summary----------')
        z_regression.show_summary()
        print ('----------------------------------------------')
        print ('--------- Linear Correlation Summary----------')
        z_correlation.show_summary()

      if False:
        for i_spot in range(RUC3.n_spots):
          fig, axs = plt.subplots(4,2)
          axs[0][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[1][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[0][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[1][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[2][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[2][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[3][0].imshow([[0, 1, 1], [0, 1, 2]]) # before
          axs[3][1].imshow([[0, 1, 1], [0, 1, 2]]) # after

          x = RUC3.store_model_Lambda[i_spot]
          y = RUC3.store_Imeas[i_spot]
          x1 = RUC2.store_model_Lambda[i_spot]
          y1 = RUC2.store_Imeas[i_spot]
          x0 = RUC.store_init_model_Lambda[i_spot]
          y0 = RUC.store_init_Imeas[i_spot]
          deltaI0 = np.abs(x0-y0)
          deltaI = np.abs(x-y)
        
          vmin = RUC3.store_vmin[i_spot]
          vmax = RUC3.store_vmax[i_spot]
          vmin1 = RUC2.store_vmin[i_spot]
          vmax1 = RUC2.store_vmax[i_spot]
          vmin0 = RUC.store_init_vmin[i_spot]
          vmax0 = RUC.store_init_vmax[i_spot]
          axs[0][0].images[0].set_data(x0)
          axs[1][0].images[0].set_data(x1)
          axs[2][0].images[0].set_data(x)
          axs[3][0].images[0].set_data(deltaI0)

          axs[0][1].images[0].set_data(y0)
          axs[1][1].images[0].set_data(y1)
          axs[2][1].images[0].set_data(y)
          axs[3][1].images[0].set_data(deltaI)

          axs[0][0].images[0].set_clim(vmin0, vmax0)
          axs[1][0].images[0].set_clim(vmin1, vmax1)
          axs[2][0].images[0].set_clim(vmin, vmax)
          # Stuff just for deltaI
          m = deltaI0[deltaI0 > 1e-9].mean()
          s = deltaI0[deltaI0 > 1e-9].std()
          dvmax = m+5*s
          dvmin = m-s
          axs[3][0].images[0].set_clim(dvmin, dvmax)
          m = deltaI[deltaI > 1e-9].mean()
          s = deltaI[deltaI > 1e-9].std()
          dvmax = m+5*s
          dvmin = m-s
          axs[3][1].images[0].set_clim(dvmin, dvmax)

          axs[0][1].images[0].set_clim(vmin0, vmax0)
          axs[1][1].images[0].set_clim(vmin1, vmax1)
          axs[2][1].images[0].set_clim(vmin, vmax)

          axs[0][0].set_title('Before Refinement: Calc')
          axs[1][0].set_title('Stage 1 Refinement: Calc')
          axs[2][0].set_title('Final Refinement: Calc')
          axs[3][0].set_title('Initial difference image')

          axs[0][1].set_title('Observation')
          axs[1][1].set_title('Observation')
          axs[2][1].set_title('Observation')
          axs[3][1].set_title('Final difference image')

          plt.suptitle("Spot number = %d at %.2f resolution"%(i_spot,RUC3.spot_resolution[i_spot]))
          if show_pixel_values:
            ds=3
            df=3
            smax=np.argmax(x0)//x0.shape[0]
            fmax=np.argmax(x0)%x0.shape[0]
            if smax-ds < 0 or smax +ds >= x0.shape[0]:
              ds=0
            if fmax-df <0 or fmax+df >= x0.shape[1]:
              df=0
            for s in range(smax-ds, smax+ds):
              for f in range(fmax-df, fmax+df):
                text=axs[0][0].text(f,s, int(x0[s,f]), ha="center", va="center", color="w")
        plt.show()

      C3 = deepcopy(C2)
      if RUC.refine_Umatrix:
        ang, ax = RUC3.get_correction_misset(as_axis_angle_deg=True)
        C3.rotate_around_origin(ax, ang)
        C3.set_B(RUC3.get_refined_Bmatrix())

####################################################################################



########################################################################################
    if True:
      print (' ================== Stage 4 refinement ======================= ')
      refined_ncells = RUC3.x[-4]
      refined_scale = RUC3.x[-1]
      refined_gain = RUC3.x[-2]
      refined_local_spotscale = RUC3.x[RUC3.n_bg:RUC3.n_bg+RUC3.n_spots]

      refined_tilt_abc=[]
      for i in range(RUC3.n_spots):
        refined_tilt_abc.append([RUC3.x[i], RUC3.x[RUC3.n_spots+i], RUC3.x[2*RUC3.n_spots+i]])
      refined_tilt_abc=np.array(refined_tilt_abc)

      Ncells_abc = refined_ncells #np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
      nbcryst = nanoBragg_crystal.nanoBragg_crystal()
      nbcryst.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
      nbcryst.mos_spread_deg = mos_spread_deg
      nbcryst.n_mos_domains = n_mos_domains
      nbcryst.thick_mm = 0.005
      nbcryst.miller_array = data["sfall"]
      nbcryst.dxtbx_crystal = C3

      nbbeam = nanoBragg_beam.nanoBragg_beam()
      nbbeam.unit_s0 = B.get_unit_s0()
      nbbeam.spectrum = data["spectrum"]

      SIM4 = SimData()
      SIM4.crystal = nbcryst
      SIM4.detector = D
      SIM4.beam = nbbeam

      SIM4.instantiate_diffBragg(adc_offset=0,
                                oversample=0,
                                interpolate=0,
                                verbose=0)
      SIM4.D.show_params()
      ucell4 = C3.get_unit_cell().parameters()
      UcellMan4 = MonoclinicManager(a=ucell4[0], b=ucell4[1], c=ucell4[2], beta=ucell4[4]*np.pi / 180.)
      RUC4 = RefineAll(
          spot_rois=data["bboxes_x1x2y1y2"],
          spot_resolution=data['resolution'],
          abc_init=refined_tilt_abc,
          img=data["data_img"],
          SimData_instance=SIM4,
          plot_images=args.plot,
          ucell_manager=UcellMan4,
          init_gain=refined_gain,
          init_scale=refined_scale,
          init_local_spotscale=refined_local_spotscale)
      RUC4.trad_conv = True
      RUC4.trad_conv_eps = 1e-5
      RUC4.max_calls = 250
      RUC4.refine_background_planes = False
      RUC4.refine_Umatrix = True
      RUC4.refine_Bmatrix = True
      RUC4.verbose=True
      RUC4.plot_stride=1
      RUC4.refine_ncells=False #
      RUC4.refine_crystal_scale = False # scale factor
      RUC4.refine_gain_fac = False
      RUC4.refine_detdist=False
      RUC4.refine_local_spotscale=False
      RUC4.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
      RUC4.use_curvatures=False # using curvatures asmit
      RUC4.calc_curvatures=True # if set to False, never uses curvatures
      RUC4.refine_with_psf=refine_with_psf 
      RUC4.run()

      print("Done.")
      print("Refined scale =%f", RUC4.x[-1])
      # Refinement analysis: use a linear correlation and linear regression model ?
      if True:
        all_x = flex.double()
        all_y = flex.double()
        for i_spot in range(RUC4.n_spots):
          x = RUC4.store_model_Lambda[i_spot]
          y = RUC4.store_Imeas[i_spot]
          all_x.extend(flex.double(x.flatten()))
          all_y.extend(flex.double(y.flatten()))
        z_regression=flex.linear_regression(all_x, all_y)
        z_correlation=flex.linear_correlation(all_x, all_y)
        print ('--------- Linear Regression Summary----------')
        z_regression.show_summary()
        print ('----------------------------------------------')
        print ('--------- Linear Correlation Summary----------')
        z_correlation.show_summary()

      if True:
        for i_spot in range(RUC4.n_spots):
          fig, axs = plt.subplots(4,2)
          axs[0][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[1][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[0][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[1][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[2][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[2][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[3][0].imshow([[0, 1, 1], [0, 1, 2]]) # before
          axs[3][1].imshow([[0, 1, 1], [0, 1, 2]]) # after

          x = RUC4.store_model_Lambda[i_spot]
          y = RUC4.store_Imeas[i_spot]
          x1 = RUC2.store_model_Lambda[i_spot]
          y1 = RUC2.store_Imeas[i_spot]
          x0 = RUC.store_init_model_Lambda[i_spot]
          y0 = RUC.store_init_Imeas[i_spot]
          deltaI0 = np.abs(x0-y0)
          deltaI = np.abs(x-y)
        
          vmin = RUC4.store_vmin[i_spot]
          vmax = RUC4.store_vmax[i_spot]
          vmin1 = RUC2.store_vmin[i_spot]
          vmax1 = RUC2.store_vmax[i_spot]
          vmin0 = RUC.store_init_vmin[i_spot]
          vmax0 = RUC.store_init_vmax[i_spot]
          axs[0][0].images[0].set_data(x0)
          axs[1][0].images[0].set_data(x1)
          axs[2][0].images[0].set_data(x)
          axs[3][0].images[0].set_data(deltaI0)

          axs[0][1].images[0].set_data(y0)
          axs[1][1].images[0].set_data(y1)
          axs[2][1].images[0].set_data(y)
          axs[3][1].images[0].set_data(deltaI)

          axs[0][0].images[0].set_clim(vmin0, vmax0)
          axs[1][0].images[0].set_clim(vmin1, vmax1)
          axs[2][0].images[0].set_clim(vmin, vmax)
          # Stuff just for deltaI
          m = deltaI0[deltaI0 > 1e-9].mean()
          s = deltaI0[deltaI0 > 1e-9].std()
          dvmax = m+5*s
          dvmin = m-s
          axs[3][0].images[0].set_clim(dvmin, dvmax)
          m = deltaI[deltaI > 1e-9].mean()
          s = deltaI[deltaI > 1e-9].std()
          dvmax = m+5*s
          dvmin = m-s
          axs[3][1].images[0].set_clim(dvmin, dvmax)

          axs[0][1].images[0].set_clim(vmin0, vmax0)
          axs[1][1].images[0].set_clim(vmin1, vmax1)
          axs[2][1].images[0].set_clim(vmin, vmax)

          axs[0][0].set_title('Before Refinement: Calc')
          axs[1][0].set_title('Stage 1 Refinement: Calc')
          axs[2][0].set_title('Final Refinement: Calc')
          axs[3][0].set_title('Initial difference image')

          axs[0][1].set_title('Observation')
          axs[1][1].set_title('Observation')
          axs[2][1].set_title('Observation')
          axs[3][1].set_title('Final difference image')

          plt.suptitle("Spot number = %d at %.2f resolution"%(i_spot,RUC4.spot_resolution[i_spot]))
          if show_pixel_values:
            ds=3
            df=3
            smax=np.argmax(x0)//x0.shape[0]
            fmax=np.argmax(x0)%x0.shape[0]
            if smax-ds < 0 or smax +ds >= x0.shape[0]:
              ds=0
            if fmax-df <0 or fmax+df >= x0.shape[1]:
              df=0
            for s in range(smax-ds, smax+ds):
              for f in range(fmax-df, fmax+df):
                text=axs[0][0].text(f,s, int(x0[s,f]), ha="center", va="center", color="w")
        plt.show()

      C4 = deepcopy(C)
      ang, ax = RUC4.get_correction_misset(as_axis_angle_deg=True)
      C4.rotate_around_origin(ax, ang)
      C4.set_B(RUC3.get_refined_Bmatrix())

####################################################################################




    # Now display prediction on jungfrau
    if False:
      D_jung = ExperimentListFactory.from_json_file('../jungfrau_scripts/rotated_plus_90.json', check_format=False)[0].detector

      nbcryst_final = nanoBragg_crystal.nanoBragg_crystal()
      nbcryst_final.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
      nbcryst_final.mos_spread_deg = mos_spread_deg
      nbcryst_final.n_mos_domains = n_mos_domains
      nbcryst_final.thick_mm = 0.005
      nbcryst_final.miller_array = data["sfall"]
      nbcryst_final.dxtbx_crystal = C2

      SIM_jung = SimData()
      SIM_jung.crystal = nbcryst_final
      SIM_jung.detector = D #D_jung
      SIM_jung.beam = nbbeam

      SIM_jung.instantiate_diffBragg(adc_offset=0,
                                oversample=0,
                                interpolate=0,
                                verbose=0)
      SIM_jung.D.show_params()
      SIM_jung.D.spot_scale = 1e6 # to guide your eye
      SIM_jung.panel_id = 0
      SIM_jung.D.add_diffBragg_spots()
      img=SIM_jung.D.raw_pixels.as_numpy_array()
      import matplotlib.pyplot as plt
      m = img[img>1.e-9].mean()
      s = img[img>1.e-9].std()
      vmin = 1
      vmax = 1000 #m+5*s
      plt.imshow(img, vmax=vmax, vmin=vmin)
      plt.show()
      plt.imshow(img, vmax=vmax, vmin=vmin)
      plt.show()
      exit()
    exit()

### Error message
## 14. Assertion error
'''
Traceback (most recent call last):
  File "load_ls49.py", line 246, in <module>
    RUC.run(setup=False) # Now it will use curvatures
  File "/Users/abhowmick/Desktop/software/dials/modules/cctbx_project/simtbx/diffBragg/refiners/pixel_refinement.py", line 224, in run
    verbose=curvature_min_verbose)
  File "/Users/abhowmick/Desktop/software/dials/modules/cctbx_project/scitbx/lbfgs/tst_curvatures.py", line 39, in lbfgs_run
    requests_diag=requests_diag)
  File "/Users/abhowmick/Desktop/software/dials/modules/cctbx_project/scitbx/lbfgs/tst_curvatures.py", line 116, in __call__
    self._verify_diag()
  File "/Users/abhowmick/Desktop/software/dials/modules/cctbx_project/scitbx/lbfgs/tst_curvatures.py", line 127, in _verify_diag
    assert self.d.select(sel).all_gt(0)
AssertionError
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

'''


