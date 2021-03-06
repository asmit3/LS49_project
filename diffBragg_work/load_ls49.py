#
# File courtesy Derek Mendez
# Copying it from cctbx_project/simtbx/diffBragg to make more local edits here instead of in cctbx
#
#
def process_reference(reference):
    """Load the reference spots.
       Code taken from dials.stills_process"""
    if reference is None:
        return None, None
    assert "miller_index" in reference
    assert "id" in reference
    mask = reference.get_flags(reference.flags.indexed)
    rubbish = reference.select(~mask)
    if mask.count(False) > 0:
        reference.del_selected(~mask)
        logger.info(" removing %d unindexed reflections" % mask.count(True))
    if len(reference) == 0:
        raise Sorry(
            """
    Invalid input for reference reflections.
    Expected > %d indexed spots, got %d
  """
            % (0, len(reference))
        )
    mask = reference["miller_index"] == (0, 0, 0)
    if mask.count(True) > 0:
        rubbish.extend(reference.select(mask))
        reference.del_selected(mask)
        logger.info(" removing %d reflections with hkl (0,0,0)" % mask.count(True))
    mask = reference["id"] < 0
    if mask.count(True) > 0:
        raise Sorry(
            """
    Invalid input for reference reflections.
    %d reference spots have an invalid experiment id
  """
            % mask.count(True)
        )
    return reference, rubbish

def integrate(experiments, indexed, stills_process_params):
    """ Code taken from dials.stills_process"""
    indexed, _ = process_reference(indexed)

    # Get the integrator from the input parameters
    from dials.algorithms.profile_model.factory import ProfileModelFactory
    from dials.algorithms.integration.integrator import IntegratorFactory
    from dials.array_family import flex
    from dxtbx.model.experiment_list import ExperimentListFactory

    # Compute the profile model
    # Predict the reflections
    # Match the predictions with the reference
    # Create the integrator
    experiments = ProfileModelFactory.create(stills_process_params, experiments, indexed)
    predicted = flex.reflection_table.from_predictions_multi(
        experiments,
        dmin=stills_process_params.prediction.d_min,
        dmax=stills_process_params.prediction.d_max,
        margin=stills_process_params.prediction.margin,
        force_static=stills_process_params.prediction.force_static,
    )
    predicted.match_with_reference(indexed)
    integrator = IntegratorFactory.create(stills_process_params, experiments, predicted)

    # Integrate the reflections
    integrated = integrator.integrate()
    return experiments, integrated

    if False:
        # Dump experiments to disk
        experiments.as_file(integrated_experiments_filename)
        # Save the reflections
        integrated.as_file(integrated_reflections_filename)
        # Now 

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
                            #mtz_file='anom_ls49_oxy_2.1_unit_pr_lorentz_fromsampleresid_m008_s0_mark0.mtz',
                            mtz_file='anom_ls49_oxy_2.1_gentle_doubleloop_pr_lorentz_m021_s0_mark0.mtz',
                            #mtz_file='anom_ls49_oxy_2.1_unit_pr_lorentz_double_primeref_m008_s0_mark0.mtz',
                            #mtz_file='anom_ls49_oxy_2.3_unit_pr_lorentz_primeref_m008_s0_mark0.mtz',
                            outlier_with_diffBragg=True,
                            ls49_data_dir=None,
                            outdir=None):
    import os, pickle, numpy as np
    from scipy.interpolate import interp1d
    import dxtbx
    from dxtbx.model.experiment_list import ExperimentListFactory
    from simtbx.diffBragg import utils
    from dials.array_family import flex
    from iotbx import mtz
    import libtbx.load_env
    from dials.util import Sorry
    from outlier_rejection import outlier_rejection_ls49

    if outdir is None:
      outdir='.'

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
    exp_list = ExperimentListFactory.from_json_file(os.path.join(ls49_data_dir,'../reindexed_rayonix','idx-%s_integrated.expt'%tstamp), check_format=False)[0:1]
    refls = flex.reflection_table.from_file(os.path.join(ls49_data_dir,'../reindexed_rayonix','idx-%s_integrated.refl'%tstamp))
    refls=refls.select(refls['id']==0)
    # Filter reflections if necessary depending on predictions from diffBragg ?
    if outlier_with_diffBragg:
      exp_list, refls=outlier_rejection_ls49(exp_list, refls,ls49_data_dir=ls49_data_dir, ts=tstamp, outdir=outdir, dump_output_files=False)

    exp = exp_list[0]
    C = exp.crystal
    B = exp.beam
    D = exp.detector

    #refls = flex.reflection_table.from_file('idx-%s_indexed.refl' % tstamp)
    Nbefore = len(refls)
    # Keep validation reflections here using inverse logic of refls selection
    # SHould only use low res reflections for validation
    bool_list = [2.5 < d < 13.5 for d in refls['d']] 
    #bool_list = [not myval for myval in bool_list]
    refls_validation = refls.select(flex.bool(bool_list))
    min_snr=3.0
    snr = refls_validation['intensity.sum.value'] / flex.sqrt(refls_validation['intensity.sum.variance'])
    order = np.argsort(snr)[::-1] 
    refls_validation=refls_validation.select(snr>min_snr) 
    # End validation  related 
    #
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


    # Special case please remove
    #refls_validation=refls

    bboxes = [list(refls['shoebox'][i].bbox)[:4] for i in range(len(refls)) ]
    resolutions = [D[0].get_resolution_at_pixel(B.get_s0(), refls[i]['xyzobs.px.value'][0:2]) for i in range(len(refls))]
    bboxes = np.array(bboxes)
    bboxes[bboxes > 960] = 959
    bboxes[bboxes < 0] = 0


    # validation stuff
    bboxes_validation = [list(refls_validation['shoebox'][i].bbox)[:4] for i in range(len(refls_validation)) ]
    bboxes_validation = np.array(bboxes_validation)
    bboxes_validation[bboxes_validation > 960] = 959
    bboxes_validation[bboxes_validation < 0] = 0
    # end validation stuff
    mill_idx = [ list(refls['miller_index'][i]) for i in range(len(refls)) ]
    R2 = flex.reflection_table.from_file(os.path.join(ls49_data_dir, '../reindexed_rayonix', 'idx-%s_indexed.refl'%tstamp))
    R2=R2.select(R2['id']==0)
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
    # validation start
    num_spots_validation=len(refls_validation)
    tilt_abc_validation = np.zeros((num_spots_validation, 3))
    tilts_validation = []
    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes_validation):
        shoebox_img = img[j1:j2, i1:i2]
        shoebox_mask = is_bg_pixel[j1:j2, i1:i2]
        tilt, bgmask, coeff, _ = utils.tilting_plane(shoebox_img, mask=shoebox_mask, zscore=2)
        tilts_validation.append(tilt)
        tilt_abc_validation[i_spot] = (coeff[1], coeff[2], coeff[0])
    # end validation

    fee_file='idx-fee_data_%s.pickle'%tstamp
    chann_lambda, channI = np.array(pickle.load(open(os.path.join(ls49_data_dir,fee_file), 'r'))[tstamp]).T
    I = interp1d(chann_lambda, channI)
    max_energy  = chann_lambda[np.argmax(channI)]
    min_energy_interpol = max(max_energy - 35, min(chann_lambda))
    max_energy_interpol = min(max_energy + 35, max(chann_lambda))
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
    sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs(+)')]
    #sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'IMEAN')]
    sfall = sfall.as_amplitude_array()
    return {'dxcrystal': C, 'dxdetector': D, 'dxbeam': B, 'mill_idx': mill_idx, 'data_img': img, 'bboxes_x1x2y1y2': bboxes, 
       'tilt_abc': tilt_abc, 'spectrum': spectrum, 'sfall': sfall,
            'mask': is_BAD_pixel, 'experiment':exp, 'indexed_reflections': R2, 'resolution': resolutions, 'cbf_imageset':cbf_imageset, 'refls_validation': refls_validation, 'bboxes_validation': bboxes_validation, 'tilt_abc_validation': tilt_abc_validation}


def run_all_refine_ls49(ts=None, ls49_data_dir=None, show_plotted_images=False, outdir=None, params=None):
    from simtbx.diffBragg.refiners import RefineAll
    from simtbx.diffBragg.sim_data import SimData
    import numpy as np
    from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
    from simtbx.diffBragg import nanoBragg_crystal, nanoBragg_beam
    from copy import deepcopy
    from dxtbx.model.experiment_list import Experiment, ExperimentList, ExperimentListFactory
    from dials.array_family import flex
    import os
    import dxtbx
    from libtbx.easy_pickle import dump
    from outlier_rejection import outlier_rejection_ls49

    if outdir is None:
      outdir='.'

    print ('Inside run_all_refine_ls49: Starting processing')
    data = process_ls49_image_real(tstamp=ts,Nstrongest=7, resmin=2.1, resmax=5.0, ls49_data_dir=ls49_data_dir, outdir=outdir)
    refine_with_psf=True
    plot_images=True # This is a lie. Mostly need this to store model_Lambda for statistics etc

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
    dump_explist.as_file(os.path.join(outdir, 'before_refinement_%s.expt'%ts))
    indexed_reflections.as_file(os.path.join(outdir, 'before_refinement_%s.refl'%ts))

    # Some global variables here for LS49
    mos_spread_deg=0.001
    n_mos_domains=1
    a,b,c,_,_,_ = C.get_unit_cell().parameters()
    Deff = 1000
    C.set_domain_size_ang(Deff)
    C.set_half_mosaicity_deg(mos_spread_deg) 
    Ncells_abc_0 = np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
    n_spots = len(data['tilt_abc'])
    init_local_spotscale = flex.double([1.0]*n_spots)

    # Define number of macrocycles and strategies
    n_macrocycles=5
    total_cycles=4*n_macrocycles

    #ncells_strategy =          [True,  False, False,  False, True,  False, False, False,  False, False, False, True,  False,  False, False, True,  False,  False,  False, True]
    #local_spotscale_strategy = [False, False, False,  False, False, False, False, False,  False, False, False, False, False,  False, False, False, False,  False,  False, False]
    #crystal_scale_strategy =   [True,  False, False,  False, True,  False, False, False,  False, False, False, True,  False,  False, False, True,  False,  False,  False, True]

    #background_strategy =      [False, False,  False, True,  False,  False, False, True,  False, False, True,  False, False,  False, True,  False,  False, False,  True,  True]
    #umat_strategy =            [False, True,   False, False, False,  True,  False, False, True,  False, False, False, True,   False, False, False,  True,  False,  False, True] 
    #bmat_strategy =            [False, False,  True,  False, False,  False, True,  False, False, True,  False, False, False,  True,  False, False,  False, True,   False, True] 
    ncells_strategy =           [True,  False,  False,   False ]*n_macrocycles  + [True,]
    local_spotscale_strategy =  [False, False,  False,   False ]*n_macrocycles  + [False,]
    crystal_scale_strategy =    [True,  False,  False,   False ]*n_macrocycles + [True,]

    background_strategy =       [False, False,  False,   True]*n_macrocycles +  [True,]
    umat_strategy =             [False, True,   False,   False]*n_macrocycles + [True,]
    bmat_strategy =             [False, False,  True,    False]*n_macrocycles + [True,]
    

    for n_cycle in range(total_cycles):
      print ('Refinement cycle %d for timestamp %s'%(n_cycle, ts))
      if n_cycle==0:
        refined_ncells = Ncells_abc_0
        refined_scale = 1.0
        refined_gain = 1.0
        refined_local_spotscale = init_local_spotscale 
        refined_tilt_abc=data['tilt_abc']

      n_refls_to_refine= len(data['bboxes_x1x2y1y2'])
      total_refls_considered= len(data['bboxes_x1x2y1y2'])
      na=0
      nb=n_refls_to_refine-1
  
      # Special case here:
      # If refining umat, then only refine top few reflections
      # Hoping these will be the high resolution ones
      is_special_umat_refine=False
      #if umat_strategy[n_cycle]\
      #   and not background_strategy[n_cycle]\
      #   and not local_spotscale_strategy[n_cycle]\
      #   and (n_cycle) % 4 == 2:
      #  is_special_umat_refine=True

      #if is_special_umat_refine:
      #  n_refls_to_refine=5
      #  na=0
      #  nb=n_refls_to_refine-1
      #  print ('Resolution range of Umat Refinement: %.2f - %.2f'%(data['resolution'][na], data['resolution'][nb]))

      is_special_ncells_scale_refine=False
      #if ncells_strategy[n_cycle]\
      #   and crystal_scale_strategy[n_cycle]\
      #   and not background_strategy[n_cycle]\
      #   and not local_spotscale_strategy[n_cycle]\
      #   and (n_cycle) % 4 == 1:
      #  is_special_ncells_scale_refine=True
      #if is_special_ncells_scale_refine:
      #  n_refls_to_refine=5
      #  na=total_refls_considered-n_refls_to_refine-1
      #  nb=total_refls_considered-1
      #  print ('Resolution range of Ncells+Scale Refinement: %.2f - %.2f'%(data['resolution'][-n_refls_to_refine], data['resolution'][-1]))


      Ncells_abc=refined_ncells # Setting it here for cycles not equal to 0
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

      ucell = C.get_unit_cell().parameters()
      UcellMan = MonoclinicManager(a=ucell[0], b=ucell[1], c=ucell[2], beta=ucell[4]*np.pi / 180.)

      RUC = RefineAll(
          spot_rois=data["bboxes_x1x2y1y2"][na:nb],
          spot_resolution=data['resolution'][na:nb],
          abc_init=refined_tilt_abc[na:nb, :],
          img=data["data_img"],
          SimData_instance=SIM,
          plot_images=plot_images,
          show_plotted_images=show_plotted_images,
          ucell_manager=UcellMan,
          init_gain=refined_gain,
          init_scale=refined_scale,
          init_local_spotscale=refined_local_spotscale[na:nb])

      RUC.trad_conv_eps = 1e-5
      RUC.use_curvatures=False #
      RUC.calc_curvatures=True # if set to False, never uses curvatures
      RUC.refine_with_psf=refine_with_psf

      RUC.refine_ncells=ncells_strategy[n_cycle]
      RUC.refine_crystal_scale=crystal_scale_strategy[n_cycle]
      RUC.refine_local_spotscale=local_spotscale_strategy[n_cycle]
      RUC.refine_background_planes = background_strategy[n_cycle]
      RUC.refine_Umatrix = umat_strategy[n_cycle]
      RUC.refine_Bmatrix = bmat_strategy[n_cycle]
      #RUC.refine_detdist = detdist_strategy[n_cycle]

      RUC.run()

      show_pixel_values=False

      refined_ncells = RUC.x[-4]
      refined_scale = RUC.x[-1]
      refined_gain = RUC.x[-2]
      if not is_special_umat_refine and not is_special_ncells_scale_refine:
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
        print ('Statistics for cycle %d and timestamp %s'%(n_cycle, ts))
        print ('--------- Linear Regression Summary----------')
        z_regression.show_summary()
        print ('----------------------------------------------')
        print ('--------- Linear Correlation Summary----------')
        z_correlation.show_summary()
        # Skewness should match up ?
        z_calc=all_x-flex.mean(all_x)
        z_obs=all_y-flex.mean(all_y)
        skew_calc=flex.mean(flex.pow(z_calc,3))/flex.mean(flex.pow(z_calc,2))**(3/2.)
        skew_obs=flex.mean(flex.pow(z_obs,3))/flex.mean(flex.pow(z_obs,2))**(3/2.)
        print ('Skewness of Raw data = ', skew_obs) 
        print ('Skewness of Modelled data = ', skew_calc) 
      
      if n_cycle==0:
        RUC0=RUC

      if show_plotted_images and n_cycle==total_cycles-1:
        import pylab as plt
        # for paper
        dump_xy_for_paper = []
        for i_spot in range(RUC.n_spots):
          fig, axs = plt.subplots(4,2)
          axs[0][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[1][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[0][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[1][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[2][0].imshow([[0, 1, 1], [0, 1, 2]])
          axs[2][1].imshow([[0, 1, 1], [0, 1, 2]])
          axs[3][0].imshow([[0, 1, 1], [0, 1, 2]]) # before
          axs[3][1].imshow([[0, 1, 1], [0, 1, 2]]) # after

          x = RUC.store_model_Lambda[i_spot]
          y = RUC.store_Imeas[i_spot]
          x1 = RUC.store_model_Lambda[i_spot]
          y1 = RUC.store_Imeas[i_spot]
          x0 = RUC0.store_init_model_Lambda[i_spot]
          y0 = RUC0.store_init_Imeas[i_spot]
          deltaI0 = np.abs(x0-y0)
          deltaI = np.abs(x-y)

        
          vmin = RUC.store_vmin[i_spot]
          vmax = RUC.store_vmax[i_spot]
          vmin1 = RUC.store_vmin[i_spot]
          vmax1 = RUC.store_vmax[i_spot]
          vmin0 = RUC0.store_init_vmin[i_spot]
          vmax0 = RUC0.store_init_vmax[i_spot]
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

          dump_xy_for_paper.append((x,y, vmin, vmax, x0, y0, vmin0, vmax0))

          plt.suptitle("Spot number = %d at %.2f resolution"%(i_spot,RUC.spot_resolution[i_spot]))
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
        from libtbx.easy_pickle import dump
        dump('Fig5a_data_%s.pickle'%ts, dump_xy_for_paper)

      C2 = deepcopy(C)
      if RUC.refine_Umatrix:
        try:
          ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True)
          C2.rotate_around_origin(ax, ang)
        except Exception as e:
          print ('RUC Rotation not possible:',str(e))
      C2.set_B(RUC.get_refined_Bmatrix())
      

      # refined unit cell parameters
      ucell_ref = C2.get_unit_cell().parameters()

      print("ground truth unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell)
      print("refined unit cell: %2.7g,%2.7g,%2.7g,%2.7g,%2.7g,%2.7g" % ucell_ref)
      print("")

      ### DANGER POINT ##
      ### Memory related ###
      SIM.D.free_all()

      C.show()
      C2.show()
      dump_exp = Experiment(imageset=data['cbf_imageset'], 
                           beam=B,
                           detector=D,
                           goniometer=exp.goniometer,
                           scan=exp.scan,
                           crystal=C2)
      dump_explist = ExperimentList([dump_exp])
      if n_cycle==total_cycles-1:
        dump_explist.as_file(os.path.join(outdir, 'after_refinement_%s_%s.expt'%(ts,n_cycle)))
      # Dump refl file as well based on prediction from refined model
      if True:
        from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
        ref_predictor = ExperimentsPredictorFactory.from_experiments(
                        dump_explist,
                        force_stills=True,
                        spherical_relp=False)
        ref_predictor(indexed_reflections)
        if n_cycle==total_cycles-1:
          indexed_reflections.as_file(os.path.join(outdir, 'after_refinement_%s_%s.refl'%(ts,n_cycle)))
          # Also do the integration here for good measure
          if True:
            print ('Integrating reflection in diffBragg for timestamp: %s'%ts)
            integrated_expts, integrated_refls=integrate(dump_explist, indexed_reflections, params)
            # Compare before and after outlier rejection
            #integrated_expts.as_file(os.path.join(outdir, 'integrated_b4_outlier_%s_%s.expt'%(ts, n_cycle)))
            #integrated_refls.as_file(os.path.join(outdir, 'integrated_b4_outlier_%s_%s.refl'%(ts, n_cycle)))
            integrated_expts, integrated_refls=outlier_rejection_ls49(integrated_expts, integrated_refls,ls49_data_dir=ls49_data_dir, ts=ts, outdir=outdir, dump_output_files=False)
            #integrated_expts.as_file(os.path.join(outdir, 'integrated_a4_outlier_%s_%s.expt'%(ts, n_cycle)))
            #integrated_refls.as_file(os.path.join(outdir, 'integrated_a4_outlier_%s_%s.refl'%(ts, n_cycle)))
            # Also write out integration pickles for cxi.merge
            from xfel.command_line.frame_extractor import ConstructFrame
            from libtbx import easy_pickle
            frame = ConstructFrame(integrated_refls, integrated_expts[0]).make_frame()
            frame["pixel_size"] = integrated_expts[0].detector[0].get_pixel_size()[0]
            fname='int-a4-outlier-'+str(n_cycle)+'-0-'+ts+'.pickle'
            outfile=os.path.join(outdir, fname)
            easy_pickle.dump(outfile, frame)
            print ('All integration files written out for timestamp: %s'%ts)
            #from IPython import embed; embed(); exit()
      C=C2 # For next round

###################################################################################
#  Rayonix validation: Predict on all shoeboxes in rayonix and then calculate statistics
      if False:
        refls_validation=data['refls_validation']

        nbcryst_validate_ray = nanoBragg_crystal.nanoBragg_crystal()
        nbcryst_validate_ray.Ncells_abc = refined_ncells, refined_ncells, refined_ncells
        nbcryst_validate_ray.mos_spread_deg = mos_spread_deg
        nbcryst_validate_ray.n_mos_domains = n_mos_domains
        nbcryst_validate_ray.thick_mm = 0.005
        nbcryst_validate_ray.miller_array = data["sfall"]
        nbcryst_validate_ray.dxtbx_crystal = C2

        SIM_validate_ray = SimData()
        SIM_validate_ray.crystal = nbcryst_validate_ray
        SIM_validate_ray.detector = D
        SIM_validate_ray.beam = nbbeam

        SIM_validate_ray.instantiate_diffBragg(adc_offset=0,
                                  oversample=0,
                                  interpolate=0,
                                  verbose=0)
        SIM_validate_ray.D.show_params()
        #SIM_jung.D.spot_scale = 1e6 # to guide your eye
        scale=refined_scale
        n_spots_validation = len(data['tilt_abc_validation'])
        show_validation_images=False
        if show_validation_images:
          import pylab as plt
          plt.clf()
          fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
          ax1.imshow([[0, 1, 1], [0, 1, 2]])
          ax2.imshow([[0, 1, 1], [0, 1, 2]])
        all_x = flex.double()
        all_y = flex.double()
        for i_spot in range(n_spots_validation):
          print ('Spot number = ', i_spot)
          bbox=data['bboxes_validation'][i_spot]
          a,b,c=data['tilt_abc_validation'][i_spot]
          x1,x2,y1,y2=bbox
          nx_lim, ny_lim = data['data_img'].shape
          if x1 < 0: x1=0
          if y1 < 0: y1=0
          if x2 >=nx_lim: x2=nx_lim-1
          if y2 >=ny_lim: y2=ny_lim-1
          yr, xr = np.indices((y2-y1+1, x2-x1+1))
          tilt_plane=xr*a+yr*b+c
          SIM_validate_ray.D.region_of_interest=((bbox[0], bbox[1]), (bbox[2], bbox[3]))
          SIM_validate_ray.D.add_diffBragg_spots()
          img2=SIM_validate_ray.D.raw_pixels.as_numpy_array()
          img2=img2[y1:y2+1, x1:x2+1]

          img=tilt_plane+scale*scale*img2

          if refine_with_psf:
            from simtbx.diffBragg.utils import convolve_with_psf
            img=convolve_with_psf(img)
          data_img=data['data_img'][y1:y2+1, x1:x2+1]
          x = img 
          y = data_img
          #from IPython import embed; embed(); exit()
          all_x.extend(flex.double(x.flatten()))
          all_y.extend(flex.double(y.flatten()))

          if show_validation_images:
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

            ax1.images[0].set_data(img)
            ax1.images[0].set_clim(vmin, vmax)
            ax2.images[0].set_data(data_img)
            ax2.images[0].set_clim(vmin0, vmax0)
            plt.suptitle('Spot Number: %d, Res = %.2f'%(i_spot, refls_validation['d'][i_spot]))
            fig.canvas.draw()
            plt.pause(1)
        z_regression=flex.linear_regression(all_x, all_y)
        z_correlation=flex.linear_correlation(all_x, all_y)
        print ('Validation for cycle %d and timestamp %s'%(n_cycle, ts))
        print ('--------- Linear Regression Summary----------')
        z_regression.show_summary()
        print ('----------------------------------------------')


###################################################################################
###################################################################################

####################################################################################

      # Now display prediction on jungfrau
      # Do it on a cycle by cycle basis
      if False:
        cbf_jungfrau_path=os.path.join(ls49_data_dir, '../', 'jungfrau_cbf','jungfrauhit_%s.cbf'%ts )
        loader = dxtbx.load(cbf_jungfrau_path)
        cbf_imageset = loader.get_imageset([cbf_jungfrau_path])
        img=cbf_imageset.get_raw_data(0)
        assembled_img=np.zeros([1024, 1024])
        assembled_img[0:256, 0:256]=img[4].as_numpy_array()[::-1, :]
        assembled_img[0:256, 256:512]=img[5].as_numpy_array()[::-1, :]
        assembled_img[0:256, 512:768]=img[6].as_numpy_array()[::-1, :]
        assembled_img[0:256, 768:1024]=img[7].as_numpy_array()[::-1, :]

        assembled_img[256:512, 0:256]=img[0].as_numpy_array()[::-1, :]
        assembled_img[256:512, 256:512]=img[1].as_numpy_array()[::-1, :]
        assembled_img[256:512, 512:768]=img[2].as_numpy_array()[::-1, :]
        assembled_img[256:512, 768:1024]=img[3].as_numpy_array()[::-1, :]

        assembled_img[512:768, 0:256]=img[12].as_numpy_array()[::-1, :]
        assembled_img[512:768, 256:512]=img[13].as_numpy_array()[::-1, :]
        assembled_img[512:768, 512:768]=img[14].as_numpy_array()[::-1, :]
        assembled_img[512:768, 768:1024]=img[15].as_numpy_array()[::-1, :]

        assembled_img[768:1024, 0:256 ]=img[8].as_numpy_array()[::-1, :]
        assembled_img[768:1024, 256:512]=img[9].as_numpy_array()[::-1, :]
        assembled_img[768:1024, 512:768]=img[10].as_numpy_array()[::-1, :]
        assembled_img[768:1024, 768:1024]=img[11].as_numpy_array()[::-1, :]

        D_jung = ExperimentListFactory.from_json_file(os.path.join(ls49_data_dir, '../', 'fake_jungfrau_from_rayonix.expt'), check_format=False)[0].detector
        #D_jung = ExperimentListFactory.from_json_file('../jungfrau_scripts/rotated_plus_90.json', check_format=False)[0].detector

        nbcryst_final = nanoBragg_crystal.nanoBragg_crystal()
        nbcryst_final.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
        nbcryst_final.mos_spread_deg = mos_spread_deg
        nbcryst_final.n_mos_domains = n_mos_domains
        nbcryst_final.thick_mm = 0.005
        nbcryst_final.miller_array = data["sfall"]
        nbcryst_final.dxtbx_crystal = C2

        SIM_jung = SimData()
        SIM_jung.crystal = nbcryst_final
        SIM_jung.detector = D_jung
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

        stuff_to_dump=[assembled_img, img]

        if show_plotted_images:
          import pylab as plt
          plt.clf()
          fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
          ax1.imshow([[0, 1, 1], [0, 1, 2]])
          ax2.imshow([[0, 1, 1], [0, 1, 2]])
          m = assembled_img[assembled_img>1.e-9].mean()
          s = assembled_img[assembled_img>1.e-9].std()
          vmin = 1
          vmax = m+10*s
          ax1.images[0].set_data(assembled_img)
          ax1.images[0].set_clim(vmin, vmax)
          ax2.images[0].set_data(img)
          ax2.images[0].set_clim(1, 1000)
          dump(os.path.join(outdir, 'jungfrau_pred_%s_%d.pickle'%(ts, n_cycle)), fig )

        dump(os.path.join(outdir, 'jungfrau_comparison_%s_%d.pickle'%(ts, n_cycle)), stuff_to_dump )
        #print ('Dumped Jungfrau Prediction for %d cycle'%n_cycle)
        #plt.show()
    # DANGER return here
    return

if __name__ == "__main__":
    from dials.command_line.stills_process import phil_scope as stills_process_phil_scope
    from libtbx.phil import parse
    import os, sys
    ''' Example phil file

        output {
          experiments_filename = None
          strong_filename = None
        }
        spotfinder {
          filter {
            min_spot_size = 2
            d_min = 2
          }
          threshold {
            dispersion {
              gain = 0.46
              global_threshold = 50
            }
          }
        }
        indexing {
          known_symmetry {
            space_group = P 1 21 1
            unit_cell = 63.6,28.8,35.6,90,106.5,90
          }
          refinement_protocol {
            d_min_start = 2
          }
        }
        indexing.basis_vector_combinations.max_refine=1
        indexing.stills.candidate_outlier_rejection=False
        indexing.index_assignment.simple.hkl_tolerance=0.30
        profile.gaussian_rs.centroid_definition=com
        spotfinder.threshold.algorithm=dispersion
        spotfinder.lookup.mask='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/input/mask_r4.pickle'
        integration.lookup.mask='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/input/mask_r4.pickle'
        integration.summation.detector_gain=0.46
        integration.debug.output=True
        integration.debug.separate_files=False'''

    args = sys.argv[1:]
    user_phil = []
    for arg in args:
      if os.path.isfile(arg):
        user_phil.append(parse(file_name=arg))
      else:
        try:
          user_phil.append(parse(arg))
        except Exception as e:
          raise Sorry("Unrecognized argument: %s"%arg)
    stills_process_params = stills_process_phil_scope.fetch(sources=user_phil).extract()

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

    #ts = timestamps_of_interest[3]
    #ts='20180501132216201'
    #ts='20180501120317142'
    #ts='20180501164500977'
    #ts='20180501151638379'
    ts='20180501165223017'
    # Keep this as rayonix_expt even though the integrated expt/refls are in a different directory reindexed_rayonix
    ls49_data_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/all_files/rayonix_expt'
    # On cori here is the path
    #ls49_data_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_grid_search_4_or_more_regression/rayonix_images_4_or_more_spots_r183_255'
    #ls49_data_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/all_files/rayonix_expt'
    # On my Macbook Pro, here is the path
    #ls49_data_dir='/Users/abhowmick/Desktop/software/dials/modules/LS49_regression/diffBragg_work/jungfrau_grid_search_4_or_more_regression/rayonix_images_4_or_more_spots_r183_255'
    #ts='20180501114703722' # Image used in blog to compare on jungfrau
    #ts='20180501120317142'
    #ts='20180501163914779' # weird , low rotation change but blows up in prediction ?? !!
    outdir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_grid_search_4_or_more_regression/temp_2'
    #ts='20180501114755146'
    #outdir=None
    #ls49_data_dir=None
    run_all_refine_ls49(ts=ts, ls49_data_dir=ls49_data_dir, outdir=outdir, show_plotted_images=True, params=stills_process_params)
