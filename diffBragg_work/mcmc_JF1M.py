#
# File courtesy Derek Mendez
# Copying it from cctbx_project/simtbx/diffBragg to make more local edits here instead of in cctbx
#
#
from load_ls49_JF1M import strong_spot_mask, process_ls49_image_real

#def process_ls49_image_real(tstamp='20180501143555114', #tstamp='20180501143559313',
#                            Nstrongest = 30,
#                            resmax=12.0, resmin=3.0,
#                            #mtz_file='5cmv_Iobs.mtz',
#                            #mtz_file='anom_ls49_oxy_2.3_t3_gentle_pr_s0_mark0.mtz',
#                            #mtz_file='anom_ls49_oxy_2.1_unit_pr_lorentz_double_primeref_m008_s0_mark0.mtz',
#                            pdb_file='Refine47_withFE.pdb',
#                            seed=0.0,
#                            ls49_data_dir=None):
#    import os, pickle, sys, numpy as np
#    from scipy.interpolate import interp1d
#    import dxtbx
#    from dxtbx.model.experiment_list import ExperimentListFactory
#    from simtbx.diffBragg import utils
#    from dials.array_family import flex
#    from iotbx import mtz
#    import libtbx.load_env
#    from dials.util import Sorry
#    import math
#
#    if ls49_data_dir is None:
#      LS49_regression = libtbx.env.find_in_repositories(
#        relative_path="LS49_regression",
#        test=os.path.isdir)
#      if LS49_regression is None:
#        raise Sorry('LS49_regression folder needs to be present or else specify ls49_data_dir')
#      ls49_data_dir = os.path.join(LS49_regression, 'diffBragg_work', 'iota_r0222_cori', 'rayonix_expt')
#
#    GAIN = 1.0
#    loader = dxtbx.load(os.path.join(ls49_data_dir,'../../all_files/jungfrau_cbf','jungfrauhit_%s.cbf'%tstamp))
#    cbf_imageset = loader.get_imageset([os.path.join(ls49_data_dir,'../../all_files/jungfrau_cbf','jungfrauhit_%s.cbf'%tstamp)])
#    
#    img = []
#    for raw_data in loader.get_raw_data():
#      img.append(raw_data.as_numpy_array()/GAIN)
#    exp_list = ExperimentListFactory.from_json_file(os.path.join(ls49_data_dir,'jungfrau_shoeboxes_%s.expt'%tstamp), check_format=False)
#    exp = exp_list[0]
#    C = exp.crystal
#    B = exp.beam
#    D = exp.detector
#    #refls = flex.reflection_table.from_file('idx-%s_indexed.refl' % tstamp)
#    refls = flex.reflection_table.from_file(os.path.join(ls49_data_dir,'jungfrau_shoeboxes_%s.refl'%tstamp))
#    Nbefore = len(refls)
#    #refls = refls.select(flex.bool([resmin < d < resmax for d in refls['d']]))
#    #print("Kept %d out of %d refls in the res range %2.2f to %2.2f"
#    #      % (len(refls), Nbefore, resmin, resmax))
#
#    snr = flex.double()
#    for refl in refls.rows():
#      i1,i2,j1,j2 = refl['bbox'][:4]
#      if i1 <0: i1=0
#      if j1<0: j1=0
#       
#      temp_img=img[refl['panel']][j1:j2, i1:i2]
#      #from IPython import embed; embed(); exit()
#      #sorted_img=np.sort(temp_img.flatten())
#      #bg=sorted_img[10:100]
#      #fg=sorted_img[-10:]
#      #bgm=np.mean(bg)
#      #fgm=np.mean(fg)
#      snr.append(np.max(temp_img)/np.std(temp_img))
#      #snr.append((fgm)/bgm)
#
#    # FIXME WHICH SNR TO USE ??
#    #snr = refls['intensity.sum.value'] / flex.sqrt(refls['intensity.sum.variance'])
#    #from IPython import embed; embed(); exit()
#    order = np.argsort(snr)[::-1]
#    min_snr=1.0
#    refls=refls.select(snr>min_snr)
#
#    # Filter out based on whether bbox is near the edge of detector
#    # Spurious refinement noticed
#    keep_refl=flex.bool(len(refls), True)
#    #from IPython import embed; embed(); exit()
#    for ii,refl in enumerate(refls.rows()):
#      #from IPython import embed; embed(); exit()
#      x1,x2,y1,y2=refl['bbox'][0:4]
#      #if refl['panel'] not in [7,3,15,11]: continue
#      if x2 > 255 and refl['panel'] in [7,3,15,11]:
#        keep_refl[ii]=False
#        continue
#      #if x2 < 1 or x2 > 255:
#      #  keep_refl[ii]=False
#      #  continue
#      #if y1 < 1 or y1 > 255:
#      #  keep_refl[ii]=False
#      #  continue
#      #if y2 < 1 or y2 > 255:
#      #  keep_refl[ii]=False
#      #  continue
#    refls=refls.select(keep_refl)  
#    # Filter here based on only strong spots. This is to avoid blank shoeboxes being processed by diffbragg
#    strong_refls=flex.reflection_table.from_file(os.path.join(ls49_data_dir, '../../out_jungfrau_spotfinding_v3/idx-jungfrauhit_%s_strong.refl'%tstamp))
#    int_px=refls['xyzobs.px.value']
#    strong_px=strong_refls['xyzobs.px.value']
#    keep_spots = flex.bool(len(int_px), False)
#  
#    critical_dist = 10 # in pixels
#    for ii, spot in enumerate(int_px):
#      x0,y0,z0  = spot
#      for jj, strong_spot in enumerate(strong_px):
#        x,y,z = strong_spot
#        dist = math.sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0))
#        if dist < critical_dist:
#          keep_spots[ii]=True
#          break
#    refls=refls.select(keep_spots)
#    #from IPython import embed; embed(); exit()
#
#    #refls = refls.select(snr > snr[order[Nstrongest]])
#
#    bboxes = [list(refls['shoebox'][i].bbox)[:4] for i in range(len(refls)) ]
#    panels = [refl['panel'] for refl in refls.rows()]
#    resolutions = [D[0].get_resolution_at_pixel(B.get_s0(), refls[i]['xyzobs.px.value'][0:2]) for i in range(len(refls))]
#    bboxes = np.array(bboxes)
#    bboxes[bboxes > 256] = 255
#    bboxes[bboxes < 0] = 0
#    mill_idx = [ list(refls['miller_index'][i]) for i in range(len(refls)) ]
#    R2 = flex.reflection_table.from_file(os.path.join(ls49_data_dir, '../out_using_diffBragg_rayonix_v2/no_duplicate_millers/select_best_jf_refined/jf_refined_%s.refl'%tstamp))
#    #R2 = flex.reflection_table.from_file(os.path.join(ls49_data_dir, '../out_using_diffBragg_rayonix_v3/jungfrau_indexed_only_%s.refl'%tstamp))
#    n_panels=16
#    strong_mask = strong_spot_mask(refls=R2, panel_size=img[0].shape, n_panels=n_panels)
#    is_bg_pixel = []
#    BAD_pixel_mask=pickle.load(open(os.path.join(ls49_data_dir,'../','new_jungfrau_mask_panel13.pickle'), 'r'))
#    is_BAD_pixel = []
#    for pid in range(n_panels):
#      panel_is_bg_pixel = np.logical_not(strong_mask[pid]) 
#      #panel_BAD_pixel_mask=BAD_pixel_mask[pid].as_numpy_array()
#      #if pid==15:
#      #  from IPython import embed; embed(); exit()
#      #panel_is_bg_pixel[panel_BAD_pixel_mask]=False
#      
#      #is_bg_pixel.append(np.logical_not(strong_mask[pid]))
#      is_bg_pixel.append(panel_is_bg_pixel)
#    is_bg_pixel = np.array(is_bg_pixel)
#    #is_bg_pixel = np.logical_not(strong_mask)
#    
#    #for pid in range(n_panels):
#    #  is_BAD_pixel.append(BAD_pixel_mask[pid].as_numpy_array())
#    #is_BAD_pixel = np.array(is_BAD_pixel)
#    #is_BAD_pixel = np.logical_not(list(pickle.load(open(os.path.join(ls49_data_dir,'../','new_jungfrau_mask_panel13.pickle'), 'r'))))
#    #is_bg_pixel[is_BAD_pixel] = False
#    num_spots = len(refls)
#    tilt_abc = np.zeros((num_spots, 3))
#    tilts = []
#    for i_spot, (i1, i2, j1, j2) in enumerate(bboxes):
#        pid=panels[i_spot]
#        if i1<0: i1=0
#        if j1<0: j1=0
#        shoebox_img = img[pid][j1:j2, i1:i2]
#        #from IPython import embed; embed(); exit()
#        shoebox_mask = is_bg_pixel[pid, j1:j2, i1:i2]
#        tilt, bgmask, coeff, _ = utils.tilting_plane(shoebox_img, mask=shoebox_mask, zscore=2.0)
#        tilts.append(tilt)
#        tilt_abc[i_spot] = (coeff[1], coeff[2], coeff[0])
#
#    fee_file='idx-fee_data_%s.pickle'%tstamp
#    chann_lambda, channI = np.array(pickle.load(open(os.path.join(ls49_data_dir,'../../all_files/rayonix_expt/',fee_file), 'r'))[tstamp]).T
#    I = interp1d(chann_lambda, channI)
#    max_energy  = chann_lambda[np.argmax(channI)]
#    min_energy_interpol = max(max_energy - 35, min(chann_lambda))
#    max_energy_interpol = min(max_energy + 35, max(chann_lambda))
#    print ('INTERPOLATION ENERGIES = ', min_energy_interpol, max_energy_interpol)
#    interp_energies = np.arange(min_energy_interpol, max_energy_interpol, 0.5)
#    interp_fluxes = I(interp_energies)
#    interp_fluxes /= interp_fluxes.sum()
#    interp_fluxes *= 1000000000000.0
#    spectrum = zip(12398.419739640716 / interp_energies, interp_fluxes)
#   
#    ########## FP/FDP stuff #################### 
#    # Generate fp/fdp as desired for Fe
#    from cctbx.eltbx import henke
#    import numpy as np
#
#    from libtbx.easy_pickle import load
#
#    # Special cases of seed
#    # -1 --> Fe0 from Fe.dat table in data
#    # -2 --> Fe2
#    # -3 --> Fe3
#    # -20 --> Fe0 translated 20 ev
#    # -10 --> Fe3 translated 10 ev
#
#    # All positive numbers > 3 are random seeds used to generate monte carlo trial curves
#
#    if seed == -3:
#      energies, all_fp, all_fdp = load(os.path.join(ls49_data_dir, '../scattering_factor_refinement/data_sherrell/oxidized_Fe.pickle'))
#    if seed == -2:
#      energies, all_fp, all_fdp = load(os.path.join(ls49_data_dir, '../scattering_factor_refinement/data_sherrell/reduced_Fe.pickle'))
#    if seed == -1:
#      energies, all_fp, all_fdp = load(os.path.join(ls49_data_dir, '../scattering_factor_refinement/data_sherrell/neutral_Fe.pickle'))
#
#
#    # Create mcmc random curves if seed > 3
#    # Reserved numbers
#    # 0 --> neutral Fe, no fp/fdp added
#    # 1 --> Use fp/fdp discretized curve from Fe.dat neutral
#    # 3 --> use fp/fdp discretized curve from Fe3 sherrell data 
#    sherrell_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_setup_before_scattering_factor/scattering_factor_refinement/data_sherrell'
#    if seed > 0:
#      from scipy.interpolate import CubicSpline
#      sherrell_files=['Fe.dat','pf-rd-ox_fftkk.out', 'pf-rd-red_fftkk.out']
#      u_ev = [7102,7107,7112,7117,7122, 7127, 7132, 7137, 7142]
#      #u_ev = [7100,7105,7110,7115,7120, 7125, 7130, 7135, 7140]
#      #u_ev = [x+2 for x in u_ev] 
#      storage = {}
#      fp0 = [] # for high and low energy remotes
#      fdp0 = [] # for high and low energy remotes
#      ev0=[]
#      for sherrell_file in sherrell_files:
#        energies = []
#        all_fp = []
#        all_fdp = []
#        with open(os.path.join(sherrell_dir, sherrell_file), 'r') as f:
#          for line in f:
#            if line !='\n':
#              ax = line.split()
#              ev = float(ax[0])
#              fp = float(ax[1])
#              fdp = float(ax[2])
#              if ev < 7200 and ev > 7050:
#                energies.append(ev)
#                all_fp.append(fp)
#                all_fdp.append(fdp)
#                if sherrell_file == sherrell_files[0]:
#                  fp0.append(fp)
#                  fdp0.append(fdp)
#                  ev0.append(ev)
#        x=(energies, all_fp, all_fdp)
#        csp=CubicSpline(energies, all_fp)
#        csdp=CubicSpline(energies, all_fdp)
#        y_fp=csp(u_ev)
#        y_fdp=csdp(u_ev)
#        storage[sherrell_file]=(y_fp, y_fdp)  
#      if seed == 1:
#        energies, all_fp, all_fdp =  u_ev, list(storage['Fe.dat'][0]), list(storage['Fe.dat'][1])
#      if seed == 3:
#        energies, all_fp, all_fdp =  u_ev, list(storage['pf-rd-ox_fftkk.out'][0]), list(storage['pf-rd-ox_fftkk.out'][1]) 
#      if seed == 2:
#        energies, all_fp, all_fdp =  u_ev, list(storage['pf-rd-red_fftkk.out'][0]), list(storage['pf-rd-red_fftkk.out'][1]) 
#
#
#      # Creat mcmc moves
#      if seed > 3:
#        all_fp = []
#        all_fdp = []
#        energies = []
#        flex.set_random_seed(seed)
#        alpha=flex.random_double()
#        print ('Random seed and aplha = %d  %5.2f'%(seed, alpha))
#        padding = 0.0
#        
#        for ii, ev in enumerate(u_ev):
#          # get Fe0 data
#          fp00=storage[sherrell_files[0]][0][ii]
#          fdp00=storage[sherrell_files[0]][1][ii]
#          fp3=storage[sherrell_files[1]][0][ii]
#          fdp3=storage[sherrell_files[1]][1][ii]
#          #from IPython import embed; embed(); exit()
#        
#          t_fp = fp00+alpha*(fp3-fp00+padding)
#          t_fdp = fdp00+alpha*(fdp3-fdp00+padding)
#          all_fp.append(t_fp)
#          all_fdp.append(t_fdp)
#          energies.append(ev)
#
#      # End monte carlo bit
#
#      # Make sure you take care of remote fp/fdp values on both sides
#      #from IPython import embed; embed(); exit()
#      energies = [ev0[0]]+energies
#      all_fp = [fp0[0]]+all_fp
#      all_fdp = [fdp0[0]]+all_fdp
#      #
#      energies = energies + [ev0[-1]]
#      all_fp = all_fp + [fp0[-1]]
#      all_fdp = all_fdp + [fdp0[-1]]
#      
# 
#
#    #if ev_offset == -20:
#    #  energies, all_fp, all_fdp = load(os.path.join(ls49_data_dir, '../scattering_factor_refinement/data_sherrell/neutral_Fe.pickle'))
#    #  energies = [x+20 for x in energies] 
##
#    #if ev_offset == -10:
#    #  energies, all_fp, all_fdp = load(os.path.join(ls49_data_dir, '../scattering_factor_refinement/data_sherrell/oxidized_Fe.pickle'))
#    #  energies = [x+10 for x in energies] 
#    
#    from scipy.interpolate import interp1d
#    
#
#    interpolator = interp1d(energies, all_fp)
#    interpolator2 = interp1d(energies, all_fdp)
#
#    interp_fp = interpolator(interp_energies)
#    interp_fdp = interpolator2(interp_energies)
#    
#    #from IPython import embed; embed(); exit()
#    if False:
#      ev3, all_fp, all_fdp = load(os.path.join(ls49_data_dir, '../scattering_factor_refinement/data_sherrell/oxidized_Fe.pickle'))
#      #ev3, all_fp, all_fdp = load(os.path.join(ls49_data_dir, '../scattering_factor_refinement/data_sherrell/neutral_Fe.pickle'))
#      import matplotlib.pyplot as plt
#      #from IPython import embed; embed(); exit()
#      plt.scatter(interp_energies, interp_fdp)
#      plt.scatter(ev3, all_fdp)
#      plt.scatter(interp_energies, interp_fp)
#      plt.scatter(ev3, all_fp)
#      plt.show()
#    fp_fdp = (flex.double(interp_fp), flex.double(interp_fdp))
#    fp_fdp_0 = (flex.double([0.0]*len(interp_fp)), flex.double([0.0]*len(interp_fdp)))
#    
#    #############################################
#
#
#    #M = mtz.object('ls49_oxy_2.5_s0_mark0.mtz')
#    #sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs')]
#
#    #M = mtz.object(os.path.join(ls49_data_dir,'../../all_files/',mtz_file))
#    #sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs')]
#    #sfall = M.as_miller_arrays_dict()[('crystal', 'dataset', 'Iobs(+)')]
#    #sfall = sfall.as_amplitude_array()
#
#
#    # Generate complex structure factors here from the PDB file and pass along
#    #  
#    sfall = utils.get_complex_fcalc_from_pdb(os.path.join(ls49_data_dir, '../../all_files/', pdb_file))
#
#    return {'dxcrystal': C, 'dxdetector': D, 'dxbeam': B, 'mill_idx': mill_idx, 'data_img': img, 'bboxes_x1x2y1y2': bboxes, 
#       'tilt_abc': tilt_abc, 'spectrum': spectrum, 'sfall': sfall,
#            'mask': is_BAD_pixel, 'experiment':exp, 'indexed_reflections': R2, 'resolution': resolutions, 'cbf_imageset':cbf_imageset, 'panels':panels, 'fp_fdp':fp_fdp, 'fp_fdp_0':fp_fdp_0}
#

def run_all_refine_ls49_JF1M(ts=None, ls49_data_dir=None, show_plotted_images=False, outdir=None, params=None, seed=0, short_circuit_dir=None):
    from argparse import ArgumentParser
    from simtbx.diffBragg.refiners import RefineAll_JF1M_MultiPanel
    from simtbx.diffBragg.sim_data import SimData
    import numpy as np
    from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
    from simtbx.diffBragg import nanoBragg_crystal, nanoBragg_beam
    from copy import deepcopy
    from dxtbx.model.experiment_list import Experiment, ExperimentList, ExperimentListFactory
    from dials.array_family import flex
    import os
    from libtbx.easy_pickle import load

    data = process_ls49_image_real(tstamp=ts,Nstrongest=10, resmin=2.1, resmax=3.5, ls49_data_dir=ls49_data_dir, seed=seed)
    #C = data["dxcrystal"]
    #D = data["dxdetector"]
    #B = data["dxbeam"]

    # Load stuff here from short_circuit_dir
    if short_circuit_dir is None:
      raise Sorry('Need to supply short-circuit dir') 

    info=load(os.path.join(short_circuit_dir, 'RUC_info_%s.pickle'%ts))
    C=info['C']
    D=info['D']
    B=info['B']

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
    Ncells_abc_0 = np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
    C.set_domain_size_ang(Deff)
    C.set_half_mosaicity_deg(mos_spread_deg)
    Ncells_abc_0 = np.power(4/3*np.pi*Deff**3/a/b/c, 1/3.)
    n_spots = len(data['tilt_abc'])


    # Define number of macrocycles and strategies
    n_macrocycles=5
    total_cycles=3*n_macrocycles
    ncells_strategy =           [True,  False, False, ]*n_macrocycles
    local_spotscale_strategy =  [False, False, False,]*n_macrocycles
    crystal_scale_strategy =    [True,  False, False,]*n_macrocycles

    background_strategy =       [False, True,  False,]*n_macrocycles
    umat_strategy =             [False, False, True, ]*n_macrocycles 
    bmat_strategy =             [False, False, False,]*n_macrocycles 

    for n_cycle in range(total_cycles):
      print ('Refinement cycle %d for timestamp %s'%(n_cycle, ts))
      if n_cycle==0:
        refined_ncells = info['refined_ncells']
        refined_scale = info['refined_scale']
        refined_gain = info['refined_gain']
        refined_tilt_abc=info['refined_tilt_abc']
        init_local_spotscale = info['refined_local_spotscale']
        refined_local_spotscale = init_local_spotscale


      Ncells_abc=refined_ncells # Setting it here for cycles not equal to 0
      nbcryst = nanoBragg_crystal.nanoBragg_crystal()
      nbcryst.Ncells_abc = Ncells_abc, Ncells_abc, Ncells_abc
      nbcryst.mos_spread_deg = mos_spread_deg
      nbcryst.n_mos_domains = n_mos_domains
      nbcryst.thick_mm = 0.005
      nbcryst.miller_array = data["sfall"]
      nbcryst.dxtbx_crystal = C

      nbcryst.fp_fdp = data['fp_fdp']

      nbbeam = nanoBragg_beam.nanoBragg_beam()
      nbbeam.unit_s0 = B.get_unit_s0()
      nbbeam.spectrum = data["spectrum"]
      

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
      #SIM.D.Ncells_abc =  10
      
      SIM.D.detector_attenuation_length_mm = 0
      SIM.D.detector_thicksteps = 0
      SIM.D.point_pixel = False  # NOTE detector distance negative bug
      SIM.D.show_params()
      #SIM.D.beam_center_mm=D[data['panels'][0]].get_beam_centre(B.get_unit_s0())
      #SIM.D.spot_scale = 1e6 # to guide your eye
      #SIM.panel_id = 0
      #SIM.D.add_diffBragg_spots()
      #img=SIM.D.raw_pixels.as_numopy_array()

      ucell = C.get_unit_cell().parameters()
      UcellMan = MonoclinicManager(a=ucell[0], b=ucell[1], c=ucell[2], beta=ucell[4]*np.pi / 180.)
      if False: #args.plot:
          import pylab as plt
          I = copy.deepcopy(data["data_img"])
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
      # DANGER FIXME
      if n_spots < 2:
        print ('Too few spots for diffbragg refinement :  %s'%ts)
        return
      #
      RUC = RefineAll_JF1M_MultiPanel(
        spot_rois=data["bboxes_x1x2y1y2"],
        spot_resolution=data['resolution'],
        abc_init=refined_tilt_abc,
        img=data["data_img"],
        SimData_instance=SIM,
        plot_images=True,
        show_plotted_images=show_plotted_images,
        ucell_manager=UcellMan,
        init_gain=refined_gain,
        init_scale=refined_scale,
        init_local_spotscale = refined_local_spotscale,
        panel_id=data['panels'][0],
        panel_ids=data['panels'],
        timestamp=ts)
      RUC.trad_conv = True
      RUC.trad_conv_eps = 1e-5
      RUC.max_calls = 250
      RUC.refine_ncells=ncells_strategy[n_cycle]
      RUC.refine_crystal_scale=crystal_scale_strategy[n_cycle]
      RUC.refine_local_spotscale=local_spotscale_strategy[n_cycle]
      RUC.refine_background_planes = background_strategy[n_cycle]
      RUC.refine_Umatrix = umat_strategy[n_cycle]
      RUC.refine_Bmatrix = bmat_strategy[n_cycle]



      RUC.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
      RUC.use_curvatures=False #
      RUC.calc_curvatures=True # if set to False, never uses curvatures
    #RUC.run()
    #if RUC.hit_break_to_use_curvatures:
    #  RUC.num_positive_curvatures=0
    #  RUC.use_curvatures=True
    #  RUC.run(setup=False) # Now it will use curvatures

      RUC.run()
      refined_ncells = RUC.x[RUC.ncells_xpos]
      refined_scale = RUC.x[-1]
      refined_gain = RUC.x[-2]
      refined_local_spotscale = RUC.x[RUC.n_bg:RUC.n_bg+RUC.n_spots]
      refined_tilt_abc=[]
      for i in range(RUC.n_spots):
        refined_tilt_abc.append([RUC.x[i], RUC.x[RUC.n_spots+i], RUC.x[2*RUC.n_spots+i]])
      refined_tilt_abc=np.array(refined_tilt_abc)


      if n_cycle==0:
        RUC0=RUC

    ##########

      show_pixel_values=False
      if show_plotted_images and n_cycle==total_cycles-1:
        import pylab as plt
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
          x1 = RUC0.store_model_Lambda[i_spot]
          y1 = RUC0.store_Imeas[i_spot]
          x0 = RUC0.store_init_model_Lambda[i_spot]
          y0 = RUC0.store_init_Imeas[i_spot]
          deltaI0 = np.abs(x0-y0)
          deltaI = np.abs(x-y)
        
          vmin = RUC.store_vmin[i_spot]
          vmax = RUC.store_vmax[i_spot]
          vmin1 = RUC0.store_vmin[i_spot]
          vmax1 = RUC0.store_vmax[i_spot]
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

      best=RUC.best_image
      C2 = deepcopy(C)
      if RUC.refine_Umatrix:
        ang, ax = RUC.get_correction_misset(as_axis_angle_deg=True)
        C2.rotate_around_origin(ax, ang)
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
      # Make sure refinement passes sanity checks
      


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
          # Also dumping important stuff for mcmc bit
          #dump_dict = {}
          #dump_dict['refined_ncells'] = refined_ncells
          #dump_dict['refined_scale'] = refined_scale
          #dump_dict['refined_gain'] = refined_gain
          #dump_dict['refined_local_spotscale'] = refined_local_spotscale
          #dump_dict['refined_tilt_abc'] = refined_tilt_abc
          #dump_dict['C'] = C
          #dump_dict['D'] = D
          #dump_dict['B'] = B
          #from libtbx.easy_pickle import dump
          #dump(os.path.join(outdir, 'RUC_info_%s.pickle'%ts), dump_dict )

      # Important here
      C=C2

    # FIXME final return here
    return RUC.f_vals[-1]

if __name__ == "__main__":
    import os, sys
    #ls49_data_dir='/Users/abhowmick/Desktop/software/dials/modules/LS49_regression/diffBragg_work/jungfrau_grid_search_4_or_more_regression/out_jungfrau_shoeboxes2'
    #ls49_data_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_work_on_refinement_07May20/out_jungfrau_shoeboxes_v3'
    ls49_data_dir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_setup_before_scattering_factor/out_jungfrau_shoeboxes_v2'
    ts='20180501114703722'
    #ts='20180501132216201'
    #ts='20180501120317142'
    # Debugging on 16 june regarding weird likelihood values
    #ts='20180501131709048' # ts with np.max(fe3-fe0) 
    #ts='20180501172057358'
    #ts='20180501172221032' 
    #ts='20180501112219788'
    #ts='20180501155759914' # 13 spots ?? --> wrong, this has diverged in dials refinement. c-axis is way off ; please remove
    #ts='20180501150413592' # near blank edge on RHS
    #ts='20180501165717967'
    outdir='/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_grid_search_4_or_more_regression/temp_2'
    #outdir=None
    #ls49_data_dir=None
    run_all_refine_ls49_JF1M(ts=ts, ls49_data_dir=ls49_data_dir, outdir=outdir, show_plotted_images=True, params=None, seed=3, short_circuit_dir=outdir)
