#
#
#
from load_ls49_JF1M import strong_spot_mask, process_ls49_image_real

def run_all_refine_ls49_JF1M(ts=None, ls49_data_dir=None, show_plotted_images=False, outdir=None, params=None, seed=0, short_circuit_dir=None, swap_spectra_timestamp=False):
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

    data = process_ls49_image_real(tstamp=ts,Nstrongest=10, resmin=2.1, resmax=3.5, ls49_data_dir=ls49_data_dir, seed=seed, swap_spectra_timestamp=swap_spectra_timestamp)
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
    run_all_refine_ls49_JF1M(ts=ts, ls49_data_dir=ls49_data_dir, outdir=outdir, show_plotted_images=True, params=None, seed=3, short_circuit_dir=outdir, swap_spectra_timestamp=False)
