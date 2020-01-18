#
# MPI version of load_ls49
#


from LS49_project.diffBragg_work.load_ls49 import strong_spot_mask, process_ls49_image_real


class Script(object):
  def __init__(self):
    """initializer """

  def run():
    

if __name__ == "__main__":
    from argparse import ArgumentParser
    from simtbx.diffBragg.refiners import RefineAll
    from simtbx.diffBragg.sim_data import SimData
    import numpy as np
    from simtbx.diffBragg.refiners.crystal_systems import MonoclinicManager
    from simtbx.diffBragg import nanoBragg_crystal, nanoBragg_beam
    from copy import deepcopy
    from dxtbx.model.experiment_list import Experiment, ExperimentList, ExperimentListFactory

    parser = ArgumentParser()
    parser.add_argument("--plot", action='store_true')
    parser.add_argument("--scaleonly", action='store_true')
    args = parser.parse_args()


    # Initial r0222 regression
    timestamps_of_interest = ['20180501143533988', # bad
                              '20180501143701853'] # Does not work

    ts = timestamps_of_interest[0]
    data = process_ls49_image_real(tstamp=ts,Nstrongest=10, resmin=2.0, resmax=13.5)

    C = data["dxcrystal"]
    D = data["dxdetector"]
    B = data["dxbeam"]
    exp=data['experiment']
    dump_exp = Experiment(imageset=exp.imageset, 
                         beam=B,
                         detector=D,
                         goniometer=exp.goniometer,
                         scan=exp.scan,
                         crystal=C)
    dump_explist = ExperimentList([exp])
    dump_explist.as_file('before_refinement_%s.expt'%ts)

    # Some global variables here for LS49
    mos_spread_deg=0.01
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
        for x1, x2, y1, y2 in data["bboxes_x1x2y1y2"]:
            patch = plt.Rectangle(
                width=x2-x1,
                height=y2-y1,
                xy=(x1, y1),
                ec='r', fc='none')
            plt.gca().add_patch(patch)
        plt.show()

    RUC = RefineAll(
        spot_rois=data["bboxes_x1x2y1y2"],
        spot_resolution=data['resolution'],
        abc_init=data["tilt_abc"],
        img=data["data_img"],
        SimData_instance=SIM,
        plot_images=args.plot,
        ucell_manager=UcellMan,
        init_gain=1,
        init_scale=1.0)
    RUC.trad_conv = True
    RUC.trad_conv_eps = 1e-5
    RUC.max_calls = 250
    RUC.refine_background_planes = False
    RUC.refine_Umatrix = False
    RUC.refine_Bmatrix = False
    RUC.verbose=True
    RUC.plot_stride=1
    RUC.refine_background_planes = False
    RUC.refine_ncells=False #
    RUC.refine_crystal_scale = False # scale factor
    RUC.refine_gain_fac = False
    RUC.refine_detdist=False
    RUC.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
    RUC.use_curvatures=False #
    RUC.calc_curvatures=True # if set to False, never uses curvatures
    RUC.refine_with_restraints=True
    #RUC.run()
    #if RUC.hit_break_to_use_curvatures:
    #  RUC.num_positive_curvatures=0
    #  RUC.use_curvatures=True
    #  RUC.run(setup=False) # Now it will use curvatures

    # First refine scale factors and ncells, then refine everything else
    RUC.refine_ncells=True
    RUC.refine_crystal_scale=True
    RUC.run()
    refined_ncells = RUC.x[7]
    refined_scale = RUC.x[10]

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
    print (' ================== now refining all the variables ======================= ')
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
        abc_init=data["tilt_abc"],
        img=data["data_img"],
        SimData_instance=SIM2,
        plot_images=args.plot,
        ucell_manager=UcellMan,
        init_gain=1,
        init_scale=refined_scale)
    RUC2.trad_conv = True
    RUC2.trad_conv_eps = 1e-5
    RUC2.max_calls = 250
    RUC2.refine_background_planes = False
    RUC2.refine_Umatrix = True
    RUC2.refine_Bmatrix = True
    RUC2.verbose=True
    RUC2.plot_stride=1
    RUC2.refine_background_planes = False
    RUC2.refine_ncells=True #
    RUC2.refine_crystal_scale = True # scale factor
    RUC2.refine_gain_fac = False
    RUC2.refine_detdist=False
    RUC2.use_curvatures_threshold=7 # Keep calculating them and after 7 times switch over to using them
    RUC2.use_curvatures=False #
    RUC2.calc_curvatures=True # if set to False, never uses curvatures
    RUC2.run()
    #from IPython import embed; embed(); exit()

    print("Done.")
    print("Refined scale =%f", RUC.x[-1])

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
    dump_exp = Experiment(imageset=exp.imageset, 
                         beam=B,
                         detector=D,
                         goniometer=exp.goniometer,
                         scan=exp.scan,
                         crystal=C2)
    dump_explist = ExperimentList([exp])
    dump_explist.as_file('after_refinement_%s.expt'%ts)

