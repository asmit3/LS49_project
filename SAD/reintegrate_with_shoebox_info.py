from __future__ import absolute_import, division, print_function

import copy
import logging
import os, glob
import sys
import tarfile
import time
import six.moves.cPickle as pickle
from six.moves import StringIO

import dials.util
from dials.util import log
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model.experiment_list import ExperimentList
from libtbx.utils import Abort, Sorry
from dials.command_line.stills_process import Processor, do_import
from dials.command_line.stills_process import control_phil_str, dials_phil_str, program_defaults_phil_str
from libtbx.phil import parse
import shutil
from libtbx import easy_pickle
from LS49_project.FEE_spec.fee_spectra import FEE_spec

logger = logging.getLogger("LS49_project.reintegrate_with_shoebox_info")


help_message = """
Input will be a text file with full path of refined.expt files
Output will be a)integrated.refl and b) integrated.expt files with shoebox information saved
Additionally saves (c) refined.expt (d) indexed.refl (e) CBF image file (f) fee_data pickle file

Useful when initial integration did not save the shoebox information

If running on NERSC, please source these 2 environment variables
Only then can you save the fee data file otherwise the xtc file cannot be read
export SIT_PSDM_DATA=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3
export SIT_DATA=/global/cscratch1/sd/asmit/LS10/psdm/psdm/data

"""

reintegrate_phil_str= """
  reintegrate {
    refined_expt_fnames = None 
      .type = str
      .help = Text file with refined expt filenames (with absolute path). Will be used to backout input indexed.refl and refined.expt filenames for integration
    run_number = 222
      .type = int
      .help = LS49 run number to dump fee spectra information
}
 
"""


phil_scope = parse(reintegrate_phil_str+control_phil_str + dials_phil_str, process_includes=True).fetch(
    parse(program_defaults_phil_str)
)


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser

        # The script usage
        usage = "usage: dials.stills_process [options] [param.phil] filenames"

        self.tag = None
        self.reference_detector = None

        # Create the parser
        self.parser = OptionParser(usage=usage, phil=phil_scope, epilog=help_message)

    def load_reference_geometry(self):
        if self.params.input.reference_geometry is None:
            return

        from dxtbx.model.experiment_list import ExperimentListFactory

        try:
            ref_experiments = ExperimentListFactory.from_json_file(
                self.params.input.reference_geometry, check_format=False
            )
        except Exception:
            try:
                import dxtbx

                img = dxtbx.load(self.params.input.reference_geometry)
            except Exception:
                raise Sorry(
                    "Couldn't load geometry file %s"
                    % self.params.input.reference_geometry
                )
            else:
                self.reference_detector = img.get_detector()
        else:
            assert len(ref_experiments.detectors()) == 1
            self.reference_detector = ref_experiments.detectors()[0]

    def run(self):
        """Execute the script."""
        from libtbx import easy_mp

        # Parse the command line
        params, options, all_paths = self.parser.parse_args(
            show_diff_phil=False, return_unhandled=True, quick_parse=True
        )

        if not all_paths and params.input.file_list is not None:
            all_paths.extend(
                [path.strip() for path in open(params.input.file_list).readlines()]
            )

        # Check we have some filenames
        # Asmit - dont need filenames here
        #if not all_paths:
        #    self.parser.print_help()
        #    return

        # Mask validation
        for mask_path in params.spotfinder.lookup.mask, params.integration.lookup.mask:
            if mask_path is not None and not os.path.isfile(mask_path):
                raise Sorry("Mask %s not found" % mask_path)

        # Save the options
        self.options = options
        self.params = params

        st = time.time()

        # Configure logging
        #log.config(verbosity=options.verbose, logfile="dials.process.log")

        bad_phils = [f for f in all_paths if os.path.splitext(f)[1] == ".phil"]
        if len(bad_phils) > 0:
            self.parser.print_help()
            logger.error(
                "Error: the following phil files were not understood: %s"
                % (", ".join(bad_phils))
            )
            return

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        for abs_params in self.params.integration.absorption_correction:
            if abs_params.apply:
                if not (
                    self.params.integration.debug.output
                    and not self.params.integration.debug.separate_files
                ):
                    raise Sorry(
                        "Shoeboxes must be saved to integration intermediates to apply an absorption correction. "
                        + "Set integration.debug.output=True, integration.debug.separate_files=False and "
                        + "integration.debug.delete_shoeboxes=True to temporarily store shoeboxes."
                    )

        # Import stuff
        logger.info("Loading files...")

        tags = []
        all_refined_expt_files = []
        all_indexed_refl_files = []

        with open(self.params.reintegrate.refined_expt_fnames) as fin:
          for line in fin:
            if line !='/n':
              refined_expt=line.strip()
              ts_index=refined_expt.index('idx-2018')
              ts = refined_expt[ts_index+4:ts_index++21]
              print (refined_expt, ts)
              indexed_refl = refined_expt.split('_refined.expt')[0]+'_indexed.refl'
              all_indexed_refl_files.append(indexed_refl)
              all_refined_expt_files.append(refined_expt)
              tag=ts
              tags.append(tag)
  
        if True:
            # Wrapper function
            def do_work(i, item_list):
                processor = Processor_reintegrate(
                    copy.deepcopy(params), composite_tag="%04d" % i, rank=i
                )
                for item in item_list:
                    processor.process_experiments(item)
                processor.finalize()

            iterable = list(zip(tags, all_refined_expt_files, all_indexed_refl_files))

        import pdb; pdb.set_trace()
        # Process the data
        if params.mp.method == "mpi":
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
            size = comm.Get_size()  # size: number of processes running in this job

            # Configure the logging
            if params.output.logging_dir is None:
                logfile = None
            else:
                log_path = os.path.join(
                    params.output.logging_dir, "log_rank%04d.out" % rank
                )
                error_path = os.path.join(
                    params.output.logging_dir, "error_rank%04d.out" % rank
                )
                print("Redirecting stdout to %s" % log_path)
                print("Redirecting stderr to %s" % error_path)
                sys.stdout = open(log_path, "a", buffering=0)
                sys.stderr = open(error_path, "a", buffering=0)
                print("Should be redirected now")

                logfile = os.path.join(
                    params.output.logging_dir, "info_rank%04d.out" % rank
                )

            log.config(verbosity=options.verbose, logfile=logfile)

            if size <= 2:  # client/server only makes sense for n>2
                subset = [
                    item for i, item in enumerate(iterable) if (i + rank) % size == 0
                ]
                do_work(rank, subset)
            else:
                if rank == 0:
                    # server process
                    for item in iterable:
                        print("Getting next available process")
                        rankreq = comm.recv(source=MPI.ANY_SOURCE)
                        print("Process %s is ready, sending %s\n" % (rankreq, item[0]))
                        comm.send(item, dest=rankreq)
                    # send a stop command to each process
                    print("MPI DONE, sending stops\n")
                    for rankreq in range(size - 1):
                        rankreq = comm.recv(source=MPI.ANY_SOURCE)
                        print("Sending stop to %d\n" % rankreq)
                        comm.send("endrun", dest=rankreq)
                    print("All stops sent.")
                else:
                    # client process
                    while True:
                        # inform the server this process is ready for an event
                        print("Rank %d getting next task" % rank)
                        comm.send(rank, dest=0)
                        print("Rank %d waiting for response" % rank)
                        item = comm.recv(source=0)
                        if item == "endrun":
                            print("Rank %d received endrun" % rank)
                            break
                        print("Rank %d beginning processing" % rank)
                        try:
                            do_work(rank, [item])
                        except Exception as e:
                            print(
                                "Rank %d unhandled exception processing event" % rank,
                                str(e),
                            )
                        print("Rank %d event processed" % rank)
        else:
            from dxtbx.command_line.image_average import splitit

            if params.mp.nproc == 1:
                do_work(0, iterable)
            else:
                result = list(
                    easy_mp.multi_core_run(
                        myfunction=do_work,
                        argstuples=list(enumerate(splitit(iterable, params.mp.nproc))),
                        nproc=params.mp.nproc,
                    )
                )
                error_list = [r[2] for r in result]
                if error_list.count(None) != len(error_list):
                    print(
                        "Some processes failed excecution. Not all images may have processed. Error messages:"
                    )
                    for error in error_list:
                        if error is None:
                            continue
                        print(error)

        # Total Time
        logger.info("")
        logger.info("Total Time Taken = %f seconds" % (time.time() - st))


class Processor_reintegrate(Processor):
    def setup_filenames(self, tag):
        # before processing, set output paths according to the templates
        if (
            self.experiments_filename_template is not None
            and "%s" in self.experiments_filename_template
        ):
            self.params.output.experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.experiments_filename_template % ("idx-" + tag),
            )
        if (
            self.strong_filename_template is not None
            and "%s" in self.strong_filename_template
        ):
            self.params.output.strong_filename = os.path.join(
                self.params.output.output_dir,
                self.strong_filename_template % ("idx-" + tag),
            )
        if (
            self.indexed_filename_template is not None
            and "%s" in self.indexed_filename_template
        ):
            self.params.output.indexed_filename = os.path.join(
                self.params.output.output_dir,
                self.indexed_filename_template % ("idx-" + tag),
            )
        if (
            self.refined_experiments_filename_template is not None
            and "%s" in self.refined_experiments_filename_template
        ):
            self.params.output.refined_experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.refined_experiments_filename_template % ("idx-" + tag),
            )
        if (
            self.integrated_filename_template is not None
            and "%s" in self.integrated_filename_template
        ):
            self.params.output.integrated_filename = os.path.join(
                self.params.output.output_dir,
                self.integrated_filename_template % ("idx-" + tag),
            )
        if (
            self.integrated_experiments_filename_template is not None
            and "%s" in self.integrated_experiments_filename_template
        ):
            self.params.output.integrated_experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.integrated_experiments_filename_template % ("idx-" + tag),
            )

    def process_experiments(self, item):

        tag, refined_expt, indexed_refl = item

        if not self.params.output.composite_output:
            self.setup_filenames(tag)
        
        self.tag = tag
        self.debug_start(tag)

        # Do the processing
        #from IPython import embed; embed(); exit()
        experiments=ExperimentListFactory.from_json_file(refined_expt)
        indexed=flex.reflection_table.from_file(indexed_refl)

        self.debug_write("refine_start")
        try:
            experiments, indexed = self.refine(experiments, indexed)
        except Exception as e:
            print("Error refining", tag, str(e))
            self.debug_write("refine_failed_%d" % len(indexed), "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        try:
            if self.params.dispatch.integrate:
                # change phil params for shoeboxes to ensure they get written out
                self.params.integration.debug.delete_shoeboxes=False 
                self.params.integration.debug.output=True
                self.params.integration.debug.separate_files=False
                self.debug_write("integrate_start")
                integrated = self.integrate(experiments, indexed)
            else:
                print("Integration turned off. Exiting")
                self.debug_write("index_ok_%d" % len(indexed), "done")
                return
        except Exception as e:
            print("Error integrating", tag, str(e))
            self.debug_write("integrate_failed_%d" % len(indexed), "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        self.debug_write("integrate_ok_%d" % len(integrated), "done")
        # Now copy CBF file from out or out_reint directory
        if 'out_reint' in refined_expt:
          refined_expt=refined_expt.replace('out_reint', 'out')
        cbf_file_path_and_name = refined_expt.replace('_refined.expt', '.cbf')
        ts_index = cbf_file_path_and_name.index('idx-2018')
        cbf_file_name=cbf_file_path_and_name[ts_index:]
        dest = shutil.copyfile(cbf_file_path_and_name, os.path.join(self.params.output.output_dir, cbf_file_name)) 
        # Now write out spectra information
        FEE=easy_pickle.load('/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/input/FEE_r143_v1.pickle')
        FEE.run=self.params.reintegrate.run_number
        FEE.save_spectra_of_select_events([tag], fname=os.path.join(self.params.output.output_dir, 'idx-fee_data_%s.pickle'%(tag))) 

if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        script = Script()
        script.run()
