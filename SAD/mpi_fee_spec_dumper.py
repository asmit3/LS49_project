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

logger = logging.getLogger("LS49_project.mpi_fee_spec_dumper")


help_message = """
Input will be a text file with full path of refined.expt files
output will be a fee spec pickle file


If running on NERSC, please source these 2 environment variables
Only then can you save the fee data file otherwise the xtc file cannot be read
export SIT_PSDM_DATA=/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3
export SIT_DATA=/global/cscratch1/sd/asmit/LS10/psdm/psdm/data

"""

fee_spec_phil_str= """
  fee_spec {
    refined_expt_fnames = None 
      .type = str
      .help = Text file with refined expt filenames (with absolute path). Will be used to backout input indexed.refl and refined.expt filenames for integration
    run_number = 222
      .type = int
      .help = LS49 run number to dump fee spectra information
}
 
"""


phil_scope = parse(fee_spec_phil_str+control_phil_str + dials_phil_str, process_includes=True).fetch(
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

        with open(self.params.fee_spec.refined_expt_fnames) as fin:
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
                processor = Processor_fee_spec(
                    copy.deepcopy(params), composite_tag="%04d" % i, rank=i
                )
                for item in item_list:
                    processor.process_experiments(item)
                processor.finalize()

            iterable = list(zip(tags, all_refined_expt_files, all_indexed_refl_files))

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


class Processor_fee_spec(Processor):
    def process_experiments(self, item):
        tag, refined_expt, indexed_refl = item
        self.tag = tag
        self.debug_start(tag)
        # Now write out spectra information
        FEE=easy_pickle.load('/global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/input/FEE_r143_v1.pickle')
        FEE.run=self.params.fee_spec.run_number
        FEE.save_spectra_of_select_events([tag], fname=os.path.join(self.params.output.output_dir, 'idx-fee_data_%s.pickle'%(tag))) 

if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        script = Script()
        script.run()
