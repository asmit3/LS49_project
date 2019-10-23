#!/bin/bash
# copy over the following after a phenix.refine run
# the PDB from the phenix.refine folder ==> name it ls49.pdb
# the MTZ from the phenix.refine folder ==> name it as ls49.mtz
# the Rfree mtz, usually called *_data.mtz in the phenix.refine folder ==> rfree.mtz

trial=ls49

libtbx.python ../../exafel_project/nks/map_height_at_atoms.py \
${trial}.pdb \
${trial}.mtz \
input.xray_data.labels=I-obs \
input.symmetry_safety_check=warning \
xray_data.r_free_flags.file_name=rfree.mtz \
xray_data.r_free_flags.label=R-free-flags \
selection="element Fe" |tee ${trial}_peakht.log
