#!/usr/bin/env python
from __future__ import print_function

import yt
from yt.analysis_modules.halo_finding.halo_objects import FOFHaloFinder 
yt.enable_parallelism()

def create_fof_yt(SnapNum, SnapDir, OutDir):
    """
    Creates a FOF catalogue using yt.

    Parameters
    ----------

    SnapNum: Integer. Required.
        Snapshot number to create the catalogue for.

    OutDir: String. Required.
        Path where the catalogue is written to. 

    Returns
    ----------

    None. The FOF catalogue is written to the directory specified by
    ``outdir``.
    """

    if SnapNum < 66: # Slightly different naming scheme depending upon the snapshot number.
        snapshot_name = "snapdir_{0:03d}/snapdir_{0:03d}/snapshot_{0:03d}.0.hdf5".format(SnapNum)
    else:
        snapshot_name = "snapdir_{0:03d}/snapshot_{0:03d}.0.hdf5".format(SnapNum)

    fname = "{0}/{1}".format(SnapDir, snapshot_name)

    unit_base = {'length' : (1.0, "Mpc/h"), 'mass' : (1.0e10, "Msun/h"), 'velocity' : (1.0, "km/s")} # Set the yt units to the GADGET ones.

    ds = yt.load(fname)
    halos = FOFHaloFinder(ds, link = -0.0049, dm_only = False, ptype = "PartType1", padding = 0.2)

    outname = "{0}/yt_fof_tab_{1}".format(OutDir, SnapNum)
    halos.write_particle_lists(outname)
    print("Successfully wrote to {0}".format(outname))

if __name__ == "__main__":

    SnapDir = "/lustre/projects/p134_swin/jseiler/simulations/1.6Gpc/means/halo_1721673/dm_gadget/data/"
    yt_outdir = "/lustre/projects/p004_swin/jseiler/yt_fofs"

    create_fof_yt(30, SnapDir, yt_outdir) 
