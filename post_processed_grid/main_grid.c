#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

char bufz0[1000];

int main(int argc, char **argv)
{

  int filenr, p, i, GridNr;
 
  read_parameter_file(argv[1]);
  init(); // Initialize all the parameters (set units, create scale factor/age arrays etc).

  for (GridNr = 0; GridNr < NGrid; ++GridNr)
  {
    init_grid(GridNr); // Initialize the grid.  
//    update_grid_diffuse(GridNr); // Read in all the diffuse gas.
    for (filenr = FirstFile; filenr < LastFile + 1; ++filenr)
    {
      if(filenr % 30 == 0)
      printf("Doing file %d.\n", filenr);

    //  load_halos(filenr); // Load the halos.
    //  update_grid_halo(totNHalos, GridNr); // Update the properties associated with halos
    //  myfree(Halo); // Don't need halos anymore.

      for (i = 0; i < 2; ++i) // i = 0 does the normal galaxies, i = 1 does the merged galaxies.
      {
        if(i == 0)      
	  sprintf(bufz0, "%s/%s_%d", GalaxiesInputDir, FileNameGalaxies, filenr);
        else       
	  sprintf(bufz0, "%s/%s_%d", GalaxiesInputDir, FileNameMergedGalaxies, filenr); 

        if ( access(bufz0, F_OK ) == -1) // Sanity check.
        {
          printf("-- input for file %s does not exist, exiting now.\n", bufz0);
          exit(0); 
        }

//        if (Verbose == 1)
//         printf("Loading galaxies for file %d, name '%s'\n", filenr, bufz0); 
        load_gals(bufz0);    

//        if (Verbose == 1)
//          printf("Gridding properties.\n");

        for(p = 0; p < NtotGals; ++p)
        {    	
          update_grid_properties(p, i, GridNr); // Go through each galaxy and read it's grid history and grid the properties.
        }

	free_gals();	

      }

//      printf("Done File %d.\n\n", filenr);
    }
  
      if (fescPrescription == 1)
	    normalize_photon(GridNr);
      if (fescPrescription == 4)
            normalize_slope_photons(GridNr);

//    update_grid_density(GridNr); // Go through the grid and convert to overdensity.
    if (PhotonPrescription == 0) 
      update_grid_nion_halo(GridNr); // If we are using a halo-based photon prescription, calculate the number of photons based on halo mass.

    if (Verbose == 1)
      count_grid_properties(GridNr); // Counts how many halos/galaxies/Photons are in the grid at each redshift.
    save_grid(GridNr); // Saves grid.
 
  }

  myfree(Grid);
  save_redshift();

  return 0;

} 
