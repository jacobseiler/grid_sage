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

// Takes index of a galaxy (p) and maps the properties (SFR, StellarMass etc) onto the grid. //

void update_grid_properties(int p, int merged, int GridNr)
{

  int i, grid_position;
  for(i = 0; i < MAXSNAPS; ++i)
  {
    grid_position = GalGrid[p].History[i]; 
    if (grid_position == -1) continue; // If the galaxy hasn't formed yet, move on to next snapshot. 
    if (i > ListOutputGrid[GridNr]) break; // If the current snapshot redshift is lower than an output redshift, proceed to next one.

    //These are instantaneous properties measured at a specific redshift. //
    if((i == ListOutputGrid[GridNr]) && (GalGrid[p].StellarMass[i] > 0.0) & (GalGrid[p].SFR[i] > 0.0)) // Only want to take those galaxies with a non-zero stellar mass (they actually exist and are evolving).
    {
	
      Grid[grid_position].Sfr += GalGrid[p].SFR[i];
      Grid[grid_position].Count += 1;
      Grid[grid_position].StellarMass += GalGrid[p].StellarMass[i];
      if (PhotonPrescription == 1)  // If our photon prescription is using the galaxy photon calculated from the STARBURST spectra.
      {	
        Grid[grid_position].Nion_HI += pow(10,GalGrid[p].Photons_HI[i])*f_esc; 
        Grid[grid_position].Nion_HeI += pow(10,GalGrid[p].Photons_HeI[i])*f_esc; 
        Grid[grid_position].Nion_HeII += pow(10,GalGrid[p].Photons_HeII[i])*f_esc; 
      }
    }
      
    // Tracking the Merging Types. //
    if(i == ListOutputGrid[GridNr] && merged == 0)
    {
      ++GalaxyCount;
    }
  }
}

// Takes the number of halos loaded in memory (totNHalos) and maps the properties onto the grid. //

void update_grid_halo(int totNHalos, int GridNr)
{
  int i;
  int x_grid, y_grid, z_grid, grid_position; 

  for(i = 0; i < totNHalos; ++i)
  {
    if(Halo[i].SnapNum == ListOutputGrid[GridNr] && Halo[i].Mvir > 0.0) // Checks to make sure the halo existed at the output redshift.
    {

      x_grid = Halo[i].Pos[0]*GridSize/BoxSize; // Convert the (x,y,z) position to a grid (x,y,z).
      y_grid = Halo[i].Pos[1]*GridSize/BoxSize;
      z_grid = Halo[i].Pos[2]*GridSize/BoxSize; 

      grid_position = (x_grid*GridSize+y_grid)*GridSize+z_grid; // Convert the grid (x,y,z) to a 1D value.

      Grid[grid_position].HaloMass += Halo[i].Mvir;
      ++Grid[grid_position].HaloCount; 
    
    }
  }
}

// 

void update_grid_diffuse(int GridNr)
{
 
  char bufz0[1000];
  FILE* load_fd = NULL;


  sprintf(bufz0, "%s/%d_Unbound.dat", DiffuseDir, ListOutputGrid[GridNr]);
 
  if(!(load_fd = fopen(bufz0, "r")))
  {
    printf("Can't open file `%s'\n", bufz0);
    exit(0);
  }
   
  int i;
  printf("Reading in diffuse particles for Snapshot %d.\n", ListOutputGrid[GridNr]); 
  for (i = 0; i < CUBE(GridSize); ++i)
    fread(&Grid[i].Diffuse, sizeof(int), 1, load_fd);

//  fread(&GridDiffuse, sizeof(int), CUBE(GridSize), load_fd); 

  fclose(load_fd);
}

void update_grid_density(int GridNr)
{ 
  int i, count = 0;
  double total_density = 0, max_density = 0;
  double tot_StellarMass = 0, tot_DiffuseMass = 0, tot_HaloMass = 0;

  for(i = 0; i < CUBE(GridSize); ++i)
  {
    tot_StellarMass += Grid[i].StellarMass;
    tot_DiffuseMass += Grid[i].Diffuse*PartMass;
    tot_HaloMass += Grid[i].HaloMass;
    total_density += Grid[i].StellarMass + Grid[i].Diffuse*PartMass + Grid[i].HaloMass;  
  }
    
  printf("=======================================\n\n");
  printf("Calculating density for Snapshot %d (z = %.3f).\n", ListOutputGrid[GridNr], ZZ[ListOutputGrid[GridNr]]);
  printf("The total mass in the grid is %.4e (10^10 MSun/h)\n", total_density);
  printf("With %.4e Stellar, %.4e Diffuse and %.4e Halo Mass\n", tot_StellarMass, tot_DiffuseMass, tot_HaloMass);

  total_density /= CUBE(GridSize);

  printf("Average density for a grid cell is %.4f (10^10 MSun/h)\n", total_density);
  printf("One particle has a mass of %.4f (10^10MSun/h)\n", PartMass);
  printf("Hence we make a cell with no particles have a density of %.4f\n", PartMass/total_density);
    
  // Convert to overdensity. //
  for(i = 0; i < CUBE(GridSize); ++i)
  {
    Grid[i].Density = (Grid[i].StellarMass + Grid[i].Diffuse*PartMass + Grid[i].HaloMass)/total_density;
    if (Grid[i].Density > max_density) // Tracks maximum overdensity.
	    max_density = Grid[i].Density;
    if(Grid[i].Density == 0.0) // We don't want a grid cell to have zero density so if it does, just populate it with a single dark matter particle.
    {
      ++count;
      Grid[i].Density = PartMass/total_density;
    }
      
  }

  printf("%d cells with zero Density\n", count);
  printf("The maximum density was %.4f\n", max_density); 
  printf("=======================================\n\n");

}

void update_grid_nion_halo(int GridNr) // Calculates number of ionizing photons using the halos. 
{
  int i;
  double totHaloMass;
  double evolve_time;

  double UnitConversion;

  printf("Calculating number of photons using halo-based prescription.\n");
  
 
  if (GridNr == NGrid - 1)
	  evolve_time = (Age[ListOutputGrid[GridNr]-1] - Age[ListOutputGrid[GridNr]]) * UnitTime_in_Megayears;
  else
	  evolve_time = (Age[ListOutputGrid[GridNr+1]] - Age[ListOutputGrid[GridNr]]) * UnitTime_in_Megayears;

  totHaloMass = 0;
  for (i = 0; i < CUBE(GridSize); ++i)
  {
    // Using Illiev 2012. 
    UnitConversion = SOLAR_MASS/0.53/PROTONMASS/SEC_PER_MEGAYEAR;
    Grid[i].Nion_HI = Grid[i].HaloMass*1e10/Hubble_h*SOLAR_MASS*SourceEfficiency*BaryonFrac/0.53/PROTONMASS/evolve_time/SEC_PER_MEGAYEAR;
  }
  
}

void count_grid_properties(int GridNr) // Count number of galaxies/halos in the grid.
{
  int i;

  int GalCount = 0, HaloCount = 0, DiffuseCount = 0;
  double totPhotons_HI = 0, totPhotons_HeI = 0, totPhotons_HeII = 0, HaloMass = 0;
  for (i = 0; i < CUBE(GridSize); ++i)
  {
    GalCount += Grid[i].Count;
    HaloCount += Grid[i].HaloCount;
    totPhotons_HI += Grid[i].Nion_HI;
    totPhotons_HeI += Grid[i].Nion_HeI;
    totPhotons_HeII += Grid[i].Nion_HeII;
    HaloMass += Grid[i].HaloMass;
    DiffuseCount += Grid[i].Diffuse;
  }

  printf("At redshift %.3f (Snapshot %d) there was %d galaxies, %d halos and %d diffuse particles with [%.4e, %.4e, %.4e] {HI, HeI, HeII}  ionizing Photons emitted per second ([%.4e %.4e %.4e] s^-1 Mpc^-3) and %.4e (Msun) halo mass.\n", ZZ[ListOutputGrid[GridNr]], ListOutputGrid[GridNr], GalCount, HaloCount, DiffuseCount, totPhotons_HI, totPhotons_HeI, totPhotons_HeII, totPhotons_HI / pow(BoxSize/Hubble_h,3), totPhotons_HeI / pow(BoxSize/Hubble_h, 3), totPhotons_HeII / pow(BoxSize/Hubble_h,3), HaloMass*1e10/Hubble_h);

}
