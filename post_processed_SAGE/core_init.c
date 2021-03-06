#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "core_allvars.h"
#include "core_proto.h"

void init(void)
{
  int i;

  count_gal = 0;
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	 // start-up seed 

  set_units();
  srand((unsigned) time(NULL));

  read_snap_list();

  for(i = 0; i < Snaplistlen; i++)
  {
    ZZ[i] = 1 / AA[i] - 1;
    Age[i] = time_to_present(ZZ[i]);

  }

  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);

  read_cooling_functions();

  count_onehalo = 0;  

  zeromass_count = 0;
  suppression_count = 0;
  previous_tree = 0;

  smallest_mass = 100000;
  lowmass_halo = 0;

  outside_box = 0;
  inside_box = 0;

  count_Mvir = 0;
  count_Len = 0;
  if (IMF == 1)
  {
    // Salpeter IMF //
    IMF_norm =0.1706;
    IMF_slope = -2.35;
    Eta_SNII = 7.432e-3; //Number fraction of stars that will end their lives as type II supernovae.
    m_SNII = 0.144; // Mass fraction of stars that will end their lives as type II supernovae.

  } else if (IMF == 2)
  {
    // Chabrier IMF //
    IMF_norm = 0.23638; 
    IMF_slope = -2.3;
    Eta_SNII = 1.1893e-2; //Number fraction of stars that will end their lives as type II supernovae.
    m_SNII = 0.23638; // Mass fraction of stars that will end their lives as type II supernovae.

  }
 
  // Testing Parameters //

  /*
  V_energy = 70.0;
  alpha_energy = 1.80; 
  beta_energy = 16.0; 

  V_mass = 80.0; 
  alpha_mass = 40.0; 
  beta_mass = 16.0; 
  */

  // Tiamat Parameters //

  V_energy = 70.0;
  alpha_energy = 0.5;
  //beta_energy = 2.0;
  beta_energy = 10.0;

  V_mass = 70.0; 
  alpha_mass = 6.0; 
  //beta_mass = 0.0; 
  beta_mass = 10.0; 

  epsilon_mass_max = 30.0;

  //

  count_firstSF = 0;
  count_notfirstSF = 0;

  if(TimeResolutionSN > 50)
  {
    fprintf(stderr, "The selected time resolution for SN feedback (TimeResolutionSN) is set too high (%d Myr).  Using TimeResolutionSN > 50Myr is the same as using the instantaneous recycling approximation; set 'IRA' to 1 instead!\n", TimeResolutionSN); 
    ABORT(EXIT_FAILURE);  
  } else if(TimeResolutionSN > 35)
  {
    fprintf(stderr, "Your selected time resolution for SN feedback (TimeResolutionSN) is quite high (%d Myr).  Beyond 50Myr the instantaneous recycling approximation is valid hence with your value it would likely be correct to set 'IRA' to 1.\n", TimeResolutionSN);
  } else
  {
    Time_SFH = 0;
    SN_Array_Len = 0;
    while(Time_SFH < 50)
    {
      Time_SFH += TimeResolutionSN;
      ++SN_Array_Len;
    }

    //SN_Array_Len = 50;
  }
  
  mergedgal_mallocs = 0;
  gal_mallocs = 0 ;

  mergedgal_frees = 0;
  gal_frees = 0;
 
  small_steps = 0;
  large_steps = 0;
  times_written = 0; 
}


int32_t init_grid()
{

  // There are two modes of operation for accounting for the photoionization feedback. //
  // In the first we read the photoionization rates and redshift of reionization for ALL redshifts. //
  // Feedback is then applied to all of redshifts.  This has extensive memory requirements as assuming a 512^3 grid to float precision, then each grid will be ~0.5GB in size. //


  int32_t i; 
  FILE *photoion, *reionredshift;
  char buf[MAXLEN];

  printf("Reading in the data for the reionization feedback using cifog.\n");
  Grid = mymalloc(sizeof(struct GRID_STRUCT));
  if (Grid == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the high level grid struct.\n");
    return EXIT_FAILURE;
  }

  Grid->GridSize = GridSize;
  Grid->NumCellsTotal = CUBE(GridSize); 

  Grid->ReionRedshift = malloc(sizeof(*(Grid->ReionRedshift)) * Grid->NumCellsTotal);
  if (Grid->ReionRedshift == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the reionization redshift grid.\n");
    return EXIT_FAILURE;
  } 

  snprintf(buf, MAXLEN, "%s/%s", PhotoionDir, ReionRedshiftName); 
  if(!(reionredshift= fopen(buf, "rb")))
  {
    fprintf(stderr, "Cannot open file %s\n", buf);
    return EXIT_FAILURE;
  }

  printf("Reading the reionization redshift grid from %s\n", buf);
  fread(Grid->ReionRedshift, sizeof(*(Grid->ReionRedshift)), Grid->NumCellsTotal, reionredshift);
  fclose(reionredshift);
   
  Grid->NumGrids = MAXSNAPS;
    
  Grid->PhotoGrid = malloc(sizeof(struct PHOTO_GRID) * Grid->NumGrids); // Allocate enough memory to hold the photoionization grids.
  if (Grid->PhotoGrid == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the photoionization grid struct\n");
    return EXIT_FAILURE;
  }

  for (i = 0; i < Grid->NumGrids; ++i)
  {
    Grid->PhotoGrid[i].PhotoRate = malloc(sizeof(*(Grid->PhotoGrid[i].PhotoRate)) * Grid->NumCellsTotal); // Then allocate memory for each photoionization grid.
    if (Grid->PhotoGrid[i].PhotoRate == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for the photoionization grid.\n");
      return EXIT_FAILURE;
    }

 
    // For some of the early snapshots we don't have photoionization grids (because ionization hasn't started yet at z=100).  
    // Let's put a flag to know whether we have any valid data for this snapshot so we don't have to create empty grids.
    if (i > LowSnap && i <= HighSnap)
    {
      Grid->PhotoGrid[i].valid_grid = 1;
    }
    else
    { 
      Grid->PhotoGrid[i].valid_grid = 0;
      printf("Snapshot %d is not a valid snapshot for reionization -- SKIPPING! --\n", i);
      continue;
    }
    snprintf(buf, MAXLEN, "%s/%s_%03d", PhotoionDir, PhotoionName, i);

    if(!(photoion = fopen(buf, "rb")))
    {
      fprintf(stderr, "Cannot open file %s\n", buf);
      return EXIT_FAILURE;
    }

    printf("Reading photoionization grid %s\n", buf);    
    fread(Grid->PhotoGrid[i].PhotoRate, sizeof(*(Grid->PhotoGrid[i].PhotoRate)), Grid->NumCellsTotal, photoion);
    fclose(photoion);

  }

  printf("All reionization feedback stuff read in successfully\n");

  return EXIT_SUCCESS;   
}

int32_t init_reion_lists(int32_t filenr)
{

  FILE *ListFile;
  char ListFile_name[MAXLEN];
  int32_t SnapNum, SnapNum_Read;
 
  if (ReionSnap == LowSnap) // This is the first iteration of the self-consistent run so there will be no ionization yet.
  {
    return EXIT_SUCCESS;
  }

  snprintf(ListFile_name, MAXLEN, "%s/reionization_modifiers/treefile_%03d", PhotoionDir, filenr);    

  ListFile = fopen(ListFile_name, "rb");
  if (ListFile == NULL)
  {
    fprintf(stderr, "Cannot open file %s\n", ListFile_name);
    return EXIT_FAILURE;
  }
 
  printf("Reading in the reionization modifier lists.\n");
 
  ReionList = malloc(sizeof(struct REIONMOD_STRUCT));
  if (ReionList == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for Reionization List struct\n");
    return EXIT_FAILURE;
  }

  ReionList->NumLists = ReionSnap; 

  ReionList->ReionMod_List = malloc(sizeof(struct REIONMOD_LIST) * ReionList->NumLists);

  for (SnapNum = 0; SnapNum < ReionList->NumLists; ++SnapNum)
  {

    fread(&SnapNum_Read, sizeof(int32_t), 1, ListFile);

    fread(&ReionList->ReionMod_List[SnapNum].NHalos_Ionized, sizeof(int32_t), 1, ListFile);
    //printf("Snapshot %d has %d Halos in the list.\n", SnapNum_Read, ReionList->ReionMod_List[SnapNum].NHalos_Ionized);

    if (SnapNum_Read != SnapNum)
    { 
      fprintf(stderr, "When attempting to read the reionization modifier lists, the read file had a snapshot number %d when we expected a number %d\n", SnapNum_Read, SnapNum);
      return EXIT_FAILURE;
    }

    ReionList->ReionMod_List[SnapNum].NHalos_Found = 0;
 
    if (ReionList->ReionMod_List[SnapNum].NHalos_Ionized == 0) // There were no halos within ionized regions for this snapshot, reionization hasn't started or is in the beginning.
    {
      continue;
    }
   
    ReionList->ReionMod_List[SnapNum].HaloID = malloc(sizeof(*(ReionList->ReionMod_List[SnapNum].HaloID)) * ReionList->ReionMod_List[SnapNum].NHalos_Ionized);
    if (ReionList->ReionMod_List[SnapNum].HaloID == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for the HaloIDs for the Reionization Modifier lists for snapshot %d\n", SnapNum);
      return EXIT_FAILURE;
    }

    fread(ReionList->ReionMod_List[SnapNum].HaloID, sizeof(*(ReionList->ReionMod_List[SnapNum].HaloID)), ReionList->ReionMod_List[SnapNum].NHalos_Ionized, ListFile);

    /*
    int32_t i;
    for (i = 0; i < ReionList->ReionMod_List[SnapNum].NHalos_Ionized; ++i)
    {
      int64_t ID;
      ID = ReionList->ReionMod_List[SnapNum].HaloID[i];
      printf("File %d: HaloID %ld is in the list, corresponding to tree %d and Halo number %d\n", filenr, ID, (int32_t)(ID >> 32), (int32_t)ID); 
    }
    */

    ReionList->ReionMod_List[SnapNum].ReionMod = malloc(sizeof(*(ReionList->ReionMod_List[SnapNum].ReionMod)) * ReionList->ReionMod_List[SnapNum].NHalos_Ionized);
    if (ReionList->ReionMod_List[SnapNum].ReionMod == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for the ReionMods for the Reionization Modifier lists for snapshot %d\n", SnapNum);
      return EXIT_FAILURE;
    }

    fread(ReionList->ReionMod_List[SnapNum].ReionMod, sizeof(*(ReionList->ReionMod_List[SnapNum].ReionMod)), ReionList->ReionMod_List[SnapNum].NHalos_Ionized, ListFile);
  
    
  }
  fclose(ListFile);
 
  return EXIT_SUCCESS;

}


void set_units(void)
{

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;

  // convert some physical input parameters to internal units 
  Hubble = HUBBLE * UnitTime_in_s;

  // compute a few quantitites 
  RhoCrit = 3 * Hubble * Hubble / (8 * M_PI * G);

}



void read_snap_list(void)
{
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithSnapList);

  if(!(fd = fopen(fname, "r")))
  {
    printf("can't read output list in file '%s'\n", fname);
    ABORT(0);
  }

  Snaplistlen = 0;
  do
  {
    if(fscanf(fd, " %lg ", &AA[Snaplistlen]) == 1)
      Snaplistlen++;
    else
      break;
  }
  while(Snaplistlen < MAXSNAPS);

  fclose(fd);

#ifdef MPI
  if(ThisTask == 0)
#endif
    printf("found %d defined times in snaplist\n", Snaplistlen);
}



double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time, result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;

  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / Hubble,
    1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / Hubble * result;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift 
  return time;
}

double integrand_time_to_present(double a, void *param)
{
  return 1 / sqrt(Omega / a + (1 - Omega - OmegaLambda) + OmegaLambda * a * a);
}



