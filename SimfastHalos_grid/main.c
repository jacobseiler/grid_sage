#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#define MAXLEN 1024

#define ABORT(sigterm)                                                  \
do {                                                                \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
  myexit(sigterm);                                                \
} while(0)

int malloc_grid(float **Grid, int GridSize)
{
  int i;
  int status = 0;

  //if (ThisTask_GridNr == 0) // Only need to malloc the grid once.
  //  Grid = mymalloc(sizeof(struct GRID)*CUBE(GridSize));

  *Grid = NULL;
  *Grid = (float *)malloc(sizeof(float)*GridSize*GridSize*GridSize);

  if(*Grid == NULL)
  {
    status = 1; 
  }

  return status; 
  
}

void init_grid(float **Grid, int GridSize)
{
  *Grid[0] = 1.000;

}

int main(int argc, char **argv)
{

  char inputdir[MAXLEN];
  snprintf(inputdir, MAXLEN, "/lustre/projects/p004_swin/jseiler/Simfast21/Halos");

  int full_ngrid = 512; // Number of cells (along an axis) of the full resolution grid.
  int smooth_ngrid = 128; // Number of cells (along an axis) for the smoothed grid.

  char halofile[MAXLEN]; 
  double z;
  double z_high = 15.000;
  double z_low = 5.000;
  double dz = 0.250;

  float *Nion;
  int count = 0;
  

  for(z = z_high; z > z_low - dz; z -= dz)
  {
    if(count == 0)
    {
      if(malloc_grid(&Nion, smooth_ngrid) == 1)
      {
        fprintf(stderr, "Malloc of Nion grid failed.\n");
	exit(0);
      }  
    }
    
    init_grid(&Nion, smooth_ngrid);
    fprintf(stderr, "Grid[0] = %.2f\n", Nion[0]);
    ++count;   
 
  }
  return 0;

}
