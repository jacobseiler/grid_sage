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

#ifdef NDEBUG
#define XASSERT(EXP, ...)                                do{} while(0)
#else
#define XASSERT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            fflush(stdout);                                             \
            exit(EXIT_FAILURE);                                         \
        } \
    } while (0)
#endif

#ifdef NDEBUG
#define XPRINT(EXP, ...)                                do{} while(0)
#else
#define XPRINT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Warning in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            fflush(stdout);                                             \
        } \
    } while (0)
#endif

// Function taken from https://phoxis.org/2013/05/04/generating-random-numbers-from-normal-distribution-in-c/ 
// Generates a random number from a normal distribution with mean mu and standard deviation sigma. 

double randn (double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow(U1, 2) + pow(U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


int malloc_grid(float **Grid, int GridSize)
{
  int status = 0;

  fprintf(stderr, "Mallocing the grid.\n");

  *Grid = NULL;
  *Grid = (float *)malloc(sizeof(float)*GridSize*GridSize*GridSize);
  fprintf(stderr, "Malloced.\n");
  if(*Grid == NULL)
  {
    status = 1; 
  }

  return status; 
  
}

void init_grid(float *Grid, int GridSize)
{
  fprintf(stderr, "Setting grid to 0\n");
  long int i;
  for (i = 0; i < GridSize*GridSize*GridSize; ++i)
  {
  
    Grid[i] = 0.000;
  }
  fprintf(stderr, "Initalized\n");
}

int main(int argc, char **argv)
{

  FILE *fd;

  char inputdir[MAXLEN];
  snprintf(inputdir, MAXLEN, "/lustre/projects/p004_swin/jseiler/Simfast21/Halos");
  char halo_name[MAXLEN];
  char outfile[MAXLEN];
  long int Nhalos;
  
  typedef struct Halo_t_
  {
    float Mass;
    float x,y,z;
  } Halo_t;

  Halo_t *halo_v;

  char mvir_fit_basedir[MAXLEN];
  snprintf(mvir_fit_basedir, MAXLEN, "/lustre/projects/p004_swin/jseiler/tiamat/");
  char line[MAXLEN];

  int full_ngrid = 512; // Number of cells (along an axis) of the full resolution grid.
  int smooth_ngrid = 512; // Number of cells (along an axis) for the smoothed grid.
  float conversion = (float) smooth_ngrid/ (float) full_ngrid;

  char halofile[MAXLEN]; 
  double z;
  double z_high = 15.000;
  double z_low = 6.000;
  double dz = 0.250;
  int Nsnaps = 37;
  long int idx; 

  char mvir_fit_file[MAXLEN];
  double mvir_fit_z[37] = {15.086, 14.658, 14.658, 14.258, 13.882, 13.882, 13.529, 13.195, 12.881, 12.881, 12.300, 12.300, 11.777, 11.777, 11.534, 11.302, 11.081, 10.869, 10.473, 10.286, 10.107, 9.770, 9.457, 9.309, 9.027, 8.765, 8.518, 8.287, 8.069, 7.765, 7.485, 7.227, 7.065, 6.764, 6.556, 6.238, 6.007}; 

  int Nbins = 60;
  double mvir_bin[60];
  double mvir_ngamma_mean[60];
  double mean;
  double mvir_ngamma_std[60];
  double std;
  int i, j;
  for (i = 0; i < Nbins; ++i)
  {
    mvir_bin[i] = 6.05 + .1 * i;
  //  fprintf(stderr, "%.2f\n", mvir_bin[i]);
  }

  float *Nion;
  int count = 0, ngamma_counter;
  double Nion_tmp; 

  for(z = z_high; z > z_low - dz; z -= dz)
  {
    if(count == 0) // Only need to allocate memory once.
    {

      if(malloc_grid(&Nion, smooth_ngrid) == 1)
      {
        fprintf(stderr, "Malloc of Nion grid failed.\n");
	exit(0);
      }  

    }
    
    init_grid(Nion, smooth_ngrid); // Sets the Nion grid to all 0s. 
 
    // Reading in the Ngamma Means // 

    snprintf(mvir_fit_file, MAXLEN, "%s/mean_mvir_ngammafesc_z%.3f.dat", mvir_fit_basedir, mvir_fit_z[count]);
    if(!(fd = fopen(mvir_fit_file, "r")))
    {
      fprintf(stderr, "Cannot find file %s\n", mvir_fit_file);
      fprintf(stderr, "count = %d\n", count); 
      exit(0);
    }

    ngamma_counter = 0;
    while (fgets(line, MAXLEN, fd))
    {
      if(*line == '#')
	continue; 
      sscanf(line, "%lf", &mean);
      mvir_ngamma_mean[ngamma_counter] = mean; 
      ngamma_counter++;
    } 
    fclose(fd);
    fprintf(stderr, "Read in the means\n"); 
   
    // Reading in the Ngamma stds // 
    
    snprintf(mvir_fit_file, MAXLEN, "%s/std_mvir_ngammafesc_z%.3f.dat", mvir_fit_basedir, mvir_fit_z[count]);
    if(!(fd = fopen(mvir_fit_file, "r")))
    {
      fprintf(stderr, "Cannot find file %s\n", mvir_fit_file);
      fprintf(stderr, "count = %d\n", count); 
      exit(0);
    }

    ngamma_counter = 0;
    while (fgets(line, MAXLEN, fd))
    {
      if(*line == '#')
	continue; 
      sscanf(line, "%lf", &std);
      mvir_ngamma_std[ngamma_counter] = std;
      ngamma_counter++;
    } 
    fclose(fd);
    fprintf(stderr, "Read in the std\n"); 
      
    // At this point I have:
    // The mvir mass bins (mvir_bin)
    // The mvir-Ngamma*fesc relationship mean (mvir_ngamma_mean) and stds (mvir_ngamma_std).
    // Now wish to open up the halos, place the halo in the correct grid cell, get its Mvir and then assign it an emissivity based upon where it falls in the mvir bin.


    /* read halo catalog */
    snprintf(halo_name, MAXLEN, "%s/halonl_z%.3f_N512_L100.0.dat.catalog", inputdir,  z);
    if((fd=fopen(halo_name,"rb"))==NULL){  
      printf("Halo file: %s does not exist.\n", halo_name);
      exit(0);
    } 
    fread(&Nhalos,sizeof(long int),1,fd);
    fprintf(stderr, "Reading %ld halos from file %s\n", Nhalos, halo_name);
    if(!(halo_v=(Halo_t *) malloc(Nhalos*sizeof(Halo_t)))) { 
      printf("Memory problem - halos...\n");
      exit(0);
    }
    fread(halo_v,sizeof(Halo_t),Nhalos,fd);  
    fclose(fd);


    for (i = 0; i < Nhalos; ++i)
    {
      
      int grid_x = halo_v[i].x*conversion;
      if(grid_x > smooth_ngrid - 1)
	grid_x = smooth_ngrid - 1;
      if(grid_x < 0)
	grid_x = 0;

      int grid_y = halo_v[i].y*conversion;
      if(grid_y > smooth_ngrid - 1)
	grid_y = smooth_ngrid - 1;
      if(grid_y < 0)
	grid_y = 0;

      int grid_z = halo_v[i].z*conversion;
      if(grid_z > smooth_ngrid - 1)
	grid_z = smooth_ngrid- 1;
      if(grid_z < 0)
	grid_z = 0;

      idx=(long int)(grid_x*smooth_ngrid*smooth_ngrid+grid_y*smooth_ngrid+grid_z);

      XPRINT(idx < smooth_ngrid*smooth_ngrid*smooth_ngrid && idx >= 0, "For redshift %.3f the grid cell index was %d for halo %d. The halo had x = %.2f = %.2f, y = %.2f = %.2f, z = %.2f = %.2f\n", z, idx, i, halo_v[i].x, grid_x, halo_v[i].y, grid_y, halo_v[i].z, grid_z);
      if (idx < 0)
        idx = 0;
      if (idx > smooth_ngrid*smooth_ngrid*smooth_ngrid)
       idx = smooth_ngrid*smooth_ngrid*smooth_ngrid; 


      for(j = 0; j < 59; ++j)
      {
        if(log10(halo_v[i].Mass) > mvir_bin[j] && log10(halo_v[i].Mass) < mvir_bin[j+1])
        {

	  mean = log10(mvir_ngamma_mean[j]);
          std = 0.484 * (mvir_ngamma_std[j]/mvir_ngamma_mean[j]);


	  if(mean > 0.0) // If our fit had no Halos in this mass bin then the input would 'Nan'. 
          {
            Nion_tmp = randn(log10(mvir_ngamma_mean[j]), 0.484 * (mvir_ngamma_std[j]/mvir_ngamma_mean[j]));
          } 
          else
          { 
            Nion_tmp = 0;
          }

	  break;	
        } 
      }

      if (Nion_tmp > 0 && Nion_tmp < 70.0) 
      {
        Nion[idx] += pow(10, Nion_tmp)/1.0e50;
      }
      Nion_tmp = 0;
    }

    snprintf(outfile, MAXLEN, "/lustre/projects/p004_swin/jseiler/Simfast21/grid_halos/z%.3f", z); 
    if((fd=fopen(outfile,"w"))==NULL){  
      printf("Output file: %s does not exist.\n", outfile);
      exit(0);
    } 

    for(i = 0; i < smooth_ngrid*smooth_ngrid*smooth_ngrid; ++i)
    {
//	fprintf(stderr, "At write %.4e\n", Nion[i]);
      fwrite(&Nion[i], sizeof(float), 1, fd);
    }
    fprintf(stderr, "Saved to file %s\n", outfile);
    fclose(fd); 
    ++count;   
 
  }
  return 0;

}

