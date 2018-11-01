#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>


int32_t read_IDs(int64_t *particle_IDs)
{

  FILE *infile;
  char fname[1024];
  uint64_t file_size;

  snprintf(fname, 1024, "%s", "/lustre/projects/p004_swin/jseiler/Rockstar_output/Britton_fsub_test/particles_ID_50.dat");
  infile = fopen(fname, "r");

  if (infile == NULL)
  {
    fprintf(stderr, "Could not open %s\n", fname);
    return EXIT_FAILURE;
  }
 
  fseek(infile, 0L, SEEK_END); // Move to the end of the file 
  file_size = ftell(infile); // Then count how many bytes we have in it.
  fseek(infile, 0, SEEK_SET);
 
  printf("The file has %ld bytes\n", file_size);


  fclose(infile);

  return EXIT_SUCCESS;

}

int main(int argc, char **argv)
{

  int64_t *particle_IDs;
  int32_t status;

  status = read_IDs(particle_IDs);
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }
  

  return 0;

}
