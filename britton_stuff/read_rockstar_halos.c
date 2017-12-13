#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#define ROCKSTAR_MAGIC (uint64_t)0xfadedacec0c0d0d0
#define BINARY_HEADER_SIZE 256
#define PARTICLE_TYPE_NONE 0
#define PARTICLE_TYPE_IDS 1
#define PARTICLE_TYPE_FULL 2

#define OUTPUT_BUFFER_SIZE 1000000
#define VERSION_MAX_SIZE 12

struct binary_output_header {
  uint64_t magic;
  int64_t snap, chunk;
  float scale, Om, Ol, h0;
  float bounds[6];
  int64_t num_halos, num_particles;
  float box_size, particle_mass;
  int64_t particle_type;
  int32_t format_revision;
  char rockstar_version[VERSION_MAX_SIZE];
  char unused[BINARY_HEADER_SIZE - (sizeof(char)*VERSION_MAX_SIZE) - (sizeof(float)*12) - sizeof(int32_t) - (sizeof(int64_t)*6)];
};

struct halo {
  int64_t id;
  float pos[6], corevel[3], bulkvel[3];
  float m, r, child_r, vmax_r, mgrav, vmax, rvmax, rs, klypin_rs, vrms,
    J[3], energy, spin, alt_m[4], Xoff, Voff, b_to_a, c_to_a, A[3],
    b_to_a2, c_to_a2, A2[3],
    bullock_spin, kin_to_pot, m_pe_b, m_pe_d;
  int64_t num_p, num_child_particles, p_start, desc, flags, n_core;
  float min_pos_err, min_vel_err, min_bulkvel_err;
};


void binary_halos(char *readin,
                          struct binary_output_header *bheader,
                          struct halo **halos,
                          int64_t **part_ids,
                          int64_t n_halos_max,
                          int64_t n_particles_max){
  int64_t i, j;
  int32_t nitems;

  FILE *input;

  input=fopen(readin,"rb");

  nitems = fread(bheader, sizeof(struct binary_output_header), 1, input);
  assert(bheader->magic == ROCKSTAR_MAGIC);
  assert(bheader->num_halos >= 0);
  assert(bheader->num_particles >= 0);
  assert(bheader->num_halos<=n_halos_max);
  assert(bheader->num_particles<=n_particles_max);
  i = fread(*halos, sizeof(struct halo), bheader->num_halos, input);
  j = fread(*part_ids, sizeof(int64_t), bheader->num_particles, input);
  if (i!=bheader->num_halos || j!=bheader->num_particles) {
    fprintf(stderr, "[Error] Truncated input file %s!\n", readin);
    exit(1);
  }
  fclose(input);
}

int32_t get_header_counts(char *filename, struct binary_output_header *bhead, int64_t *N_halos, int64_t *N_part)
{

  FILE *infile;
  int32_t nitems;

  infile = fopen(filename, "r");
  if (infile == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", filename);
    return EXIT_FAILURE;
  }

  nitems = fread(bhead, sizeof(struct binary_output_header), 1, infile);
  *N_halos = bhead->num_halos;
  *N_part = bhead->num_particles;
 
  fclose(infile); 
  return EXIT_SUCCESS;

} 

int main(int argc, char **argv)
{

  int64_t N_halos, N_part;
  int32_t status;
  struct binary_output_header *bhead;
  bhead = malloc(sizeof(struct binary_output_header));
  char filename[1024];
  
  snprintf(filename, 1024, "%s", "/lustre/projects/p004_swin/jseiler/Rockstar_output/Britton_Halos_final_noLL/halos_20.0.bin"); 

  status = get_header_counts(filename, bhead, &N_halos, &N_part);
  printf("This file contains %ld halos consisting of %ld particles\n", N_halos, N_part);
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}
