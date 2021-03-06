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

#define XASSERT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            fflush(stdout);                                             \
            exit(EXIT_FAILURE);                                         \
        } \
    } while (0)

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

int cmpfunc (const void * a, const void * b) {
   return ( *(int64_t*)a - *(int64_t*)b );
}

int32_t load_binary_halos(char *readin,
                          struct binary_output_header *bheader,
                          struct halo **halos,
                          int64_t **part_ids,
                          int64_t n_halos_max,
                          int64_t n_particles_max){
  int64_t i, j;
  int32_t nitems;

  FILE *input;

  input=fopen(readin,"rb");
  if (input == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", readin);
    return(EXIT_FAILURE);
  }

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

  return EXIT_SUCCESS;

}

int32_t get_header_counts(char *filename, struct binary_output_header *bhead, int64_t *N_halos, int64_t *N_part)
{

  FILE *infile;
  int32_t nitems;

  infile = fopen(filename, "rb");
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

int32_t sort_particles(int64_t *part_id, int64_t N_part)
{

  printf("Performing a bubble sort.\n");
  int64_t i, j, tmp;

  for (i = 0; i < N_part; ++i)
  {
    if (i % 10000 == 0)
    {
      printf("%.4f%% done\n", (float) i / N_part * 100.0);
    } 
    for (j = i+1; j < N_part; ++j)
    {

      if (part_id[i] > part_id[j])
      {

        tmp = part_id[i];
        part_id[i] = part_id[j];
        part_id[j] = tmp;

      }
    }
  }

  return EXIT_SUCCESS;

}

int64_t check_duplicate_particles(int64_t *part_id, int64_t N_part)
{

  int64_t part_idx, compare_idx;
  int64_t duplicate_count = 0;
  for (part_idx = 0; part_idx < N_part - 1; ++part_idx)
  {
    if (part_id[part_idx] == part_id[part_idx + 1])
    {
     ++duplicate_count;
    }
    
  }
  printf("There were %ld duplicates\n", duplicate_count);
  return duplicate_count; 

}

int32_t slice_particles_into_global(int64_t *part_id_local, int64_t *part_id_global, int64_t offset, int64_t N_part_local)
{

  int64_t part_idx;

  for (part_idx = 0; part_idx < N_part_local; ++part_idx)
  {
    part_id_global[offset+part_idx] = part_id_local[part_idx];
  }
   
  return EXIT_SUCCESS; 

} 

int32_t create_unique_list(int64_t *part_id_global, int64_t N_part_global, int64_t *part_id_unique, int64_t duplicate_count)
{

  int64_t part_idx, unique_count = 1;

  part_id_unique[0] = part_id_global[0];


  for (part_idx = 0; part_idx < N_part_global; ++part_idx)
  {
    if (part_id_global[part_idx] != part_id_global[part_idx -1])
    {
      part_id_unique[unique_count] = part_id_global[part_idx];
      ++unique_count;
    } 
  }

  XASSERT(unique_count - 1 == N_part_global - duplicate_count, "Unique count = %ld \t N_part_global = %ld \t duplicate_count = %ld\nShould have unique count - 1 == N_part_global - duplicate_count (%ld)\n", unique_count, N_part_global, duplicate_count, N_part_global - duplicate_count);

  return EXIT_SUCCESS; 

}

int32_t save_unique_list(int64_t *part_id_unique, int64_t N_part_unique, char *fbase, int32_t N_files)
{

  FILE *outfile;
  char buf[1024];
  int32_t file_idx;

  int64_t N_ids_per_file = N_part_unique / N_files; // Evenly spread the IDs over each file.
  int64_t N_final_file_extra = N_part_unique - (N_ids_per_file * N_files); // However the last file will need a few more if the number isn't perfectly divisible.
  int64_t N_ids_thisfile;
  int32_t N_ids_highword = N_part_unique >> 32;
  int64_t write_idx, offset;
 
  for (file_idx = 0; file_idx < N_files; ++file_idx)
  {
    if (file_idx == N_files - 1)
    {
      N_ids_thisfile = N_ids_per_file + N_final_file_extra; 
    }
    else
    {
      N_ids_thisfile = N_ids_per_file;
    }
 
    snprintf(buf, 1024, "%s.%d", fbase, file_idx);

    outfile = fopen(buf, "wb");
    if (outfile == NULL)
    {
      fprintf(stderr, "Could not open file %s for writing.\n", buf);
      return EXIT_FAILURE;
    }

    // Header information //
    fwrite(&file_idx, sizeof(int32_t), 1, outfile);
    fwrite(&N_files, sizeof(int32_t), 1, outfile);
    fwrite(&N_ids_thisfile, sizeof(int64_t), 1, outfile);
    fwrite(&N_part_unique, sizeof(int64_t), 1, outfile);
    fwrite(&N_ids_highword, sizeof(int32_t), 1, outfile); 

    for (write_idx = 0; write_idx < N_ids_thisfile; ++write_idx)
    {
      fwrite(&part_id_unique[write_idx + offset], sizeof(int64_t), 1, outfile);
    }

    offset += N_ids_thisfile;

    printf("Successfully wrote %ld IDs to %s\n", N_ids_thisfile, buf);
  }
  
  return EXIT_SUCCESS;

  

} 

int main(int argc, char **argv)
{

  int64_t N_halos_local, N_part_local;
  int64_t N_halos_global = 0, N_part_global = 0, N_part_unique;
  int32_t status;
  char filename[1024], outname[1024];

  struct binary_output_header *bhead;
  struct halo *halos_global;
  struct halo *halos_local;
  int64_t *part_id_global; 
  int64_t *part_id_local; 

  int32_t n_files = 64, i_file;
  int64_t offset = 0;

  int64_t duplicate_count;
  int64_t *part_id_unique;

  int32_t N_files = 64;

  if (argc != 2)
  { 
    fprintf(stderr, "Usage : rockstar_halos <snapshot>\n");
    exit(EXIT_FAILURE); 
  } 

  int32_t snapshot_idx = atoi(argv[1]);

  for (i_file = 0; i_file < n_files; ++i_file)
  {
    bhead = malloc(sizeof(struct binary_output_header));
    if (bhead == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the header\n");
      exit(EXIT_FAILURE);
    } 
   
    snprintf(filename, 1024, "%s_%d.%d.bin", "/lustre/projects/p004_swin/jseiler/Rockstar_output/Britton_Halos_final_noLL/halos", snapshot_idx,i_file); 

//    printf("Reading file %s\n", filename);
    status = get_header_counts(filename, bhead, &N_halos_local, &N_part_local);
    printf("This file contains %ld halos consisting of %ld particles\n", N_halos_local, N_part_local);
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }

    N_halos_global += N_halos_local;
    N_part_global += N_part_local;

    free(bhead);
  }

  printf("In total there are %ld halos and %ld particles for this snapshot\n", N_halos_global, N_part_global); 
  printf("This means we need a total of %.4f MB for the halos and %.4f MB for the particles.\n", (float) sizeof(struct halo) * N_halos_global / 1024 / 1024, (float) sizeof(int64_t) * N_part_global / 1024 / 1024);

  halos_global = malloc(sizeof(struct halo) * N_halos_global);
  if (halos_global == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the global halos\n");
    exit(EXIT_FAILURE);
  } 
  
  part_id_global = malloc(sizeof(int64_t) * N_part_global);   
  if (part_id_global == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the global particle IDs\n");
    exit(EXIT_FAILURE);
  } 

  for (i_file = 0; i_file < n_files; ++i_file)
  {
    bhead = malloc(sizeof(struct binary_output_header));
    if (bhead == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the header\n");
      exit(EXIT_FAILURE);
    } 
   
    snprintf(filename, 1024, "%s_%d.%d.bin", "/lustre/projects/p004_swin/jseiler/Rockstar_output/Britton_Halos_final_noLL/halos", snapshot_idx,i_file); 
    
//    printf("Reading file %s\n", filename);
    status = get_header_counts(filename, bhead, &N_halos_local, &N_part_local);
    printf("This file contains %ld halos consisting of %ld particles\n", N_halos_local, N_part_local);
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }

    halos_local = malloc(sizeof(struct halo) * N_halos_local);
    if (halos_local == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the local halos\n");
      exit(EXIT_FAILURE);
    } 
    
    part_id_local = malloc(sizeof(int64_t) * N_part_local);   
    if (part_id_local == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the local particle IDs\n");
      exit(EXIT_FAILURE);
    } 

    status = load_binary_halos(filename, bhead, &halos_local, &part_id_local, N_halos_local, N_part_local);
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }

    // We now have the halos and particles for this file loaded.
    // Need to slice them into the global array.

    slice_particles_into_global(part_id_local, part_id_global, offset, N_part_local); 
    offset += N_part_local;

    free(part_id_local);
    free(halos_local);
    free(bhead);

  }

  printf("Sorting Particles\n");     
  qsort(part_id_global, N_part_global, sizeof(int64_t), cmpfunc);    
  duplicate_count = check_duplicate_particles(part_id_global, N_part_global);
  N_part_unique = N_part_global - duplicate_count; 
  part_id_unique = malloc(sizeof(int64_t) * N_part_unique);
  if (part_id_unique == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the unique particle IDs\n");
    exit(EXIT_FAILURE);
  } 
 
  create_unique_list(part_id_global, N_part_global, part_id_unique, duplicate_count);

  snprintf(outname, 1024, "/lustre/projects/p134_swin/jseiler/rockstar_particles/snap%03d", snapshot_idx);
  status = save_unique_list(part_id_unique, N_part_unique, outname, N_files);
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }


  free(part_id_unique);
  free(part_id_global);
  free(halos_global);
     
  return EXIT_SUCCESS;
}
