#include "core_allvars_grid.h"

void myexit(int signum);

int32_t init(void);
int32_t init_grid(struct GRID_STRUCT *grid);
int32_t free_grid(void);
void set_units(void);
void read_snap_list(void);
void determine_fesc_constants(void);

void read_parameter_file(char *fname);

int32_t load_gals(char *fname);
void free_gals(void);

double time_to_present(double z);
double integrand_time_to_present(double a, void *param);
void estimate_grid_memory(void);
void estimate_gal_memory(int NtotGals);

int32_t update_grid_properties(int32_t filenr);
int32_t update_quasar_tracking(int64_t gal_idx, int32_t snapshot_idx);
void count_grid_properties(struct GRID_STRUCT *count_grid);
void normalize_photon(int GridNr); 
void normalize_slope_photons(int GridNr);
void calculate_photons(float SFR, float Z, float *Ngamma_HI, float *Ngamma_HeI, float *Ngamma_HeII);
int32_t calculate_fesc(int p, int i, int filenr, float *fesc_local);
double get_metallicity(double gas, double metals);

const char* getfield(char* line, int num);

int32_t save_fesc_properties(int filenr, int32_t merged_gal_flag); 
int32_t save_grid(struct GRID_STRUCT *save_grid);
void save_redshift(void);

void calculate_halomass(void);

void *mymalloc(size_t n);
void myfree(void *p);

#ifdef MPI
struct GRID_STRUCT *MPI_sum_grids(void);
#endif
