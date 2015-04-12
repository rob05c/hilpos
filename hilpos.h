#ifndef hilposh
#define hilposh
#include <stdint.h>
#include <vector>
#include <string>

struct star {
  int    id;
  double x;
  double y;
  double z;
  uint64_t curve_position;
};

struct cuda_star {
  int    id;
  double x;
  double y;
  double z;
};

bool              is_outside_bounds(const star& s);

const double      NORMALIZE_MULTIPLIER = 2000;
double            denormalize_distance(const double distance);
std::vector<star> normalize(std::vector<star>& stars);

double            get_mid(double min, double max);
void              set_new_dimension(const double val, const double mid, double& min, double& max);
uint64_t          get_hilbert_position(double x, double y, double z, 
                                       double xmin, double ymin, double zmin,
                                       double xmax, double ymax, double zmax);
void              print_curveposition(const uint64_t pos);
void              print_star_header();
bool              star_comparator(const star& a, const star& b);

std::vector<star> get_stars(const std::string& path);
void              solve_tsp(const std::string& path);
void              solve_tsp_cuda(const std::string& path);

uint64_t* create_hilbert_codes(cuda_star* stars, size_t len,
                               double xmin, double ymin, double zmin,
                               double xmax, double ymax, double zmax);
void      cuda_sort(cuda_star* stars, uint64_t* codes, size_t len);

#endif
