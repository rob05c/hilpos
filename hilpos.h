#ifndef hilposh
#define hilposh
#include <stdint.h>
#include <vector>
#include <string>

struct star {
  int id;
  double x;
  double y;
  double z;
  uint64_t curve_position;
};


bool              is_outside_bounds(const star& s);

const double      NORMALIZE_MULTIPLIER = 2000;
star              denormalize(star star);
double            denormalize_distance(const double distance);
std::vector<star> normalize(std::vector<star>& stars);

double            get_mid(double min, double max);
void              set_new_dimension(const double val, const double mid, double& min, double& max);
uint64_t          get_hilbert_position(double x, double y, double z, 
                                       double xmin, double ymin, double zmin,
                                       double xmax, double ymax, double zmax);
double            distance(const star& a, const star& b);

void              print_curveposition(const uint64_t pos);
void              print_star_header();
void              print_star(const star& s, const double distance);
bool              star_comparator(const star& a, const star& b);

std::vector<star> get_stars(const std::string& path);
void              solve_tsp(const std::string& path);

#endif
