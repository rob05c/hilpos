#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <iostream>
//#include <chrono>
#include <vector>
#include <utility>
#include <algorithm>
#include <numeric>
#include <iomanip>
#include <cstdio>
//#include "tbb/tbb.h"
#include "sqlite3.h"

#include "hilpos.h"

using std::vector;
using std::pair;
using std::make_pair;
using std::cout;
using std::endl;
using std::string;

vector<star> get_stars(const string& path) {
  vector<star> stars;
  sqlite3 *db;
  char *zErrMsg = 0;
  int rc = sqlite3_open(path.c_str(), &db);
  if(rc){
    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
    sqlite3_close(db);
    return stars;
  }
  const string query = "select starid, x, y, z from hygxyz";
  sqlite3_stmt* res;
  rc = sqlite3_prepare_v2(db, query.c_str(), -1, &res, 0);
  if(rc != SQLITE_OK){
    fprintf(stderr, "SQL error: %s\n", zErrMsg);
    sqlite3_free(zErrMsg);
  }


  while(sqlite3_step(res) == SQLITE_ROW) {
    star s;
    s.id = sqlite3_column_int(res, 0);
    s.x  = sqlite3_column_double(res, 1);
    s.y  = sqlite3_column_double(res, 2);
    s.z  = sqlite3_column_double(res, 3);
    s.curve_position = 0;
    stars.push_back(s);
  }
  sqlite3_close(db);
  return stars;
}

/// This should be a lambda, but we're not using C++11 so we can compile with Hydra's old NVCC.
bool is_outside_bounds(const star& s) {
   return s.x < -1000.0 || s.x > 1000.0 
     || s.y < -1000.0 || s.y > 1000.0
     || s.z < -1000.0 || s.z > 1000.0;
}

static inline star denormalize(star star) {
  star.x = star.x * 2000.0 - 1000.0;
  star.y = star.y * 2000.0 - 1000.0;
  star.z = star.z * 2000.0 - 1000.0;
  return star;
}

double denormalize_distance(const double distance) {
  return distance * 2000.0;
}

static const double NORMALIZE_MULTIPLIER = 2000;
/// Removes stars which have coordinates greater than 1000 or less than -1000
/// and normalizes remaining stars to be between 0 and 1.
/// \param[in] stars
vector<star> normalize(vector<star>& stars) {
  stars.erase(std::remove_if(stars.begin(), stars.end(), is_outside_bounds), stars.end());
  for(vector<star>::iterator i = stars.begin(), end = stars.end(); i != end; ++i) {
    star& s = *i;
    s.x = (s.x + 1000.0) / 2000.0;
    s.y = (s.y + 1000.0) / 2000.0;
    s.z = (s.z + 1000.0) / 2000.0;
  }
  return stars;
}

static inline double get_mid(double min, double max) {
  return min + (max - min) / 2;
}

/// \param[inout] min
/// \param[inout] max
static inline void set_new_dimension(const double val, const double mid, double& min, double& max) {
    if(val < mid) {
      max = mid;
    } else {
      min = mid;
    }
}

/// Returns the hilbert curve position for the given 3D point 
///   to the 21st degree (floort(64/3) = 21).
/// Assumes the space is normalized to 0..1
uint64_t get_hilbert_position(double x, double y, double z, 
                             double xmin, double ymin, double zmin,
                             double xmax, double ymax, double zmax) {
  const int THREE_INTO_SIXTY_FOUR = 21;
  uint64_t code = 0;
  for(int i = 0, end = THREE_INTO_SIXTY_FOUR; i != end; ++i) {
    code = code << 3;

    const double xmid = get_mid(xmin, xmax);
    const double ymid = get_mid(ymin, ymax);
    const double zmid = get_mid(zmin, zmax);

    if(x <= xmid && y <= ymid && z <= zmid)
      code = code | 0;
    else if(x >= xmid && y <= ymid && z <= zmid)
      code = code | 1;
    else if(x >= xmid && y <= ymid && z >= zmid)
      code = code | 1;
    else if(x >= xmid && y <= ymid && z >= zmid)
      code = code | 2;
    else if(x <= xmid && y <= ymid && z >= zmid)
      code = code | 3;
    else if(x <= xmid && y >= ymid && z >= zmid)
      code = code | 4;
    else if(x >= xmid && y >= ymid && z >= zmid)
      code = code | 5;
    else if(x >= xmid && y >= ymid && z <= zmid)
      code = code | 6;
    else if(x <= xmid && y >= ymid && z <= zmid)
      code = code | 7;

    set_new_dimension(x, xmid, xmin, xmax);
    set_new_dimension(y, ymid, ymin, ymax);
    set_new_dimension(z, zmid, zmin, zmax);
  }
  return code;
}

double distance(const star& a, const star& b) {
  return hypot(hypot(a.x - b.x, a.y - b.y), a.z - b.z);

//  return fabs(a.x - b.x) + fabs(a.y - b.y) + fabs(a.z - b.z);
}

/// \todo make this a stream operator?
void print_curveposition(const uint64_t pos) {
  for(int i = 0, end = sizeof(pos) * CHAR_BIT; i != end; ++i) {
    cout << ((pos & (1 << i)) ? '1' : '0');
  }
}

void print_star_header() {
  printf("%6s, %20s, %20s, %20s, %20s\n", "starid", "x", "y", "z", "distance to previous");
}

/// \todo make this a stream operator?
void print_star(const star& s, const double distance) {
  printf("%6d, %20.15f, %20.15f, %20.15f, %20.15f\n", s.id, s.x, s.y, s.z, distance);
//  cout << std::fixed << std::setprecision(15) << "{" << s.id << "," << s.x << "," << s.y << "," << s.z << ",";
//  cout << distance;
//  print_curveposition(s.curve_position);
//  cout << "}" << endl;
}

bool star_comparator(const star& a, const star& b) {
  return a.curve_position < b.curve_position;
}

int main(const int argc, const char** argv) {
  srand(time(NULL));
  const string db = "hyg.sqlite"; /// \todo make parameter
  vector<star> stars = get_stars(db);
  cout << "num stars: " << stars.size() << endl;
  stars = normalize(stars);
  cout << "num normalised: " << stars.size() << endl;


  const double xmin = 0.0;
  const double ymin = 0.0;
  const double zmin = 0.0;
  const double xmax = 1.0;
  const double ymax = 1.0;
  const double zmax = 1.0;
  /// \todo parallelize

  for(vector<star>::iterator i = stars.begin(), end = stars.end(); i != end; ++i) {
    star& s = *i;
    s.curve_position = get_hilbert_position(s.x, s.y, s.z, 
                                               xmin, ymin, zmin, 
                                               xmax, ymax, zmax);
  }

  std::random_shuffle(stars.begin(), stars.end()); // randomize stars, in case database was sorted;
  double random_distance_sum;
  double random_distance_average;
  {
    vector<double> distances;
    for(vector<star>::iterator i = stars.begin() + 1, end = stars.end(); i != end; ++i) {
      distances.push_back(distance(*i, *(i - 1)));
    }
    random_distance_sum = std::accumulate(distances.begin(), distances.end(), 0.0);
    random_distance_average = random_distance_sum / (double)distances.size();
  }

  std::sort(stars.begin(), stars.end(), star_comparator);

  print_star_header();

  if(!stars.empty())
    print_star(*stars.begin(), 0.0);
  vector<double> distances;
  for(vector<star>::iterator i = stars.begin() + 1, end = stars.end(); i != end; ++i) {
    distances.push_back(distance(*i, *(i - 1)));
    print_star(denormalize(*i), denormalize_distance(distance(*i, *(i - 1))));
  }

  const double distance_sum = std::accumulate(distances.begin(), distances.end(), 0.0);
  printf("path length: %.15f\n", denormalize_distance(distance_sum));
  const double distance_average = distance_sum / (double)distances.size();
  printf("average distance: %.15f\n", denormalize_distance(distance_average));
  
  printf("random path length: %.15f\n", denormalize_distance(random_distance_sum));
  printf("random average distance: %.15f\n", denormalize_distance(random_distance_average));

  return 0;
}
