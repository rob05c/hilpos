#ifndef hilposhh
#define hilposhh
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

void solve_tsp(const std::string& path);
void solve_tsp_cuda(const std::string& path);

#endif
