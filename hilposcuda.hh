#ifndef hilposcudahh
#define hilposcudahh

#include <stdint.h>

struct cuda_star {
  int    id;
  double x;
  double y;
  double z;
};

uint64_t* create_hilbert_codes(cuda_star* stars, size_t len,
                               double xmin, double ymin, double zmin,
                               double xmax, double ymax, double zmax);
void      cuda_sort(cuda_star* stars, uint64_t* codes, size_t len);

#endif
