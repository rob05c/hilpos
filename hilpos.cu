#define CUB_STDERR

#include "hilpos.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>    
#include <vector>
#include <utility>

#include <linux/cuda.h>
#include <cub/cub.cuh>
#include <cub/util_allocator.cuh>
#include <cub/device/device_radix_sort.cuh>

using std::vector;
using std::pair;

using namespace cub; // debug

/// \todo fix to not be global
CachingDeviceAllocator g_allocator(true); // CUB caching allocator for device memory

/// \returns the device totalGlobalMem
inline size_t GetDeviceMemory() {
  cudaDeviceProp properties;
  int deviceNum;
  CubDebugExit(cudaGetDevice(&deviceNum));
  CubDebugExit(cudaGetDeviceProperties(&properties, deviceNum));
  return properties.totalGlobalMem;
}

#ifndef MIN
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#endif

#define GET_MID(MIN, MAX) (MIN + (MAX - MIN) / 2)

#define NEW_DIMENSION_MIN(VAL, MID, MIN, MAX) (MID * (VAL >= MID) + MIN * (VAL < MID))
#define NEW_DIMENSION_MAX(VAL, MID, MIN, MAX) (MID * (VAL < MID) + MAX * (VAL >= MID))

__global__ void k_create_hilbert_codes(cuda_star* stars, uint64_t* codes, size_t len,
                                          double xmin, double ymin, double zmin,
                                          double xmax, double ymax, double zmax) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i >= len)
    return; // skip the final block remainder

  const cuda_star& star = stars[i];
  const double x = star.x;
  const double y = star.y;
  const double z = star.z;

  uint64_t code = 0;
  #pragma unroll
  for(int i = 0, end = 21; i != end; ++i) {
    code = code << 3;

    const double xmid = GET_MID(xmin, xmax);
    const double ymid = GET_MID(ymin, ymax);
    const double zmid = GET_MID(zmin, zmax);

    const char bit0 = (y <= ymid);
    const char bit1 = (y <= ymid && z >= zmid) 
      || (y >= ymid && z <= zmid);
    const char bit2 = (x >= xmid && y <= ymid && z <= zmid)
      || (x <= xmid && y <= ymid && z >= zmid)
      || (x >= xmid && y >= ymid && z >= zmid)
      || (x <= xmid && y >= ymid && z <= zmid);
    const char bits = (bit0 << 2) | (bit1 << 1) | (bit2 << 0);
    code = code | bits;

    xmin = NEW_DIMENSION_MIN(x, xmid, xmin, xmax);
    xmax = NEW_DIMENSION_MAX(x, xmid, xmin, xmax);
    ymin = NEW_DIMENSION_MIN(y, ymid, ymin, ymax);
    ymax = NEW_DIMENSION_MAX(y, ymid, ymin, ymax);
    zmin = NEW_DIMENSION_MIN(z, zmid, zmin, zmax);
    zmax = NEW_DIMENSION_MAX(z, zmid, zmin, zmax);
  }
  codes[i] = code;
}

/// \return array of calculated hilbert codes of length len. Caller takes ownership.
uint64_t* create_hilbert_codes(cuda_star* stars, size_t len,
                               double xmin, double ymin, double zmin,
                               double xmax, double ymax, double zmax) {
  const unsigned int THREADS_PER_BLOCK = 512;
  cuda_star* cuda_stars;
  uint64_t*  cuda_codes;

  cudaMalloc((void**)&cuda_stars, len * sizeof(cuda_star));
  cudaMalloc((void**)&cuda_codes, len * sizeof(uint64_t));
  cudaMemcpy(cuda_stars, stars, len * sizeof(cuda_star), cudaMemcpyHostToDevice);
  k_create_hilbert_codes<<<(len + (THREADS_PER_BLOCK - 1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK>>>(cuda_stars, cuda_codes, len,
                                                                                                        xmin, ymin, zmin,
                                                                                                        xmax, ymax, zmax);
  uint64_t* codes = new uint64_t[len];
  cudaMemcpy(codes, cuda_codes, len * sizeof(uint64_t), cudaMemcpyDeviceToHost);
  cudaFree(cuda_stars);
  cudaFree(cuda_codes);
  return codes;
}

/// note: sorts the stars, according to the codes. Does NOT sort the codes, for efficiency. Easily could.
void cuda_sort(cuda_star* stars, uint64_t* codes, size_t len) {
  DoubleBuffer<uint64_t> d_keys;
  DoubleBuffer<cuda_star> d_values;
  CubDebugExit( g_allocator.DeviceAllocate((void**)&d_keys.d_buffers[0], sizeof(uint64_t) * len));
  CubDebugExit( g_allocator.DeviceAllocate((void**)&d_keys.d_buffers[1], sizeof(uint64_t) * len));
  CubDebugExit( g_allocator.DeviceAllocate((void**)&d_values.d_buffers[0], sizeof(cuda_star) * len));
  CubDebugExit( g_allocator.DeviceAllocate((void**)&d_values.d_buffers[1], sizeof(cuda_star) * len));

  CubDebugExit( cudaMemcpy(d_keys.d_buffers[0], codes, sizeof(uint64_t) * len, cudaMemcpyHostToDevice));
  CubDebugExit( cudaMemcpy(d_values.d_buffers[0], stars, sizeof(cuda_star) * len, cudaMemcpyHostToDevice));


  size_t temp_storage_bytes = 0;
  void* d_temp_storage = NULL;
  CubDebugExit( DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, len));
  CubDebugExit( g_allocator.DeviceAllocate(&d_temp_storage, temp_storage_bytes));

  CubDebugExit( DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys, d_values, len));

//  CubDebugExit( cudaMemcpy(codes, d_keys.Current(), len * sizeof(uint64_t), cudaMemcpyDeviceToHost));
  CubDebugExit( cudaMemcpy(stars, d_values.Current(), len * sizeof(cuda_star), cudaMemcpyDeviceToHost));

  CubDebugExit( g_allocator.DeviceFree(d_keys.d_buffers[0]));
  CubDebugExit( g_allocator.DeviceFree(d_keys.d_buffers[1]));
  CubDebugExit( g_allocator.DeviceFree(d_values.d_buffers[0]));
  CubDebugExit( g_allocator.DeviceFree(d_values.d_buffers[1]));
  CubDebugExit( g_allocator.DeviceFree(d_temp_storage));
}
