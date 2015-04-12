#include "hilpos.hh"
#include <cstdlib>
#include <string>

using std::string;

int main(const int argc, const char** argv) {
  srand(time(NULL));
  const string dbpath = "hyg.sqlite"; /// \todo make parameter
  solve_tsp_cuda(dbpath);
  return 0;
}
