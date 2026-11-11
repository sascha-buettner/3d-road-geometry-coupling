#include <cmath>
#include <unistd.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

#include "tum_road_geometry_coupling_cpp/ribbon_height_map.hpp"
using tam::road_geometry_coupling::RibbonHeightMap;
template <class T>
void put(std::ofstream & f, T x)
{
  f.write(reinterpret_cast<const char *>(&x), sizeof(x));
}
int main()
{
  auto path =
    std::filesystem::temp_directory_path() / ("bounds_" + std::to_string(getpid()) + ".rhm");
  try {
    std::mt19937 rng(5112);
    std::uniform_int_distribution<int> h(-700, 700);
    {
      std::ofstream f(path, std::ios::binary);
      f.write("TAMRHG2\0", 8);
      put(f, uint32_t(65));
      put(f, uint32_t(0));
      for (double v : {.025, .025, 0., .00005}) put(f, v);
      put(f, uint64_t(65 * 53));
      put(f, uint64_t(0));
      for (int i = 0; i < 65; ++i) {
        put(f, int32_t(-20));
        put(f, uint32_t(53));
        put(f, uint64_t(i * 53));
      }
      for (int i = 0; i < 65; ++i)
        for (int j = 0; j < 53; ++j) put(f, int16_t((i == 35 && j == 20) ? -32768 : h(rng)));
    }
    RibbonHeightMap map(path);
    std::uniform_real_distribution<double> station(-5., 6.), lateral(-.35, .55), angle(-3.14, 3.14);
    double err = 0., grad = 0.;
    int valid = 0;
    for (int i = 0; i < 2000; ++i) {
      double s = station(rng), n = lateral(rng), a = angle(rng), ds = std::cos(a), dn = std::sin(a);
      auto fast = map.envelope(s, n, ds, dn, .31, true),
           reference = map.envelope(s, n, ds, dn, .31, false);
      if (fast.valid != reference.valid) throw std::runtime_error("Bound changed coverage");
      if (!fast.valid) continue;
      ++valid;
      err = std::max(err, std::abs(fast.height_m - reference.height_m));
      grad = std::max(
        {grad, std::abs(fast.gradient_s - reference.gradient_s),
         std::abs(fast.gradient_n - reference.gradient_n)});
      if (err > 1e-10 || grad > 1e-8)
        throw std::runtime_error("Conservative bounds changed circle maximum/gradient");
    }
    std::cout << valid
              << " periodic/negative-lateral random rough-road comparisons; max height difference "
              << err << " m, gradient difference " << grad << "; resident storage "
              << map.storage_bytes() << " bytes\n";
  } catch (...) {
    std::filesystem::remove(path);
    throw;
  }
  std::filesystem::remove(path);
}
