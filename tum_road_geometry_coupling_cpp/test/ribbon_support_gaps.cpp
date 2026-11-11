#include "tum_road_geometry_coupling_cpp/ribbon_height_map.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <unistd.h>
using tam::road_geometry_coupling::RibbonHeightMap;
template <class T> void put(std::ofstream &f, T v) {
  f.write(reinterpret_cast<const char *>(&v), sizeof(v));
}
void make_map(const std::filesystem::path &path, int kind) {
  std::ofstream f(path, std::ios::binary);
  f.write("TAMRHG2\0", 8);
  put(f, std::uint32_t(128));
  put(f, std::uint32_t(0));
  for (double v : {.025, .025, -.5, .00005})
    put(f, v);
  put(f, std::uint64_t(128 * 41));
  put(f, std::uint64_t(0));
  for (int i = 0; i < 128; i++) {
    put(f, std::int32_t(0));
    put(f, std::uint32_t(41));
    put(f, std::uint64_t(i * 41));
  }
  for (int i = 0; i < 128; i++)
    for (int j = 0; j < 41; j++) {
      bool missing = kind == 0   ? (i == 64 || i == 58)
                     : kind == 1 ? (i >= 58 && i <= 72)
                                 : (i >= 64);
      put(f, std::int16_t(missing ? -32768 : 0));
    }
}
int main() {
  auto root = std::filesystem::temp_directory_path() /
              ("ribbon_support_" + std::to_string(getpid()));
  std::filesystem::create_directory(root);
  try {
    auto path = root / "map.rhm";
    make_map(path, 0);
    RibbonHeightMap map(path);
    double s = 1.6125, R = .31;
    if (map.sample(s, 0.).valid || map.sample(s - .5 * R, 0.).valid)
      throw std::runtime_error(
          "fixture must contain centre and separate bracket holes");
    auto q = map.envelope(s, 0., 1., 0., R);
    if (!q.valid || !q.bridged_missing_centre)
      throw std::runtime_error("closer two-sided measured support rejected");
    double brute = -1e9;
    for (int i = 0; i <= 200000; i++) {
      double u = -.995 * R + 1.99 * R * i / 200000.;
      auto a = map.sample(s + u, 0.);
      if (a.valid)
        brute = std::max(brute, a.height_m - R + std::sqrt(R * R - u * u));
    }
    if (std::abs(q.height_m - brute) > 5e-5)
      throw std::runtime_error(
          "bridge envelope differs from measured-only dense search");
    if (map.sample(s, 0.).valid)
      throw std::runtime_error("bridge filled a missing map cell");
    for (int kind : {1, 2}) {
      make_map(path, kind);
      RibbonHeightMap bad(path);
      if (bad.envelope(s, 0., 1., 0., R).valid)
        throw std::runtime_error("wide or one-sided unknown coverage accepted");
    }
    std::cout << "separate endpoint hole: closer support accepted; "
                 "measured-only envelope error "
              << std::abs(q.height_m - brute)
              << " m; wide and one-sided gaps rejected; map cells unchanged\n";
  } catch (...) {
    std::filesystem::remove_all(root);
    throw;
  }
  std::filesystem::remove_all(root);
}
