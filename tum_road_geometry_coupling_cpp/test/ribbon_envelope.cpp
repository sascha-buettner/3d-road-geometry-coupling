#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>

#include "tum_road_geometry_coupling_cpp/ribbon_height_map.hpp"
using tam::road_geometry_coupling::RibbonHeightMap;
template <class T> void put(std::ofstream &f, T x) {
  f.write(reinterpret_cast<const char *>(&x), sizeof(x));
}
int main() {
  const auto path = std::filesystem::temp_directory_path() /
                    ("ribbon_envelope_" + std::to_string(getpid()) + ".rhm");
  try {
    {
      std::ofstream f(path, std::ios::binary);
      f.write("TAMRHG2\0", 8);
      put(f, std::uint32_t(128));
      put(f, std::uint32_t(0));
      for (double v : {.025, .025, -.5, .00005})
        put(f, v);
      put(f, std::uint64_t(128 * 41));
      put(f, std::uint64_t(0));
      for (int i = 0; i < 128; ++i) {
        put(f, std::int32_t(0));
        put(f, std::uint32_t(41));
        put(f, std::uint64_t(i * 41));
      }
      for (int i = 0; i < 128; ++i)
        for (int j = 0; j < 41; ++j)
          put(f, std::int16_t((i >= 64 ? 1600 : 0) + (j - 20) * 5));
    }
    RibbonHeightMap map(path);
    const double radius = .31;
    auto q = map.envelope(1.4, 0., 1., 0., radius);
    const double root = std::sqrt(radius * radius - .2 * .2);
    if (!q.valid || std::abs(q.height_m - (.08 - radius + root)) > 1e-10 ||
        std::abs(q.gradient_s - .2 / root) > 1e-9 ||
        std::abs(q.longitudinal_offset_m - .2) > 1e-10)
      throw std::runtime_error("round-tire step geometry mismatch");
    std::mt19937 rng(712);
    std::uniform_real_distribution<double> station(1.3, 1.8), lateral(-.2, .2),
        heading(-.6, .6);
    double max_height = 0., max_gradient = 0.;
    for (int i = 0; i < 250; ++i) {
      const double s = station(rng), n = lateral(rng), h = heading(rng),
                   c = std::cos(h), d = std::sin(h);
      auto e = map.envelope(s, n, c, d, radius);
      double brute = -1e9;
      for (int j = 0; j <= 20000; ++j) {
        double u = -.995 * radius + 1.99 * radius * j / 20000.;
        auto a = map.sample(s + c * u, n + d * u);
        if (a.valid)
          brute = std::max(brute, a.height_m - radius +
                                      std::sqrt(radius * radius - u * u));
      }
      max_height = std::max(max_height, std::abs(brute - e.height_m));
      const double eps = 1e-6;
      double gs = (map.envelope(s + eps, n, c, d, radius).height_m -
                   map.envelope(s - eps, n, c, d, radius).height_m) /
                  (2 * eps);
      double gn = (map.envelope(s, n + eps, c, d, radius).height_m -
                   map.envelope(s, n - eps, c, d, radius).height_m) /
                  (2 * eps);
      max_gradient = std::max({max_gradient, std::abs(gs - e.gradient_s),
                               std::abs(gn - e.gradient_n)});
      if (!e.valid || std::abs(brute - e.height_m) > 5e-5 ||
          std::abs(gs - e.gradient_s) > 1e-3 ||
          std::abs(gn - e.gradient_n) > 1e-3) {
        std::cerr << "CASE " << i << " s " << s << " n " << n << " heading "
                  << h << " height " << e.height_m << " brute " << brute
                  << " gradient " << e.gradient_s << "," << e.gradient_n
                  << " FD " << gs << "," << gn << " offset "
                  << e.longitudinal_offset_m << std::endl;
        throw std::runtime_error("continuous envelope/gradient mismatch");
      }
    }
    std::cout << "250 random step/heading cases; dense-grid height error "
              << max_height << " m, finite-difference gradient error "
              << max_gradient << "; analytic leading-edge geometry passed\n";
  } catch (...) {
    std::filesystem::remove(path);
    throw;
  }
  std::filesystem::remove(path);
}
