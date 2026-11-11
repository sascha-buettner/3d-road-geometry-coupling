#include "tum_road_geometry_coupling_cpp/ribbon_height_map.hpp"

#include <unistd.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
using tam::road_geometry_coupling::RibbonHeightMap;
template <class T>
void put(std::ofstream & f, T value)
{
  f.write(reinterpret_cast<const char *>(&value), sizeof(value));
}
void near(double a, double b)
{
  if (std::abs(a - b) > 1e-12) throw std::runtime_error("bilinear grid mismatch");
}
int main()
{
  auto root =
    std::filesystem::temp_directory_path() / ("ribbon_map_test_" + std::to_string(getpid()));
  std::filesystem::create_directory(root);
  try {
    const auto path = root / "map.rhm", ref = root / "ref.csv";
    {
      std::ofstream f(ref);
      f << "hello";
    }
    {
      std::ofstream f(path, std::ios::binary);
      f.write("TAMRHG2\0", 8);
      put(f, std::uint32_t(4));
      put(f, std::uint32_t(0));
      for (double x : {1., .5, -.5, .0001}) put(f, x);
      put(f, std::uint64_t(12));
      put(f, RibbonHeightMap::fingerprint_file(ref));
      for (int i = 0; i < 4; ++i) {
        put(f, std::int32_t(0));
        put(f, std::uint32_t(3));
        put(f, std::uint64_t(i * 3));
      }
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
          put(f, std::int16_t(i == 3 && j == 2 ? -32768 : i * 100 + j * 10));
    }
    RibbonHeightMap map(path);
    auto a = map.sample(.3, -.2);
    if (!a.valid) throw std::runtime_error("valid sample rejected");
    near(a.height_m, .0036);
    near(a.gradient_s, .01);
    near(a.gradient_n, .002);
    near(map.sample(4.3, -.2).height_m, a.height_m);
    near(map.sample(-3.7, -.2).height_m, a.height_m);
    if (map.sample(3.1, .4).valid || map.sample(1., -1.).valid || map.sample(NAN, 0).valid)
      throw std::runtime_error("missing/outside accepted");
    near(double(map.reference_fingerprint() == 0xa430d84680aabd0bULL), 1.);
    std::filesystem::resize_file(path, std::filesystem::file_size(path) - 1);
    bool rejected = false;
    try {
      RibbonHeightMap broken(path);
    } catch (const std::runtime_error &) {
      rejected = true;
    }
    if (!rejected) throw std::runtime_error("truncated payload accepted");
    std::cout << "map interpolation, derivatives, periodicity, holes, reference fingerprint and "
                 "truncation passed\n";
  } catch (...) {
    std::filesystem::remove_all(root);
    throw;
  }
  std::filesystem::remove_all(root);
}
