// Copyright 2025 Simon Sagmeister
#pragma once
#include <cmath>

#include <cstdint>
namespace tam::road_geometry_coupling {
inline std::size_t wrap_index(std::size_t const &index,
                              std::size_t const &max_index) {
  // Wrap the index to be within the range of the track size
  return index % max_index;
}
inline double wrap_index(double const &index, std::size_t const &max_index) {
  if (max_index == 0)
    return 0.0;
  double wrapped = std::fmod(index, static_cast<double>(max_index));
  if (wrapped < 0)
    wrapped += max_index;
  return wrapped;
}
} // namespace tam::road_geometry_coupling
