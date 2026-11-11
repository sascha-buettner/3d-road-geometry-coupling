// Copyright 2026 TUM Autonomous Motorsport
#pragma once
#include <cmath>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
namespace tam::road_geometry_coupling
{
// Packed strip grid; rows vary in lateral extent. SI units, periodic s, explicit
// missing cells. File byte order is little endian. No trees or allocation at query.
class RibbonHeightMap
{
public:
  struct Sample
  {
    double height_m{0}, gradient_s{0}, gradient_n{0};
    bool valid{false};
    double longitudinal_offset_m{0.};
    double cross_derivative{0.};
    bool bridged_missing_centre{false};
    double cell_upper_height_m{-std::numeric_limits<double>::infinity()};
  };

private:
  struct Row
  {
    std::int32_t first_n;
    std::uint32_t count;
    std::uint64_t offset;
  };
  std::uint32_t ns_{0};
  std::uint64_t reference_fingerprint_{0};
  double ds_{0}, dn_{0}, n_min_{0}, quantum_{0};
  std::vector<Row> rows_;
  std::vector<std::int16_t> heights_;
  // Conservative 8x8 vertex maxima are built once in memory, never serialized.
  // They only reject impossible contact candidates; the raw bilinear surface and
  // exact circle/root calculations still determine every accepted maximum.
  std::vector<std::int16_t> upper_tiles_;
  std::int64_t upper_first_n_{0};
  std::uint64_t upper_columns_{0};
  static std::int64_t floor_eighth(std::int64_t n) { return n >= 0 ? n / 8 : -((7 - n) / 8); }
  void build_upper_tiles()
  {
    std::int64_t first = std::numeric_limits<std::int64_t>::max(), last = -first;
    for (const auto & row : rows_)
      if (row.count) {
        first = std::min(first, floor_eighth(row.first_n));
        last = std::max(last, floor_eighth(std::int64_t(row.first_n) + row.count - 1));
      }
    if (last < first) return;
    const std::uint64_t columns = last - first + 1, tile_rows = (std::uint64_t(ns_) + 7) / 8;
    // Unusual sparse files with a huge lateral span must not cause a large cache
    // allocation. They retain the tight per-cell upper-bound path below.
    if (columns > 64000000 / tile_rows || columns * tile_rows > heights_.size()) return;
    upper_first_n_ = first;
    upper_columns_ = columns;
    upper_tiles_.assign(columns * tile_rows, std::numeric_limits<std::int16_t>::min());
    for (std::uint32_t i = 0; i < ns_; ++i) {
      const auto & row = rows_[i];
      const auto base = std::uint64_t(i / 8) * columns;
      for (std::uint32_t j = 0; j < row.count; ++j) {
        const auto column = floor_eighth(std::int64_t(row.first_n) + j) - first;
        auto & upper = upper_tiles_[base + column];
        upper = std::max(upper, heights_[row.offset + j]);
      }
    }
  }
  double rectangle_upper_height(double s0, double s1, double n0, double n1) const
  {
    const auto first_n =
      std::max(floor_eighth(std::int64_t(std::floor((n0 - n_min_) / dn_))), upper_first_n_);
    const auto last_n = std::min(
      floor_eighth(std::int64_t(std::floor((n1 - n_min_) / dn_)) + 1),
      upper_first_n_ + std::int64_t(upper_columns_) - 1);
    const auto last_s = std::int64_t(std::floor(s1 / ds_)) + 1;
    auto upper = std::numeric_limits<std::int16_t>::min();
    for (auto row = std::int64_t(std::floor(s0 / ds_)); row <= last_s;) {
      const auto wrapped = (row % std::int64_t(ns_) + ns_) % ns_;
      const auto base = std::uint64_t(wrapped / 8) * upper_columns_;
      for (auto n = first_n; n <= last_n; ++n)
        upper = std::max(upper, upper_tiles_[base + n - upper_first_n_]);
      row += std::min<std::int64_t>(8 - wrapped % 8, ns_ - wrapped);
    }
    return upper == std::numeric_limits<std::int16_t>::min()
             ? -std::numeric_limits<double>::infinity()
             : upper * quantum_;
  }
  double tile_upper_height(double s, double n) const
  {
    double ss = s;
    if (ss < 0. || ss >= length_m()) {
      ss = std::fmod(ss, length_m());
      if (ss < 0.) ss += length_m();
    }
    auto si = std::uint32_t(std::floor(ss / ds_));
    if (si >= ns_) si = 0;
    const auto ni = std::int64_t(std::floor((n - n_min_) / dn_));
    const auto left = floor_eighth(ni), right = floor_eighth(ni + 1);
    auto upper = std::numeric_limits<std::int16_t>::min();
    for (const auto row : {si / 8, ((si + 1) % ns_) / 8})
      for (const auto column : {left, right}) {
        const auto local = column - upper_first_n_;
        if (local >= 0 && std::uint64_t(local) < upper_columns_)
          upper = std::max(upper, upper_tiles_[std::uint64_t(row) * upper_columns_ + local]);
      }
    return upper == std::numeric_limits<std::int16_t>::min()
             ? -std::numeric_limits<double>::infinity()
             : upper * quantum_;
  }
  template <class T>
  static T read(std::ifstream & f)
  {
    T value{};
    f.read(reinterpret_cast<char *>(&value), sizeof(value));
    if (!f) throw std::runtime_error("Truncated ribbon height map");
    return value;
  }
  bool at(std::uint32_t s, std::int64_t n, double & value) const
  {
    const auto & row = rows_[s];
    if (n < row.first_n || n >= std::int64_t(row.first_n) + row.count) return false;
    const auto h = heights_[row.offset + std::uint64_t(n - row.first_n)];
    if (h == std::numeric_limits<std::int16_t>::min()) return false;
    value = h * quantum_;
    return true;
  }

public:
  explicit RibbonHeightMap(const std::string & path)
  {
    const std::uint16_t endian = 1;
    if (*reinterpret_cast<const unsigned char *>(&endian) != 1)
      throw std::runtime_error("Ribbon map loader requires little-endian host");
    std::ifstream f(path, std::ios::binary);
    std::array<char, 8> magic{};
    f.read(magic.data(), 8);
    if (std::memcmp(magic.data(), "TAMRHG2\0", 8) != 0)
      throw std::runtime_error("Invalid ribbon map magic");
    ns_ = read<std::uint32_t>(f);
    const auto reserved = read<std::uint32_t>(f);
    ds_ = read<double>(f);
    dn_ = read<double>(f);
    n_min_ = read<double>(f);
    quantum_ = read<double>(f);
    const auto samples = read<std::uint64_t>(f);
    reference_fingerprint_ = read<std::uint64_t>(f);
    if (
      reserved != 0 || ns_ < 2 || ns_ > 10000000 || samples > 1000000000 || !std::isfinite(ds_) ||
      !std::isfinite(dn_) || !std::isfinite(n_min_) || !std::isfinite(quantum_) || ds_ <= 0 ||
      dn_ <= 0 || quantum_ <= 0)
      throw std::runtime_error("Invalid ribbon map dimensions");
    const auto header_end = f.tellg();
    f.seekg(0, std::ios::end);
    const auto file_bytes = f.tellg();
    const auto expected_bytes = std::uint64_t(header_end) + std::uint64_t(ns_) * 16 + samples * 2;
    if (file_bytes < 0 || std::uint64_t(file_bytes) != expected_bytes)
      throw std::runtime_error("Invalid ribbon map file size");
    f.seekg(header_end);
    rows_.resize(ns_);
    for (auto & row : rows_) {
      row.first_n = read<std::int32_t>(f);
      row.count = read<std::uint32_t>(f);
      row.offset = read<std::uint64_t>(f);
      if (row.offset > samples || row.count > samples - row.offset)
        throw std::runtime_error("Invalid ribbon row");
    }
    heights_.resize(samples);
    f.read(reinterpret_cast<char *>(heights_.data()), samples * sizeof(std::int16_t));
    if (!f || f.peek() != std::char_traits<char>::eof())
      throw std::runtime_error("Invalid ribbon payload length");
    build_upper_tiles();
  }
  std::uint64_t reference_fingerprint() const { return reference_fingerprint_; }
  // FNV-1a binds the exact small reference file to this map at load time.
  // SHA-256 provenance remains in the asset manifest for archival verification.
  static std::uint64_t fingerprint_file(const std::string & path)
  {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("Cannot read map reference");
    std::array<char, 16384> bytes{};
    std::uint64_t hash = 14695981039346656037ULL;
    while (stream) {
      stream.read(bytes.data(), bytes.size());
      for (std::streamsize i = 0; i < stream.gcount(); ++i) {
        hash ^= static_cast<unsigned char>(bytes[i]);
        hash *= 1099511628211ULL;
      }
    }
    if (!stream.eof()) throw std::runtime_error("Cannot hash map reference");
    return hash;
  }
  double length_m() const { return ns_ * ds_; }
  std::size_t storage_bytes() const
  {
    return rows_.size() * sizeof(Row) + heights_.size() * 2 + upper_tiles_.size() * 2;
  }
  // Circular longitudinal support over the continuous bilinear surface. The
  // grid is only an index: candidates include every crossed cell edge and all
  // stationary points inside each cell. No online mesh/tree search or allocation.
  // ds_du and dn_du map one metre along the wheel heading into ribbon coordinates.
  Sample envelope(
    double s, double n, double ds_du, double dn_du, double radius, bool allow_bounds = true) const
  {
    if (!(std::isfinite(radius) && radius > 0.) || !std::isfinite(ds_du) || !std::isfinite(dn_du))
      return {};
    if (!std::isfinite(s) || !std::isfinite(n)) return {};
    if (s < 0. || s >= length_m()) {
      s = std::fmod(s, length_m());
      if (s < 0.) s += length_m();
    }
    // Reject impossible query extents before converting grid indices to integers.
    if (
      !std::isfinite(radius * (std::abs(ds_du) / ds_ + std::abs(dn_du) / dn_)) ||
      radius * (std::abs(ds_du) / ds_ + std::abs(dn_du) / dn_) > 1000000.)
      return {};
    const auto centre = sample(s, n);
    // A finite tire can bridge a compact unobserved patch on measured support.
    // Require support on both sides within half the wheel radius. No height is
    // invented in the gap: only measured cells enter the circular maximum. This
    // treats unknown space as unsupported, so it cannot certify a hidden bump.
    // Exported mesh holes remain unchanged; bridges are explicitly flagged.
    // Wider unknown regions and one-sided map boundaries remain invalid.
    bool left_support = true, right_support = true;
    const double bracket = .5 * radius;
    if (!centre.valid) {
      left_support = sample(s - ds_du * bracket, n - dn_du * bracket).valid;
      right_support = sample(s + ds_du * bracket, n + dn_du * bracket).valid;
    }
    struct Edge
    {
      double u;
      int axis;
    };
    std::array<Edge, 512> edges;  // Only initialized entries [0,count) are ever read.
    std::size_t count = 0;
    double limit = radius * .995;
    if (allow_bounds && centre.valid && !upper_tiles_.empty()) {
      const double max_height = rectangle_upper_height(
        s - std::abs(ds_du) * limit, s + std::abs(ds_du) * limit, n - std::abs(dn_du) * limit,
        n + std::abs(dn_du) * limit);
      const double rise = max_height - centre.height_m + 1e-9;
      // No point outside this interval can exceed the centre: even the highest
      // measured vertex is below the circle there. This is an exact support
      // exclusion, not a shortened physical tire or a smoothed terrain map.
      if (rise >= 0. && rise < radius)
        limit = std::min(limit, std::sqrt(std::max(0., 2. * radius * rise - rise * rise)) + 1e-9);
    }
    edges[count++] = {-limit, 0};
    edges[count++] = {limit, 0};
    edges[count++] = {0., 0};
    auto add_edges = [&](
                       double coordinate, double direction, double spacing, double origin,
                       int axis) {
      if (std::abs(direction) < 1e-12) return true;
      const double a = coordinate - direction * limit, b = coordinate + direction * limit;
      const auto first = static_cast<std::int64_t>(std::ceil((std::min(a, b) - origin) / spacing));
      const auto last = static_cast<std::int64_t>(std::floor((std::max(a, b) - origin) / spacing));
      if (last - first > std::int64_t(edges.size() - count - 1)) return false;
      for (auto k = first; k <= last; ++k)
        edges[count++] = {(origin + k * spacing - coordinate) / direction, axis};
      return true;
    };
    if (!add_edges(s, ds_du, ds_, 0., 1) || !add_edges(n, dn_du, dn_, n_min_, 2)) return {};
    std::sort(edges.begin(), edges.begin() + count, [](auto a, auto b) { return a.u < b.u; });
    std::size_t unique = 0;
    for (std::size_t i = 0; i < count; ++i) {
      if (unique && std::abs(edges[i].u - edges[unique - 1].u) < 1e-10) {
        if (edges[i].axis) edges[unique - 1].axis = edges[i].axis;
      } else
        edges[unique++] = edges[i];
    }
    count = unique;
    if (!centre.valid && !(left_support && right_support)) {
      // A fixed bracket endpoint can itself lie in a separate small scan hole.
      // Find actual measured support on each side within the SAME half-radius
      // bound. Cell interiors have constant coverage, so one clipped midpoint
      // per crossed bilinear cell is sufficient. No missing height is filled.
      for (std::size_t i = 0; i + 1 < count && !(left_support && right_support); ++i) {
        const double lo = std::max(edges[i].u, -bracket);
        const double hi = std::min(edges[i + 1].u, bracket);
        if (hi - lo < 1e-10) continue;
        const double mid = .5 * (lo + hi);
        if ((mid < 0. && left_support) || (mid > 0. && right_support)) continue;
        if (sample(s + ds_du * mid, n + dn_du * mid).valid) {
          if (mid < 0.) left_support = true;
          if (mid > 0.) right_support = true;
        }
      }
      if (!(left_support && right_support)) return {};
    }
    Sample best = centre;
    if (!best.valid) best.height_m = -std::numeric_limits<double>::infinity();
    best.longitudinal_offset_m = 0.;
    auto consider = [&](double u, int axis) {
      auto q = sample(s + ds_du * u, n + dn_du * u);
      if (!q.valid) return;
      const double root = std::sqrt(radius * radius - u * u);
      const double value = q.height_m - radius + root;
      if (value <= best.height_m) return;
      q.height_m = value;
      q.longitudinal_offset_m = u;
      // At a grid edge the optimum moves with that edge. Its derivative is
      // the circle slope, not the steep face slope on either adjacent cell.
      if (axis == 1 && std::abs(ds_du) > 1e-12)
        q.gradient_s = (u / root - q.gradient_n * dn_du) / ds_du;
      if (axis == 2 && std::abs(dn_du) > 1e-12)
        q.gradient_n = (u / root - q.gradient_s * ds_du) / dn_du;
      best = q;
    };
    for (std::size_t i = 0; i + 1 < count; ++i) {
      const double lo = edges[i].u, hi = edges[i + 1].u, width = hi - lo;
      if (width < 1e-10) continue;
      // The mixed derivative gives the quadratic coefficient along this line.
      // One interior lookup is sufficient for the whole crossed bilinear cell.
      const double mid = .5 * (lo + hi);
      const double nearest = lo > 0. ? lo : (hi < 0. ? hi : 0.);
      // Circle sag >= u^2/(2R). Combining that lower sag bound with a
      // conservative height maximum gives an upper bound without square roots.
      const double minimum_sag = nearest * nearest / (2. * radius);
      if (
        allow_bounds && !upper_tiles_.empty() &&
        tile_upper_height(s + ds_du * mid, n + dn_du * mid) - minimum_sag + 1e-10 < best.height_m)
        continue;
      const auto q = sample(s + ds_du * mid, n + dn_du * mid);
      if (allow_bounds && q.valid && q.cell_upper_height_m - minimum_sag + 1e-10 < best.height_m)
        continue;
      consider(lo, edges[i].axis);
      consider(hi, edges[i + 1].axis);
      if (!q.valid) continue;
      const double b = 2. * q.cross_derivative * ds_du * dn_du;
      const double a = q.gradient_s * ds_du + q.gradient_n * dn_du - b * mid;
      if (std::abs(b) < 1e-10) {
        const double u = radius * a / std::sqrt(1. + a * a);
        if (u > lo && u < hi) consider(u, 0);
        continue;
      }
      auto derivative = [&](double u) {
        return a + b * u - u / std::sqrt(radius * radius - u * u);
      };
      std::array<double, 4> splits{lo, hi, hi, hi};
      std::size_t size = 2;
      if (b > 1. / radius) {
        const double critical =
          std::sqrt(std::max(0., radius * radius - std::pow(radius * radius / b, 2. / 3.)));
        for (double u : {-critical, critical})
          if (u > lo && u < hi) splits[size++] = u;
      }
      std::sort(splits.begin(), splits.begin() + size);
      for (std::size_t k = 0; k + 1 < size; ++k) {
        double left = splits[k], right = splits[k + 1];
        if (derivative(left) < 0. || derivative(right) > 0.) continue;
        // Analytic linear-road solution is also a strong Newton start.
        double u = std::clamp(radius * a / std::sqrt(1. + a * a), left, right);
        for (int iteration = 0; iteration < 8; ++iteration) {
          const double root = std::sqrt(radius * radius - u * u), f = derivative(u);
          if (f > 0.)
            left = u;
          else
            right = u;
          if (std::abs(f) < 1e-10) break;
          const double slope = b - radius * radius / (root * root * root);
          const double candidate = u - f / slope;
          u = std::isfinite(candidate) && candidate > left && candidate < right
                ? candidate
                : .5 * (left + right);
        }
        consider(u, 0);
      }
    }
    best.bridged_missing_centre = best.valid && !centre.valid;
    return best;
  }
  Sample sample(double s_m, double n_m) const
  {
    if (!std::isfinite(s_m) || !std::isfinite(n_m)) return {};
    double ss = s_m;
    if (ss < 0. || ss >= length_m()) {
      ss = std::fmod(ss, length_m());
      if (ss < 0.) ss += length_m();
    }
    const double u = ss / ds_, v = (n_m - n_min_) / dn_;
    if (
      v < double(std::numeric_limits<std::int32_t>::min()) ||
      v > double(std::numeric_limits<std::int32_t>::max()))
      return {};
    auto si = static_cast<std::uint32_t>(std::floor(u));
    if (si >= ns_) si = 0;
    const auto ni = static_cast<std::int64_t>(std::floor(v));
    const double a = u - std::floor(u), b = v - std::floor(v);
    double h00, h01, h10, h11;
    if (
      !at(si, ni, h00) || !at(si, ni + 1, h01) || !at((si + 1) % ns_, ni, h10) ||
      !at((si + 1) % ns_, ni + 1, h11))
      return {};
    return {
      (1 - a) * ((1 - b) * h00 + b * h01) + a * ((1 - b) * h10 + b * h11),
      ((1 - b) * (h10 - h00) + b * (h11 - h01)) / ds_,
      ((1 - a) * (h01 - h00) + a * (h11 - h10)) / dn_,
      true,
      0.,
      (h11 - h10 - h01 + h00) / (ds_ * dn_),
      false,
      std::max({h00, h01, h10, h11})};
  }
};
}  // namespace tam::road_geometry_coupling
