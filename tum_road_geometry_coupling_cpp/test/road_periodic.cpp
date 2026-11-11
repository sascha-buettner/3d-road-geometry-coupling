#include <cmath>

#include <iostream>
#include <stdexcept>

#include "tum_road_geometry_coupling_cpp/road_plane_transformer.hpp"
using namespace tam::road_geometry_coupling;
int main() {
  constexpr double R = 100., L = 2 * M_PI * R;
  constexpr int N = 2000;
  TrackData d;
  for (int i = 0; i <= N; ++i) {
    const double a = 2 * M_PI * i / N;
    d.s_m.push_back(R * a);
    d.ref_line_x_m.push_back(R * std::sin(a));
    d.ref_line_y_m.push_back(R * (1 - std::cos(a)));
    d.ref_line_z_m.push_back(0.);
    d.theta_rad.push_back(std::atan2(std::sin(a), std::cos(a)));
    d.mu_rad.push_back(0.);
    d.phi_rad.push_back(0.);
    d.omega_x_radpm.push_back(0.);
    d.omega_y_radpm.push_back(0.);
    d.omega_z_radpm.push_back(1 / R);
    d.track_width_left_m.push_back(5.);
    d.track_width_right_m.push_back(-5.);
  }
  Track track(std::move(d));
  RoadPlaneGeometryTransformer t(&track, {0., 0., 0., 0.}, {0., 0., 0.}, 10);
  double worst_s = 0., worst_n = 0.;
  int count = 0;
  for (double s = 0; s < 3 * L; s += .1) {
    const double a = s / R;
    t.transform_to_global_curvilinear_pose(
        {R * std::sin(a), R * (1 - std::cos(a)),
         std::atan2(std::sin(a), std::cos(a))});
    for (double offset : {0., -2.4, 1.8})
      for (double n : {-1., 0., 1.}) {
        const double q = (s + offset) / R;
        auto p = t.project_point(
            {(R - n) * std::sin(q), R - (R - n) * std::cos(q), q});
        const double error_s = std::remainder(p.pose.s_m - (s + offset), L);
        worst_s = std::max(worst_s, std::abs(error_s));
        worst_n = std::max(worst_n, std::abs(p.pose.n_m - n));
        ++count;
      }
  }
  std::cout << "three laps queries " << count << " max station error_m "
            << worst_s << " max lateral error_m " << worst_n << '\n';
  if (worst_s > .0002 || worst_n > .0002)
    throw std::runtime_error("Periodic front/rear projection failed");
}
