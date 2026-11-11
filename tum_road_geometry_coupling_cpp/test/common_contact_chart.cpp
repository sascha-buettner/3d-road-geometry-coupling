#include "tum_road_geometry_coupling_cpp/road_geometry_coupler.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace tam::road_geometry_coupling;
int main() {
  constexpr double R = 100., L = 2 * M_PI * R;
  constexpr int N = 2000;
  TrackData d;
  for (int i = 0; i <= N; ++i) {
    double a = 2 * M_PI * i / N;
    d.s_m.push_back(R * a);
    d.ref_line_x_m.push_back(R * sin(a));
    d.ref_line_y_m.push_back(R * (1 - cos(a)));
    d.ref_line_z_m.push_back(0.);
    d.theta_rad.push_back(atan2(sin(a), cos(a)));
    d.mu_rad.push_back(0.);
    d.phi_rad.push_back(0.);
    d.omega_x_radpm.push_back(0.);
    d.omega_y_radpm.push_back(0.);
    d.omega_z_radpm.push_back(1 / R);
    d.track_width_left_m.push_back(5.);
    d.track_width_right_m.push_back(-5.);
  }
  Track track(d);
  auto pit = d;
  for (auto &z : pit.ref_line_z_m)
    z = 1.;
  for (auto &v : pit.track_width_left_m)
    v = 10.;
  RoadGeometryCoupler model(d, pit);
  auto pmg = model.get_param_manager();
  pmg->set_value("road_plane.allow_pitlane_switching", false);
  pmg->set_value("road_plane.lookbehind_points", std::int64_t{10});
  pmg->set_value("road_plane.initial_pose.enabled", true);
  pmg->set_value("road_plane.initial_pose.x_m", 0.);
  pmg->set_value("road_plane.initial_pose.y_m", 0.);
  pmg->set_value("road_plane.initial_pose.yaw_rad", 0.);
  RoadPlaneGeometryTransformer contact(&track, {0., 0., 0., 0.}, {0., 0., 0.},
                                       10);
  double maximum = 0.;
  for (double s = .123456; s < 3 * L; s += .2) {
    double a = s / R, n = s > 20. ? 6. : 0.;
    CartesianPose2D p{(R - n) * sin(a), R - (R - n) * cos(a),
                      atan2(sin(a), cos(a))};
    auto q = contact.transform_to_global_curvilinear_pose(p).pose;
    tam::types::control::Odometry od;
    od.position_m = {p.x_m, p.y_m, 0.};
    od.orientation_rad = {0., 0., p.yaw_rad};
    od.velocity_mps = {40., 0., 0.};
    tam::types::control::AccelerationwithCovariances ac;
    model.set_odometry(od);
    model.set_acceleration(ac);
    model.step();
    auto got = model.get_transformed_odometry();
    double x = track.ref_line_x(q.s_m) + q.n_m * track.normal_x(q.s_m),
           y = track.ref_line_y(q.s_m) + q.n_m * track.normal_y(q.s_m);
    maximum = std::max(maximum,
                       std::hypot(x - got.position_m.x, y - got.position_m.y));
    if (std::abs(got.position_m.z) > 1e-9)
      throw std::runtime_error(
          "Closed-lap body jumped to auxiliary pit surface");
  }
  std::cout << "three-lap common contact/body chart maximum XY mismatch_m "
            << maximum << '\n';
  if (maximum > 1e-8)
    throw std::runtime_error("Body/contact initialization differs");
}
