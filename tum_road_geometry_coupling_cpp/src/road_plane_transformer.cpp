// Copyright 2025 Simon Sagmeister

#include "tum_road_geometry_coupling_cpp/road_plane_transformer.hpp"
namespace tam::road_geometry_coupling
{
RoadPlaneGeometryTransformer::RoadPlaneSegment
RoadPlaneGeometryTransformer::calculate_road_plane_segment(
  CartesianPose2D const & road_plane_starting_pose, std::size_t const & root_idx,
  std::size_t const & segment_length) const
{
  // Init return value
  RoadPlaneSegment segment;
  segment.root_idx = root_idx;
  segment.num_points = segment_length;

  // Localize cosy segment starting pos on the track
  auto const & s_track = track_->s_coord();
  auto const & omega_z_track = track_->omega_z();

  // Build the reference curve in the vehicle frame
  Eigen::VectorXd rpl_x(segment_length), rpl_y(segment_length), rpl_z(segment_length),
    rpl_theta(segment_length), rpl_sigma(segment_length);
  Eigen::VectorXd rpl_tx(segment_length), rpl_ty(segment_length), rpl_tz(segment_length);

  // Init starting point in the vectors
  rpl_x(0) = road_plane_starting_pose.x_m;
  rpl_y(0) = road_plane_starting_pose.y_m;
  rpl_z(0) = 0.0;  // equal to 0
  rpl_theta(0) = road_plane_starting_pose.yaw_rad;
  rpl_tx(0) = std::cos(rpl_theta(0));
  rpl_ty(0) = std::sin(rpl_theta(0));
  rpl_tz(0) = 0.0;     // equal to 0 due to projection to 2D plane
  rpl_sigma(0) = 0.0;  // sigma is not used
  for (size_t i = 1; i < segment_length; i++) {
    size_t idx{(root_idx + i) % (s_track.size()-1)};
    size_t prev_idx = (idx == 0) ? 0 : idx - 1;  // Safe access of index
    double s_dist{0.0};
    if (idx == 0) {
      prev_idx = s_track.size() - 2;
      s_dist = s_track.back() - s_track[prev_idx];
    } else {
      s_dist = s_track[idx] - s_track[prev_idx];
    }
    // Euler integration of the reference curve pose
    double delta_lat{0.0}, delta_lon{s_dist};
    double omega_z_median =
      (omega_z_track[prev_idx] + omega_z_track[idx]) * 0.5;  // Use the median curvature
    // Do a circular arc integration if the curvature is not zero
    // This avoid the problems of forward euler without
    // requiring microsteps
    if (std::abs(omega_z_median) > 1e-6) {
      delta_lat = (1 - std::cos(omega_z_median * s_dist)) / omega_z_median;
      delta_lon = std::sin(omega_z_median * s_dist) / omega_z_median;
    }
    rpl_x(i) = rpl_x(i - 1) + delta_lon * std::cos(rpl_theta(i - 1)) -
               delta_lat * std::sin(rpl_theta(i - 1));
    rpl_y(i) = rpl_y(i - 1) + delta_lon * std::sin(rpl_theta(i - 1)) +
               delta_lat * std::cos(rpl_theta(i - 1));
    rpl_theta(i) =
      tam::helpers::geometry::normalize_angle(rpl_theta(i - 1) + s_dist * omega_z_median);
    rpl_z(i) = 0.0;  // equal to 0
    // compute tangent vector based on reference curve heading
    rpl_tx(i) = cos(rpl_theta(i));
    rpl_ty(i) = sin(rpl_theta(i));
    rpl_tz(i) = 0.0;  // equal to 0 due to projection to 2D plane
    rpl_sigma(i) = rpl_sigma(i - 1) + s_dist;
  }
  segment.ref_line = tam::helpers::cosy::CurvilinearCosy::create(rpl_x, rpl_y, rpl_z)
                       ->set_tangent(rpl_tx, rpl_ty, rpl_tz)
                       ->set_s(rpl_sigma)
                       ->build();
  segment.sigma_m = rpl_sigma;

  return segment;
}
void RoadPlaneGeometryTransformer::extend_road_plane_segment(
  CartesianPose2D const & current_pose_road_plane)
{
  // Match onto the road plane ref line to get n and the index
  auto matching_result = current_segment_.ref_line->convert_to_sn_and_get_idx_global_2d(
    current_pose_road_plane.x_m, current_pose_road_plane.y_m, current_pose_road_plane.yaw_rad);
  size_t local_idx = std::get<1>(matching_result);
  local_idx = local_idx > lookbehind_points_ ? local_idx - lookbehind_points_ : 0;
  double sigma_start = current_segment_.sigma_m[local_idx];

  CartesianPose2D root_pose;
  auto start_pose = current_segment_.ref_line->convert_to_cartesian(sigma_start, 0, 0);
  root_pose.x_m = start_pose[0];
  root_pose.y_m = start_pose[1];
  root_pose.yaw_rad = start_pose[2];

  current_segment_ = calculate_road_plane_segment(
    root_pose, wrap_index(current_segment_.root_idx + local_idx, track_max_idx_));
}
void RoadPlaneGeometryTransformer::calculate_curvilinear_pose(
  CartesianPose2D const & pose, CurvilinearPoseWithDiagnostics & result)
{
  result = project_point(pose);
  if (result.diagnostics.current_local_index > REBUILD_RATIO * DEFAULT_SEGMENT_LENGTH) {
    extend_road_plane_segment(pose);
    result = project_point(pose);
  }
}
RoadPlaneGeometryTransformer::CurvilinearPoseWithDiagnostics
RoadPlaneGeometryTransformer::project_point(CartesianPose2D const & pose) const
{
  CurvilinearPoseWithDiagnostics result;
  const CurvilinearPose local(current_segment_.ref_line->convert_to_sn_and_get_idx_global(
    pose.x_m, pose.y_m, 0., pose.yaw_rad));
  result.pose.n_m = local.n_m;
  result.pose.chi_rad = local.chi_rad;
  result.pose.idx = wrap_index(local.idx + current_segment_.root_idx, track_max_idx_);
  // Preserve the double-precision local station. The helper's fractional
  // index is float and would quantize a 10,000-point track by about 0.25 mm.
  const double length = track_->s_coord().back();
  result.pose.s_m = std::fmod(track_->s_coord()[current_segment_.root_idx] + local.s_m, length);
  if (result.pose.s_m < 0.) result.pose.s_m += length;
  result.diagnostics.current_local_index = local.idx;
  result.diagnostics.current_sigma = local.s_m;
  result.diagnostics.current_n = local.n_m;
  result.diagnostics.segment_root_idx = current_segment_.root_idx;
  const auto xy = current_segment_.ref_line->convert_to_cartesian(local.s_m, 0., 0.);
  result.diagnostics.ref_line_pose = {xy[0], xy[1], xy[2]};
  return result;
}
RoadPlaneGeometryTransformer::RoadPlaneGeometryTransformer(
  Track const * track, CurvilinearPose const & curvilinear_pose_on_track,
  CartesianPose2D const & current_pose_road_plane, std::size_t lookbehind_points)
: track_{track},
  track_max_idx_{static_cast<std::size_t>(track_->s_coord().size() - 1)},
  lookbehind_points_{std::min(lookbehind_points, DEFAULT_SEGMENT_LENGTH / 4)}
{
  // Calculate the starting point of the roadplane segment so that it gives
  // the same curvilinear pose as on the track
  std::size_t root_index = static_cast<std::size_t>(curvilinear_pose_on_track.idx) % track_max_idx_;

  CartesianPose2D root_pose_road_plane;

  // Find the offset between the current s and the last point in the raceline file
  // (Start the segment directly on a point of the ref line to align the indices)
  const double fraction = curvilinear_pose_on_track.idx - std::floor(curvilinear_pose_on_track.idx);
  const double next_s =
    root_index + 1 < track_max_idx_ ? track_->s_coord()[root_index + 1] : track_->s_coord().back();
  double delta_s = fraction * (next_s - track_->s_coord()[root_index]);

  // Calculate the root pose of the road plane segment
  root_pose_road_plane.yaw_rad =
    current_pose_road_plane.yaw_rad - curvilinear_pose_on_track.chi_rad;
  root_pose_road_plane.x_m = current_pose_road_plane.x_m -
                             delta_s * std::cos(root_pose_road_plane.yaw_rad) +
                             curvilinear_pose_on_track.n_m * std::sin(root_pose_road_plane.yaw_rad);
  root_pose_road_plane.y_m = current_pose_road_plane.y_m -
                             delta_s * std::sin(root_pose_road_plane.yaw_rad) -
                             curvilinear_pose_on_track.n_m * std::cos(root_pose_road_plane.yaw_rad);

  // Integrate backwards with the same circular arcs used in forward segment
  // construction. This reserves physical space behind the vehicle at the seam.
  for (std::size_t k = 0; k < lookbehind_points_; ++k) {
    const auto previous = (root_index + track_max_idx_ - 1) % track_max_idx_;
    const double ds = root_index == 0 ? track_->s_coord().back() - track_->s_coord()[previous]
                                      : track_->s_coord()[root_index] - track_->s_coord()[previous];
    const double curvature = .5 * (track_->omega_z()[root_index] + track_->omega_z()[previous]);
    const double angle = root_pose_road_plane.yaw_rad - ds * curvature;
    double lon = ds, lat = 0.;
    if (std::abs(curvature) > 1e-6) {
      lon = std::sin(curvature * ds) / curvature;
      lat = (1 - std::cos(curvature * ds)) / curvature;
    }
    root_pose_road_plane.x_m -= lon * std::cos(angle) - lat * std::sin(angle);
    root_pose_road_plane.y_m -= lon * std::sin(angle) + lat * std::cos(angle);
    root_pose_road_plane.yaw_rad = angle;
    root_index = previous;
  }
  // Now calculate the segment starting from the corresponding root pose
  current_segment_ = calculate_road_plane_segment(root_pose_road_plane, root_index);
}
RoadPlaneGeometryTransformer::CurvilinearPoseWithDiagnostics
RoadPlaneGeometryTransformer::transform_to_global_curvilinear_pose(
  CartesianPose2D const & road_plane_pose)
{
  CurvilinearPoseWithDiagnostics result;

  // Calculate the curvilinear pose in the segment
  calculate_curvilinear_pose(road_plane_pose, result);

  return result;
}
}  // namespace tam::road_geometry_coupling
