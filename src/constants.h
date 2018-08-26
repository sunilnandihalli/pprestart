#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <limits>
static const double jmax = 9.5;
static const double amin = -9.5;
static const double amax = 9.5;
static const double max_speed = 20;
static const double dt = 0.02;
static const double lane_width = 4.0;
static const double max_s = 6945.554;
static const int num_waypoints = 50;
static const double safe_dist = 2*lane_width;
static const double sentinel = std::numeric_limits<double>::max(); 
static double loop(double s) {
  while (s > max_s) s -= max_s;
  while (s < 0) s += max_s;
  return s;
};
#endif
