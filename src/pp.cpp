#include <cmath>
#include <functional>
#include <utility>
#include <vector>
#include "spline.h"
const double lane_width = 4.0;
std::pair<double, double> toRefFrame(double x, double y, double theta,
                                     double x0, double y0) {
  double dx(x - x0);
  double dy(y - y0);
  double x1 = dx * cos(theta) + dy * sin(theta);
  double y1 = -dx * sin(theta) + dy * cos(theta);
  return std::make_pair(x1, y1);
}

std::pair<std::vector<double>, std::vector<double>> toRefFrame(
    std::vector<double> xs, std::vector<double> ys, double theta, double x0,
    double y0) {
  std::vector<double> rx1s, ry1s;
  for (int i = 0; i < xs.size(); i++) {
    double x1, y1;
    std::tie(x1, y1) = toRefFrame(xs[i], ys[i], theta, x0, y0);
    rx1s.push_back(x1);
    ry1s.push_back(y1);
  }
  return std::make_pair(rx1s, ry1s);
}

std::pair<double, double> fromRefFrame(double x1, double y1, double theta,
                                       double x0, double y0) {
  double dx = x1 * cos(theta) - y1 * sin(theta);
  double dy = x1 * sin(theta) + y1 * cos(theta);
  double x = dx + x0;
  double y = dy + y0;
  return std::make_pair(x, y);
}

std::pair<std::vector<double>, std::vector<double>> fromRefFrame(
    std::vector<double> x1s, std::vector<double> y1s, double theta, double x0,
    double y0) {
  std::vector<double> rxs, rys;
  for (int i = 0; i < x1s.size(); i++) {
    double x, y;
    std::tie(x, y) = fromRefFrame(x1s[i], y1s[i], theta, x0, y0);
    rxs.push_back(x);
    rys.push_back(y);
  }
  return std::make_pair(rxs, rys);
}

double length(tk::spline& s, double x0, double x1) {
  const double tol = 0.1;
  const double distPerTimeStep = 0.2;
  std::function<double(double, double, double, double, double, double)> f =
      [&s, &f, &distPerTimeStep](double x0, double y0, double x1, double y1,
                                 double l, double tol) {
        if (l < distPerTimeStep) return l;
        double x = (x0 + x1) * 0.5;
        double y = s(x);
        double dx0(x - x0), dx1(x1 - x), dy0(y - y0), dy1(y1 - y);
        double lf1 = sqrt(dx0 * dx0 + dy0 * dy0);
        double lf2 = sqrt(dx1 * dx1 + dy1 * dy1);
        double lf = lf1 + lf2;
        if (fabs(l - lf) < tol)
          return lf;
        else
          return f(x0, y0, x, y, lf1, tol * 0.5) +
                 f(x, y, x1, y1, lf2, tol * 0.5);
      };
  double y0(s(x0)), y1(s(x1));
  double dx(x1 - x0), dy(y1 - y0);
  double l = sqrt(dx * dx + dy * dy);
  return f(x0, y0, x1, y1, l, tol);
}

double yaw(double x1, double y1, double x2, double y2) {
  std::atan2(y2 - y1, x2 - x1);
}
struct carData {
  double s, d, x, y, v;
  carData(double _s, double _d, double _x, double _y, double _v)
      : s(_s), d(_d), x(_x), y(_y), v(_v) {}
};

bool willCollideDuringLaneChange(carData ego,carData* carBehind,carData* carFront,std::function<std::vector<double>(double, double)>& xy) {
  
}

double cost(


// target d, target v
std::pair<std::vector<double>, std::vector<double>> choseLane(
    const carData& ego, const std::vector<carData>& cars,
    std::function<std::vector<double>(double, double)>& xy) {
  
}

std::vector<double> path(double a0, double v0, double v1, double delta_d);
// assuming the car is oriented along the spline and wants to achieve a velocity
// of v1 and zero acceleration covering a distance of d in time t
std::pair<std::vector<double>, std::vector<double>> jerkFreePoints(
    tk::spline& s, double v0, double a0, double v1, double d) {
  std::vector<double> xs, ys;
  double x0(0), y0(s(x0));
  double s0 = 0;
  auto s_vals = path(a0, v0, v1, d);
  for (int i = 0; i < s_vals.size(); i++) {
    double ds = (s_vals[i] - s0);
    double x1p = x0 + ds;
    double y1p = s(x1p);
    double dx1p = x1p - x0;
    double dy1p = y1p - y0;
    double ds1p = sqrt(dx1p * dx1p + dy1p * dy1p);
    double x1 = x0 + dx1p * ds / ds1p;
    double y1 = s(x1);
    xs.push_back(x1);
    ys.push_back(y1);
    s0 = s_vals[i];
    x0 = x1;
    y0 = y1;
  }
  return std::make_pair(xs, ys);
}

std::pair<std::vector<double>, std::vector<double>> smoothPath(
    double x, double y, double s, double d, double speed,
    std::function<std::vector<double>(double, double)> xy) {}

std::pair<std::vector<double>, std::vector<double>> path_plan(
    double max_s, double car_x, double car_y, double car_s, double car_d,
    double car_yaw, double car_speed, std::vector<double>& prev_x,
    std::vector<double>& prev_y, double end_path_s, double end_path_d,
    std::vector<std::vector<double>> cars_data,
    std::function<std::vector<double>(double, double)> xy) {
  std::pair<std::vector<double>, std::vector<double>> ret;
  return ret;
}
