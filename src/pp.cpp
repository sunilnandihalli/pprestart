#include <cmath>
#include <functional>
#include <utility>
#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include "spline.h"
#include "constants.h"
#include "jmt.h"

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


std::tuple<std::vector<double>, std::vector<double>>
discretizeSpline(tk::spline& s, std::vector<double>& s_vals) {
  std::vector<double> xs, ys;
  double sp(0), x0(0), y0(s(x0));
  for (int i = 0; i < s_vals.size(); i++) {
    double ds = s_vals[i] - sp;
    double x1p = x0 + ds;
    double y1p = s(x1p);
    double dx1p = x1p - x0;
    double dy1p = y1p - y0;
    double ds1p = sqrt(dx1p * dx1p + dy1p * dy1p);
    double x1 = x0 + dx1p * ds / ds1p;
    double y1 = s(x1);
    xs.push_back(x1);
    ys.push_back(y1);
    sp = s_vals[i];
    x0 = x1;
    y0 = y1;
  }
  return std::make_tuple(xs, ys);
};


// assuming the car is oriented along the spline and wants to achieve a velocity
// of v1 and zero acceleration covering a distance of d in time t
// pts-x,pts-y,accs,vels
std::tuple<std::vector<double>,
           std::vector<double>,
           std::vector<double>,
           std::vector<double>,
           bool>
jerkFreePoints(
    tk::spline& s, double a0, double v0, double v1, double d) {
  bool goalAchieved;
  std::vector<double> xs, ys;
  double x0(0), y0(s(x0));
  double s0 = 0;
  std::vector<double> s_vals;
  std::vector<double> v_vals, a_vals;
  std::tie(s_vals, a_vals, v_vals, goalAchieved) = path(a0, v0, v1, d);
  std::tie(xs, ys) = discretizeSpline(s, s_vals);
  return std::make_tuple(xs, ys, a_vals, v_vals, goalAchieved);
}

// spline,refx,refy,theta,lastx in reference coordinates
std::tuple<tk::spline, double, double, double, double> getSpline(double x,
                                                                 double y,
                                                                 double yaw,
                                                                 double sf,
                                                                 double df,
                                                                 std::function<
                                                                     std::vector<
                                                                         std::tuple<
                                                                             double,
                                                                             double>>(
                                                                         std::vector<
                                                                             double>,
                                                                         std::vector<
                                                                             double>)> xy) {
  std::cout << " getSpline x : " << x << " y : " << y << " yaw : " << yaw
            << " sf : " << sf << " df : " << df << std::endl;
  double theta = yaw; // transform angle
  double r = 0.5;
  auto pts = xy({sf, sf - 6 * r}, {df, df});
  double xl, yl, xl_1, yl_1;
  std::tie(xl, yl) = pts[0];
  std::tie(xl_1, yl_1) = pts[1];
  std::vector<double> xs({x, x + r * cos(yaw), xl_1, xl}),
      ys({y, y + r * sin(yaw), yl_1, yl}), x1s, y1s;
  std::tie(x1s, y1s) = toRefFrame(xs, ys, theta, x, y);
  tk::spline spline;
  std::cout << " spline points : ";
  for (int i = 0; i < 4; i++)
    std::cout << "(" << x1s[i] << "  " << y1s[i] << ")" << " ";
  std::cout << std::endl;
  spline.set_points(x1s, y1s);
  return std::make_tuple(spline, x, y, theta, x1s.back());
}
// xpts,ypts,accs,vels
std::tuple<std::vector<double>,
           std::vector<double>,
           std::vector<double>,
           std::vector<double>, bool> smoothPath(
    double x,
    double y,
    double yaw,
    double a0,
    double v0,
    double vf,
    double sf,
    double df,
    std::function<std::vector<std::tuple<double, double>>(std::vector<double>,
                                                          std::vector<double>)> xy) {
  tk::spline spline;
  double xref, yref, theta, spline_x_max;
  std::tie(spline, xref, yref, theta, spline_x_max) =
      getSpline(x, y, yaw, sf, df, xy);
  double totalDist = length(spline, 0, spline_x_max);
  std::vector<double> x1s, y1s;
  std::vector<double> a_vals, v_vals;
  bool goalAchieved;
  std::tie(x1s, y1s, a_vals, v_vals, goalAchieved) =
      jerkFreePoints(spline, a0, v0, vf, totalDist);
  std::vector<double> xs, ys;
  std::tie(xs, ys) = fromRefFrame(x1s, y1s, theta, xref, yref);
  return std::make_tuple(xs, ys, a_vals, v_vals, goalAchieved);
}

struct car {
  double s, d, v;
  car(double _s, double _d, double _v) : s(_s), d(_d), v(_v) {}
};

// s , d , v
car transformCarData(const std::vector<double>& c/*carData*/, double delta_t) {
  //id,x,y,vx,vy,s,d
  double vx = c[3];
  double vy = c[4];
  double v = sqrt(vx * vx + vy * vy);
  double s = c[5];
  double d = c[6];
  return car(s + v * delta_t, d, v);
};


//double s,double d,double sdot
std::vector<car>
transformCarsData(std::vector<std::vector<double>>& cars_data, double delta_t) {
  std::vector<car> ret;
  for (std::vector<double>& cdata:cars_data) {
    ret.push_back(transformCarData(cdata, delta_t));
    if (ret.back().d < 0 || ret.back().d > 12)
      ret.pop_back();
  }
  return ret;
};

struct lane_car_data {
  std::vector<int> lane_choices;
  std::vector<car*> front;
  std::vector<car*> back;
  lane_car_data(std::vector<int> _lane_choices,
                std::vector<car*> _front,
                std::vector<car*> _back)
      : lane_choices(_lane_choices), front(_front), back(_back) {}
};

double sdist(double from, double to) {
  if (to >= from) {
    double d1 = to - from;
    double d2 = max_s - (to - from);
    if (d1 < d2)
      return d1;
    else
      return -d2;
  } else
    return -sdist(to, from);
}
int lane_id(double d) {
  double lane_id;
  modf(d / lane_width, &lane_id);
  return int(lane_id);
};

lane_car_data possibleLanes(car ego, std::vector<car>& cardata) {
  int ego_lane_id = lane_id(ego.d);
  lane_car_data
      lc({}, {nullptr, nullptr, nullptr}, {nullptr, nullptr, nullptr});
  for (car& c:cardata) {
    int lid = lane_id(c.d);
    double ego_dist = sdist(ego.s, c.s);
    if (ego_dist > 0 && (lc.front[lid] == nullptr
        || sdist(ego.s, lc.front[lid]->s) > ego_dist))
      lc.front[lid] = &c;
    if (ego_dist < 0
        && (lc.back[lid] == nullptr
            || sdist(ego.s, lc.back[lid]->s) < ego_dist))
      lc.back[lid] = &c;
  }
  for (int lid = std::max(0, ego_lane_id - 1);
       lid <= std::min(2, ego_lane_id + 1); lid++) {
    car* f(lc.front[lid]), * b(lc.back[lid]);
    if (lid == ego_lane_id ||
        ((f == nullptr || sdist(ego.s, f->s) > lane_width) &&
            (b == nullptr
                || sdist(b->s, ego.s) > (lane_width + (b->v - ego.v) * 3))))
      lc.lane_choices.push_back(lid);
  }
  return lc;
}

struct lane_change_path {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> a;
  std::vector<double> v;
  double sf, lane_id;
};

lane_change_path get_path(int tlid/*target_lane_id*/,
                          int elid/*ego_lane_id*/,
                          lane_car_data& lc,
                          double a0,
                          double v0,
                          double x0,
                          double y0,
                          double yaw0,
                          double s0,
                          double d0,
                          double lane_yaw,
                          std::function<std::vector<std::tuple<double, double>>(
                              std::vector<double>,
                              std::vector<double>)> xy) {
  car* b(lc.back[tlid]), * f(lc.front[tlid]);
  lane_change_path ret;
  double lctime = 3.0;
  double target_d = (tlid + 0.5) * lane_width;
  double target_v, target_s;
  bool goalAchieved;
  ret.lane_id = tlid;
  if (f && b) {

  } else if (f) {
    target_v = max_speed;
    target_s = loop(s0 + (target_v + v0) * 0.5 * lctime);
    tk::spline s;
    double refx, refy, theta, spline_x_max;
    std::tie(s, refx, refy, theta, spline_x_max) =
        getSpline(x0, y0, yaw0, target_s, target_d, xy);
    double dist = length(s, 0, spline_x_max);
    std::function<double(double)> jfn;
    double timeTaken;
    std::tie(jfn, timeTaken, goalAchieved) =
        achieveTargetVelocityAndDistanceInShortestTime(a0,
                                                       v0,
                                                       target_v,
                                                       0,
                                                       max_speed,
                                                       dist);
    double distanceBetweenEgoAndOtherCar =
        sdist(target_s, loop(f->s + timeTaken * f->v));
    if (distanceBetweenEgoAndOtherCar > (safe_dist
        + std::get<0>(distDuringVelocityChange(0, target_v, f->v)))) {
      std::vector<double> s_vals;
      std::tie(s_vals, ret.a, ret.v) =
          genCompletePath(jfn, a0, v0, 0.0, timeTaken);
      std::vector<double> x1s, y1s;
      std::tie(x1s, y1s) = discretizeSpline(s, s_vals);
      std::tie(ret.x, ret.y) = fromRefFrame(x1s, y1s, theta, refx, refy);
    } else {
      target_v = f->v;
      target_s = loop(s0 + (target_v + v0) * 0.5 * lctime);
      std::tie(s, refx, refy, theta, spline_x_max) =
          getSpline(x0, y0, yaw0, target_s, target_d, xy);
      double spline_length = length(s, 0, spline_x_max);
      double es0 = loop(s0 - (spline_length - sdist(s0, target_s)));
      double effective_distance = sdist(es0, f->s);
      double changeInEffectiveDistance = effective_distance - lane_width;
      std::tie(jfn,timeTaken,goalAchieved) = achieveTargetVelocityAndDistanceInShortestTime(a0,
                                                     v0 - target_v,
                                                     0,
                                                     -target_v,
                                                     max_speed - target_v,
                                                     changeInEffectiveDistance);
      std::vector<double> s_vals;
      std::tie(s_vals,ret.a,ret.v) = genCompletePath(jfn,a0,v0,0.0,timeTaken);
      while(s_vals.back()<spline_length) {
        double sl = s_vals.back();
        double vl = ret.v.back();
        s_vals.push_back(sl+vl*dt);
        ret.a.push_back(0.0);
        ret.v.push_back(vl);
      }
    }
  } else if (b) {

  } else {
    target_v = max_speed;
    target_s = s0 + target_v * 3;
    std::tie(ret.x, ret.y, ret.a, ret.v, goalAchieved) =
        smoothPath(x0, y0, yaw0, a0, v0, target_v, target_s, target_d, xy);
    ret.sf = target_s;
  }
  return ret;
}


std::pair<std::vector<double>, std::vector<double>> path_plan(double car_x,
                                                              double car_y,
                                                              double car_s,
                                                              double car_d,
                                                              double car_yaw,
                                                              double car_speed,
                                                              double lane_yaw,
                                                              std::vector<double>& prev_x,
                                                              std::vector<double>& prev_y,
                                                              double end_path_s,
                                                              double end_path_d,
                                                              std::vector<std::vector<
                                                                  double>> cars_data,
                                                              std::function<std::vector<
                                                                  std::tuple<
                                                                      double,
                                                                      double>>(
                                                                  std::vector<
                                                                      double>,
                                                                  std::vector<
                                                                      double>)> xy,
                                                              std::function<std::vector<
                                                                  double>(double,
                                                                          double,
                                                                          double)> sd) {
  std::vector<double> xs, ys;
  unsigned long np = prev_x.size();
  assert(prev_x.size() == prev_y.size());
  double ego_x(np ? prev_x.back() : car_x), ego_y(np ? prev_y.back() : car_y);
  double ego_s(np ? end_path_s : car_s), ego_d(np ? end_path_d : car_d);
  double ego_yaw(np > 1 ? yaw(prev_x[np - 2],
                              prev_y[np - 2],
                              prev_x[np - 1],
                              prev_y[np - 1]) : car_yaw);
  int ego_lane_id = lane_id(ego_d);
  static double ego_v, ego_a;
  static bool first_call(true);
  static int target_lane_id;
  if (first_call) {
    ego_v = car_speed;
    ego_a = 0.0;
    target_lane_id = ego_lane_id;
  }
  double target_d = (target_lane_id + 0.5) * lane_width;
  bool disable_lane_change_choice = fabs(ego_d - target_d) > 0.5;
  std::vector<car> transformedCarsData =
      transformCarsData(cars_data, np * dt);
  lane_car_data
      lcdata = possibleLanes(car(ego_s, ego_d, ego_v), transformedCarsData);
  lane_change_path change_path;
  if (disable_lane_change_choice) {
    change_path = get_path(target_lane_id,
                           ego_lane_id,
                           lcdata,
                           ego_a,
                           ego_v,
                           ego_x,
                           ego_y,
                           ego_yaw,
                           ego_s,
                           ego_d,
                           lane_yaw,
                           xy);
  } else {
    std::map<int, lane_change_path> options;
    int max_points = 0;

    for (int lid:lcdata.lane_choices) {
      options[lid] = get_path(lid,
                              ego_lane_id,
                              lcdata,
                              ego_a,
                              ego_v,
                              ego_x,
                              ego_y,
                              ego_yaw,
                              ego_s,
                              ego_d,
                              lane_yaw, xy);
      if (options[lid].x.size() > max_points)
        max_points = options[lid].x.size();
    }
    std::function<double(int)>
        option_dist = [&options, ego_s, max_points](int lid) {
      const lane_change_path& lc = options[lid];
      return sdist(ego_s,
                   loop(lc.sf + lc.v.back() * (max_points - lc.x.size()) * dt));
    };
    int best_lid = ego_lane_id;
    double best_dist = option_dist(ego_lane_id);

    for (auto lid_lc : options) {
      if (lid_lc.first == ego_lane_id)
        continue;
      double cur_dist = option_dist(lid_lc.first);
      if (cur_dist > best_dist) {
        best_lid = lid_lc.first;
        best_dist = cur_dist;
      }
    }
    std::cout << "best_lane_id : " << best_lid << " best_dist : " << best_dist
              << " options : " << options.size() << std::endl;
    if (options.count(best_lid) != 1) {
      throw ("best_lane_id not found");
    } else {
      change_path = options[best_lid];
    }
  }

  xs = prev_x;
  ys = prev_y;
  for (int i = prev_x.size(), j = 0; i < num_waypoints; i++, j++) {
    xs.push_back(change_path.x[j]);
    ys.push_back(change_path.y[j]);
    ego_a = change_path.a[j];
    ego_v = change_path.v[j];
  }
  first_call = false;
  return std::make_pair(xs, ys);
}
