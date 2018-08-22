#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "constants.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

int ClosestWaypoint(double x,
                    double y,
                    const vector<double>& maps_x,
                    const vector<double>& maps_y) {

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x,
                 double y,
                 double theta,
                 const vector<double>& maps_x,
                 const vector<double>& maps_y) {

  int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y - y), (map_x - x));

  double angle = fabs(theta - heading);
  angle = min(2 * pi() - angle, angle);

  if (angle > pi() / 4) {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x,
                         double y,
                         double theta,
                         const vector<double>& maps_x,
                         const vector<double>& maps_y) {
  int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = maps_x.size() - 1;
  }

  double n_x = maps_x[next_wp] - maps_x[prev_wp];
  double n_y = maps_y[next_wp] - maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
  double proj_x = proj_norm * n_x;
  double proj_y = proj_norm * n_y;

  double frenet_d = distance(x_x, x_y, proj_x, proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000 - maps_x[prev_wp];
  double center_y = 2000 - maps_y[prev_wp];
  double centerToPos = distance(center_x, center_y, x_x, x_y);
  double centerToRef = distance(center_x, center_y, proj_x, proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s,
                     double d,
                     const vector<double>& maps_s,
                     const vector<double>& maps_x,
                     const vector<double>& maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
    prev_wp++;
  }

  int wp2 = (prev_wp + 1) % maps_x.size();

  double heading =
      atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s - maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

  double perp_heading = heading - pi() / 2;

  double x = seg_x + d * cos(perp_heading);
  double y = seg_y + d * sin(perp_heading);

  return {x, y};

}
std::vector<std::tuple<double, double>> getXY(std::vector<double> s,
                                              std::vector<double> d,
                                              const vector<double>& maps_s,
                                              const vector<double>& maps_x,
                                              const vector<double>& maps_y) {
  std::vector<int> id;
  for (int i = 0; i < s.size(); i++) {
    id.push_back(i);
  }
  std::sort(id.begin(),
            id.end(),
            [&s, &d](int a, int b) {
              if (s[a] < s[b]) return true;
              else if (s[a] > s[b]) return false;
              else if (d[a] < d[b]) return true;
              else if (d[a] > d[b]) return false;
              else return true;
            });

  std::vector<std::tuple<double, double>> ret(s.size());


  int wp1, wp2;
  wp2 = std::lower_bound(maps_s.begin(),
                         maps_s.end(),
                         s[id[0]]) - maps_s.begin();
  if(wp2==maps_s.size())
    wp2 = 0;
  wp1 = wp2 > 0 ? wp2 - 1 : maps_s.size() - 1;

  double
      heading = atan2((maps_y[wp2] - maps_y[wp1]), (maps_x[wp2] - maps_x[wp1]));
  double perp_heading = heading - pi() / 2;
  double sin_h = sin(heading);
  double cos_h = cos(heading);
  double sin_ph = sin(perp_heading);
  double cos_ph = cos(perp_heading);
  double s1(maps_s[wp1]), s2(maps_s[wp2]);
  int i = 0;
  while (i < id.size()) {
    int cid = id[i];
    double cs = s[cid];
    double cd = d[cid];
    while (wp1!=(maps_s.size()-1) && cs > s2) {
      wp1 = wp2;
      wp2++;
      if (wp2 == maps_s.size())
        wp2 = 0;
      heading = atan2((maps_y[wp2] - maps_y[wp1]), (maps_x[wp2] - maps_x[wp1]));
      perp_heading = heading - pi() / 2;
      sin_h = sin(heading);
      cos_h = cos(heading);
      sin_ph = sin(perp_heading);
      cos_ph = cos(perp_heading);
      s1 = s2;
      s2 = maps_s[wp2];
    }

    double seg_s = cs - s1;
    double seg_x = maps_x[wp1] + seg_s * cos_h;
    double seg_y = maps_y[wp1] + seg_s * sin_h;
    double x = seg_x + cd * cos_ph;
    double y = seg_y + cd * sin_ph;
    ret[cid] = std::make_tuple(x, y);
    i++;
  }
  return ret;
};
#include <random>

void testGetXY(const vector<double>& ms,const vector<double>& mx,const vector<double>& my) {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  std::vector<double> S,D,errors;
  std::vector<std::tuple<double,double>> ret;
  const int num_tests = 100000;
  for(int i=0;i<num_tests;i++) {
    double s = distribution(generator)*max_s;
    S.push_back(s);
    double d = distribution(generator)*12.0;
    D.push_back(d);
    auto pt = getXY(s,d,ms,mx,my);
    ret.push_back(std::make_tuple(pt[0],pt[1]));
  }
  auto ret1 = getXY(S,D,ms,mx,my);
  for(int i=0;i<ret.size();i++){
    double x1,y1,x2,y2,s,d;
    std::tie(x1,y1) = ret[i];
    std::tie(x2,y2) = ret1[i];
    s=S[i];
    d=D[i];
    double dx(x2-x1),dy(y2-y1);
    double error = sqrt(dx*dx+dy*dy);
    if(error>1e-5)
      std::cout<<"problem"<<std::endl;
    errors.push_back(error);
  }
  std::vector<double> errorStats(errors);
  std::sort(errorStats.begin(),errorStats.end());
  std::cout<<"error 50th percentile : "<<errorStats[num_tests/2]
           <<" 75th : "<<errorStats[3*(num_tests/4)]
           <<" 90th : "<<errorStats[9*(num_tests/10)]
           <<" 95th : "<<errorStats[19*(num_tests/20)]
           <<" 99th : "<<errorStats[99*(num_tests/100)]
           <<" 100th : "<<errorStats[num_tests-1]<<std::endl;
}
double lane_yaw(double s,std::vector<double>& ms,std::vector<double>& mx,std::vector<double>& my){
  int wp2 = std::lower_bound(ms.begin(),ms.end(),s)-ms.begin();
  if(wp2==ms.size())
    wp2 = 0;
  int wp1 = wp2>0?wp2-1:ms.size()-1;
  return atan2(my[wp2]-my[wp1],mx[wp2]-mx[wp1]);
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
                                                                          double)> sd);

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](
      uWS::WebSocket<uWS::SERVER> ws,
      char* data,
      size_t length,
      uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];

          std::vector<double>
              prev_x(previous_path_x.begin(), previous_path_x.end());
          std::vector<double>
              prev_y(previous_path_y.begin(), previous_path_y.end());

          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];
          std::vector<std::vector<double>> cars_data;
          for (auto sf:sensor_fusion)
            cars_data.push_back(std::vector<double>(sf.begin(), sf.end()));
          std::function<std::vector<std::tuple<double, double>>(std::vector<
              double>, std::vector<double>)> xy =
              [&map_waypoints_s, &map_waypoints_x, &map_waypoints_y](std::vector<
                  double> s,
                                                                     std::vector<
                                                                         double> d) {
                return getXY(s,
                             d,
                             map_waypoints_s,
                             map_waypoints_x,
                             map_waypoints_y);
              };

          std::function<std::vector<double>(double, double, double)> sd =
              [&map_waypoints_x, &map_waypoints_y](double x,
                                                   double y,
                                                   double yaw) {
                return getFrenet(x, y, yaw, map_waypoints_x, map_waypoints_y);
              };

          json msgJson;
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          std::tie(next_x_vals, next_y_vals) = path_plan(car_x,
                                                         car_y,
                                                         car_s,
                                                         car_d,
                                                         car_yaw,
                                                         car_speed,
                                                         lane_yaw(end_path_s,map_waypoints_s,map_waypoints_x,map_waypoints_y),
                                                         prev_x,
                                                         prev_y,
                                                         end_path_s,
                                                         end_path_d,
                                                         cars_data,
                                                         xy, sd);

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\"," + msgJson.dump() + "]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse* res, uWS::HttpRequest req, char* data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char* message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
