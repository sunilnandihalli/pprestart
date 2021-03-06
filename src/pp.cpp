#include <cmath>
#include <functional>
#include <utility>
#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include "spline.h"
#include "constants.h"
#include "jmt.h"
#include <string>

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

// spline,refx,refy,theta,x of passed-in s_vals reference coordinates
std::tuple<tk::spline, double, double, double, double, std::vector<double>> getSpline(double x,
                                                                 double y,
									      double px,
									      double py,
                                                                 double yaw,
                                                                 std::vector<double> sf,
                                                                 std::vector<double> df,
                                                                 std::function<
                                                                     std::vector<
                                                                         std::tuple<
                                                                             double,
                                                                             double>>(
                                                                         std::vector<
                                                                             double>,
                                                                         std::vector<
                                                                             double>)> xy) {
  double theta; // transform angle
  double r = 0.5;
  auto pts = xy(sf, df);
 
  std::vector<double> xs,
      ys, x1s, y1s;
  double refx,refy;
  if(px!=sentinel && py!=sentinel) {
    xs = {px,x};
    ys = {py,y};
    theta = atan2(y-py,x-px);
    refx = px;
    refy = py;
  } else {
    xs = {x, x + r * cos(yaw)};
    ys = {y, y + r * sin(yaw)};
    theta = yaw;
    refx = x;
    refy = y;
  }
  for(auto pt:pts) {
    double x,y;
    std::tie(x,y) = pt;
    xs.push_back(x);
    ys.push_back(y);
  }
  double startx;
  std::tie(x1s, y1s) = toRefFrame(xs, ys, theta, refx, refy);
  if(px!=sentinel && py!=sentinel) {
    startx = x1s[1];
  } else {
    startx = x1s[0];
  }
  tk::spline spline;
  spline.set_points(x1s, y1s);
  return std::make_tuple(spline, refx, refy, theta, startx, std::vector<double>(x1s.begin()+2,x1s.end()));
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
  return lc;
}

struct lane_change_path {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> a;
  std::vector<double> v;
  double sf;
  int lane_id;
};

lane_change_path get_path_h(int elid/*ego_lane_id*/,
                          lane_car_data& lc,
                          double a0,
                          double v0,
                          double x0,
                          double y0,
			  double px,
			  double py,
			  double yaw0,
                          double s0,
                          double d0,
                          double lane_yaw,
                          std::function<std::vector<std::tuple<double, double>>(
                              std::vector<double>,
                              std::vector<double>)> xy) {
  car *ef(lc.front[elid]);
  lane_change_path ret;
  double lctime = 3.0;
  
  double target_v, target_s;
  bool goalAchieved;
  //  (spline,refx,refy,theta,jfn,total_time,total_distance,final_velocity,lane_change_time,target_lane_id,start_spline_x)
  typedef   std::tuple<tk::spline,double,double,double,std::function<double(double)>,double,double,double,double,int,double> option_type;
  option_type best_option;

  {
    // vector of tuple of option_type
    std::vector<option_type> options;
    std::map<int,std::string> rejection_reason;
    double max_time = 0;
    for(int tlid : {elid-1,elid,elid+1}) {
      if(tlid<0||tlid>2)
	continue;
      if(v0<5.0 && tlid!=elid) { // caution when removing this condition. 
	rejection_reason[tlid]="too slow to change lane";
	continue;
      }

      double target_d = (tlid + 0.5) * lane_width;
      car* b(lc.back[tlid]), * f(lc.front[tlid]);
      std::vector<double> spline_s_vals;
      std::vector<double> spline_d_vals;
      int lane_change_marker = std::numeric_limits<int>::max();
      {
	double vspln = std::max(v0,10.0);
	if(elid==tlid) {
	  //spline_s_vals.push_back(s0+vspln*0.25*lctime);
	  //spline_d_vals.push_back(target_d);
	  spline_s_vals.push_back(s0+vspln*0.5*lctime);
	  spline_d_vals.push_back(target_d);
	  spline_s_vals.push_back(s0+vspln*0.75*lctime);
	  spline_d_vals.push_back(target_d);
	}
	lane_change_marker = spline_s_vals.size();;

	spline_s_vals.push_back(s0+vspln*lctime);
	spline_d_vals.push_back(target_d);
	spline_s_vals.push_back(s0+vspln*1.25*lctime);
	spline_d_vals.push_back(target_d);
	spline_s_vals.push_back(s0+vspln*1.5*lctime);
	spline_d_vals.push_back(target_d);
      }
      tk::spline spln;
      double refx,refy,theta,start_spline_x;
      std::vector<double> spline_xs;
      std::tie(spln,refx,refy,theta,start_spline_x,spline_xs) = getSpline(x0,y0,px,py,yaw0,spline_s_vals,spline_d_vals,xy);
      std::vector<double> s_vals,a_vals,v_vals;
      std::function<double(double)> jfn;
      double total_time,total_dist,final_v,lane_change_time(0.0);
      if(tlid!=elid) {
	double lane_change_dist = length(spln,start_spline_x,spline_xs[lane_change_marker]);
	double changed_lane_time;
	double changed_lane_dist;
	std::function<double(double)> jfn1,jfn2;
	bool goalAchieved;
	double target_v;
	std::tie(jfn1,lane_change_time,goalAchieved) = achieveTargetVelocityAndDistanceInShortestTime
	  (a0,v0,v0,std::max(0.0,v0+a0*fabs(a0/(2*jmax))-0.1),std::min(max_speed,v0+a0*fabs(a0/(2*jmax))+0.1),lane_change_dist);
	double lane_changed_s = spline_s_vals[lane_change_marker];
	// should not collide with the vehicle in the back in the next lane
	if(b==nullptr || sdist(b->s,lane_changed_s)>safe_dist+(lane_change_time*b->v)) {
	  if(f==nullptr || (sdist(lane_changed_s,f->s)>safe_dist+lane_change_time*f->v+150.0)) {
	    target_v = max_speed;
	    std::tie(changed_lane_dist,jfn2,changed_lane_time) = distDuringVelocityChange(0,v0,target_v);
	  } else if (sdist(s0,f->s)>safe_dist) {
	    double delta_d = sdist(lane_changed_s,f->s+lane_change_time*f->v)-safe_dist;
	    target_v = f->v;
	    std::tie(jfn2,changed_lane_time,goalAchieved) = achieveTargetVelocityAndDistanceInShortestTime
	      (0,v0-target_v,0,-target_v,max_speed-target_v,delta_d);
	    changed_lane_dist = delta_d+target_v*changed_lane_time;
	  } else {
	    rejection_reason[tlid] = "car in the front is too close";
	    continue;
	  }
	} else {
	  rejection_reason[tlid]="car behind may hit you";
	  continue;
	}
	total_time = lane_change_time+changed_lane_time;
	total_dist = sdist(s0,lane_changed_s)+changed_lane_dist;
	jfn = [jfn1,jfn2,lane_change_time,total_time](double t) {
	  if(t<lane_change_time) return jfn1(t);
	  else if (t<total_time) return jfn2(t-lane_change_time);
	  else return 0.0;
	};
	final_v = target_v;
      } else {
	if(f) {
	  double delta_d = sdist(s0,f->s)-safe_dist;
	  bool goalAchieved;
	  std::tie(jfn,total_time,goalAchieved) = achieveTargetVelocityAndDistanceInShortestTime(a0,v0-f->v,0,-f->v,max_speed-f->v,delta_d);
	  final_v = f->v;
	  total_dist = delta_d+total_time*f->v;
	} else {
	  std::tie(total_dist,jfn,total_time) = distDuringVelocityChange(a0,v0,max_speed);
	  final_v = max_speed;
	}
      }
      if(total_time>max_time)
	max_time = total_time;
      options.push_back(std::make_tuple(spln,refx,refy,theta,jfn,total_time,total_dist,final_v,lane_change_time,tlid,start_spline_x));
    }
    double max_dist_at_max_time=0;
    for(option_type& option:options) {
      double t,d,v;
      std::tie(std::ignore,std::ignore,std::ignore,std::ignore,std::ignore,t,d,v,std::ignore,std::ignore,std::ignore) = option;
      double dist_at_max_time = d+(v*(max_time-t))+(std::get<9>(option)==elid?2*max_time:0.0);
      if(dist_at_max_time>max_dist_at_max_time) {
	max_dist_at_max_time = dist_at_max_time;
	best_option = option;
      }
    }
    {
      std::cout<<"num options : "<<options.size()<<std::endl;
      std::map<int,option_type*> option_map;
      for(option_type& x:options) {
	option_map[std::get<9>(x)] = &x;
      }
      for(int i=0;i<3;i++) {
	auto f(lc.front[i]),b(lc.back[i]);
	std::cout<<"   "<<i;
	std::cout<<" b : ";
	if(b) {
	  std::cout<<"("<<std::setw(4)<<int(sdist(s0,b->s))<<","<<std::setw(4)<<int(b->v)<<")  ";
	} else {
	  std::cout<<"(    ,    )  ";
	}
	std::cout<<" f : ";
	if(f) {
	  std::cout<<"("<<std::setw(4)<<int(sdist(s0,f->s))<<","<<std::setw(4)<<int(f->v)<<")  ";
	} else {
	  std::cout<<"(    ,    )  ";
	}
	std::cout<<(option_map[i]?"yes":" no")<<"  "<<rejection_reason[i]<<" ";
	if(elid==i)
	  std::cout<<"ego";
	std::cout<<std::endl;
      }
      std::cout<<std::flush;
    }
    {
      tk::spline s;
      std::function<double(double)> jfn;
      double total_time,total_distance,lane_change_time;
      double refx,refy,theta,start_spline_x;
      std::tie(s,refx,refy,theta,jfn,total_time,total_distance,std::ignore,lane_change_time,ret.lane_id,start_spline_x) = best_option;
      std::vector<double> s_vals;
      std::tie(s_vals,ret.a,ret.v) = genCompletePath(jfn,a0,v0,start_spline_x,std::max(lane_change_time,1.1));
      std::vector<double> x1s,y1s;
      std::tie(x1s,y1s) = discretizeSpline(s,s_vals);
      std::tie(ret.x,ret.y) = fromRefFrame(x1s,y1s,theta,refx,refy);
    }
  }
  return ret;
}
lane_change_path get_path(int elid/*ego_lane_id*/,
                          lane_car_data& lc,
                          double a0,
                          double v0,
                          double x0,
                          double y0,
			  double px,
			  double py,
			  double yaw0,
                          double s0,
                          double d0,
                          double lane_yaw,
                          std::function<std::vector<std::tuple<double, double>>(
                              std::vector<double>,
                              std::vector<double>)> xy) {
  static int call_id = 0;
  static bool first = true;
  static char fname[100];
  if(first) {
    std::sprintf(fname,"get_path_calls_and_responses_%d.csv",std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count() );
    first = false;
  }
  static std::ofstream fout(fname);
  lane_change_path lcp = get_path_h(elid,lc,a0,v0,x0,y0,px,py,yaw0,s0,d0,lane_yaw,xy);
  return lcp;
}

void update(double& s,double v) {
  if(s!=sentinel) {
    s+=v*dt;
  }
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
  static int call_id(0);
  static int point_id(0);
  static double ego_v, ego_a;
  static bool first_call(true);
  static int target_lane_id;
  static std::ofstream fout("rundata.csv");
  const char* header = "id,callid,end_path_s,end_path_d,x,y,v,a,np,f0s,f0v,b0s,b0v,f1s,f1v,b1s,b1v,f2s,f2v,b2s,b2v,target_lane_id";
  call_id++;
  std::vector<double> xs, ys;
  unsigned long np = prev_x.size();
  assert(prev_x.size() == prev_y.size());
  double ego_x(np ? prev_x.back() : car_x), ego_y(np ? prev_y.back() : car_y);
  double ego_px(np>1?prev_x[np-2]:sentinel),ego_py(np>1?prev_y[np-2]:sentinel);
  double ego_s(np ? end_path_s : car_s), ego_d(np ? end_path_d : car_d);
  double ego_yaw(np > 1 ? yaw(prev_x[np - 2],
                              prev_y[np - 2],
                              prev_x[np - 1],
                              prev_y[np - 1]) : car_yaw);
  int ego_lane_id = lane_id(ego_d);
  std::cout<<"ego_v : "<<ego_v<<" ego_a : "<<ego_a<<std::endl;
  if (first_call) {
    ego_v = car_speed;
    ego_a = 0.0;
    target_lane_id = ego_lane_id;
    fout<<header<<std::endl;
  }
  xs = prev_x;
  ys = prev_y;
  if(np<50) {
    std::vector<car> transformedCarsData =
      transformCarsData(cars_data, np * dt);
    lane_car_data
      lcdata = possibleLanes(car(ego_s, ego_d, ego_v), transformedCarsData);
    lane_change_path change_path;
    change_path = get_path(ego_lane_id,
			   lcdata,
			   ego_a,
			   ego_v,
			   ego_x,
			   ego_y,
			   ego_px,
			   ego_py,
			   ego_yaw,
			   ego_s,
			   ego_d,
			   lane_yaw,
			   xy);
    auto f0(lcdata.front[0]),f1(lcdata.front[1]),f2(lcdata.front[2]);
    auto b0(lcdata.back[0]),b1(lcdata.back[1]),b2(lcdata.back[2]);
    double f0s(f0?f0->s:sentinel),f0v(f0?f0->v:sentinel),b0s(b0?b0->s:sentinel),b0v(b0?b0->v:sentinel),
      f1s(f1?f1->s:sentinel),f1v(f1?f1->v:sentinel),b1s(b1?b1->s:sentinel),b1v(b1?b1->v:sentinel),
      f2s(f2?f2->s:sentinel),f2v(f2?f2->v:sentinel),b2s(b2?b2->s:sentinel),b2v(b2?b2->v:sentinel);
    for (int i = 0;i<change_path.x.size(); i++) {
      xs.push_back(change_path.x[i]);
      ys.push_back(change_path.y[i]);
      ego_a = change_path.a[i];
      ego_v = change_path.v[i];
      point_id++;
      fout<<point_id<<","<<call_id<<","<<end_path_s<<","<<end_path_d<<","<<xs.back()<<","<<ys.back()<<","<<ego_v<<","<<ego_a<<","<<np;
      update(f0s,f0v);
      update(b0s,b0v);
      update(f1s,f1v);
      update(b1s,b1v);
      update(f2s,f2v);
      update(b2s,b2v);
      for(double x : {f0s,f0v,b0s,b0v,f1s,f1v,b1s,b1v,f2s,f2v,b2s,b2v}) {
	fout<<",";
	if(x!=sentinel)
	  fout<<x;
      }
      fout <<","<<change_path.lane_id<<std::endl;
    }
    fout<<std::flush;
  }
  if(fabs(ego_a)>amax*1.05 || fabs(ego_v)>max_speed*1.05) {
    std::cout<<"error : ego_a : "<<ego_a<<" ego_v : "<<ego_v<<std::endl<<" num_points : "<<xs.size()<<" prev_size : "<<np<<std::endl ;
  }
  first_call = false;
  return std::make_pair(xs, ys);
}
