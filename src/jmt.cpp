#include <iostream>
#include <functional>
#include <cmath>
#include <fstream>
#include <vector>

#include "constants.h"

double secant(const std::function<double(double)>& f,double min,double max,double tol,const char* call_name) {
  double x0(min),x1(max);
  double f0(f(x0)),f1(f(x1));
  while(true) {
    //std::cout<<" x0 : "<<x0<<" f0 : "<<f0<<" x1 : "<<x1<<" f1 : "<<f1<<std::endl;
    double xs = (f0*x1-f1*x0)/(f0-f1);
    double fs = f(xs);
    if(fabs(fs)<tol) {
      std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return xs "<<std::endl;
      return xs;
    } else {
      if(f0*fs<0) {
        f1 = fs;
        x1 = xs;
      } else if (f1*fs<0) {
        f0 = fs;
        x0 = xs;
      } else {
        if(fabs(f0)<fabs(f1)) {
	  std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return x0 "<<std::endl;
	  return x0;
        } else {
	  std::cout<<call_name<<" xs "<<xs<<" fs "<<fs<<" x0 "<<x0<<" f0 "<<f0<<" x1 "<<x1<<" f1 "<<f1<<" return x1 "<<std::endl;
          return x1;
	}
      }
    }
  }
}

// assuming maximum jerk is constant
// minimum velocity change when changing acceleration
double deltaVelocity(double a0,double a1) {
  return (a0+a1)*0.5*fabs(a1-a0)/jmax;
}

// acceleration and velocity are initial values
double dist(double v,double a,double j,double t) {
  return v*t+a*t*t*0.5+j*t*t*t/6.0;
}

// returns <distChange,velChange,jerk_fn,time>
std::tuple<double,double,std::function<double(double)>,double> distDuringAccelChangeH(double v0,double a0,double a1) {
  double j = (a1>a0?1:-1)*jmax;
  double t1 = (a1-a0)/j;
  double d = dist(v0,a0,j,t1);
  double dv = deltaVelocity(a0,a1);
  auto jfn = [t1,j](double t) { if(t>=0 && t<=t1) return j;};
  return std::make_tuple(d,dv,jfn,t1);
}

std::tuple<double,double,std::function<double(double)>,double> distDuringAccelChange(double v0,double a0,double a1) {
      auto x = distDuringAccelChangeH(v0,a0,a1);
      static char buffer[1000];
      sprintf(buffer,"distDuringAccelChange : v0 %f a0 %f a1 %f returns dist %f velChange %f deltaT %f",v0,a0,a1,std::get<0>(x),
	      std::get<1>(x),std::get<3>(x));
      std::cout<<buffer<<std::endl;
      return x;
}
    

// dist,jfn,totaltime
std::tuple<double,std::function<double(double)>,double> distDuringVelocityChangeH(double a0,double v0,double v1) {
  double dv = v1-v0;
  // change in velocities
  auto direct = distDuringAccelChange(v0,a0,0);
  auto toAmax = distDuringAccelChange(v0,a0,amax);
  auto fromAmax = distDuringAccelChange(v0+std::get<1>(toAmax),amax,0);
  auto toAmin = distDuringAccelChange(v0,a0,amin);
  auto fromAmin = distDuringAccelChange(v0+std::get<1>(toAmin),amin,0);
  auto h = [v0,v1](std::tuple<double,double,std::function<double(double)>,double> to,
						       std::tuple<double,double,std::function<double(double)>,double> from,
											      double apeak) {
    double t1(std::get<3>(to)),t3(std::get<3>(from)),d1(std::get<0>(to)),d3(std::get<0>(from));
    double u1(v0+std::get<1>(to)),u3(v1-std::get<1>(from));
    double t2 = (u3-u1)/apeak;
    double d2 = (u1+u3)*0.5*t2;
    auto jfn1 = std::get<2>(to);
    auto jfn3 = std::get<2>(from);
    auto jfn = [jfn1,jfn3,t1,t2,t3](double t) {
      if(t>=0 && t<t1) {
	return jfn1(t);
      } else if (t<=t1+t2) {
	return 0.0;
      } else if (t<=t1+t2+t3) {
	return jfn3(t-t1-t2);
      }
    };
    if(t2>=0) {
      return std::make_tuple(d1+d2+d3,jfn,t1+t2+t3);
    } else {
      throw("error");
    }
  };
  if(fabs(dv-std::get<1>(direct))<1e-6) {
    return std::make_tuple(std::get<0>(direct),std::get<2>(direct),std::get<3>(direct));
  } else if (std::get<1>(toAmax)+std::get<1>(fromAmax)<dv) {
    return h(toAmax,fromAmax,amax);
  } else if (std::get<1>(toAmin)+std::get<1>(fromAmin)>dv) {
    return h(toAmin,fromAmin,amin);
  } else {
    auto f = [a0,dv](double apeak) {
      return deltaVelocity(a0,apeak)+deltaVelocity(apeak,0)-dv;
    };
    static char buffer[1000];
    sprintf(buffer,"apeak : a0 %f v0 %f ",a0,v0);
    double apeak = secant(f,amin,amax,1e-6,buffer);
    auto toApeak = distDuringAccelChange(v0,a0,apeak);
    auto fromApeak = distDuringAccelChange(v0+std::get<1>(toApeak),apeak,0);
    return h(toApeak,fromApeak,apeak);
  }
}
// dist,jfn,totaltime
std::tuple<double,std::function<double(double)>,double> distDuringVelocityChange(double a0,double v0,double v1) {
      auto x = distDuringVelocityChangeH(a0,v0,v1);
      static char buffer[1000];
      sprintf(buffer,"distDuringVelocityChange : a0 %f v0 %f v1 %f returns dist %f totaltime %f",a0,v0,v1,std::get<0>(x),std::get<2>(x));
      std::cout<<buffer<<std::endl;
      return x;
    }
std::tuple<std::function<double(double)>,double> achieveZeroAcelAndVel(double a0,double v0,double vmin,double vmax,double delta_d) {
  auto direct = distDuringVelocityChange(a0,v0,0);
  auto toVmin = distDuringVelocityChange(a0,v0,vmin);
  auto fromVmin = distDuringVelocityChange(0,vmin,0);
  auto toVmax = distDuringVelocityChange(a0,v0,vmax);
  auto fromVmax = distDuringVelocityChange(0,vmax,0);
  double directDist = std::get<0>(direct);
  double vminDist = std::get<0>(toVmin)+std::get<0>(fromVmin);
  double vmaxDist = std::get<0>(toVmax)+std::get<0>(fromVmax);
  auto h = [delta_d](std::tuple<double,std::function<double(double)>,double> to,std::tuple<double,std::function<double(double)>,double> from,double vpeak) {
    double t3 = std::get<2>(from);
    double d3 = std::get<0>(from);
    double t1 = std::get<2>(to);
    double d1 = std::get<0>(to);
    double d2 = delta_d-d1-d3;
    double t2 = d2/vpeak;
    auto jfn3 = std::get<1>(from);
    auto jfn1 = std::get<1>(to);
    auto jfn = [t1,t2,t3,jfn1,jfn3](double t) {
      if(t>=0 && t < t1) {
	return jfn1(t);
      } else if (t<=t1+t2) {
	return 0.0;
      } else if (t<=t1+t2+t3) {
	return jfn3(t-t1-t2);
      }
    };
    return std::make_tuple(jfn,t1+t2+t3);
  };
  if(fabs(delta_d-directDist) < 1e-5 || (delta_d<vminDist && !(vmin<0)) || (delta_d>vmaxDist && !(vmax>0))) {
    return std::make_tuple(std::get<1>(direct),std::get<2>(direct));
  } else if(delta_d < vminDist) {
    if(vmin<0) {
      return h(toVmin,fromVmin,vmin);
    } else {
      throw("should not come here");
    }
  } else if(delta_d > vmaxDist) {
    if(vmax>0) {
      return h(toVmax,fromVmax,vmax);
    } else {
      throw("should not come here");
    }    
  } else {
    auto f = [a0,v0,delta_d](double vpeak) {
      return std::get<0>(distDuringVelocityChange(a0,v0,vpeak))+std::get<0>(distDuringVelocityChange(0,vpeak,0))-delta_d;
    };
    static char buffer[1000];
    sprintf(buffer,"peakVel a0 %f v0 %f delta_d %f",a0,v0,delta_d);
    double vpeak = secant(f,vmin,vmax,1e-5,buffer);
    auto toVpeak = distDuringVelocityChange(a0,v0,vpeak);
    auto fromVpeak = distDuringVelocityChange(0,vpeak,0);
    return h(toVpeak,fromVpeak,vpeak);
  }
}

std::vector<double> path(double a0,double v0,double v1,double delta_d) {
  double totalTimeEstimate = delta_d/((v0+v1)*0.5);
  std::function<double(double)> jfn;
  double timeTaken;
  do {
    std::tie(jfn,timeTaken) = achieveZeroAcelAndVel(a0,v0-v1,-v1,max_speed-v1,delta_d-v1*totalTimeEstimate);
    totalTimeEstimate = (timeTaken+totalTimeEstimate)*0.5;
  } while(fabs(timeTaken-totalTimeEstimate)*v1 > max_speed*dt);

  double ap(a0),vp(v0),sp(0);
  std::vector<double> ret;
  for(double t = dt; t<=totalTimeEstimate;t+=dt) {
    double ac,vc,sc,jc;
    jc = jfn(t);
    ac = ap+jc*dt;
    vc = vp+(ap+ac)*dt*0.5;
    sc = sp+(vp+vc)*dt*0.5;
    ret.push_back(sc);
    ap = ac;
    vp = vc;
    sp = sc;
  }
  return ret;
}

void path(double a0,double v0,double vmin,double vmax,double delta_d,const char* fname) {
  auto jfn_t = achieveZeroAcelAndVel(a0,v0,vmin,vmax,delta_d);
  auto jfn = std::get<0>(jfn_t);
  auto total_t = std::get<1>(jfn_t);
  std::cout<<"fname : "<<fname<<" total_t : "<<total_t<<std::endl;
  double ap(a0),vp(v0),sp(0);
  double dt = 0.02;
  std::ofstream fout(fname);
  fout<<"t,j,a,v,s"<<std::endl;
  fout<<"0,0,"<<a0<<","<<v0<<",0"<<std::endl;
  for(double t = dt; t<=total_t;t+=dt) {
    double ac,vc,sc,jc;
    jc = jfn(t);
    ac = ap+jc*dt;
    vc = vp+(ap+ac)*dt*0.5;
    sc = sp+(vp+vc)*dt*0.5;
    fout<<t<<","<<jc<<","<<ac<<","<<vc<<","<<sc<<std::endl;
    ap = ac;
    vp = vc;
    sp = sc;
  }
}

int main_jmt() {

  double a0(0),dv(5),tol(1e-7);
  //void path(double a0,double v0,double vmin,double vmax,double delta_d,const char* fname) 
  path(0.0,0.0, 0,23, 100,"trial1.csv");
  path(5.0,0.0,0,23,100,"trial2.csv");
  path(10.0,2.0,0,23,100,"trial3.csv");

  path(10.0,10.0,0,23,100,"trial4.csv");
  path(0.0,23.0,0,23,50,"trial5.csv");
  path(0.0,0.0,-10,10,10,"trial6.csv");
  path(0.0,0.0,-10,10,-10,"trial7.csv");
  path(10,5,-10,10,-20,"trial8.csv");
  path(-10,-5,-10,10,20,"trial9.csv");
  
  
}
