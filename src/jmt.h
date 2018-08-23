//
// Created by sunilsn on 8/22/18.
//

#ifndef PATH_PLANNING_JMT_H
#define PATH_PLANNING_JMT_H
//path points,final acceleration,final velocity
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, bool>
path(double a0, double v0, double v1, double delta_d);

// dist,jfn,total_time
std::tuple<double, std::function<double(double)>, double>
distDuringVelocityChange(double a0, double v0, double v1);

//jfn,total_time,goalAchieved
std::tuple<std::function<double(double)>, double, bool>
achieveTargetVelocityAndDistanceInShortestTime(double a0,
                                               double v0,
                                               double v1,
                                               double vmin,
                                               double vmax,
                                               double delta_d);

// s_vals,a_vals,v_vals
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
genCompletePath(std::function<double(double)> jfn,
                double a0, double v0, double s0, double time);
#endif //PATH_PLANNING_JMT_H
