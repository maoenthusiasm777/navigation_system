#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <algorithm>

const double omega = 7.292115 * 1E-5;   // rad / s
const double R = 6378137.0;             // meters
const double e = 0.0818191908426;
const double Lat0= 0.492538083125361;   // rad 
const double Long0 = 1.972094934421829; // rad
const double H0 = 50.385756076313555;
const double RN = R*(1-e*e)/std::pow(1-e*e*(sin(Lat0))*(sin(Lat0)),3/2); // meters
const double RE = R/std::sqrt((1-e*e*(sin(Lat0))*(sin(Lat0))));   // meters
const double g=9.7803267711905*(1+0.00193185138639*(sin(Lat0))*(sin(Lat0)))
               / sqrt(1-0.00669437999031*(sin(Lat0))*(sin(Lat0))) / (1 + H0 / R) / (1 + H0 / R); // m * s^-2
const double c = 2.99792458e8;
//  NED coordinate system.
extern Eigen::Vector3d gn;
