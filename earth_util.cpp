#include <Eigen/Core>
#include <Eigen/Geometry>

#include "constant.h"
#include "earth_util.h"

Vector3d LLHToXYZ(double latitude,double longitude,double height) {
  Vector3d pos_xyz;
  const double re = R/std::sqrt((1-e*e*(sin(latitude))*(sin(latitude))));
  pos_xyz[0] = (re+height)*cos(latitude)*cos(longitude);
  pos_xyz[1] = (re+height)*cos(latitude)*sin(longitude);
  pos_xyz[2] = ((1-e * e)*re+height)*sin(latitude);
  return pos_xyz;
}

Matrix3d LLHToXYZDiff(double latitude,double longitude,double height) { 
  Matrix3d transform;
  transform << 
  -(RE+height)*sin(latitude)*cos(longitude),-(RE+height)*cos(latitude)*sin(longitude),cos(latitude)*cos(longitude),
  -(RE+height)*sin(latitude)*sin(longitude),(RE+height)*cos(latitude)*cos(longitude),cos(latitude)*sin(longitude),
  ((1-e * e)*RE+height)*cos(latitude),0,sin(latitude);
  return transform;
}

Matrix3d EcefToNed(double latitude,double longitude) {
  Matrix3d transform;
  transform << 
  -sin(latitude)*cos(longitude),-sin(latitude)*sin(longitude),cos(latitude),
  -sin(longitude),cos(longitude),0,
  -cos(latitude)*cos(longitude),-cos(latitude)*sin(longitude),-sin(latitude);
  return transform;
}

Matrix3d NedToEcef (double latitude,double longitude) {
  const Matrix3d transform = EcefToNed(latitude,longitude);
  return transform.transpose();
}
