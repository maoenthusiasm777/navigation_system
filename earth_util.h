#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace Eigen;

Vector3d LLHToXYZ(double latitude,double longitude,double height);

Matrix3d LLHToXYZDiff(double latitude,double longitude,double height);

Matrix3d EcefToNed(double latitude,double longitude);

Matrix3d NedToEcef (double latitude,double longitude);
