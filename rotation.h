#include <Eigen/Core>
#include <Eigen/Geometry>

Matrix3d CrossMat(const Vector3d& vec);

Matrix3d EulerAngleToRotion (const Vector3d& euler_angle);

Vector3d DcmToEulerAngle (const Matrix3d& dcm);

Quaterniond DcmToQuater (const Matrix3d& dcm);

Matrix3d QuaToDcm(const Quaterniond& q);