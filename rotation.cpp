#include "rotation.h"

Matrix3d CrossMat(const Vector3d& vec) {
    Matrix3d cross_vec;
    cross_vec << 0,-vec.z(),vec.y(),
                 vec.z(),0,-vec.x(),
                 -vec.y(),vec.x(),0;
    return cross_vec;
}

Matrix3d EulerAngleToRotion (const Vector3d& euler_angle) {
    Matrix3d a;
    a<< 1,0,0,
        0,cos(euler_angle(2)*M_PI/180),sin(euler_angle(2)*M_PI/180),
        0,-sin(euler_angle(2)*M_PI/180),cos(euler_angle(2)*M_PI/180);
    Matrix3d b;
    b<< cos(euler_angle(1)*M_PI/180),0,-sin(euler_angle(1)*M_PI/180),
        0,1,0,
        sin(euler_angle(1)*M_PI/180),0,cos(euler_angle(1)*M_PI/180);
    Matrix3d c;
    c << cos(euler_angle(0)*M_PI/180),sin(euler_angle(0)*M_PI/180),0,
        -sin(euler_angle(0)*M_PI/180),cos(euler_angle(0)*M_PI/180),0,
        0,0,1;
    return a*b*c;    
}

Vector3d DcmToEulerAngle (const Matrix3d& dcm) { // cnb
    Vector3d angle;
    angle << atan2(dcm(1,2),dcm(2,2)) * 180.0 /M_PI,-asin(dcm(0,2))* 180.0 /M_PI,atan2(dcm(0,1),dcm(0,0))* 180.0 /M_PI;
    return angle; // roll ,pitch ,yaw.
}

Quaterniond DcmToQuater (const Matrix3d& dcm) {
    double q0 = 0.5*std::sqrt(1 + dcm(0,0)+dcm(1,1) + dcm(2,2));
    return Quaterniond(q0,
                       0.25/q0 * (dcm(2,1) - dcm(1,2)),
                       0.25/q0 * (dcm(0,2) - dcm(2,0)),
                       0.25/q0 * (dcm(1,0) - dcm(0,1)));
}

Matrix3d QuaToDcm(const Quaterniond& q) {
    Matrix3d dcm;
    dcm << q.w() * q.w() + q.x() * q.x() - q.y() * q.y() - q.z() * q.z(),2*(q.x()*q.y()-q.w()*q.z()),2*(q.x()*q.z() + q.w()*q.y()),
          2*(q.x()*q.y()+q.w()*q.z()),q.w() * q.w() - q.x() * q.x() + q.y() * q.y() - q.z() * q.z(),2*(q.y()*q.z() - q.w()*q.x()),
          2*(q.x()*q.z() - q.w()*q.y()),2*(q.y()*q.z() + q.w()*q.x()),q.w() * q.w() - q.x() * q.x() - q.y() * q.y() + q.z() * q.z();
}
