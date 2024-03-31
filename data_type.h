#include <Eigen/Core>
#include <Eigen/Geometry>

struct IMURawData {
    double time;
    double frame_number;
    double angle_increment_x;
    double angle_increment_y;
    double angle_increment_z;
    double velocity_increment_x;
    double velocity_increment_y;
    double velocity_increment_z;
    double odo_vel;
};

struct GNSSResultData {
    double lat;     // rad
    double longitude;
    double height;  // m
    double vn;
    double ve;
    double vd;
};

class NavigationState {       // as a node in optimization
  public: 
    NavigationState(const Matrix3d& cbn0, double lat0, double longitude0,double h0): 
                    lat(lat0), longitude(longitude0), h(h0) { cbn = cbn0; }
    void ResetVelocity() {vn = 0.0;ve = 0.0; vd = 0.0;}
    void ResetPose() {xn.setZero();}
    void SetAngle(const Matrix3d& rotaion) {
        cbn = rotaion;
        yaw = cbn.eulerAngles(2,1,0)[0] * 180 / M_PI;
        pitch = cbn.eulerAngles(2,1,0)[1]* 180 / M_PI;
        roll = cbn.eulerAngles(2,1,0)[2]* 180 / M_PI;
        qbn = Quaterniond(cbn);
        qbn.normalized();
    }
    void SetAngle (const Quaterniond& rotaion) {
       qbn = rotaion;
       qbn.normalized();
       cbn = qbn.toRotationMatrix();
       yaw = cbn.transpose().eulerAngles(2,1,0)[0]* 180 / M_PI;
       pitch = cbn.transpose().eulerAngles(2,1,0)[1]* 180 / M_PI;
       roll = cbn.transpose().eulerAngles(2,1,0)[2]* 180 / M_PI;
    }
    double roll;
    double pitch;
    double yaw;

    double vn = 0.0;
    double ve = 0.0;
    double vd = 0.0;

    double lat;
    double longitude;
    double h;

    Matrix3d cbn;
    Vector3d xn;  //meters
    Quaterniond qbn;
};