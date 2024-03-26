#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iomanip>

using namespace Eigen;

static double omega = 7.292115 * 1E-5;   // rad / s
static double R = 6378137.0;             // meters
static double e = 0.0818191908426;
static double Lat0= 0.492538083125361;   // rad 
static double Long0 = 1.972094934421829; // rad
static double H0 = 50.385756076313555;
static double RN = R*(1-e*e)/std::pow(1-e*e*(sin(Lat0))*(sin(Lat0)),3/2); // meters
static double RE = R/std::sqrt((1-e*e*(sin(Lat0))*(sin(Lat0))));   // meters
static double g=9.7803267711905*(1+0.00193185138639*(sin(Lat0))*(sin(Lat0)))
               / sqrt(1-0.00669437999031*(sin(Lat0))*(sin(Lat0))) / (1 + H0 / R) / (1 + H0 / R); // m * s^-2

//  ENU coordinate system.
Vector3d gn(0,0,g);
// From odometry frame to body system.
Matrix3d cmb;

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

    double vn;
    double ve;
    double vd;

    double lat;
    double longitude;
    double h;

    Matrix3d cbn;
    Vector3d xn;  //meters
    Quaterniond qbn;
};

class FliterState {
    public:
      void SetInitialP() {
          P(0,0) = 0.0007 * 0.0007;
          P(1,1) = 0.0007 * 0.0007;
          P(2,2) = 0.0007 * 0.0007;
          P(3,3) = 0.01 * 0.01;
          P(4,4) = 0.01 * 0.01;
          P(5,5) = 0.01 * 0.01;
          P(6,6) = 1.5749e-6 * 1.5749e-6;
          P(7,7) = 1.5749e-6 * 1.5749e-6;
          P(8,8) = 10 * 10;
          P(9,9) = 2.909e-6 * 2.909e-6;
          P(10,10) = 2.909e-6 * 2.909e-6;
          P(11,11) = 2.909e-6 * 2.909e-6;
          P(12,12) = 0.0005 * 0.0005;
          P(13,13) = 0.0005 * 0.0005;
          P(14,14) = 0.0005 * 0.0005;
      }
      void reset() {x.setZero();}
      Matrix<double,15,1> x = MatrixXd::Zero(15,1);
      Matrix<double,15,15> P; // initial val
};

class FliterState9 {
    public:
      void SetInitialP() {
          P(0,0) = 0.0007 * 0.0007;
          P(1,1) = 0.0007 * 0.0007;
          P(2,2) = 0.0007 * 0.0007;
          P(3,3) = 0.01 * 0.01;
          P(4,4) = 0.01 * 0.01;
          P(5,5) = 0.01 * 0.01;
          P(6,6) = 8.5749e-6 * 8.5749e-6;
          P(7,7) = 4.5749e-6 * 4.5749e-6;
          P(8,8) = 10 * 10;
      }
      void reset() {x.setZero();}
      Matrix<double,9,1> x = MatrixXd::Zero(9,1);
      Matrix<double,9,9> P; // initial val
};

Matrix3d CoarseAlignment (const std::vector<IMURawData>& imu_data) {
    const double t = 120.0;
    double w_b_ie_sum_x = 0.0;
    double w_b_ie_sum_y = 0.0;
    double w_b_ie_sum_z = 0.0;
    double f_b_sum_x = 0.0;
    double f_b_sum_y = 0.0;
    double f_b_sum_z = 0.0;
    for (int i = 0;i < imu_data.size();++i) {
      w_b_ie_sum_x += imu_data[i].angle_increment_x;
      w_b_ie_sum_y += imu_data[i].angle_increment_y;
      w_b_ie_sum_z += imu_data[i].angle_increment_z;
      f_b_sum_x += imu_data[i].velocity_increment_x;
      f_b_sum_y += imu_data[i].velocity_increment_y;
      f_b_sum_z += imu_data[i].velocity_increment_z;
    }
    Vector3d w_b_ie(w_b_ie_sum_x / t,w_b_ie_sum_z / t,-w_b_ie_sum_y/t);
    Vector3d f_b(f_b_sum_x / t,f_b_sum_z / t,-f_b_sum_y /t);
    Vector3d w_f = CrossMat(w_b_ie) * f_b;
    Matrix3d mat_a;
    mat_a.col(0) << w_b_ie;
    mat_a.col(1) << -f_b;
    mat_a.col(2) << -w_f;

    Matrix3d mat_b;
    Vector3d wie_n;
    wie_n << omega*cos(Lat0),0,-omega*sin(Lat0);
    mat_b.col(0) << wie_n;
    mat_b.col(1) << gn;
    mat_b.col(2) << (CrossMat(wie_n)*gn);
    Matrix3d cnb;
    cnb = mat_a*mat_b.inverse();
    
    return cnb.transpose();
}

void SINSCalculate(const std::vector<IMURawData>& imu_data,NavigationState* state) {
    const double t = 0.02; 
    Vector3d v_old(state->vn,state->ve,state->vd);
    const Matrix3d cbn_old = state->cbn;

    Vector3d wie_n;
    wie_n << omega*cos(state->lat),0,-omega*sin(state->lat);
    Vector3d w_en_n;
    w_en_n << state->ve/RE,-state->vn/RN,-state->ve*tan(state->lat)/RE;
    Vector3d w_in_n;
    w_in_n = wie_n + w_en_n;
    Vector3d w_in_b =  state->cbn.transpose() * w_in_n;
    Vector3d det_theta1(imu_data[0].angle_increment_x,imu_data[0].angle_increment_y,imu_data[0].angle_increment_z);
    det_theta1 -= w_in_b * t / 2;
    Vector3d det_theta2(imu_data[1].angle_increment_x,imu_data[1].angle_increment_y,imu_data[1].angle_increment_z);
    det_theta2 -= w_in_b * t / 2;
    Vector3d det_theta = det_theta1 + det_theta2;
    
    const Vector3d det_vk1(imu_data[0].velocity_increment_x,imu_data[0].velocity_increment_y,imu_data[0].velocity_increment_z);
    const Vector3d det_vk2(imu_data[1].velocity_increment_x,imu_data[1].velocity_increment_y,imu_data[1].velocity_increment_z);  
    const Vector3d det_vk = det_vk1+ det_vk2;

    const Vector3d det_rot = 0.5 * det_theta.cross(det_vk);
    const Vector3d det_vscul = 2 / 3.0 * (det_theta1.cross(det_vk2) + det_vk1.cross(det_theta2));
    const Vector3d det_vsfk = cbn_old * (det_vk + det_rot + det_vscul);
    
    const Vector3d det_vgcork = t * (gn - (2 * wie_n + w_en_n).cross(v_old));
    //update velocity
    Vector3d v_new = v_old + det_vsfk + det_vgcork;
    state->vn = v_new(0);
    state->ve = v_new(1);
    state->vd = 0.0;

    //update pose
    state->lat += 0.5 * (v_old(0) + v_new(0)) / RN * t;
    state->longitude += 0.5 * (v_old(1) + v_new(1)) / (RE*cos(state->lat))*t;
    state->xn += 0.5*(v_old+v_new) * t;

    //update rotation
    wie_n << omega*cos(state->lat),0,-omega*sin(state->lat);
    w_en_n << state->ve / RE,-state->vn / RN,-state->ve*tan(state->lat) / RE;
    w_in_b = state->cbn.transpose() * (wie_n + w_en_n);
    det_theta1 = Vector3d(imu_data[0].angle_increment_x,imu_data[0].angle_increment_y,imu_data[0].angle_increment_z)
                  - w_in_b * t / 2;
    det_theta2 = Vector3d(imu_data[1].angle_increment_x,imu_data[1].angle_increment_y,imu_data[1].angle_increment_z)
                  - w_in_b * t / 2;
    Vector3d rot_vec =  det_theta1 + det_theta2 + 2 / 3.0 * det_theta1.cross(det_theta2);
    AngleAxisd angle_axis(rot_vec.norm(),rot_vec / rot_vec.norm());
    Matrix3d mat = state->cbn * angle_axis.toRotationMatrix();
    mat.normalized();
    state -> cbn = mat;
}

void DR(const std::vector<IMURawData>& imu_data,NavigationState* state) {
    const double t = 0.02; 
    Vector3d v_odo_1(imu_data[0].odo_vel,0,0);
    Vector3d v_odo_2(imu_data[1].odo_vel,0,0);
    Vector3d vn1 = state->cbn * cmb * v_odo_1;
    Vector3d vn2 = state->cbn * cmb * v_odo_2; 
    Vector3d v_new = (vn1 + vn2) / 2;
    //update velocity
    state->vn = v_new(0);
    state->ve = v_new(1);
    state->vd = v_new(2);
    //update position
    state->lat +=  v_new(0) / RN * t;
    state->longitude += v_new(1) / (RE*cos(state->lat))*t;
    state->xn += v_new * t;

    //update rotation
    Vector3d wie_n(omega*cos(state->lat),0,-omega*sin(state->lat));
    Vector3d w_en_n(state->ve / RE,-state->vn / RN,-state->ve*tan(state->lat) / RE);
    Vector3d w_in_b = state->cbn.transpose() * (wie_n + w_en_n);
    Vector3d det_theta1 = Vector3d(imu_data[0].angle_increment_x,imu_data[0].angle_increment_y,imu_data[0].angle_increment_z)
                  - w_in_b * t / 2;
    Vector3d det_theta2 = Vector3d(imu_data[1].angle_increment_x,imu_data[1].angle_increment_y,imu_data[1].angle_increment_z)
                  - w_in_b * t / 2;
    Vector3d rot_vec =  det_theta1 + det_theta2 + 2 / 3.0 * det_theta1.cross(det_theta2);
    AngleAxisd angle_axis(rot_vec.norm(),rot_vec / rot_vec.norm());
    // Quaterniond q = state -> qbn * Quaterniond(angle_axis);
    // q = q.normalized().coeffs();
    // state->qbn = q;
    //std::cout<< "rot vec" <<rot_vec.norm() << std::endl;
     Matrix3d mat = state->cbn * angle_axis.toRotationMatrix();
     mat.normalized();
     state -> cbn = mat;
    //state -> SetAngle(q);     
}



void STExtendKalmanFilter(const std::vector<IMURawData>& imu_data,
                          const GNSSResultData& gnss_data,
                          NavigationState *state,
                          FliterState *error_state) {
    
    const double t = 0.02;
    SINSCalculate(imu_data,state);
    const Vector3d det_vk1(imu_data[0].velocity_increment_x,imu_data[0].velocity_increment_y,imu_data[0].velocity_increment_z);
    const Vector3d det_vk2(imu_data[1].velocity_increment_x,imu_data[1].velocity_increment_y,imu_data[1].velocity_increment_z); 
    const Vector3d fn = state->cbn*((det_vk1 + det_vk2) / t);
    const Matrix3d cross_vn = CrossMat(Vector3d(state->vn,state->ve,state->vd));
    Matrix<double,3,3> A12;
    A12 << 0,1.0/(RE+state->h),0,
           - 1.0 / (RN+state->h),0,0,
           0,-tan(state->lat) / (RE+state->h),0;

    Matrix<double,3,3> A11;
    A11 << 0, -(omega * sin(state->lat) + state->ve / (RE+state->h) * tan(state->lat)),state->vn / (RN+state->h),
    (omega * sin(state->lat) + state->ve / (RE+state->h) * tan(state->lat)),0,(omega * cos(state->lat) + state->ve / (RE+state->h)),
    - state->vn / (RN+state->h), -(omega * cos(state->lat) + state->ve / (RE+state->h)),0;
    A11 += A12 * cross_vn;
    
    Matrix<double,3,3> A13;
    A13 << -omega * sin(state->lat),0, -state->ve / (RE+state->h) / (RE+state->h),
           0,0,state->vn / (RN+state->h) / (RN+state->h),
           -(omega * cos(state->lat) + state->ve / (RE+state->h)) /cos(state->lat) /cos(state->lat),0,state->ve * tan(state->lat) / (RE+state->h) / (RE+state->h);
    
    Matrix<double,3,3> A21;
    A21 = -(CrossMat(gn) + cross_vn * CrossMat(Vector3d(omega*cos(state->lat),0,-omega*sin(state->lat))));

    Matrix<double,3,3> A22;
    A22 = -(2.0 * CrossMat(Vector3d(omega*cos(state->lat),0,-omega*sin(state->lat))) + CrossMat(Vector3d(state->ve / (RE+state->h),-state->vn / (RN+state->h),-state->ve * tan(state->lat) / (RE+state->h))));

    Matrix<double,3,3>  A23;
    A23 << -omega * sin(state->lat),0,0,
            0,0,0,
            -omega * cos(state->lat),0,0;
    A23 = cross_vn * A23;

    Matrix<double,3,3>  A32;
    A32 << 1.0 / (RN+state->h),0,0,
           0,1.0/(RE+state->h) / cos(state->lat),0,
           0,0,-1.0;
    Matrix<double,3,3>  A33;
    A33 << 0,0,-state->vn / (RN+state->h) / (RN+state->h),
           state->ve * tan(state->lat) / (RE+state->h) / cos(state->lat), 0 ,-state->ve / (RE+state->h) / (RE+state->h) / cos(state->lat),
           0,0,0;
    Matrix<double,3,3>  A31;
    A31 = A32 * cross_vn;

    Matrix<double,15,15> A = MatrixXd::Zero(15,15);
    A.block(0,0,3,3) = A11;
    A.block(0,3,3,3) = A12;
    A.block(0,6,3,3) = A13;
    A.block(0,9,3,3) = -state->cbn;

    A.block(3,0,3,3) = A21;
    A.block(3,3,3,3) = A22;
    A.block(3,6,3,3) = A23;
    A.block(3,9,3,3) = cross_vn * state->cbn;
    A.block(3,12,3,3) = state->cbn;

    A.block(6,0,3,3) = A31;
    A.block(6,3,3,3) = A32;
    A.block(6,6,3,3) = A33;

    Matrix<double,15,6> G = MatrixXd::Zero(15,6);
    G.block(0,0,3,3) = -state->cbn;
    G.block(3,0,3,3) = cross_vn * state->cbn;
    G.block(3,3,3,3) = state->cbn;

    Matrix<double,15,15> F = MatrixXd::Identity(15,15) + A * t + 0.5 * A * A * t *t;
    Matrix<double,6,6> Q0;
    Q0(0,0) = 2.3639e-12;
    Q0(1,1) = 2.3639e-12;
    Q0(2,2) = 2.3639e-12;
    Q0(3,3) = 2.42e-7;
    Q0(4,4) = 2.42e-7;
    Q0(5,5) = 2.42e-7;
    Matrix<double,15,15> M1 = G * Q0 * G.transpose();
    Matrix<double,15,15> M2 = F * M1 + (F * M1).transpose();
    Matrix<double,15,15> M3 = F * M2 + (F * M2).transpose();
    Matrix<double,15,15> M4 = F * M3 + (F * M3).transpose();
    Matrix<double,15,15> Q=t*M1+(t*t/2.0)*M2+(t*t*t / 6.0)*M3+(t*t*t*t / 24.0)*M4;

    Matrix<double,6,6> RK = MatrixXd::Zero(6,6);
    RK(0,0) = 0.01 * 0.01;
    RK(1,1) = 0.01 * 0.01;
    RK(2,2) = 0.01 * 0.01;
    RK(3,3) = 1.5749e-6 * 1.5749e-6;
    RK(4,4) = 1.5749e-6 * 1.5749e-6;
    RK(5,5) = 10 * 10;
    Matrix<double,6,15> H = MatrixXd::Zero(6,15);
    H.block(0,0,3,3) = cross_vn;
    H.block(0,3,6,6) = MatrixXd::Identity(6,6);
    
    //filter
    error_state->reset();  //initial error is zero
    Matrix<double,15,1> xk = F * error_state->x;
    Matrix<double,15,15> pk = F * error_state->P * F.transpose() + Q;
    Matrix<double,15,6> kalman_gain = pk * H.transpose()*(H*pk*H.transpose()+RK).inverse();
    Matrix<double,6,1> Z;
    Z << state->vn - gnss_data.vn,
          state->ve - gnss_data.ve,
          state->vd - gnss_data.vd,
          state->lat - gnss_data.lat;
          state->longitude - gnss_data.longitude;
          state->h - gnss_data.height;

    error_state->x = xk + kalman_gain * (Z - H * xk);
    error_state->P = (MatrixXd::Identity(15,15) - kalman_gain * H) * pk * (MatrixXd::Identity(15,15) - kalman_gain * H).transpose()
                     + kalman_gain * RK *kalman_gain.transpose();
    // update navigation state
    const Vector3d det_vn  = error_state->x.block(3,0,3,1) + cross_vn * error_state->x.block(0,0,3,1);
    state->vn -= det_vn(0);
    state->ve -=det_vn(1);
    state->lat -= error_state->x(6,0);
    state->longitude -= error_state->x(7,0);
    Matrix3d mat = (MatrixXd::Identity(3,3) + CrossMat(error_state->x.block(0,0,3,1))) * state->cbn;
    mat.normalized();
    state -> cbn = mat;
    //std::cout <<std::setprecision(10) << "error y : " << error_state->x(7,0) * R << std::endl;
    //std::cout <<std::setprecision(10) << "error lat : " << (state->lat - gnss_data.lat) * R << std::endl;
    std::cout <<std::setprecision(10) << "error long : " << (state->longitude - gnss_data.longitude) * R << std::endl;
    std::cout<< std::setprecision(10) << "state long : " <<  state->longitude << std::endl;
}



int main() {
    //read raw data
    std::ifstream infile;
    infile.open("/home/mdk/inti_navi/exp3/imu_data.txt");
    std::string line;
    std::vector<IMURawData> imu_data_coarse;
    std::vector<IMURawData> imu_data_precise;
    std::vector<IMURawData> imu_data_sins;
    int i =1;
    std::cout << "<<read date<<<< "  <<std::endl;
    while(getline(infile,line)) {
        IMURawData data;
        std::vector<double> line_data;
        std::istringstream str(line);
        std::string out;
        while (str >> out) {
            line_data.push_back(std::stod(out));
        }
        data.time = line_data[0];
        data.frame_number = line_data[1];
        data.angle_increment_x = line_data[2];
        data.angle_increment_y = line_data[3];
        data.angle_increment_z = line_data[4];
        data.velocity_increment_x = line_data[5];
        data.velocity_increment_y = line_data[6];
        data.velocity_increment_z = line_data[7];
        data.odo_vel = line_data[8];
        if(i > 0 && i < 12001) {
          imu_data_coarse.push_back(data);
          ++i;
        } else if(i < 60001) {
            imu_data_precise.push_back(data);
            ++i;
        } else {
            imu_data_sins.push_back(data);
            ++i;
        }
        
    }
    infile.close();

    std::ifstream infile_gnss;
    infile_gnss.open("/home/mdk/inti_navi/exp3/gps_data.txt");
    std::string line_gnss;
    std::vector<GNSSResultData> gnss_result_data;
    while(getline(infile_gnss,line_gnss)) {
        GNSSResultData data;
        std::vector<double> line_data;
        std::istringstream str(line_gnss);
        std::string out;
        while (str >> out) {
            line_data.push_back(std::stod(out));
        }
        data.lat = line_data[7];
        data.longitude = line_data[9];
        data.height = line_data[8];
        data.vn = line_data[5];
        data.ve = line_data[4];
        data.vd = line_data[6];
        gnss_result_data.push_back(data);
    }
    infile_gnss.close();
    std::cout << "size " << imu_data_sins.size() << gnss_result_data.size() << std::endl;
    // set cmb
    cmb << 1,7e-3,0,
        -7e-3,1,0,
        0,0,1;
    cmb = -1.007 * cmb;
    //initial euler
    NavigationState state;
    state.ResetPose();
    state.ResetVelocity();
    state.lat = Lat0;
    state.longitude = Long0;
    state.h = H0;

    //TODO () use coarse and precise Alignment to get initial angle
    Vector3d euler_angle(23.7879,0.2580,0.5128);
    Matrix3d rot = EulerAngleToRotion(euler_angle); //cnb
    state.cbn = rot.transpose();
    std::ofstream outfile;
    outfile.open("/home/mdk/inti_navi/exp3/zuhe_pos_vel_yaw1.txt",std::ios::out | std::ios::trunc);

    //integrated navigation
    FliterState error_state;
    error_state.reset();
    error_state.SetInitialP();
    int kk = 600;
    for (int i = 0; i < imu_data_sins.size()-1; i = i+2) {
      if(i % 100 == 0 && kk < 15522) {
          STExtendKalmanFilter({imu_data_sins[i],imu_data_sins[i+1]},gnss_result_data[kk],&state,&error_state);
          ++kk;
          outfile << std::setprecision(10) << state.lat << " " << state.longitude << " " 
                  << state.vn << " " << state.ve << " " << state.vd << " " <<
                   DcmToEulerAngle(state.cbn.transpose())(2)<< std::endl;
      } else {
          DR({imu_data_sins[i],imu_data_sins[i+1]},&state);
      }
    }
     outfile.close();
    return 0;
}