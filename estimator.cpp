#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iomanip>

#include "rotation.h"
#include "estimator.h"

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

void ExtendKalmanFilter::Predict(
                const Matrix<double,15,15>& F,
                const Matrix<double,15,15>& Q) {
    error_state_ = F * error_state_;
    p_ = F * p_ * F.transpose() + Q;
    std::ofstream outfile;
    outfile.open("/home/mdk/inti_navi/st_ekf_proj/zuhe_pos_vel_yaw2.txt",std::ios::out | std::ios::app);
    outfile << p_ << std::endl;
    outfile << "......." << std::endl;
    outfile << error_state_ << std::endl;
    outfile.close();
}

void ExtendKalmanFilter::Update(
         const Matrix<double,6,1>& mesurement,
         const Matrix<double,6,15>& H,
         const Matrix<double,6,6>& Rk) {
   const Matrix<double,15,6> kalman_gain = p_ * H.transpose()*(H*p_*H.transpose()+Rk).inverse();
   error_state_ = error_state_ + kalman_gain * (mesurement - H * error_state_);
   p_ = (MatrixXd::Identity(15,15) - kalman_gain * H) * p_ * (MatrixXd::Identity(15,15) - kalman_gain * H).transpose()
                     + kalman_gain * Rk *kalman_gain.transpose();
}

void ExtendKalmanFilter::Process(
                const Matrix<double,15,15>& F,
                const Matrix<double,15,15>& Q,
                const Matrix<double,6,1>& mesurement,
                const Matrix<double,6,15>& H,
                const Matrix<double,6,6>& Rk) {
    //filter
    reset();  //initial error state is zero
    Predict(F,Q);
    Update(mesurement, H, Rk);
}

double IntegratedNavigation::t_ = 0.02;

void IntegratedNavigation::Process(
                const std::vector<IMURawData>& imu_data_sins,
                const std::vector<GNSSResultData>& gnss_result_data) {
  std::ofstream outfile;
  outfile.open("/home/mdk/inti_navi/st_ekf_proj/zuhe_pos_vel_yaw1.txt",std::ios::out | std::ios::trunc);
  int kk = 600;
  std::cout << "lat : " << navigation_state_.lat << "," << navigation_state_.longitude 
               << std::endl;
  for (int i = 0; i < imu_data_sins.size()-1; i = i+2) {
      if(i % 100 == 0 && kk < 15522) {
          STExtendKalmanFilter({imu_data_sins[i],imu_data_sins[i+1]},gnss_result_data[kk]);
          ++kk;
          outfile << std::setprecision(10) << navigation_state_.lat << " " << navigation_state_.longitude << " " 
                  << navigation_state_.vn << " " << navigation_state_.ve << " " << navigation_state_.vd << " " <<
                   DcmToEulerAngle(navigation_state_.cbn.transpose())(2)<< std::endl;
      } else {
          DR({imu_data_sins[i],imu_data_sins[i+1]});
      }
  }
  outfile.close();
}

void IntegratedNavigation::STExtendKalmanFilter(
                        const std::vector<IMURawData>& imu_data,
                        const GNSSResultData& gnss_data) {
    SINSCalculate(imu_data);
    // const Vector3d det_vk1(imu_data[0].velocity_increment_x,imu_data[0].velocity_increment_y,imu_data[0].velocity_increment_z);
    // const Vector3d det_vk2(imu_data[1].velocity_increment_x,imu_data[1].velocity_increment_y,imu_data[1].velocity_increment_z); 
    // const Vector3d fn = state->cbn*((det_vk1 + det_vk2) / t);
    const Matrix3d cross_vn = CrossMat(Vector3d(navigation_state_.vn,navigation_state_.ve,navigation_state_.vd));
    
    const Matrix<double,15,15> A = ConstructAMat(cross_vn);
    const Matrix<double,15,15> F = MatrixXd::Identity(15,15) + A * t_ + 0.5 * A * A * t_ * t_;
    const Matrix<double,15,15> Q = ConstructQMat(cross_vn,F);
    const Matrix<double,6,6> Rk = ConstructRkMat();
    const Matrix<double,6,15> H = ConstructMeasurementMat(cross_vn);
    Matrix<double,6,1> Z;
    Z << navigation_state_.vn - gnss_data.vn,
         navigation_state_.ve - gnss_data.ve,
         navigation_state_.vd - gnss_data.vd,
         navigation_state_.lat - gnss_data.lat;
         navigation_state_.longitude - gnss_data.longitude;
         navigation_state_.h - gnss_data.height;
    
    // stekf
    ekf_.Process(F,Q,Z,H,Rk);

    // update navigation state
   UpdateNavigationState(cross_vn,ekf_.GetErrorState()); 
}

void IntegratedNavigation::UpdateNavigationState(
                         const Matrix3d& cross_vn,
                         const Matrix<double,15,1>& error_state) {
   const Vector3d det_vn  = error_state.block(3,0,3,1) + cross_vn * error_state.block(0,0,3,1);
    navigation_state_.vn -= det_vn(0);
    navigation_state_.ve -=det_vn(1);
    navigation_state_.lat -= error_state(6,0);
    navigation_state_.longitude -= error_state(7,0);
    Matrix3d mat = (MatrixXd::Identity(3,3) + CrossMat(error_state.block(0,0,3,1))) * navigation_state_.cbn;
    mat.normalized();
    navigation_state_.cbn = mat;
    // std::cout <<std::setprecision(10) << "error y : " << error_state(7,0) * R << std::endl;
    //std::cout <<std::setprecision(10) << "error lat : " << (state->lat - gnss_data.lat) * R << std::endl;
    // std::cout <<std::setprecision(10) << "error long : " << (state->longitude - gnss_data.longitude) * R << std::endl;
    // std::cout<< std::setprecision(10) << "state long : " <<  navigation_state_.longitude << std::endl;
}

Matrix3d IntegratedNavigation::CoarseAlignment (
                    const std::vector<IMURawData>& imu_data) {
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

void IntegratedNavigation::SINSCalculate(
                   const std::vector<IMURawData>& imu_data) {
    Vector3d v_old(navigation_state_.vn,navigation_state_.ve,navigation_state_.vd);
    const Matrix3d cbn_old = navigation_state_.cbn;

    Vector3d wie_n;
    wie_n << omega*cos(navigation_state_.lat),0,-omega*sin(navigation_state_.lat);
    Vector3d w_en_n;
    w_en_n << navigation_state_.ve/RE,-navigation_state_.vn/RN,-navigation_state_.ve*tan(navigation_state_.lat)/RE;
    Vector3d w_in_n;
    w_in_n = wie_n + w_en_n;
    Vector3d w_in_b =  navigation_state_.cbn.transpose() * w_in_n;
    Vector3d det_theta1(imu_data[0].angle_increment_x,imu_data[0].angle_increment_y,imu_data[0].angle_increment_z);
    det_theta1 -= w_in_b * t_ / 2;
    Vector3d det_theta2(imu_data[1].angle_increment_x,imu_data[1].angle_increment_y,imu_data[1].angle_increment_z);
    det_theta2 -= w_in_b * t_ / 2;
    Vector3d det_theta = det_theta1 + det_theta2;
    
    const Vector3d det_vk1(imu_data[0].velocity_increment_x,imu_data[0].velocity_increment_y,imu_data[0].velocity_increment_z);
    const Vector3d det_vk2(imu_data[1].velocity_increment_x,imu_data[1].velocity_increment_y,imu_data[1].velocity_increment_z);  
    const Vector3d det_vk = det_vk1+ det_vk2;

    const Vector3d det_rot = 0.5 * det_theta.cross(det_vk);
    const Vector3d det_vscul = 2 / 3.0 * (det_theta1.cross(det_vk2) + det_vk1.cross(det_theta2));
    const Vector3d det_vsfk = cbn_old * (det_vk + det_rot + det_vscul);
    
    const Vector3d det_vgcork = t_ * (gn - (2 * wie_n + w_en_n).cross(v_old));
    //update velocity
    Vector3d v_new = v_old + det_vsfk + det_vgcork;
    navigation_state_.vn = v_new(0);
    navigation_state_.ve = v_new(1);
    navigation_state_.vd = 0.0;

    //update pose
    navigation_state_.lat += 0.5 * (v_old(0) + v_new(0)) / RN * t_;
    navigation_state_.longitude += 0.5 * (v_old(1) + v_new(1)) / (RE*cos(navigation_state_.lat))*t_;
    navigation_state_.xn += 0.5*(v_old+v_new) * t_;

    //update rotation
    wie_n << omega*cos(navigation_state_.lat),0,-omega*sin(navigation_state_.lat);
    w_en_n << navigation_state_.ve / RE,-navigation_state_.vn / RN,-navigation_state_.ve*tan(navigation_state_.lat) / RE;
    w_in_b = navigation_state_.cbn.transpose() * (wie_n + w_en_n);
    det_theta1 = Vector3d(imu_data[0].angle_increment_x,imu_data[0].angle_increment_y,imu_data[0].angle_increment_z)
                  - w_in_b * t_ / 2;
    det_theta2 = Vector3d(imu_data[1].angle_increment_x,imu_data[1].angle_increment_y,imu_data[1].angle_increment_z)
                  - w_in_b * t_ / 2;
    Vector3d rot_vec =  det_theta1 + det_theta2 + 2 / 3.0 * det_theta1.cross(det_theta2);
    AngleAxisd angle_axis(rot_vec.norm(),rot_vec / rot_vec.norm());
    Matrix3d mat = navigation_state_.cbn * angle_axis.toRotationMatrix();
    mat.normalized();
    navigation_state_.cbn = mat;
}

void IntegratedNavigation::DR(
           const std::vector<IMURawData>& imu_data) {
    Vector3d v_odo_1(imu_data[0].odo_vel,0,0);
    Vector3d v_odo_2(imu_data[1].odo_vel,0,0);
    Vector3d vn1 = navigation_state_.cbn * cmb * v_odo_1;
    Vector3d vn2 = navigation_state_.cbn * cmb * v_odo_2; 
    Vector3d v_new = (vn1 + vn2) / 2;
    //update velocity
    navigation_state_.vn = v_new(0);
    navigation_state_.ve = v_new(1);
    navigation_state_.vd = v_new(2);
    //update position
    navigation_state_.lat +=  v_new(0) / RN * t_;
    navigation_state_.longitude += v_new(1) / (RE*cos(navigation_state_.lat)) * t_;
    navigation_state_.xn += v_new * t_;

    //update rotation
    Vector3d wie_n(omega*cos(navigation_state_.lat),0,-omega*sin(navigation_state_.lat));
    Vector3d w_en_n(navigation_state_.ve / RE,-navigation_state_.vn / RN,-navigation_state_.ve*tan(navigation_state_.lat) / RE);
    Vector3d w_in_b = navigation_state_.cbn.transpose() * (wie_n + w_en_n);
    Vector3d det_theta1 = Vector3d(imu_data[0].angle_increment_x,imu_data[0].angle_increment_y,imu_data[0].angle_increment_z)
                  - w_in_b * t_ / 2;
    Vector3d det_theta2 = Vector3d(imu_data[1].angle_increment_x,imu_data[1].angle_increment_y,imu_data[1].angle_increment_z)
                  - w_in_b * t_ / 2;
    Vector3d rot_vec =  det_theta1 + det_theta2 + 2 / 3.0 * det_theta1.cross(det_theta2);
    AngleAxisd angle_axis(rot_vec.norm(),rot_vec / rot_vec.norm());
    // Quaterniond q = state -> qbn * Quaterniond(angle_axis);
    // q = q.normalized().coeffs();
    // state->qbn = q;
    //std::cout<< "rot vec" <<rot_vec.norm() << std::endl;
     Matrix3d mat = navigation_state_.cbn * angle_axis.toRotationMatrix();
     mat.normalized();
     navigation_state_.cbn = mat;
    //state -> SetAngle(q);     
}

Matrix<double,6,6> IntegratedNavigation::ConstructRkMat() {
    // TODO: Set once
    Matrix<double,6,6> RK = MatrixXd::Zero(6,6);
    RK(0,0) = 0.01 * 0.01;
    RK(1,1) = 0.01 * 0.01;
    RK(2,2) = 0.01 * 0.01;
    RK(3,3) = 1.5749e-6 * 1.5749e-6;
    RK(4,4) = 1.5749e-6 * 1.5749e-6;
    RK(5,5) = 10 * 10;
    return RK;
}

Matrix<double,6,15> IntegratedNavigation::ConstructMeasurementMat(
                                        const Matrix3d& cross_vn) {
    Matrix<double,6,15> H = MatrixXd::Zero(6,15);
    H.block(0,0,3,3) = cross_vn;
    H.block(0,3,6,6) = MatrixXd::Identity(6,6);
    return H;
}

Matrix<double,15,15> IntegratedNavigation::ConstructAMat(const Matrix3d& cross_vn) {
    const double& re = RE+navigation_state_.h;
    const double& rn = RN+navigation_state_.h;
    const double& latitude = navigation_state_.lat;
    const double& ve = navigation_state_.ve;
    const double& vn = navigation_state_.vn;

    Matrix<double,3,3> A12;
    A12 << 0,1.0/re,0,
           - 1.0 / rn,0,0,
           0,-tan(latitude) / re,0;

    Matrix<double,3,3> A11;
    A11 << 0, 
    -(omega * sin(latitude) + ve / re * tan(latitude)),
    vn / rn,
    (omega * sin(latitude) + ve / re * tan(latitude)),
    0,
    (omega * cos(latitude) + ve / re),
    - vn / rn, 
    -(omega * cos(latitude) + ve / re),
    0;
    A11 += A12 * cross_vn;
    
    Matrix<double,3,3> A13;
    A13 << -omega * sin(latitude),0, -ve / re / re,
           0,0,vn / rn / rn,
           -(omega * cos(latitude) + ve / re) /cos(latitude) /cos(latitude),
           0,ve * tan(latitude) / re / re;
    
    Matrix<double,3,3> A21;
    A21 = -(CrossMat(gn) + cross_vn * CrossMat(Vector3d(omega*cos(latitude),0,-omega*sin(latitude))));

    Matrix<double,3,3> A22;
    A22 = -(2.0 * CrossMat(Vector3d(omega*cos(latitude),0,-omega*sin(latitude))) + CrossMat(Vector3d(ve / re,-vn / rn,-ve * tan(latitude) / re)));

    Matrix<double,3,3>  A23;
    A23 << -omega * sin(latitude),0,0,
            0,0,0,
            -omega * cos(latitude),0,0;
    A23 = cross_vn * A23;

    Matrix<double,3,3>  A32;
    A32 << 1.0 / rn,0,0,
           0,1.0/ re / cos(latitude),0,
           0,0,-1.0;
    Matrix<double,3,3>  A33;
    A33 << 0,0,-vn / rn / rn,
           ve * tan(latitude) / re / cos(latitude), 0 ,-ve / re / re / cos(latitude),
           0,0,0;
    Matrix<double,3,3>  A31;
    A31 = A32 * cross_vn;

    Matrix<double,15,15> A = MatrixXd::Zero(15,15);
    A.block(0,0,3,3) = A11;
    A.block(0,3,3,3) = A12;
    A.block(0,6,3,3) = A13;
    A.block(0,9,3,3) = -navigation_state_.cbn;

    A.block(3,0,3,3) = A21;
    A.block(3,3,3,3) = A22;
    A.block(3,6,3,3) = A23;
    A.block(3,9,3,3) = cross_vn * navigation_state_.cbn;
    A.block(3,12,3,3) = navigation_state_.cbn;

    A.block(6,0,3,3) = A31;
    A.block(6,3,3,3) = A32;
    A.block(6,6,3,3) = A33;

    return A;
}

Matrix<double,15,15> IntegratedNavigation::ConstructQMat(
                                  const Matrix3d& cross_vn,
                                  const Matrix<double,15,15>& F) {
    Matrix<double,15,6> G = MatrixXd::Zero(15,6);
    G.block(0,0,3,3) = -navigation_state_.cbn;
    G.block(3,0,3,3) = cross_vn * navigation_state_.cbn;
    G.block(3,3,3,3) = navigation_state_.cbn;

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
    Matrix<double,15,15> Q=t_*M1+(t_*t_/2.0)*M2+(t_*t_*t_ / 6.0)*M3+(t_*t_*t_*t_ / 24.0)*M4;
    return Q;
}

void ReadData(std::vector<IMURawData>* imu_data_coarse,
    std::vector<IMURawData>* imu_data_precise,
    std::vector<IMURawData>* imu_data_sins,
    std::vector<GNSSResultData>* gnss_result_data) {
    std::ifstream infile;
    infile.open("/home/mdk/inti_navi/st_ekf_proj/imu_data.txt",std::ios::in);
    if (!infile.is_open()) {
       std::cout << "Connot read imu data" << std::endl;
       return;
    }
    std::string line;
    int i =1;
    std::cout << "<<read data<<<< "  <<std::endl;
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
          imu_data_coarse->push_back(data);
          ++i;
        } else if(i < 60001) {
            imu_data_precise->push_back(data);
            ++i;
        } else {
            imu_data_sins->push_back(data);
            ++i;
        }
        
    }
    infile.close();

    std::ifstream infile_gnss;
    infile_gnss.open("/home/mdk/inti_navi/st_ekf_proj/gps_data.txt");
    if (!infile_gnss.is_open()) {
        std::cout << "Connot read gnss data" << std::endl;
        return;
    }
    std::string line_gnss;
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
        gnss_result_data->push_back(data);
    }
    infile_gnss.close();
    std::cout << "size " << imu_data_sins->size() << gnss_result_data->size() << std::endl;
}

void SetCmb() {
    // set cmb
    cmb << 1,7e-3,0,
        -7e-3,1,0,
        0,0,1;
    cmb = -1.007 * cmb;
}

void Estimator(
    const std::vector<IMURawData>& imu_data_coarse,
    const std::vector<IMURawData>& imu_data_precise,
    const std::vector<IMURawData>& imu_data_sins,
    const std::vector<GNSSResultData>& gnss_result_data) {
    //initial euler
    //TODO (): use coarse and precise Alignment to get initial angle
    Vector3d euler_angle(23.7879,0.2580,0.5128);
    Matrix3d rot = EulerAngleToRotion(euler_angle); //cnb

    //integrated navigation
    IntegratedNavigation integrated_navigation(rot.transpose(), Lat0, Long0, H0);
    integrated_navigation.Process(imu_data_sins, gnss_result_data);
}

int main() {
    //read raw data
    std::vector<IMURawData> imu_data_coarse;
    std::vector<IMURawData> imu_data_precise;
    std::vector<IMURawData> imu_data_sins;
    std::vector<GNSSResultData> gnss_result_data;
    ReadData(&imu_data_coarse, &imu_data_precise, &imu_data_sins, &gnss_result_data);
    SetCmb();
    Estimator(imu_data_coarse, imu_data_precise,imu_data_sins,gnss_result_data);
    return 0;
}