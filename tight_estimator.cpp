#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iomanip>

#include "rotation.h"
#include "tight_estimator.h"
#include "constant.h"
#include "earth_util.h"

Vector3d gn = Vector3d(0,0,g);

void ExtendKalmanFilter::Predict(
                const MatrixXd& F,
                const MatrixXd& Q) {
    error_state_ = F * error_state_;
    p_ = F * p_ * F.transpose() + Q;
}

void ExtendKalmanFilter::Update(
         const MatrixXd& mesurement,
         const MatrixXd& H,
         const MatrixXd& Rk) {
   const MatrixXd kalman_gain = p_ * H.transpose()*(H*p_*H.transpose()+Rk).inverse();
   error_state_ = error_state_ + kalman_gain * (mesurement - H * error_state_);
//    p_ = (MatrixXd::Identity(17,17) - kalman_gain * H) * p_ * (MatrixXd::Identity(17,17) - kalman_gain * H).transpose()
//                      + kalman_gain * Rk *kalman_gain.transpose();
   p_ = (MatrixXd::Identity(17,17) - kalman_gain * H) * p_;
   p_ = (p_ + p_.transpose()) * 0.5;
   std::ofstream outfile;
}

void ExtendKalmanFilter::Process(
                const MatrixXd& F,
                const MatrixXd& Q,
                const MatrixXd& mesurement,
                const MatrixXd& H,
                const MatrixXd& Rk) {
    //filter
    reset();  //initial error state is zero
    Predict(F,Q);
    Update(mesurement, H, Rk);
}

void TightlyIntegratedNavigation::SINSCalculate(
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

void TightlyIntegratedNavigation::UpdateNavigationState(
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
}

Matrix<double,17,17> TightlyIntegratedNavigation::ConstructFMat(const Matrix3d& cross_vn) {
    const double& re = RE+navigation_state_.h;
    const double& rn = RN+navigation_state_.h;
    const double& latitude = navigation_state_.lat;
    const double& ve = navigation_state_.ve;
    const double& vn = navigation_state_.vn;

    Matrix<double,3,3> A12;
    A12 << 0.0,1.0/re,0.0,
           - 1.0 / rn,0.0,0.0,
           0.0,-tan(latitude) / re,0.0;

    Matrix<double,3,3> A11;
    A11 << 0.0, 
    -(omega * sin(latitude) + ve / re * tan(latitude)),
    vn / rn,
    (omega * sin(latitude) + ve / re * tan(latitude)),
    0.0,
    (omega * cos(latitude) + ve / re),
    - vn / rn, 
    -(omega * cos(latitude) + ve / re),
    0.0;
    A11 += A12 * cross_vn;
    
    Matrix<double,3,3> A13;
    A13 << -omega * sin(latitude),0.0, -ve / re / re,
           0.0,0.0,vn / rn / rn,
           -(omega * cos(latitude) + ve / re /cos(latitude) /cos(latitude)),
           0.0,ve * tan(latitude) / re / re;
    
    Matrix<double,3,3> A21;
    A21 = -(CrossMat(gn) + cross_vn * CrossMat(Vector3d(omega*cos(latitude),0,-omega*sin(latitude))));

    Matrix<double,3,3> A22;
    A22 = -(2.0 * CrossMat(Vector3d(omega*cos(latitude),0,-omega*sin(latitude))) + CrossMat(Vector3d(ve / re,-vn / rn,-ve * tan(latitude) / re)));

    Matrix<double,3,3>  A23;
    A23 << -omega * sin(latitude),0.0,0.0,
            0.0,0.0,0.0,
            -omega * cos(latitude),0.0,0.0;
    A23 = cross_vn * A23;

    Matrix<double,3,3>  A32;
    A32 << 1.0 / rn,0.0,0.0,
           0,1.0/ re / cos(latitude),0.0,
           0.0,0.0,-1.0;
    Matrix<double,3,3>  A33;
    A33 << 0.0,0.0,-vn / rn / rn,
           ve * tan(latitude) / re / cos(latitude), 0.0 ,-ve / re / re / cos(latitude),
           0.0,0.0,0.0;
    Matrix<double,3,3>  A31;
    A31 = A32 * cross_vn;

    Matrix<double,17,17> A = MatrixXd::Zero(15,15);
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
    A(15,16) = 1;
    
    const Matrix<double,17,17> F = MatrixXd::Identity(17,17) + A * t_ + 0.5 * A * A * t_ * t_;
    return F;
}

Matrix<double,17,17> TightlyIntegratedNavigation::ConstructQMat(
                                  const Matrix3d& cross_vn,
                                  const Matrix<double,17,17>& F) {
    Matrix<double,17,6> G = MatrixXd::Zero(17,6);
    G.block(0,0,3,3) = -navigation_state_.cbn;
    G.block(3,0,3,3) = cross_vn * navigation_state_.cbn;
    G.block(3,3,3,3) = navigation_state_.cbn;

    Matrix<double,6,6> Q0 =  MatrixXd::Zero(6,6);
    Q0(0,0) = 2.3639e-12;
    Q0(1,1) = 2.3639e-12;
    Q0(2,2) = 2.3639e-12;
    Q0(3,3) = 2.42e-6;
    Q0(4,4) = 2.42e-6;
    Q0(5,5) = 2.42e-6;
    Matrix<double,17,17> M1 = G * Q0 * G.transpose();
    Matrix<double,17,17> M2 = F * M1 + (F * M1).transpose();
    Matrix<double,17,17> M3 = F * M2 + (F * M2).transpose();
    Matrix<double,17,17> M4 = F * M3 + (F * M3).transpose();
    Matrix<double,17,17> Q=t_*M1+(t_*t_/2.0)*M2+(t_*t_*t_ / 6.0)*M3+(t_*t_*t_*t_ / 24.0)*M4;
    return Q;
}

MatrixXd TightlyIntegratedNavigation::ConstructRkMat(int sat_nums) {
    int mat_size = 2 * sat_nums;
    MatrixXd RK = MatrixXd::Zero(mat_size,mat_size);
    for (int i = 0; i < sat_nums; ++i) {
      RK(i,i) = 0.1 * 0.1;  
      RK(sat_nums + i, sat_nums + i) = 0.001 * 0.001;
    }
    return RK;
}

MatrixXd TightlyIntegratedNavigation::ConstructMeasurementMat(const Vector3d& pos,
                                 const Matrix3d& cross_vn,
                                 const Matrix3d& ned_to_ecef,
                                 const Matrix3d& llh_to_xyz_diff,
                                 const GNSSRawData& gnss_raw_data) {   
  const int sat_num = gnss_raw_data.sat_datas.size();
  MatrixXd H = MatrixXd::Zero(2 * sat_num,17);
  for (int i = 0; i < sat_num; ++i) {
    const Vector3d& sat_pos = gnss_raw_data.sat_datas[i].sat_pos;
    const Vector3d& sat_vel = gnss_raw_data.sat_datas[i].sat_vel;
    double r = (pos - sat_pos).norm();
    const Vector3d& error = (pos - sat_pos) / r;
    H.block(i,6,1,3) = error.transpose()* llh_to_xyz_diff; 
    H(i,15) = -1;
    H.block(i + sat_num,0,1,3) = error.transpose() * ned_to_ecef * cross_vn;
    H.block(i + sat_num,3,1,3) = error.transpose() * ned_to_ecef;
    H(i + sat_num,16) = -1;
  }
  return H;
}

MatrixXd TightlyIntegratedNavigation::BuildMeasurements(
                           const Vector3d& pos,
                           const Vector3d& vn,
                           const Matrix3d& ned_to_ecef,
                           const GNSSRawData& gnss_raw_data) {
  int sat_num = gnss_raw_data.sat_datas.size();
  const Vector3d ve = ned_to_ecef * vn;
  MatrixXd mesurement =  MatrixXd::Zero(2 * sat_num,1);
  for (int i = 0; i < sat_num; ++i) {
    const Vector3d& sat_pos = gnss_raw_data.sat_datas[i].sat_pos;
    const Vector3d& sat_vel = gnss_raw_data.sat_datas[i].sat_vel;
    const double r = (pos - sat_pos).norm(); 
    mesurement(i,0) =  r - 
                           gnss_raw_data.sat_datas[i].pseudorange - 
                           c * gnss_raw_data.sat_datas[i].sat_dt;
    const double error_vel = (pos - sat_pos).dot((ve - sat_vel)) / r;
    mesurement(i+sat_num,0) = gnss_raw_data.sat_datas[i].doppler_shift * c / 1575.42e6 + 
                                  error_vel  - 
                                  c * gnss_raw_data.sat_datas[i].sat_dtt;
  }
  return mesurement;
}

void TightlyIntegratedNavigation::STExtendKalmanFilter(
                        const std::vector<IMURawData>& imu_data,
                        const GNSSRawData& gnss_data) {
    SINSCalculate(imu_data);

    const Vector3d vn = Vector3d(navigation_state_.vn,navigation_state_.ve,navigation_state_.vd);
    const Matrix3d cross_vn = CrossMat(vn);
    
    const Matrix<double,17,17> F = ConstructFMat(cross_vn);
    const Matrix<double,17,17> Q = ConstructQMat(cross_vn,F);

    const MatrixXd Rk = ConstructRkMat(gnss_data.sat_datas.size());

    const Vector3d pos = 
            LLHToXYZ(navigation_state_.lat, navigation_state_.longitude, navigation_state_.h);
    const Matrix3d ned_to_ecef = NedToEcef(navigation_state_.lat, navigation_state_.longitude);
    const Matrix3d llh_to_xyz_diff = 
            LLHToXYZDiff(navigation_state_.lat, navigation_state_.longitude, navigation_state_.h);
    const MatrixXd H = ConstructMeasurementMat(pos,cross_vn,ned_to_ecef,llh_to_xyz_diff,gnss_data);

    const MatrixXd Z = BuildMeasurements(pos,vn,ned_to_ecef,gnss_data);
    
    // stekf
    ekf_.Process(F,Q,Z,H,Rk);

    // update navigation state
   UpdateNavigationState(cross_vn,ekf_.GetErrorState()); 
}

double TightlyIntegratedNavigation::t_ = 0.02;

void TightlyIntegratedNavigation::Process(
                const std::vector<IMURawData>& imu_data_sins,
                const std::vector<GNSSRawData>& gnss_result_data) {
  std::ofstream outfile;
  outfile.open("/home/mdk/inti_navi/st_ekf_proj/tight_estimator_result.txt",std::ios::out | std::ios::trunc);
  int kk = 0;
  std::cout << "llh : " << navigation_state_.lat << "," << navigation_state_.longitude 
               << "," << navigation_state_.h << std::endl;
  for (int i = 0; i < imu_data_sins.size()-1; i = i+2) {  
      if(i % 100 == 0 && kk < 809) {
          STExtendKalmanFilter({imu_data_sins[i],imu_data_sins[i+1]},gnss_result_data[kk]);
          ++kk;
          outfile << std::setprecision(10) << navigation_state_.lat << " " << navigation_state_.longitude << " " 
                  << navigation_state_.vn << " " << navigation_state_.ve << " " << navigation_state_.vd << " " <<
                   DcmToEulerAngle(navigation_state_.cbn.transpose())(2) << " " <<navigation_state_.xn(0)
                   << " " << navigation_state_.xn(1) << std::endl;
      } else {
          SINSCalculate({imu_data_sins[i],imu_data_sins[i+1]});
      }     
  }
  outfile.close();
}

void Estimator(
    const std::vector<IMURawData>& imu_data_coarse,
    const std::vector<IMURawData>& imu_data_precise,
    const std::vector<IMURawData>& imu_data_sins,
    const std::vector<GNSSRawData>& gnss_result_data) {
    //initial euler
    //TODO (): use coarse and precise Alignment to get initial angle
    Vector3d euler_angle(92.2743,-0.3395,-0.6511);  
    Matrix3d rot = EulerAngleToRotion(euler_angle); //cnb
    
    double init_lat = 0.493040846947;
    double init_long = 1.97419090135;
    double init_h = 41.69;
    //integrated navigation
    TightlyIntegratedNavigation tightly_integrated_navigation(rot.transpose(), init_lat, init_long, init_h);
    tightly_integrated_navigation.Process(imu_data_sins, gnss_result_data);
}

void ReadData(std::vector<IMURawData>* imu_data_sins,
              std::vector<GNSSRawData>* gnss_result_data) {
    std::ifstream infile;
    infile.open("/home/mdk/inti_navi/st_ekf_proj/imu_data.txt");
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
        data.angle_increment_x = line_data[0];
        data.angle_increment_y = line_data[1];
        data.angle_increment_z = line_data[2];
        data.velocity_increment_x = line_data[3];
        data.velocity_increment_y = line_data[4];
        data.velocity_increment_z = line_data[5];
        imu_data_sins->push_back(data);  
    }
    infile.close();

    std::ifstream infile_gnss;
    infile_gnss.open("/home/mdk/inti_navi/st_ekf_proj/gnss_data_1.txt");
    if (!infile_gnss.is_open()) {
        std::cout << "Connot read gnss data" << std::endl;
        return;
    }
    std::string line_gnss;
    while(getline(infile_gnss,line_gnss)) {
        GNSSRawData data;
        std::vector<double> line_data;
        std::istringstream str(line_gnss);
        std::string out;
        while (str >> out) {
            line_data.push_back(std::stod(out));
        }
        data.time_stamp = line_data[0];
        data.sat_nums = line_data[1];        
        gnss_result_data->push_back(data);
    }
    infile_gnss.close();
    
    infile_gnss.open("/home/mdk/inti_navi/st_ekf_proj/gnss_data_2.txt");
    if (!infile_gnss.is_open()) {
           std::cout << "Connot read gnss data2" << std::endl;
           return;
    }
    std::vector<GNSSRawData::SatelliteData> sat_datas;
    while(getline(infile_gnss,line_gnss)) {
            std::vector<double> line_data;
            std::istringstream str(line_gnss);
            std::string out;
            while (str >> out) {
              line_data.push_back(std::stod(out));
            }
            GNSSRawData::SatelliteData satllite_data;
            satllite_data.id = line_data[0];
            satllite_data.pseudorange = line_data[1];
            satllite_data.phase = line_data[2];
            satllite_data.doppler_shift = line_data[3];
            satllite_data.sat_pos = Vector3d(line_data[4],line_data[5],line_data[6]);
            satllite_data.sat_vel = Vector3d(line_data[7],line_data[8],line_data[9]);
            satllite_data.sat_dt = line_data[10];
            satllite_data.sat_dtt = line_data[11];
            sat_datas.push_back(satllite_data);
    }
    infile_gnss.close();
    
    int current_size = 0;
    for (int i = 0; i < gnss_result_data->size();++i) {
      const int sat_nums = (*gnss_result_data)[i].sat_nums;
      (*gnss_result_data)[i].sat_datas.reserve(sat_nums);
      for (int j = 0; j < sat_nums; ++j) {
        (*gnss_result_data)[i].sat_datas.push_back(sat_datas[j + current_size]);
      }
      current_size += sat_nums;
    }
    std::cout << "gnss data size : " << gnss_result_data->size() << std::endl;
}

int main() {
    //read raw data
    std::vector<IMURawData> imu_data_coarse;
    std::vector<IMURawData> imu_data_precise;
    std::vector<IMURawData> imu_data_sins;
    std::vector<GNSSRawData> gnss_raw_data;
    ReadData(&imu_data_sins, &gnss_raw_data);
    Estimator(imu_data_coarse, imu_data_precise,imu_data_sins,gnss_raw_data);
    return 0;
}
