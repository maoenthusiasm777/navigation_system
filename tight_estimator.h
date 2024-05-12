#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>

#include "data_type.h"

using namespace Eigen;

class ExtendKalmanFilter {
 public:
  ExtendKalmanFilter() {
      // Set initial val
      p_(0,0) = 0.0007 * 0.0007;
      p_(1,1) = 0.0007 * 0.0007;
      p_(2,2) = 0.0007 * 0.0007;
      p_(3,3) = 0.01 * 0.01;
      p_(4,4) = 0.01 * 0.01;
      p_(5,5) = 0.01 * 0.01;
      p_(6,6) = 5.5749e-6 * 5.5749e-6;
      p_(7,7) = 1.5749e-5 * 1.5749e-5;
      p_(8,8) = 10 * 10;
      p_(9,9) = 2.909e-5 * 2.909e-5;
      p_(10,10) = 2.909e-5 * 2.909e-5;
      p_(11,11) = 2.909e-5 * 2.909e-5;
      p_(12,12) = 0.0005 * 0.0005;
      p_(13,13) = 0.0005 * 0.0005;
      p_(14,14) = 0.0005 * 0.0005;
      p_(15,15) = 1 * 1;
      p_(16,16) = 0.1 * 0.1;
  }
  virtual ~ExtendKalmanFilter() = default;
  void Process(const MatrixXd& F,
               const MatrixXd& Q,
               const MatrixXd& mesurement,
               const MatrixXd& H,
               const MatrixXd& Rk);
  MatrixXd GetErrorState() const { return error_state_; }

 private:
  void Predict(const MatrixXd& F,
               const MatrixXd& Q);
  void Update(const MatrixXd& mesurement,
              const MatrixXd& H,
              const MatrixXd& Rk);
  void reset() {error_state_.setZero();}
  Matrix<double,17,1> error_state_ = MatrixXd::Zero(17,1);
  Matrix<double,17,17> p_ = MatrixXd::Zero(17,17); // Covariance
};

class TightlyIntegratedNavigation {
 public:
  TightlyIntegratedNavigation(const Matrix3d& cbn0,double lat0,double long0,double h0): 
                      navigation_state_(cbn0, lat0, long0, h0) {}
  
  virtual ~TightlyIntegratedNavigation() = default;
  void Process(const std::vector<IMURawData>& imu_data_sins,
               const std::vector<GNSSRawData>& gnss_result_data);

 private:
  static double t_;
  NavigationState navigation_state_;
  void SINSCalculate(const std::vector<IMURawData>& imu_data);
  void STExtendKalmanFilter(const std::vector<IMURawData>& imu_data,
                            const GNSSRawData& gnss_data);

  void UpdateNavigationState(const Matrix3d& cross_vn,
                             const Matrix<double,15,1>& error_state);

  MatrixXd ConstructRkMat(int sat_nums);
  MatrixXd ConstructMeasurementMat(const Vector3d& pos,
                                   const Matrix3d& cross_vn,
                                   const Matrix3d& ned_to_ecef,
                                   const Matrix3d& llh_to_xyz_diff,
                                   const GNSSRawData& gnss_raw_data);
  MatrixXd BuildMeasurements(const Vector3d& pos,
                             const Vector3d& vn,
                             const Matrix3d& ned_to_ecef,
                             const GNSSRawData& gnss_raw_data);
  Matrix<double,17,17> ConstructFMat(const Matrix3d& cross_vn);
  Matrix<double,17,17> ConstructQMat(const Matrix3d& cross_vn,
                                  const Matrix<double,17,17>& F);
  ExtendKalmanFilter ekf_;
};
