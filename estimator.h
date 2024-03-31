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
      p_(6,6) = 1.5749e-6 * 1.5749e-6;
      p_(7,7) = 1.5749e-6 * 1.5749e-6;
      p_(8,8) = 10 * 10;
      p_(9,9) = 2.909e-6 * 2.909e-6;
      p_(10,10) = 2.909e-6 * 2.909e-6;
      p_(11,11) = 2.909e-6 * 2.909e-6;
      p_(12,12) = 0.0005 * 0.0005;
      p_(13,13) = 0.0005 * 0.0005;
      p_(14,14) = 0.0005 * 0.0005;
  }
  virtual ~ExtendKalmanFilter() = default;
  void Process(const Matrix<double,15,15>& F,
               const Matrix<double,15,15>& Q,
               const Matrix<double,6,1>& mesurement,
               const Matrix<double,6,15>& H,
               const Matrix<double,6,6>& Rk);
  Matrix<double,15,1> GetErrorState() { return error_state_; }

 private:
  void Predict(const Matrix<double,15,15>& F,
               const Matrix<double,15,15>& Q);
  void Update(const Matrix<double,6,1>& mesurement,
              const Matrix<double,6,15>& H,
              const Matrix<double,6,6>& Rk);
  void reset() {error_state_.setZero();}
  Matrix<double,15,1> error_state_ = MatrixXd::Zero(15,1);
  Matrix<double,15,15> p_ = MatrixXd::Zero(15,15); // Covariance
};

class IntegratedNavigation {
 public:
  IntegratedNavigation(const Matrix3d& cbn0,double lat0,double long0,double h0): 
                      navigation_state_(cbn0, lat0, long0, h0) {}
  
  virtual ~IntegratedNavigation() = default;
  void Process(const std::vector<IMURawData>& imu_data_sins,
               const std::vector<GNSSResultData>& gnss_result_data);

 private:
  static double t_;
  NavigationState navigation_state_;
  Matrix3d CoarseAlignment(const std::vector<IMURawData>& imu_data);
  void DR(const std::vector<IMURawData>& imu_data);
  void SINSCalculate(const std::vector<IMURawData>& imu_data);
  void STExtendKalmanFilter(const std::vector<IMURawData>& imu_data,
                            const GNSSResultData& gnss_data);
  void UpdateNavigationState(const Matrix3d& cross_vn,
                             const Matrix<double,15,1>& error_state);

  Matrix<double,6,6> ConstructRkMat();
  Matrix<double,6,15> ConstructMeasurementMat(
                                        const Matrix3d& cross_vn);
  Matrix<double,15,15> ConstructAMat(const Matrix3d& cross_vn);
  Matrix<double,15,15> ConstructQMat(
                                  const Matrix3d& cross_vn,
                                  const Matrix<double,15,15>& F);
  ExtendKalmanFilter ekf_;
};
