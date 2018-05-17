#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  
  ///* augmented state vector - not sre if i need this but its here for now
  VectorXd x_aug_;

  ///* state covariance matrix
  MatrixXd P_;
  
  ///*augmented covarience matrix - again not sure if i need this but its here for now
  MatrixXd P_aug_;
  
  ///predicted augmented sigma points - again not sure if i need this but is here for now
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;
  long previous_t_;
  double delta_t;
  
  ///*set radar measurement dims
  int n_z_;
  
  ///*matrix for sigma points in radar measurement space
  MatrixXd Zsig;
  
  ///*Matrix for predicted sigma points in measurement space
  VectorXd z_pred;
  
  ///*S Matrix
  MatrixXd S;
  
  ///* R Matrix
  MatrixXd Rrad;
  ///* R Matrix
  MatrixXd Rlas;
  // H matrix for laser
  MatrixXd H_;
  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;
  
  ///*tuning parameters
  double NIS_laser;
  double NIS_radar;  
  float radar_iteration;
  float laser_iteration;
  float radar_NIS_cnt;
  float laser_NIS_cnt;
  float radar_chi2;
  float laser_chi2;

  ///*

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
};

#endif /* UKF_H */
