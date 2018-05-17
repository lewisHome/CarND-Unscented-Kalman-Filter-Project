#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // flag to indicate intialisation state of filter
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  
  //initial state dimension
  n_x_ = 5;
  
  //initial augmented state dimensions
  n_aug_ = 7;
  
  //initial lambda
  lambda_=3-n_aug_;  

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
    
  //initial predicted sigma point matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  //weights vector
  weights_ = VectorXd(2*n_aug_+1);
   
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  //set RADAR noise matrix
  Rrad = MatrixXd(3,3);
     //R measurement noise matrix
   Rrad << std_radr_*std_radr_,0,0,
        0,std_radphi_*std_radphi_,0,
        0,0,std_radrd_*std_radrd_;

  //set Laser noise matrix
  Rlas=MatrixXd(2,2);
  Rlas <<std_laspx_*std_laspx_,0,
          0,std_laspy_*std_laspy_;
          
  //
  NIS_radar = 0.0;
  //
  NIS_laser = 0.0;
  
  //NIS output variables
  radar_iteration = 0.0;
  laser_iteration = 0.0;
  radar_NIS_cnt = 0.0;
  laser_NIS_cnt = 0.0;
  radar_chi2 = 0.0;
  laser_chi2 = 0.0;
  
  //time variable
  delta_t = 0;
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (! is_initialized_)
  {
    //Initialise x_, P_, previous time and any other shizzle I need here
    x_ = VectorXd(5);
    x_.fill(0.0);
    if (use_radar_ == true && use_laser_ == true)
    {
      cout<<"Radar and Laser used for positioning"<<endl;
    }
    else if (use_laser_ == true && not use_radar_)
    {
      cout<<"Only Laser used for positioning"<<endl;
    }
    else if (use_radar_ == true && not use_laser_)
    {
      cout<<"Only Radar used for positioning"<<endl;
    }
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR  && use_radar_==true)
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //read radar data
      float ro = meas_package.raw_measurements_(0);
      float theta = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);

      //convert polar coordinates to cartesian
      float x=ro*cos(theta);
      float y=ro*sin(theta);
      
      //fill state matrix
      x_ << x,y,ro_dot,theta,0.0;    
      P_ << std_radr_*std_radr_,0,0,0,0,
            0,std_radr_*std_radr_,0,0,0,
            0,0,0.1,0,0,
            0,0,0,std_radphi_,0,
            0,0,0,0,std_radphi_;
      
      cout<<"Radar used to initialise"<<endl;
      is_initialized_ = true;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER  && use_laser_ == true)
    {
      float x=meas_package.raw_measurements_(0);
      float y=meas_package.raw_measurements_(1);

      //fill state matrix - only fill first two values as lidar has no associated      
      //veloicty data
      x_(0)=x;
      x_(1)=y;

      P_ << std_laspx_*std_laspx_,0,0,0,0,
            0,std_laspy_*std_laspy_,0,0,0,
            0,0,0.1,0,0,
            0,0,0,0.1,0,
            0,0,0,0,0.1;

      cout<<"Laser used to initialise"<<endl;
      is_initialized_ = true;        
    }

  
    //Set weights
    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i = 1; i<2*n_aug_ +1; i++)
    {
      weights_(i) = 0.5/(n_aug_ + lambda_);
    }

    previous_t_ = meas_package.timestamp_;
    return;
  }
  
  else
  {
    delta_t = (meas_package.timestamp_ - previous_t_)/1000000.0;
    previous_t_ = meas_package.timestamp_;

    Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::LASER  && use_laser_ == true)
    {
      UpdateLidar(meas_package);
    }

    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true)
    {
      UpdateRadar(meas_package);
    }
    return;
  }
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  //Prediction step 1. generate augmented sigma points lesson 7 section 18
  // initial augmented state vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  
  // initial augmented covariance matrix
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);
  
  // initial predicted sigma points matrix - need dims
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ +1);

  lambda_=3-n_aug_;  
  //create augmented mean state
  x_aug_.head(5)=x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;
  
  //create augmeneted covarience matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
  //create square root matrix
    MatrixXd L = P_aug_.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i=0; i<n_aug_; i++)
  {
    Xsig_aug_.col(i+1)=x_aug_+sqrt(lambda_+n_aug_)*L.col(i);
    Xsig_aug_.col(i+1+n_aug_)=x_aug_-sqrt(lambda_+n_aug_)*L.col(i);
   }
   
   
  //Prediction step 2. predict sigma points lesson 7 section 20
    for (int i=0; i<2*n_aug_+1;i++) 
    {
      double p_x = Xsig_aug_(0,i);
      double p_y = Xsig_aug_(1,i);
      double v = Xsig_aug_(2,i);
      double yaw = Xsig_aug_(3,i);
      double yawd = Xsig_aug_(4,i);
      double nu_a = Xsig_aug_(5,i);
      double nu_yawdd = Xsig_aug_(6,i);
    
      //predicted state values
      double px_p, py_p;
    
      //avoid division by zero
      if (fabs(yawd)> 0.000001)
      {
        px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t));
      }
      else 
      {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
      }
    
      double v_p =v;
      double yaw_p = yaw+ yawd*delta_t;
      double yawd_p = yawd;
      double delta_t2 = delta_t*delta_t;

      //add noise
      px_p = px_p+0.5*nu_a*delta_t2*cos(yaw);
      py_p = py_p+0.5*nu_a*delta_t2*sin(yaw);
      v_p = v_p + nu_a*delta_t;

      yaw_p = yaw_p + 0.5*nu_yawdd*delta_t2;
      yawd_p = yawd_p + nu_yawdd*delta_t;
    
      Xsig_pred_(0,i) = px_p;
      Xsig_pred_(1,i) = py_p;
      Xsig_pred_(2,i) = v_p;
      Xsig_pred_(3,i) = yaw_p;
      Xsig_pred_(4,i) = yawd_p;
    }
 
  //Prediction step 3. predict mean and covarience Lesson 7 Section 23
  //Step 3.2 predict state mean
  x_.fill(0.0);
  for(int i = 0; i<2*n_aug_+1; i++){
      x_=x_+weights_(i) * Xsig_pred_.col(i);
     } 
  //Step 3.3 predict state covarience matrix
  P_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++){
  //state difference
  VectorXd x_diff = Xsig_pred_.col(i) - x_;
  //angle normalisation
  while(x_diff(3)>M_PI) x_diff(3)-=2.*M_PI;
  while(x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
  P_=P_+weights_(i) * x_diff * x_diff.transpose();  
  }
//  cout<<"Predict update x \n"<<x_<<endl;
//  cout<<"Predict update P \n"<<P_<<endl;

}
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

 	//measurement matrix
  int n_z_ =2;

  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_+1);
  
  VectorXd z_pred = VectorXd(n_z_);
  
  MatrixXd S = MatrixXd(n_z_,n_z_);

  VectorXd z = VectorXd(n_z_);

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  
  //initialise variables
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);
  Tc.fill(0.0);

  double las_px = meas_package.raw_measurements_(0);
  double las_py = meas_package.raw_measurements_(1);

  z<< las_px,
      las_py;
  
  for (int i=0; i<2*n_aug_ +1; i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
                 
    z_pred +=weights_(i)*Zsig.col(i);
  }
    
  for (int i=0; i< 2*n_aug_ +1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
    
  S = S+Rlas;
       
  //calc cross correlation
  for(int i = 0; i<2*n_aug_ +1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //normalise angles
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //normalise angles
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  
  //residual
  VectorXd z_diff = z-z_pred;

  //normalise angles
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
  //Kalman gain
  MatrixXd K = Tc * S.inverse();

  //Calcualte NIS
  NIS_laser = z_diff.transpose() * S.inverse() * z_diff;

  laser_iteration+=1;
  if (NIS_laser < 7.815)
  {
    laser_NIS_cnt+=1;
  }

  laser_chi2 = 100 * (laser_NIS_cnt/laser_iteration);
  cout<<"Laser chi squared = "<<laser_chi2<<endl;

  //save NIS
  ofstream NIS_laser_f_("laser_NIS.csv", ios::app);
  NIS_laser_f_<<NIS_laser<<"\n";
  NIS_laser_f_.close();
  //cout<<"NIS laser: "<<NIS_laser<<endl;
	//new estimate	
	x_ = x_ + (K * z_diff);
  //cout<<"Lidar update x \n"<<x_<<endl;
	P_ = P_ - K*S*K.transpose();
  //cout<<"Lidar update P \n"<<P_<<endl;
	return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //Step 1 Transform sigma points into measurment space- leasson 7 section 27
  // set measurement dimensions for radar r,phi and r_dot
  n_z_ = 3;
  
  MatrixXd Zsig = MatrixXd(n_z_, 2*n_aug_+1);
  
  VectorXd z_pred = VectorXd(n_z_);
  
  MatrixXd S = MatrixXd(n_z_,n_z_);

  VectorXd z_ = VectorXd(n_z_);

  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  
  //initialise variables
  Zsig.fill(0.0);
  z_pred.fill(0.0);
  S.fill(0.0);
  Tc.fill(0.0);
  
  for (int i=0; i<2*n_aug_+1; i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
  
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
  
    // measuremetn model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);    //r
    Zsig(1,i) = atan2(p_y,p_x);             //phi
    //avoid division by zero
    if (Zsig(0,i)<0.000001)
    {
     Zsig(2,i) = (p_x*v1 + p_y*v2)/0.0000001;
    }
    else
    { 
      Zsig(2,i) = (p_x*v1 + p_y*v2)/Zsig(0,i);//r_dot
    }
  }

  //Calc mean predicted measurment
  for(int i=0; i<2*n_aug_ + 1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //Calc innovation covariance matrix S
  for (int i=0; i<2*n_aug_+1; i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //normalise angles
    while(z_diff(1)>M_PI) z_diff(1)-=2.*M_PI;
    while(z_diff(1)<-M_PI) z_diff(1) +=2.*M_PI;
    //measurement covarience matrix S
    S = S+weights_(i) * z_diff * z_diff.transpose();
   }
   
   //add noise   
   S=S+Rrad;
  
  //load radar measurement
  double meas_rho = meas_package.raw_measurements_(0);
  double meas_phi = meas_package.raw_measurements_(1);
  double meas_rhod = meas_package.raw_measurements_(2);
  
  z_ << meas_rho,
       meas_phi,
       meas_rhod;
     
  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //calculate NIS
  NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
  
  radar_iteration+=1;
  if (NIS_radar < 7.815)
  {
    radar_NIS_cnt+=1;
  }

  radar_chi2 = 100*(radar_NIS_cnt/radar_iteration);
  cout<<"Radar chi squared = "<<radar_chi2<<endl;
  //save NIS
  ofstream NIS_radar_f_("Radar_NIS.csv", ios::app);
  NIS_radar_f_<<NIS_radar<<"\n";
  NIS_radar_f_.close();
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;  
  //cout<<"Radar update x \n"<<x_<<endl;
  P_ = P_ - K*S*K.transpose();
  //cout<<"Radar update P \n"<<P_<<endl;
  return;
}


