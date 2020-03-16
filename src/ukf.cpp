#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
 

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // std_a_ = 30;
   std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // std_yawdd_ = 30;
   std_yawdd_ = 1.5;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

  // state dimension
  n_x_ = 5;

  // augmanted state dimension
  n_aug_ = 7;
  
  time_us_ = 0;

  // lambda design parameter, Sigma point spreading parameter
  lambda_ = 3 - n_x_;
  
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Intialize the weights
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  double weight = 0.5/(lambda_+n_aug_);
  weights_(0) = weight_0;

  for (int i=1; i<2*n_aug_+1; ++i) {  
    weights_(i) = weight;
  }
  
  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  
  
  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;
} 

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if(!is_initialized_)
  { 
    if(meas_package.sensor_type_ == meas_package.LASER && use_laser_){
      x_ << meas_package.raw_measurements_[0], 
            meas_package.raw_measurements_[1], 
            0,
            0, 
            0;
    }else if(meas_package.sensor_type_ == meas_package.RADAR && use_radar_){
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      x_ << rho * cos(phi),
            rho * sin(phi),
            0,
            0,
            0;
    }
    
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // Calculate delta time between k and k+1 and update time_us_
  double dt = (meas_package.timestamp_ - time_us_) /1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if(meas_package.sensor_type_ == meas_package.LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }  
  else if(meas_package.sensor_type_ == meas_package.RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }

  
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  
  /*
  *
  *       Generate Augmanted Sigma Points
  */ 
  VectorXd x_aug = VectorXd(n_aug_);                      // create augmented mean vector
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);              // create augmented state covariance
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);   // create sigma point matrix

  // augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  
  /*
  *
  *     Predict Augmanted Sigma Points   
  */ 

  // predict sigma points
  for(int i=0; i<2*n_aug_+1; i++){
    
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double vk = Xsig_aug(2, i);
    double yawk = Xsig_aug(3, i);
    double yawdotk = Xsig_aug(4, i);
    double miua = Xsig_aug(5, i);
    double miuyaw = Xsig_aug(6, i);

    if(fabs(yawdotk) <0.001){
      Xsig_pred_(0, i) = px + (vk * cos(yawk) * delta_t) + (0.5 * delta_t*delta_t*cos(yawk)*miua);  
      Xsig_pred_(1, i) = py + (vk * sin(yawk) * delta_t) + (0.5 * delta_t*delta_t*sin(yawk)*miua);  
    }else{
      Xsig_pred_(0, i) = px + vk/yawdotk * ( sin (yawk + yawdotk*delta_t) - sin(yawk)) + (0.5 * delta_t*delta_t*cos(yawk)*miua);  
      Xsig_pred_(1, i) = py + vk/yawdotk * ( cos(yawk) - cos(yawk+yawdotk*delta_t) ) + (0.5 * delta_t*delta_t*sin(yawk)*miua);   
    }
    Xsig_pred_(2, i) = vk + delta_t * miua;
    Xsig_pred_(3, i) = yawk + 0.5 * delta_t * delta_t * miuyaw;
    Xsig_pred_(4, i) = yawdotk + delta_t * miuyaw;  
  }


  /**
   * Predict mean anf covariance
   */
  x_.fill(0.0);
  
  
  // iterate over sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    // Angle Normalizatuin
    while (x_diff(3) > M_PI)    x_diff(3) -= 2.*M_PI;
    while (x_diff(3) <-M_PI)    x_diff(3) += 2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    
  }  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  
  int z_n = 2;    // measrurment parameters
  
  VectorXd Z = meas_package.raw_measurements_;	  // laser measurment vector
  
  MatrixXd Zsig   = MatrixXd(z_n, 2*n_aug_+1);    // Measurment Sigma Points Matrix
  VectorXd Z_pred = VectorXd(z_n);                // Predicted Measurment Mean Vector
  MatrixXd S      = MatrixXd(z_n, z_n);           // Predicted Measurment Covariance Matrix

  // Transform Sigma points into measurment Space
  for (int i=0; i<2*n_aug_+1; i++)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  // predicted measurment mean and covariance
  Z_pred.fill(0.0);  
  for (int i = 0; i < 2*n_aug_; i++)
  {
    Z_pred = Z_pred + weights_(i) * Zsig.col(i);
  }
  
  MatrixXd R = MatrixXd(2, 2);
  R.fill(0);
  R <<std_laspx_*std_laspx_, 0,
      0, std_laspy_*std_laspy_;
  
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - Z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  // add noise
  S = S + R;

  // Calculate Cross-Correlation between sigma points and measurment in state space
  MatrixXd Tc = MatrixXd(n_x_, z_n);
  Tc.fill(0.0);
  for (int i=0; i< 2*n_aug_+1; i++)
  {
    VectorXd Z_diff = Zsig.col(i) - Z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * Z_diff.transpose();
  }



  VectorXd z_diff = Z - Z_pred;
  // Kalman Gain
  MatrixXd K = Tc * S.inverse();
  // State Update
  x_ = x_ + K * z_diff;
  // Covariance Matrix Update
  P_ = P_ - K * S * K.transpose();
  
  //calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  
  int z_n = 3;                                    // measrurment parameters  
  VectorXd Z = meas_package.raw_measurements_;    // laser measurment vector
  MatrixXd Zsig   = MatrixXd(z_n, 2*n_aug_+1);    // Measurment Sigma Points Matrix
  VectorXd Z_pred = VectorXd(z_n);                // Predicted Measurment Mean Vector
  MatrixXd S      = MatrixXd(z_n, z_n);           // Predicted Measurment Covariance Matrix

  // Transform Sigma points into measurment Space
  double p_x, p_y, v, yaw, yawd;
  double v1, v2; 
  for (int i=0; i<2*n_aug_+1; i++)
  {
    p_x   = Xsig_pred_(0, i);
    p_y   = Xsig_pred_(1, i);
    v     = Xsig_pred_(2, i);
    yaw   = Xsig_pred_(3, i);
    yawd  = Xsig_pred_(4, i);

    v1 = v*cos(yaw);
    v2 = v*sin(yaw);

    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                         //r
    Zsig(1, i) = atan2(p_y, p_x);                                 //phi
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);     //rdot
  }

  // predicted measurment mean and covariance
  Z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    Z_pred = Z_pred + weights_(i) * Zsig.col(i);
  }
  
  MatrixXd R = MatrixXd(3, 3);
  R.fill(0);
  S.fill(0);
  R(0, 0) = std_radr_*std_radr_;
  R(1, 1) = std_radphi_*std_radphi_;
  R(2, 2) = std_radrd_*std_radrd_;
  for(int i=0; i<2 * n_aug_ +1; i++){
    VectorXd z_diff = Zsig.col(i) - Z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S = S + weights_(i) * z_diff * z_diff.transpose() ;
  }
  S = S + R;

  // Calculate Cross-Correlation between sigma points and measurment in state space
  MatrixXd Tc = MatrixXd(n_x_, z_n);
  Tc.fill(0.0);
  for (int i=0; i< 2*n_aug_+1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - Z_pred;
    // Angle Normalization
    while (z_diff(1) > M_PI)    z_diff(1) -= 2.*M_PI;
    while (z_diff(1) <-M_PI)    z_diff(1) += 2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // Angle Normalization
    while (x_diff(3) > M_PI)    x_diff(3) -= 2.*M_PI;
    while (x_diff(3) <-M_PI)    x_diff(3) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }


  VectorXd z_diff = Z - Z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  // Kalman Gain
  MatrixXd K = Tc * S.inverse();
  // State Update
  x_ = x_ + K * z_diff;
  // Covariance Matrix Update
  P_ = P_ - K * S * K.transpose();  
  
  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}