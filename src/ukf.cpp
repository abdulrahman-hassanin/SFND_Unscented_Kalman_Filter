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
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 0.0225, 0,
        0, 0, 0, 0, 0.0225;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // std_a_ = 30;
   std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // std_yawdd_ = 30;
   std_yawdd_ = 1.0;
  
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
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  // state dimension
  n_x_ = 5;

  // augmanted state dimension
  n_aug_ = 7;
  
  time_us_ = 0;

  // lambda design parameter, Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  // Initailize weights
  weights_(0)   = lambda_ / (lambda_+n_aug_);
  double weight = 0.5 * (lambda_+n_aug_);
  for(int i=1; i<2*n_aug_+1; i++)
    weights_(i) = weight;  

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
  
  // Initialize measurement noise covairance matix
  radarNoise = MatrixXd(3, 3);
  radarNoise << std_radr_*std_radr_, 0, 0,
                0, std_radphi_*std_radphi_, 0,
                0, 0, std_radrd_*std_radrd_;

  lidarNoise = MatrixXd(2,2);
  lidarNoise << std_laspx_*std_laspx_, 0,
                0, std_laspy_*std_laspy_;
  
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
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        
        double rho     = meas_package.raw_measurements_[0];     
        double phi     = meas_package.raw_measurements_[1];     
        double rho_dot = meas_package.raw_measurements_[2];     
        
        double x       = rho * cos(phi);
        double y       = rho * sin(phi);
        double vx      = rho_dot * cos(phi);
        double vy      = rho_dot * sin(phi);
        double v       = sqrt(vx * vx + vy * vy);

        x_ << x , y, v, rho, rho_dot;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      x_<< meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;        
    }
    
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // Calculate delta time between k and k+1 and update time_us_
  double dt = (meas_package.timestamp_ - time_us_) * 1e-6;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
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
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
  }
  
  
  /*
  *
  *     Predict Augmanted Sigma Points   
  */ 

  // predict sigma points
  for (int i = 0; i< 2*n_aug_+1; ++i) {
    // extract values for better readability
    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double yawd     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    double px_p, py_p, v_p, yaw_p, yawd_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (-cos (yaw + yawd*delta_t) + cos(yaw));
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    v_p = v;
    yaw_p = yaw + yawd*delta_t;
    yawd_p = yawd;

    // add noise
    px_p   = px_p   + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p   = py_p   + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p    = v_p    + nu_a*delta_t;
    yaw_p  = yaw_p  + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
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
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_; i++)
  {
    Z_pred = Z_pred + weights_(i) * Zsig.col(i);
  }
  for (int i = 0; i < 2*n_aug_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - Z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  // add noise
  S = S + lidarNoise;

  // Calculate Cross-Correlation between sigma points and measurment in state space
  MatrixXd Tc = MatrixXd(n_x_, z_n);
  Tc.fill(0.0);
  for (int i=0; i< 2*n_aug_+1; i++)
  {
    VectorXd Z_diff = Zsig.col(i) - Z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * Z_diff.transpose();
  }

  // laser measurment vector
  VectorXd Z = meas_package.raw_measurements_;

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
    v1 = v*sin(yaw);

    Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);                         //r
    Zsig(1, i) = atan2(p_y, p_x);                                 //phi
    Zsig(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);     //rdot
  }

  // predicted measurment mean and covariance
  Z_pred.fill(0.0);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_; i++)
  {
    Z_pred = Z_pred + weights_(i) * Zsig.col(i);
  }
  for (int i = 0; i < 2*n_aug_; i++)
  {
    VectorXd z_diff = Zsig.col(i) - Z_pred;

    // Angle Normalization
    while (z_diff(1) > M_PI)    z_diff(1) -= 2.*M_PI;
    while (z_diff(1) <-M_PI)    z_diff(1) += 2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  
  // Add Measurment covarince Matrix
  S = S + radarNoise;

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
    while (x_diff(2) > M_PI)    x_diff(2) -= 2.*M_PI;
    while (x_diff(2) <-M_PI)    x_diff(2) += 2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // laser measurment vector
  VectorXd Z = meas_package.raw_measurements_;

  VectorXd z_diff = Z - Z_pred;
  // Kalman Gain
  MatrixXd K = Tc * S.inverse();
  // State Update
  x_ = x_ + K * z_diff;
  // Covariance Matrix Update
  P_ = P_ - K * S * K.transpose();  
  
  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
}