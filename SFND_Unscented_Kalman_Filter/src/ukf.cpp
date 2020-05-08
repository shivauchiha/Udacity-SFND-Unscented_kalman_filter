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
  std_a_ = 1.0;//30

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;//30
  
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
  n_x_ = 5;

  n_aug_= 7;

  lambda_= 3 - n_aug_;

  P_<<1,0,0,0,0,
      0,1,0,0,0,
      0,0,1,0,0,
      0,0,0,1,0,
      0,0,0,0,1;

  weights_ = VectorXd(2*n_aug_+1);
  weights_(0)= lambda_/(lambda_+n_aug_);
  for(int i=1;i<(2*n_aug_+1);i++)
  {
     weights_(i)=0.5*(1/(lambda_+n_aug_));
  }
 
 Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);
 Xsig_pred_.fill(double(0));

 //noise matrix
 noise_lidar = MatrixXd(2,2);
noise_lidar << std_laspx_*std_laspx_,0,
               0,std_laspy_*std_laspy_;

 noise_radar = MatrixXd(3,3);
 noise_radar<<std_radr_*std_radr_,0,0,
              0,std_radphi_*std_radphi_,0,
              0,0, std_radrd_*std_radrd_;
              


}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   * 
   */
  
    
    if(is_initialized_!=true)
    { double Px,Py,v ;
      Px=Py=v=0;
      if(meas_package.sensor_type_ == MeasurementPackage::LASER)
      {
          Px=meas_package.raw_measurements_[0];
          Py=meas_package.raw_measurements_[1];
          P_(0,0)=std_laspx_*std_laspx_;
          P_(1,1)=std_laspy_*std_laspy_;
          
      }

      if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
      {
        double range = meas_package.raw_measurements_[0];
        double angle= meas_package.raw_measurements_[1];
        //we wont be taking velocity information as that velocity is with respective to ego entity
        Px=range*cos(angle);
        Py=range*sin(angle);
      }
      x_<<Px,Py,v,0,0;
      time_us_=meas_package.timestamp_;
      is_initialized_=true;
      

    }
  
    else
    {
      double delt=(meas_package.timestamp_-time_us_)/1000000.0;
      time_us_=meas_package.timestamp_;
      Prediction(delt);
      if(meas_package.sensor_type_ == MeasurementPackage::RADAR &&use_radar_)
      {
         UpdateRadar(meas_package);
      }
      else if(meas_package.sensor_type_ == MeasurementPackage::LASER &&use_laser_)
      {
         UpdateLidar(meas_package);
      }
      
    }
    
  
  
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  
  
   Eigen::VectorXd x_aug = VectorXd(7); //augmented mean vector
   Eigen::MatrixXd P_aug = MatrixXd(7,7);//augmented state vector
   Eigen::MatrixXd Xsig_aug = MatrixXd(n_aug_,2*n_aug_+1);
   x_aug.head(5) = x_;
   x_aug[5]=0;
   x_aug[6]=0;
   P_aug.fill(0.0);
   P_aug.topLeftCorner(5,5)=P_;
   P_aug(5,5)=std_a_*std_a_;
   P_aug(6,6)=std_yawdd_*std_yawdd_;
   Eigen::MatrixXd A = P_aug.llt().matrixL();
  //Augmented sigma points
   Xsig_aug.col(0) = x_aug;
   for (int i = 0; i< n_aug_; ++i) {
     Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * A.col(i);
     Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i);
   }
  
 
   
  //extractions 
  for(int i =0;i<2*n_aug_+1;i++)
  {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double psi = Xsig_aug(3,i);
    double psid = Xsig_aug(4,i);
    double na = Xsig_aug(5,i);
    double nv = Xsig_aug(6,i);
    
    double px_p,py_p;
    
    //avoid division by zero
    if (psid !=0)
    {
      px_p=p_x + v/psid * ( sin (psi + psid*delta_t) - sin(psi));
      py_p=p_y + v/psid * ( cos(psi) - cos(psi+psid*delta_t) );
    }
    else
    {
      px_p = p_x + v*delta_t*cos(psi);
      py_p = p_y + v*delta_t*sin(psi);
    
    }
    double v_p = v;
    double psi_p = psi + psid*delta_t;
    double psid_p = psid;
    
    //lets add some noise
    px_p = px_p + 0.5*na*delta_t*delta_t * cos(psi);
    py_p = py_p + 0.5*na*delta_t*delta_t * sin(psi);
    v_p = v_p + na*delta_t;
    psi_p = psi_p + 0.5*nv*delta_t*delta_t;
    psid_p = psid_p + nv*delta_t;
    
    // xsig pred
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = psi_p;
    Xsig_pred_(4,i) = psid_p;
  }
    
    x_.fill(0.0);
    //predict mean
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
    {  
    x_ = x_ +  Xsig_pred_.col(i)*weights_(i);
    }

  // predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
    {  
      Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
     while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

     P_ = P_ +   (x_diff * x_diff.transpose())*weights_(i);
    }
  
  
  
  

  
  
  
  
  
  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  
  int n_z=2;//measurement dimension of radar
  Eigen::MatrixXd Zsig = MatrixXd(n_z,2*n_aug_+1);
  Eigen::VectorXd Z_pred = VectorXd(n_z);
  Eigen::MatrixXd S = MatrixXd(n_z,n_z);
  
  for(int i =0 ; i<2*n_aug_+1;i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);
    double psid = Xsig_pred_(4, i);
    
    Zsig(0,i) = p_x;                       
    Zsig(1,i) = p_y;                                
    
    
  } 
  
  Z_pred.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    Z_pred = Z_pred + weights_(i)*Zsig.col(i);
  }
  
  S.fill(0.0);
  for (int i=0; i< 2*n_aug_+1;++i)
  {
    VectorXd Z_diff = Zsig.col(i) - Z_pred;
    

    S = S + weights_(i) * Z_diff * Z_diff.transpose();
  }
  
  S = S + noise_lidar;
  //extract sensor data
  Eigen::VectorXd Z = meas_package.raw_measurements_;
  //matrix for cross correlation Tc
  Eigen::MatrixXd Tc = MatrixXd(n_x_,n_z);
  
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    Eigen::VectorXd Z_diff = Zsig.col(i) - Z_pred;
    

   
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    

    Tc = Tc + weights_(i) * x_diff * Z_diff.transpose();
  }

  // Kalman gain K;
  Eigen::MatrixXd K = Tc * S.inverse();

  
  Eigen::VectorXd z_diff = Z - Z_pred;

  // angle normalization
 

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  
  
  
  
  
  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  
  int n_z=3;//measurement dimension of radar
  Eigen::MatrixXd Zsig = MatrixXd(n_z,2*n_aug_+1);
  Eigen::VectorXd Z_pred = VectorXd(n_z);
  Eigen::MatrixXd S = MatrixXd(n_z,n_z);
  
  for(int i =0 ; i<2*n_aug_+1;i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double psi = Xsig_pred_(3,i);
    double psid = Xsig_pred_(4, i);
    
    double v1 = cos(psi)*v;
    double v2 = sin(psi)*v;
    
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       
    Zsig(1,i) = atan2(p_y,p_x);                                
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
    
  } 
  
  Z_pred.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++)
  {
    Z_pred = Z_pred + weights_(i)*Zsig.col(i);
  }
  
  S.fill(0.0);
  for (int i=0; i< 2*n_aug_+1;++i)
  {
    VectorXd Z_diff = Zsig.col(i) - Z_pred;
    while (Z_diff(1)> M_PI) Z_diff(1)-=2.*M_PI;
    while (Z_diff(1)<-M_PI) Z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * Z_diff * Z_diff.transpose();
  }
  
  S = S + noise_radar;
  //extract sensor data
  Eigen::VectorXd Z = meas_package.raw_measurements_;
  //matrix for cross correlation Tc
  Eigen::MatrixXd Tc = MatrixXd(n_x_,n_z);
  
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    Eigen::VectorXd Z_diff = Zsig.col(i) - Z_pred;
    while (Z_diff(1)> M_PI) Z_diff(1)-=2.*M_PI;
    while (Z_diff(1)<-M_PI) Z_diff(1)+=2.*M_PI;

   
    Eigen::VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * Z_diff.transpose();
  }

  // Kalman gain K;
  Eigen::MatrixXd K = Tc * S.inverse();

  
  Eigen::VectorXd z_diff = Z - Z_pred;

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
 
  

  
}