#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {

  x_ = F_*x_;
  P_ = F_*P_*F_.transpose()+Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // LASER
  //update the state by using Kalman Filter equations
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
//std::cout << "\nLASER y:" << y << std::endl;
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;  
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // RADAR
  // map the state to the measurement coordinates
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);
  
  float rho = sqrt(px*px + py*py);
  float theta = atan2(py,px);

  if (rho < 0.01) { // prevent divide by zero
    rho = 0.01;
  }
  float rho_dot = (px*vx+py*vy)/rho;
  
  VectorXd z_pred = VectorXd(3);
  z_pred << rho, theta, rho_dot;
   
  VectorXd y = z - z_pred;

  //normalize the theta error, y(1), to be in range (-pi, pi)
  if (y(1) > M_PI) {
    y(1) = -2.*M_PI + y(1);
  }
  if (y(1) < -M_PI) {
    y(1) = 2.*M_PI + y(1);
  }

//  std::cout << "meas_theta:" << z(1) << "\npred_theta:" << theta << std::endl;
//std::cout << "\nRADAR y:" << y << std::endl;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;  
}
