#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in; // state
  P_ = P_in; // posterior covariance
  F_ = F_in; // state transfer function
  H_ = H_in; // observation model
  R_ = R_in; // measurement noise covariance
  Q_ = Q_in; // process noise covariance
}

void KalmanFilter::Predict() {
  /**
  DONE:
    * predict the state
  */

  x_ = F_*x_; // +u was in lesson T2.L5.6
	P_ = F_*P_*F_.transpose()+Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  DONE:
    * update the state by using Kalman Filter equations
  */
  MatrixXd I;
  I.setIdentity(P_.rows(), P_.cols());
  VectorXd y = z - H_*x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  x_ = x_ + K * y;
  P_ = (I-K*H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //Update(z);
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  MatrixXd I;
  I.setIdentity(P_.rows(), P_.cols());

  // map state space to measure space (lesson segment 19)
  auto px = x_(0);
  auto py = x_(1);
  auto vx = x_(2);
  auto vy = x_(3);

  // based on review, avoid update if the position
  // may cause math errors
  if ( px < 0.0001 && py < 0.0001 )
    return; // abort the update

  auto hx = VectorXd(3);
  auto r =  sqrt(px*px+py*py);
  hx[0] =r;
  
  hx[1] = atan2(py,px);
  hx[2] = fabs(r)>0?(px*vx+py*vy)/r:0;

  VectorXd y = z - hx;

  // normalize angle between -pi and pi
  while(y(1)>M_PI) {
    y(1) = y(1) - 2 * M_PI;
  }
  while(y(1)<-M_PI) {
    y(1) = y(1) + 2 * M_PI;
  }

  MatrixXd PHt = P_ * H_.transpose(); 
  MatrixXd S = H_ * PHt + R_;
  MatrixXd K = PHt * S.inverse();

  x_ = x_ + K * y;
  P_ = (I-K*H_)*P_;

}
