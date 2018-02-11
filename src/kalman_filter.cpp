#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

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
  /**
  TODO:
    * predict the state
  */

  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}


void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */

  /******************************
  z_pred =      H_     *  x_

                        / px \
  / px \   /1, 0, 0, 0\ | py |
  \ py / = \0, 1, 0, 0/ | vx |
                        \ vy /
  *******************************/
  VectorXd z_pred = H_ * x_;
  // `z`: raw measurement
  // `y`: the difference between raw measurement and prediction
  VectorXd y = z - z_pred;
  Estimate(y);
}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  double px = x_(0);
  double py = x_(1);
  double px2 = px * px;
  double py2 = py * py;
  double vx = x_(2);
  double vy = x_(3);

  double rho = sqrt(px2 + py2);
  double phi = atan2(py, px);
  double rho_dot;

  if ((px2 + py2) >= 0.0001) {
    rho_dot = (px*vx + py*vy) / sqrt(px2 + py2);
  } else {
    cout << endl << "Divided by zero !!!!!!!!!!!!!!" << endl << endl;
    rho_dot = 0;
  }

  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;
  Estimate(y);
}


//void KalmanFilter::Estimate(const VectorXd &y) {
void KalmanFilter::Estimate(VectorXd &y) {
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_* H_.transpose() * S.inverse();

  //double before_y1 = y(1);

  while (y(1) < -M_PI || M_PI < y(1)) {
    if (y(1) < -M_PI) {
      cout << "y(1) = phi is less than M_PI" << endl;
      y(1) += 2 * M_PI;
    } else if (M_PI < y(1)) {
      cout << "y(1) = phi is more than M_PI" << endl;
      y(1) -= 2 * M_PI;
    }
  }

  /*
  if (abs(before_y1-y(1)) > 0.0001) {
    cout << "[before] y(1) = phi is " << before_y1/M_PI << "*M_PI" << endl;
    cout << "[after] y(1) = phi is " << y(1)/M_PI << "*M_PI" << endl <<endl;
  }
  */

  // new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}
