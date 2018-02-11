#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() == 0) {
    cout << "The estimations vector size is zero!" << endl;
    return rmse;
  } else if (estimations.size() != ground_truth.size()) {
    cout << "The estimations verctor size in not equal to the ground truth one!" << endl;
    return rmse;
  }

  for (unsigned int i=0; i<estimations.size(); ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  rmse = rmse / estimations.size();

  return rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

  MatrixXd Hj(3, 4);
  Hj << 0, 0, 0, 0,
	      0, 0, 0, 0,
	      0, 0, 0, 0;

  float px = x_state(0);
	float py = x_state(1);
	float px2 = px * px;
	float py2 = py * py;
	float vx = x_state(2);
	float vy = x_state(3);

  //compute the Jacobian matrix
  Hj(0,0) = px / sqrt(px2+py2);
  Hj(0,1) = py / sqrt(px2+py2);
  Hj(1,0) = -py / (px2+py2);
  Hj(1,1) = px / (px2+py2);
  Hj(2,0) = py*(vx*py-vy*px) / pow((px2+py2),(3.0/2.0));
  Hj(2,1) = px*(vy*px-vx*py) / pow((px2+py2),(3.0/2.0));
  Hj(2,2) = Hj(0,0);
  Hj(2,3) = Hj(0,1);

  return Hj;
}
