#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
   * TODO: predict the state
   */
	x_ = F_ * x_; //predict new state
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_; //predict new state covariance
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Ht * Si;

	cout << "z_pred: " << z_pred << endl;
	cout << "y: " << y << endl;
	cout << "H: " << H_ << endl;
	cout << "Si: " << Si << endl;
	cout << "K: " << K << endl;


	//new state
	x_ = x_ + K * y;
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
	//recover state parameters
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);


	//convert to polar coordinates i.e. calculate h(x') = [rho, phi, rho_dot].transposed
	//1. rho + check and correct for very tiny values so that it will not cause problem at the division later
	//float r = sqrt(px*px+py*py);
	float rho = tools.CheckNearZero(sqrt(px*px+py*py));


	//2. phi
	float phi = atan2(py, px);
	//float phi = 0.0;

	if (px == 0) {
		cout << "px is zero" << px << endl;
	}

	//adjust phi value between pi and -pi
	const float pi = 3.14159265;
	if(phi < -1.0*pi){
		phi += 2*pi;
	}

	else if (phi > pi){
		phi -= 2*pi;
	}


	//3. rho_dot
	float rho_dot = (px*vx + py*vy)/rho;

	//this is the predicted data
	VectorXd z_pred = VectorXd(3);

	z_pred << rho, phi, rho_dot;

	//this is the error between the actual and the predicted data
	VectorXd y = z - z_pred;

	//calculate the Jacobian
	Hj_ = tools.CalculateJacobian(x_);

	MatrixXd Hjt = Hj_.transpose();
	MatrixXd S = Hj_ * P_ * Hjt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = P_ * Hjt * Si;


	//if (t_step == 274){
	//	cout << "x_: " << x_ << endl;
	//	cout << "rho: " << rho << endl;
	//	cout << "phi: " << phi << endl;
	//	cout << "rho_dot: " << rho_dot << endl;
		cout << "z: " << z << endl;
		cout << "py: " << py << endl;
		cout << "px: " << px << endl;
		cout << "z_pred: " << z_pred << endl;
		cout << "y: " << y << endl;
		cout << "Hj: " << Hj_ << endl;
		cout << "Si: " << Si << endl;
		cout << "K: " << K << endl;
	//}

	//update prediction with measured data
	x_ = x_ + K * y;
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;


}
