#include "tools.h"
#include <iostream>
#include <vector>
#include <cmath>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;
using std::abs;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	//check validity of inputs
	//estimations vector size is not zero
	//estimations vector size equals ground_truth vector size

	if(estimations.size() != ground_truth.size() || estimations.size() == 0){
		cout << "Invalid estimation or ground truth data" << endl;
		return rmse;
	}

	//calculate RMSE
	//squared residuals
	for (int i=0; i < estimations.size(); i++){
		VectorXd residual = estimations[i] - ground_truth[i];

		//square coefficient-wise
		residual = residual.array()*residual.array();

		//sum up
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//take the square root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//helper computations to shrink formulae + correction for very tiny values
	//double px = Tools::CheckNearZero(px_r);
	//double py = Tools::CheckNearZero(py_r);
	float c1 = Tools::CheckNearZero(px*px+py*py);
	float c2 = Tools::CheckNearZero(sqrt(c1));
	float c3 = Tools::CheckNearZero(c1*c2);
	float d1 = vx*py - vy*px;
	float d2 = vy*px - vx*py;


	//calculate the Jacobian
	Hj << (px/c2), (py/c2), 0, 0,
		  (-1.*py/c1), (px/c1), 0, 0,
		  py*d1/c3, px*d2/c3, px/c2, py/c2;

	return Hj;

}

//function to check for very close to zero values
float Tools::CheckNearZero(const float &value) {

	//threshold for comparison: values lower than this are corrected to this threshold
	const float threshold = 0.0001;
	float val = abs (value);

	if(val < threshold) {
		return threshold;
	} else {
		return value;
	}
}

