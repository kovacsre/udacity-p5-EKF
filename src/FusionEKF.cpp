#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  //Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix - laser
  ekf_.H_ = MatrixXd(2, 4);
  ekf_.H_ << 1, 0, 0, 0,
		  	 0, 1, 0, 0;
  //measurement matrix - radar
  ekf_.Hj_= MatrixXd(3, 4);

  //initial state transition matrix F
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
		     0, 1, 0, 1,
			 0, 0, 1, 0,
			 0, 0, 0, 1;

  //state covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
		     0, 1, 0, 0,
			 0, 0, 1000, 0,
			 0, 0, 0, 1000;

  //process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      //
    	float rho = measurement_pack.raw_measurements_[0];
    	float phi = measurement_pack.raw_measurements_[1];
    	//as (vx, vy) cannot be obtained from rho_dot, keep initialized velocity values in ekf_.x_
    	//position values (px, py) are calculated from rho and phi
    	ekf_.x_(0) = rho * cos (phi); //x in cartesian coordinates
    	ekf_.x_(1) = rho * sin (phi); //y in cartesian coordinates
    	cout << "state: " << ekf_.x_ << endl;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
    	//update ekf_.x_ with measured position values (px, py) initial velocity (vx, vy) values are kept
    	float px = measurement_pack.raw_measurements_[0];
    	float py = measurement_pack.raw_measurements_[1];

    	ekf_.x_(0) = px;
		ekf_.x_(1) = py;
    	cout << "state: " << ekf_.x_ << endl;

    }
    //store timestamp for dt calculation later
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  //noise values
  short noise_ax = 9;
  short noise_ay = 9;

  //calculate elapsed time in seconds (deltaT)
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  //save current timestamp to use in the next cycle
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt2 = dt*dt;
  float dt3 = dt2*dt;
  float dt4 = dt3*dt;



  //update F_ according to the elapsed time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //update process covariance matrix Q with the elapsed time
  ekf_.Q_ << dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
		     0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
			 dt3/2*noise_ax, 0, dt2*noise_ax, 0,
			 0, dt3/2*noise_ay, 0, dt2*noise_ay;


  ekf_.Predict();

  /**
   * Update
   */


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
