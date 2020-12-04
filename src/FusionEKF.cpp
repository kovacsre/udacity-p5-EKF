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

  timestep = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

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
  //ekf_.Q_ << 0, 0, 0, 0,
	//       0, 0, 0, 0,
	//		 0, 0, 0, 0,
	//		 0, 0, 0, 0;
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
    ekf_.x_ << 1, 1, 1, 1; //tweak velocity data for better RMSE values at start

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      //
    	float rho = measurement_pack.raw_measurements_[0];
    	float phi = measurement_pack.raw_measurements_[1];
    	//float rho_dot = measurement_pack.raw_measurements_[2];
    	//as (vx, vy) cannot be obtained from rho_dot, keep initialized velocity values in ekf_.x_
    	//position values (px, py) are calculated from rho and phi
    	ekf_.x_(0) = rho * cos (phi); //x in cartesian coordinates
    	ekf_.x_(1) = rho * sin (phi); //y in cartesian coordinates
    	cout << "state: " << ekf_.x_ << endl;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
    	//update ekf_.x_ with measured position values (px, py) initial velocity (vx, vy) values are kept
    	ekf_.x_(0) = measurement_pack.raw_measurements_[0];
    	ekf_.x_(1) = measurement_pack.raw_measurements_[1];
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

  //calculate elapsed time in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt2 = dt*dt;
  float dt3 = (dt2*dt) / 2;
  float dt4 = (dt2*dt2) / 4;



  //update F_ according to the elapsed time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //update process covariance matrix Q with the elapsed time
  ekf_.Q_ << dt4*noise_ax, 0, dt3*noise_ax, 0,
		     0, dt4*noise_ay, 0, dt3*noise_ay,
			 dt3*noise_ax, 0, dt2*noise_ax, 0,
			 0, dt3*noise_ay, 0, dt2*noise_ay;


  ekf_.Predict();

  /**
   * Update
   */


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
	  ekf_.R_ = R_radar_;
	  cout << "ekf_.x_: " << ekf_.x_ << endl;
	  cout << "ekf_.P_: " << ekf_.P_ << endl;
	  cout << "ekf_.F_: " << ekf_.F_ << endl;
	  cout << "R_radar_: " << ekf_.R_ << endl;
	  cout << "ekf_.Q_: " << ekf_.Q_ << endl;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // TODO: Laser updates
	  ekf_.R_ = R_laser_;
	  cout << "R_laser_: " << ekf_.R_ << endl;
	  cout << "ekf_.x_: " << ekf_.x_ << endl;
	  cout << "ekf_.P_: " << ekf_.P_ << endl;
	  cout << "ekf_.F_: " << ekf_.F_ << endl;
	  cout << "R_radar_: " << ekf_.R_ << endl;
	  cout << "ekf_.Q_: " << ekf_.Q_ << endl;
	  ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
