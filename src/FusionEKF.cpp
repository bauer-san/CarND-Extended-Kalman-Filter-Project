#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,      0,    0,
                 0, 0.0009,    0,
                 0,      0, 0.09;

  //measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
 
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Initialize state using the initial RADAR measurements converted from polar to cartesian coordinates
	  ekf_.x_ << measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]),
				 measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]),
                 measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]),
                 measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);
;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state using the initial LASER measurements
	  ekf_.x_ << measurement_pack.raw_measurements_[0],
				 measurement_pack.raw_measurements_[1], 
				 0,
				 0;
    }

    //Initialize the state covariance matrix, P
	ekf_.P_ = MatrixXd(4,4);
	ekf_.P_ << 1., 0., 0., 0., // We are *SO* sure of the initial state
	           0., 1., 0., 0.,
	           0., 0., 1., 0.,
	           0., 0., 0., 1.;
	
    noise_ax = 9;
    noise_ay = 9;
	
	// initialize the previous timestamp to this initial measurement
	previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
	
  //1. Modify the F matrix so that the time is integrated
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, dt, 0,
			 0, 1, 0, dt,
			 0, 0, 1, 0,
			 0, 0, 0, 1;
  //2. Set the process covariance matrix Q
  float dt2 = dt*dt;    //dt^2
  float dt3 = dt2*dt/2; //dt^3/2
  float dt4 = dt3*dt/2; //dt^4/4

  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << dt4*noise_ax,     0,              dt3*noise_ax,    0,
		     0,                dt4*noise_ay,   0,               dt3*noise_ay,
			 dt3*noise_ax,     0,              dt2*noise_ax,    0,
			 0,                dt3*noise_ay,   0,               dt2*noise_ay;
  
//  if (dt > 0.01) {
    ekf_.Predict();
//  }

//  cout << "x_ after prediction\n" << ekf_.x_ << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	ekf_.R_ = R_radar_;	
	
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
	
  } else {
    // Laser updates
	ekf_.H_ = H_laser_;
	ekf_.R_ = R_laser_;
	
	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
//  cout << "x_ = " << ekf_.x_ << endl;
//  cout << "P_ = " << ekf_.P_ << endl;
}
