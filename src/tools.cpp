#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the inputs:
	if (estimations.size()==0) {        //  * the estimation vector size should not be zero
	    cout << "CalculateRMSE() - ERROR - zero length estimation vector";
	}
	if (estimations.size() != ground_truth.size()) {//  * the estimation vector size should equal ground truth vector size
	    cout << "CalculateRMSE() - ERROR - estimation and ground_truth vectors are different size";
	}

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
	    VectorXd residual = estimations[i]-ground_truth[i];
        residual = residual.array()*residual.array();
        rmse += residual;
	}

	//calculate the mean
    rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

    //check division by zero
	if (abs(px) < 0.01) {
	  px = 0.01;
	}
	if (abs(py) < 0.01) {
	  py = 0.01;
	}
	//compute the Jacobian matrix
      Hj(0,0) = px / (sqrt(pow(px,2)+pow(py,2)));
      Hj(0,1) = py / (sqrt(pow(px,2)+pow(py,2)));
      Hj(0,2) = 0.;
      Hj(0,3) = 0.;
    
      Hj(1,0) = -py / (pow(px,2)+pow(py,2));
      Hj(1,1) = px / (pow(px,2)+pow(py,2));
      Hj(1,2) = 0;
      Hj(1,3) = 0;
    
      Hj(2,0) = py*(vx*py-vy*px) / pow((pow(px,2)+pow(py,2)), 3/2);
      Hj(2,1) = px*(vy*px-vx*py) / pow((pow(px,2)+pow(py,2)), 3/2);
      Hj(2,2) = px / (sqrt(pow(px,2)+pow(py,2)));
      Hj(2,3) = py / (sqrt(pow(px,2)+pow(py,2)));
	
	return Hj;
}
