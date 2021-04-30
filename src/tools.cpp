#include "tools.h"
#include <iostream>

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  unsigned int estimation_size = estimations.size();
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimation_size != ground_truth.size()
      || estimation_size == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }
  
  //accumulate squared residuals
  for (unsigned int i=0; i < estimation_size; ++i) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  
// calculate the mean
rmse = rmse/estimation_size;
//calculate the squared root  
rmse = rmse.array().sqrt();
  
return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  //To avoid divide by zero.
  if (px == 0 && py == 0)
  {
      cout<<"Divion by zero";
  } else{
  //precompute few terms to avoid repeatations of calculations.
    float denominator = (px * px) + (py * py);
    float sqrt_denominator = sqrt(denominator);
    float denominator_three_by_two = denominator * sqrt_denominator;
    float der_x = px / sqrt_denominator;
    float der_y = py / sqrt_denominator;
    float der_rho_x = (px * ((vy*px) - (vx*py)))/(denominator_three_by_two);
    float der_rho_y = (py * ((vx*py) - (vy*px)))/(denominator_three_by_two);
    
    //Calculating Jacobian Matric
    Hj << der_x , der_y, 0, 0,
           -py/denominator, px/denominator,0,0,
           der_rho_y, der_rho_x, der_x, der_y;
  }
  
  return Hj;
}
