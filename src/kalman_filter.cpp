#include "kalman_filter.h"

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
    * Predict the state and covariance matrix
    */
    x_ = F_ * x_ ;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
    * update the state by using Kalman Filter equations
    */
    VectorXd y = z - H_ * x_;
    Update_(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
    * update the state by using Extended Kalman Filter equations
    */
    double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
    double theta = atan2(x_(1),  x_(0));
    double rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;

    VectorXd h = VectorXd(3); // h(x_)
    h << rho, theta, rho_dot;

    VectorXd y = z - h;
    y(1) = NormalizeAngle_(y(1));

    Update_(y);
}

void KalmanFilter::Update_(const VectorXd &y){
    // Precompute matrices involved
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;

    // Update the state
    x_ = x_ + (K * y);
    int x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

#define PI 3.141592653589793
double KalmanFilter::NormalizeAngle_(double theta){

    if(theta > PI) theta -= 2*PI;
    else if(theta < -PI) theta += 2*PI;

    return theta;
}
