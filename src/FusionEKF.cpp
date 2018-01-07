#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPS 0.0000001

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

    Hj_ = MatrixXd(3, 4);

    // measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    // measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

    // The initial value of the process covariance.
    // In the first frame, we are more certain of the positions than we are
    // about velocities. So we initialize the velocity covariances to high values
    P_initial_ = MatrixXd(4, 4);
    P_initial_ << 1, 0, 0   , 0,
                  0, 1, 0   , 0,
                  0, 0, 1000, 0,
                  0, 0, 0, 1000;

    noise_ax_ = 9.0;
    noise_ay_ = 9.0;
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
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        ekf_.x_ = VectorXd(4);

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             * Convert radar from polar to cartesian coordinates and initialize state.
             */
            float rho = measurement_pack.raw_measurements_[0]; // range
            float phi = measurement_pack.raw_measurements_[1]; // bearing
            float rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho
            // Coordinates convertion from polar to cartesian
            float x = rho * cos(phi);
            float y = rho * sin(phi);
            float vx = rho_dot * cos(phi);
            float vy = rho_dot * sin(phi);
            ekf_.x_ << x, y, vx , vy;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             * Initialize state.
             */
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        }
        // Deal with the special case initialisation problems
        if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){
            ekf_.x_(0) = EPS;
            ekf_.x_(1) = EPS;
        }
        // Initial covariance matrix
        ekf_.P_ = P_initial_;

        // Initialize the state transition matrix. We initialize dt to 1 and update it later
        ekf_.F_ = MatrixXd(4, 4);
        ekf_.F_ << 1, 0, 1, 0,
                   0, 1, 0, 1,
                   0, 0, 1, 0,
                   0, 0, 0, 1;

        // Save the initiall timestamp for dt calculation
        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_);
    dt /= 1000000.0; // convert micros to s
    previous_timestamp_ = measurement_pack.timestamp_;

    // State transition matrix update
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;

    // Noise covariance matrix computation
    float dt_2 = dt * dt; //dt^2
    float dt_3 = dt_2 * dt; //dt^3
    float dt_4 = dt_3 * dt; //dt^4

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4 * noise_ax_/4, 0, dt_3 * noise_ax_/2, 0,
               0, dt_4 * noise_ay_/4, 0, dt_3 * noise_ay_/2,
               dt_3 * noise_ax_/2, 0, dt_2 * noise_ax_, 0,
               0, dt_3 * noise_ay_/2, 0, dt_2 * noise_ay_;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        // Use Jacobian instead of H
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
    //cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
}

