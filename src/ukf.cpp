#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setIdentity();
  P_ *= 0.5;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  // Radar measurement noise covariance
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
	          0, std_radphi_*std_radphi_, 0,
	          0, 0, std_radrd_*std_radrd_;
  
  // Laser measurement noise covariance
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
	          0, std_laspy_*std_laspy_;

  is_initialized_ = false;
  // CTRV model
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //initialize weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_[0] = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i<2 * n_aug_ + 1; i++) {
	  weights_[i] = 1 / (2 * (lambda_ + n_aug_));
  }

  // Process Noise Covariance Matrix
  Q_ = MatrixXd(2, 2);
  Q_ << std_a_ * std_a_, 0,
	    0, std_yawdd_ * std_yawdd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  // this function was largely taken from Udacity EKF lesson solutions and my EKF solutions

  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
	if (!is_initialized_) {
		/**
		* Initialize the state ukf_.x_ with the first measurement.
		* Create the covariance matrix.
		* Remember: you'll need to convert radar from polar to cartesian coordinates.
		*/
		// first measurement
		cout << "UKF Initialization " << endl;
		x_ << 1, 1, 1, 1, 1;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = meas_package.raw_measurements_[0];
			float phi = meas_package.raw_measurements_[1];
			x_ << rho * cos(phi), rho*sin(phi), 0, 0, 0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
		time_us_ = meas_package.timestamp_;
		return;
	}

	float delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
	time_us_ = meas_package.timestamp_;

	Prediction(delta_t);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLidar(meas_package);
	}

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  This function was largely taken from the UKF lesson solutions.

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  // Augmented Sigma Point Generation
	//create augmented mean vector
	VectorXd x_aug = VectorXd(7);
	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(7, 7);
	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

	//create augmented mean state
	x_aug.fill(0.0);
	x_aug.head(n_x_) = x_;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(n_x_, n_x_) = P_;
	P_aug.bottomRightCorner(2, 2) = Q_;

	//create square root matrix
	MatrixXd A = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0) = x_aug;
	for (int i = 0; i<n_aug_; i++) {
		Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_)*A.col(i);
		Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_)*A.col(i);
	}

	//cout << "x_:" << endl << x_ << endl;
	//cout << "Xsig_aug:" << endl << Xsig_aug << endl;
	//cout << "P_aug:" << endl << P_aug << endl;

  // Sigma Point Prediction
  //create matrix with predicted sigma points as columns
	float dt2 = delta_t * delta_t;
	Xsig_pred_.fill(0.0);
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		float v = Xsig_aug(2, i);
		float psi = Xsig_aug(3, i);
		float psi_dot = Xsig_aug(4, i);
		float cos_psi = cos(psi);
		float sin_psi = sin(psi);

		if (psi_dot > 0.001 || psi_dot < -0.001) {
			Xsig_pred_.col(i) << v / psi_dot * (sin(psi + psi_dot * delta_t) - sin_psi),
				v / psi_dot * (-1 * cos(psi + psi_dot * delta_t) + cos_psi),
				0,
				psi_dot*delta_t,
				0;
		}
		else {
			Xsig_pred_.col(i) << v * cos_psi*delta_t,
				v*sin_psi*delta_t,
				0,
				0,
				0;
		}
		MatrixXd noise(n_x_, 1);
		float nu_a = Xsig_aug(5, i);
		float nu_psi = Xsig_aug(6, i);
		noise << 0.5*dt2*cos_psi*nu_a,
			     0.5*dt2*sin_psi*nu_a,
			     delta_t*nu_a,
			     0.5*dt2*nu_psi,
			     delta_t*nu_psi;
		Xsig_pred_.col(i) += Xsig_aug.col(i).head(n_x_) + noise;
	}
	//cout << "Xsig_pred_:" << endl << Xsig_pred_ << endl;

  // Predict mean and covariance
	x_.fill(0.0);
  //predict state mean
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		x_ += weights_[i] * Xsig_pred_.col(i);
	}
	P_.fill(0.0);
	//predict state covariance matrix
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		MatrixXd error = Xsig_pred_.col(i) - x_;
		NormalizeAngle(error(3));

		P_ += weights_[i] * error*error.transpose();
	}

	cout << "x_:" << endl << x_ << endl;
	cout << "P_:" << endl << P_ << endl;


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Taken largely from UKF lesson solutions

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 2;
	VectorXd z = meas_package.raw_measurements_;
	//cout << "z:" << endl << z << endl;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	//transform sigma points into measurement space
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		Zsig(0, i) = Xsig_pred_(0, i);
		Zsig(1, i) = Xsig_pred_(1, i);
	}

	//calculate mean predicted measurement
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		z_pred += weights_(i)*Zsig.col(i);
	}

	VectorXd z_diff = VectorXd(n_z);
	//calculate innovation covariance matrix S
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		z_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_diff(1));
		S += weights_(i)*z_diff*z_diff.transpose();
	}

	S += R_laser_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

	//calculate cross correlation matrix
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		NormalizeAngle(x_diff(3));
		z_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_diff(1));
		Tc += weights_(i)*x_diff*z_diff.transpose();
	}
	//calculate Kalman gain K;
	MatrixXd S_I = S.inverse();
	MatrixXd K = Tc * S_I;
	//update state mean and covariance matrix
	z_diff = z - z_pred;
	NormalizeAngle(z_diff(1));

	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	float NIS = z_diff.transpose() * S_I * z_diff;
	cout << "Laser NIS: " << NIS << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Taken largely from UKF lesson solutions

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
    //set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;
	VectorXd z = meas_package.raw_measurements_;
	//cout << "z:" << endl << z << endl;

	//create matrix for sigma points in measurement space
	MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z, n_z);
	S.fill(0.0);

	//transform sigma points into measurement space
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		float px = Xsig_pred_(0, i);
		float py = Xsig_pred_(1, i);
		float v = Xsig_pred_(2, i);
		float psi = Xsig_pred_(3, i);
		float r = sqrt(px * px + py * py);
		Zsig(0, i) = r;
		Zsig(1, i) = atan2(py, px);
		Zsig(2, i) = (px*cos(psi)*v + py * sin(psi)*v) / r;
	}
	//calculate mean predicted measurement
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		z_pred += weights_(i)*Zsig.col(i);
	}

	VectorXd z_diff = VectorXd(n_z);
	//calculate innovation covariance matrix S
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		z_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_diff(1));
		S += weights_(i)*z_diff*z_diff.transpose();
	}

	S += R_radar_;

	//create matrix for cross correlation Tc
	MatrixXd Tc = MatrixXd(n_x_, n_z);
	Tc.fill(0.0);

	//calculate cross correlation matrix
	for (int i = 0; i<2 * n_aug_ + 1; i++) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		NormalizeAngle(x_diff(3));
		z_diff = Zsig.col(i) - z_pred;
		NormalizeAngle(z_diff(1));
		Tc += weights_(i)*x_diff*z_diff.transpose();
	}
	//calculate Kalman gain K;
	MatrixXd S_I = S.inverse();
	MatrixXd K = Tc * S_I;
	//update state mean and covariance matrix
	z_diff = z - z_pred;
	NormalizeAngle(z_diff(1));
	x_ += K * z_diff;
	P_ -= K * S * K.transpose();

	float NIS = z_diff.transpose() * S_I * z_diff;
	cout << "Radar NIS: " << NIS << endl;
}

void UKF::NormalizeAngle(double& phi)
{
	// Function recommended by Udacity EKF reviewer and copied from their feedback
	phi = atan2(sin(phi), cos(phi));
}
