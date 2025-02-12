/**
* @Function: Pure INS navigation without external sensor aiding
*
* @Author  : Saurav
* @Email   : 
*
**/
#pragma once

#include "gici/imu/imu_estimator_base.h"
#include "gici/gnss/gnss_types.h"

namespace gici {

struct InsEstimatorOptions {
	// Max window length
	size_t max_window_length = 10;
};

// IMU initialization options
struct InsInitializerOptions {
	// Initial pose time (seconds) (gps_t)
	double t0 = -1;

	// Initial attitude (deg) (Euler ZYX rotation order)
	Eigen::Vector3d rpy_0 = Eigen::Vector3d::Zero();

	// Initial velocity (m/s) (ENU)
	Eigen::Vector3d vel_0 = Eigen::Vector3d::Zero();
	
	// Initial position (deg/m) (lat/long/height)
	Eigen::Vector3d lla_0 = Eigen::Vector3d::Zero();
	
	// Initial accelerometer biases (m/s^2) (IMU frame: x(up), y(right), z(forward))
	Eigen::Vector3d accel_bias_0 = Eigen::Vector3d::Zero();
	
	// Initial gyroscope biases (rad/s) (IMU frame: x(up), y(right), z(forward))
	Eigen::Vector3d gyro_bias_0 = Eigen::Vector3d::Zero();
};


class InsEstimator:	public ImuEstimatorBase {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	InsEstimator(const InsEstimatorOptions& options,
			const InsInitializerOptions& init_options, 
			const ImuEstimatorBaseOptions& imu_base_options,
			const EstimatorBaseOptions& base_options);

	~InsEstimator();

	// Add measurement
	bool addMeasurement(const EstimatorDataCluster& measurement) override;

	// Estimate current graph
  bool estimate() override;



protected:
	// Add measurement and state
	bool addState(const ImuMeasurement& imu_measurement);

	// Marginalization
  bool marginalization();
	
	// Add initialization from YAML file
	void setInitialization(const InsInitializerOptions &init_options); 


	// Shift memory for states and measurements
  inline void shiftMemory() {
    states_.push_back(State());
    while (states_.size() > ins_options_.max_window_length) {
      states_.pop_front();
    }
  } 
	// Options 
	InsEstimatorOptions ins_options_;
	InsInitializerOptions init_options_;
	
	// Last update to Ins (NHC/HMC/&/or/ZUPT)
	double last_update_;
	double update_intrvl_;
	int32_t update_cnt_;

};

}