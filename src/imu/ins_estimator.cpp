/**
* @Function: Pure INS navigation without external sensor aiding
*
* @Author  : Saurav
* @Email   : 
*
* 
**/
#include "gici/estimate/estimator_base.h"
#include "gici/gnss/geodetic_coordinate.h"
#include "gici/imu/imu_estimator_base.h"
#include "gici/imu/ins_estimator.h"

#include "gici/estimate/speed_and_bias_parameter_block.h"
#include "gici/imu/imu_error.h"
#include "gici/imu/speed_and_bias_error.h"
#include "gici/estimate/pose_error.h"
#include "gici/imu/yaw_error.h"
#include "gici/imu/roll_and_pitch_error.h"
#include "gici/imu/hmc_error.h"
#include "gici/imu/nhc_error.h"
#include "gici/utility/transform.h"

namespace gici {

// default constructor
InsEstimator::InsEstimator(const InsEstimatorOptions& options,
      const InsInitializerOptions &init_options, 
			const ImuEstimatorBaseOptions& imu_base_options,
			const EstimatorBaseOptions& base_options):
	last_update_(-1.0),
	update_intrvl_(1.0),
	update_cnt_(0),
	ins_options_(options), 
  init_options_(init_options),
	ImuEstimatorBase(imu_base_options, base_options),
	EstimatorBase(base_options)
{
	type_ = EstimatorType::Ins;
	shiftMemory();
  setInitialization(init_options_);
}

// defautl destructor
InsEstimator::~InsEstimator() {}		

bool InsEstimator::addMeasurement(const EstimatorDataCluster& measurement) {
	// Initialization
  if (coordinate_ == nullptr || !gravity_setted_) return false;

  double timestamp = measurement.imu->timestamp;
  
  if (measurement.imu && measurement.imu_role == ImuRole::Major) {
    addImuMeasurement(*measurement.imu);
    // if(last_update_==-1){ 
    //   last_update_=timestamp;
    // } else if ((timestamp - last_update_) > update_intrvl_) {
    //     last_update_ = timestamp;
    //     return addState(*measurement.imu);
    // }
  }
  
	return true;
}


bool InsEstimator::addState(const ImuMeasurement& imu_measurement){

	double timestamp = imu_measurement.timestamp;

	const int32_t bundle_id = update_cnt_++;
	BackendId pose_id = createGnssPoseId(bundle_id);
	size_t index = insertImuState(timestamp, pose_id);
	CHECK(index == states_.size() - 1);
	curState().status = GnssSolutionStatus::DeadReckoning;
	// ZUPT
 // addZUPTResidualBlock(curState());
	
  // Car motion
  if (imu_base_options_.car_motion) {
    // heading measurement constraint
    //addHMCResidualBlock(curState());
    // non-holonomic constraint
    //addNHCResidualBlock(curState());
  }
	return true;
}

bool InsEstimator::estimate(){

	optimize();

	// Log information
  if (base_options_.verbose_output) {
    LOG(INFO) << estimatorTypeToString(type_) << ": " 
      << "Iterations: " << graph_->summary.iterations.size() << ", "
      << std::scientific << std::setprecision(3) 
      << "Initial cost: " << graph_->summary.initial_cost << ", "
      << "Final cost: " << graph_->summary.final_cost;
  }

	// Apply marginalization
  bool ret_ = marginalization();

  // Shift memory for states and measurements
  shiftMemory();

  return ret_;

}

// Marginalization
bool InsEstimator::marginalization()
{
  // Check if we need marginalization
  if (states_.size() < ins_options_.max_window_length) {
    return true;
  }

  // Erase old marginalization item
  if (!eraseOldMarginalization()) return false;

  // Add marginalization items
  // IMU states and residuals
  addImuStateMarginBlockWithResiduals(oldestState());

  // Apply marginalization and add the item into graph
  return applyMarginalization();
}


void InsEstimator::setInitialization(const InsInitializerOptions &init_options){

  Eigen::Vector3d lla = init_options_.lla_0.eval();
  Eigen::Vector3d lla_rad = GeoCoordinate::degToRad(lla);
  setGravity(earthGravity(lla_rad));
  auto llaPtr = std::make_shared<GeoCoordinate>(lla_rad,
            GeoType::LLA);

  setCoordinate(llaPtr); 

  // Clear old graph
  states_.clear();
  Graph::ParameterBlockCollection parameters = graph_->parameters();
  for (auto& parameter : parameters) {
    graph_->removeParameterBlock(parameter.first);
  }

  double timestamp = init_options.t0;

  SpeedAndBias speed_and_bias_0;
  
  speed_and_bias_0.head<3>() = init_options.vel_0;
  speed_and_bias_0.segment<3>(3) = init_options.gyro_bias_0;
  speed_and_bias_0.segment<3>(6) = init_options.accel_bias_0;

  // Add new pose
  Eigen::Quaterniond q_WS = eulerAngleToQuaternion(init_options.rpy_0);

  Transformation T_WS = Transformation(Eigen::Vector3d::Zero(), q_WS);
  

  // pose prior
  Eigen::Vector3d cur_position = coordinate_->convert(
        lla_rad, GeoType::LLA, GeoType::ENU);

  T_WS.getPosition() = cur_position;

  const int32_t bundle_id = update_cnt_++;
  BackendId pose_id = createGnssPoseId(bundle_id);
  size_t index = insertImuState(timestamp, pose_id, T_WS, speed_and_bias_0, true);
  states_[index].status = GnssSolutionStatus::Fixed;

  // Initial errors
    if (isFirstEpoch()) {
      // bias error
      double speed_std = speed_and_bias_0.head<3>().norm() * 2.0;
      if (speed_std < 1.0) speed_std = 1.0;
      addImuSpeedAndBiasResidualBlock(states_[index], speed_and_bias_0, 
        speed_std,  
        imu_base_options_.imu_parameters.sigma_bg,
        imu_base_options_.imu_parameters.sigma_ba);
      // yaw error
      if (imu_base_options_.car_motion) {
        double std_yaw = sqrt(square(imu_base_options_.body_to_imu_rotation_std * D2R) + 
          square(0.1 / speed_and_bias_0.head<2>().norm()));
        addPoseResidualBlock(states_[index], T_WS, 0.1, std_yaw);
      }
    }

    // Car motion
    if (imu_base_options_.car_motion) {
      if (!isFirstEpoch()) addHMCResidualBlock(states_[index]);
      addNHCResidualBlock(states_[index]);
    } 
  return;
}
}