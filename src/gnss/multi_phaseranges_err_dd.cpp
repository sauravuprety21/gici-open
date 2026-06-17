/**
 * @Function: Double-differenced phaserange residual block with differencing
 * operator covariance propagation for ceres backend.
 *
 * @Author  : Saurav Uprety
 * @Email   :
 *
 * Copyright (C) 2023 by Cheng Chi, All rights reserved.
 **/
#include "gici/gnss/multi_phaseranges_err_dd.h"

#include "gici/estimate/pose_local_parameterization.h"
#include "gici/gnss/gnss_common.h"
#include "gici/utility/transform.h"

namespace gici {

// Construct with measurement and information matrix
template <int... Ns>
MultiPhaserangesErrorDD<Ns...>::MultiPhaserangesErrorDD(
        const GnssMeasurement &measurement_rov,
        const GnssMeasurement &measurement_ref,
        const GnssMeasurementDDIndexPairs &index_pairs,
        const GnssMeasurementIndex index_rov_base,
        const GnssMeasurementIndex index_ref_base,
        const GnssErrorParameter &error_parameter) :
    index_pairs_(index_pairs), index_rov_base_(index_rov_base),
    index_ref_base_(index_ref_base) {
  
  CHECK(!checkZero(measurement_ref.position))
  << "The position of reference station is not setted!";

  // Check parameter block types
  // Group 1
  if (dims_.kNumParameterBlocks == 2 && dims_.GetDim(0) == 3 &&
      dims_.GetDim(1) == 1) {
    is_estimate_body_ = false;
    is_estimate_atmosphere_ = false;
    parameter_block_group_ = 1;
  }
  // Group 2
  else if (dims_.kNumParameterBlocks == 3 && dims_.GetDim(0) == 7 &&
           dims_.GetDim(1) == 3 && dims_.GetDim(2) == 1) {
    is_estimate_body_ = true;
    is_estimate_atmosphere_ = false;
    parameter_block_group_ = 2;
  } else {
    LOG(FATAL) << "MultiPhaserangesErrorDD parameter blocks setup invalid!";
  }

  set_num_residuals(index_pairs_.size());

  *mutable_parameter_block_sizes() = std::vector<int32_t>{Ns...};

  // Ambiguity parameters 
  for(size_t i  = 0; i < num_residuals(); i++ ){
    mutable_parameter_block_sizes()->push_back(1);
  }

  setMeasurement(measurement_rov, measurement_ref);

  satellite_rov_base_ = measurement_rov_.getSat(index_rov_base);
  satellite_ref_base_ = measurement_ref_.getSat(index_ref_base);
  observation_rov_base_ = measurement_rov_.getObs(index_rov_base);
  observation_ref_base_ = measurement_ref_.getObs(index_ref_base);

  system_ = satellite_rov_base_.getSystem();

  setInformation(error_parameter);
}

// Set the information.
template <int... Ns>
void MultiPhaserangesErrorDD<Ns ...>::setInformation(
    const GnssErrorParameter &error_parameter) {

  // General error parameters
  error_parameter_ = error_parameter;
  Eigen::Vector3d factor;
  for (size_t i = 0; i < 3; i++)
    factor(i) = error_parameter_.phase_error_factor[i];
  double system_factor =
      square(error_parameter_.system_error_ratio.at(system_));

  // Base satellite
  double elev_rov_base, var_sd_base;

  elev_rov_base = gnss_common::satelliteElevation(
      satellite_rov_base_.sat_position, measurement_rov_.position);
  var_sd_base =
      (square(factor(0)) + square(factor(1) / sin(elev_rov_base))) * 2.0;
  var_sd_base *= system_factor;

  // Covaraince construction
  covariance_.resize(num_residuals(), num_residuals());
  covariance_.setConstant(var_sd_base);

  for (size_t i = 0; i < num_residuals(); i++) {
    // Add diagonal entries for non-base satellites
    double elevation_rov, var_sd;
    Satellite satellite_rov;

    const auto &index_pair = index_pairs_[i];

    satellite_rov = measurement_rov_.getSat(index_pair.rov);

    elevation_rov = gnss_common::satelliteElevation(satellite_rov.sat_position,
                                                    measurement_rov_.position);
    var_sd = (square(factor(0)) + square(factor(1) / sin(elevation_rov))) * 2.0;
    var_sd *= system_factor;
    covariance_(i, i) += var_sd;
  }
  Eigen::LLT<covariance_t> lltOfCovariance(covariance_);
  CHECK(lltOfCovariance.info() == Eigen::Success)
      << "DD multi-phaseranges covariance - LLT decomp. fail";
  whitening_cholesky_factor_ = lltOfCovariance.matrixL();
}

template <int... Ns>
bool MultiPhaserangesErrorDD<Ns ...>::Evaluate(double const *const *parameters,
                                           double *residuals,
                                           double **jacobians) const {
  return EvaluateWithMinimalJacobians(parameters, residuals, jacobians,
                                      nullptr);
}

// This evaluates the error term and additionally computes
// the Jacobians in the minimal internal representation.
template <int... Ns>
bool MultiPhaserangesErrorDD<Ns...>::EvaluateWithMinimalJacobians(
    double const *const *parameters, double *residuals, double **jacobians,
    double **jacobians_minimal) const {
  Eigen::Vector3d t_WR_ECEF, t_WS_W, t_SR_S;
  Eigen::Quaterniond q_WS;

  // Position and ambiguities
  if (!is_estimate_body_) {
    t_WR_ECEF = Eigen::Map<const Eigen::Vector3d>(parameters[0]);
  } else {
    // pose in ENU frame
    t_WS_W = Eigen::Map<const Eigen::Vector3d>(&parameters[0][0]);
    q_WS = Eigen::Map<const Eigen::Quaterniond>(&parameters[0][3]);

    // relative position
    t_SR_S = Eigen::Map<const Eigen::Vector3d>(parameters[1]);

    // receiver position
    Eigen::Vector3d t_WR_W = t_WS_W + q_WS * t_SR_S;

    if (!coordinate_) {
      LOG(FATAL) << "Coordinate not set!";
    }
    if (!coordinate_->isZeroSetted()) {
      LOG(FATAL) << "Coordinate zero not set!";
    }
    t_WR_ECEF = coordinate_->convert(t_WR_W, GeoType::ENU, GeoType::ECEF);
  }

  // base satellite ambiguity
  double dambiguity_base = parameters[dims_.kNumParameterBlocks - 1][0];
  
  /* Jacobians */
  Eigen::Matrix<double, Eigen::Dynamic, 3> J_t_ECEF =
      Eigen::Matrix<double, Eigen::Dynamic, 3>::Zero(num_residuals(), 3);

  Eigen::Matrix<double, Eigen::Dynamic, 6> J_T_WS =
      Eigen::Matrix<double, Eigen::Dynamic, 6>::Zero(num_residuals(), 6);

  Eigen::Matrix<double, Eigen::Dynamic, 6> J0_minimal =
      Eigen::Matrix<double, Eigen::Dynamic, 6>::Zero(num_residuals(), 6);

  Eigen::Matrix<double, Eigen::Dynamic, 3> J_t_SR_S =
      Eigen::Matrix<double, Eigen::Dynamic, 3>::Zero(num_residuals(), 3);

  Eigen::Matrix<double, 6, 7, Eigen::RowMajor> J_lift =
      Eigen::Matrix<double, 6, 7, Eigen::RowMajor>::Zero();

  Eigen::VectorXd J_dambiguity_base = Eigen::VectorXd::Constant(num_residuals(), -1.0);

  Eigen::MatrixXd J_dambiguities(num_residuals(), num_residuals());
  // Rest satellites ambiguities
  J_dambiguities =
      whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
          Eigen::MatrixXd::Identity(num_residuals(), num_residuals()));

  /* Residuals */
  Eigen::VectorXd Residuals = Eigen::VectorXd::Zero(num_residuals());

  double rho_rov_base = gnss_common::satelliteToReceiverDistance(
      satellite_rov_base_.sat_position, t_WR_ECEF);
  double rho_ref_base = gnss_common::satelliteToReceiverDistance(
      satellite_ref_base_.sat_position, measurement_ref_.position);

  for (size_t i = 0; i < num_residuals(); i++) {
    const auto &index_pair = index_pairs_[i];

    const Satellite &satellite_rov = measurement_rov_.getSat(index_pair.rov);
    const Satellite &satellite_ref = measurement_ref_.getSat(index_pair.ref);

    const Observation &observation_rov =
        measurement_rov_.getObs(index_pair.rov);
    const Observation &observation_ref =
        measurement_ref_.getObs(index_pair.ref);

    double rho_rov, rho_ref;
    rho_rov = gnss_common::satelliteToReceiverDistance(
        satellite_rov.sat_position, t_WR_ECEF);
    rho_ref = gnss_common::satelliteToReceiverDistance(
        satellite_ref.sat_position, measurement_ref_.position);
    
    if (jacobians != nullptr) {
      // Receiver position in ECEF
      if (jacobians[0] != nullptr) {
        J_t_ECEF.block(i, 0, 1, 3) =
            (-((t_WR_ECEF - satellite_rov.sat_position) / rho_rov) +
             (t_WR_ECEF - satellite_rov_base_.sat_position) / rho_rov_base)
                .transpose();
      }
      // Poses
      if (is_estimate_body_) {
        if (jacobians[0] != nullptr) {
          // Body position in ENU
          J_T_WS.block(i, 0, 1, 3) =
              J_t_ECEF.block(i, 0, 1, 3) *
              coordinate_->rotationMatrix(GeoType::ENU, GeoType::ECEF);
          // Body rotation in ENU
          J_T_WS.block(i, 3, 1, 3) =
              J_T_WS.block(i, 0, 1, 3) *
              -skewSymmetric(q_WS.toRotationMatrix() * t_SR_S);
        }

        if (jacobians[1] != nullptr) {
          // Relative position
          J_t_SR_S.block(i, 0, 1, 3) =
              J_T_WS.block(i, 0, 1, 3) * q_WS.toRotationMatrix();
        }
      }
    }
    double dambiguity = is_estimate_body_ ? parameters[3 + i][0] : parameters[2 + i][0];
    double dphaserange_estimate = rho_rov - rho_ref - rho_rov_base +
                                  rho_ref_base + dambiguity - dambiguity_base;

    double dphaserange =
        observation_rov.phaserange - observation_ref.phaserange -
        observation_rov_base_.phaserange + observation_ref_base_.phaserange;

    Residuals(i) = dphaserange - dphaserange_estimate;
  }

  Eigen::Map<Eigen::VectorXd> weighted_residuals(residuals, num_residuals());
  weighted_residuals =
      whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
          Residuals);

  // Group 1
  if (parameter_block_group_ == 1 && jacobians != nullptr) {
    // Position
    if (jacobians[0] != nullptr) {
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> J0(
          jacobians[0], num_residuals(), 3);
      J0 = whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
          J_t_ECEF);
    }

    // Base satellite ambiguity
    if (jacobians[1] != nullptr) {
      Eigen::Map<Eigen::VectorXd> J1(jacobians[1], num_residuals());
      J1 = whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
              J_dambiguity_base);
    }
    
    for(size_t i = 0; i < num_residuals(); i++){
      if (jacobians[2 + i] != nullptr){
        Eigen::Map<Eigen::VectorXd> J2(jacobians[2 + i], num_residuals());
        J2 = J_dambiguities.block(0, i, num_residuals(), 1).eval();
      }
    }
  }

  // Group 2
  if (parameter_block_group_ == 2 && jacobians != nullptr) {
    // Pose
    if (jacobians[0] != nullptr ) {
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 7, Eigen::RowMajor>> J0(
          jacobians[0], num_residuals(), 7);

      J0_minimal =
          whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
              J_T_WS);

      // pseudo inverse of the local parametrization Jacobian:
      PoseLocalParameterization::liftJacobian(parameters[0], J_lift.data());

      J0 = J0_minimal * J_lift;
    }

    // Relative position
    if (jacobians[1] != nullptr) {
      Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>> J1(
          jacobians[1], num_residuals(), 3);
      J1 = whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
          J_t_SR_S);
    }
    // Base satellite ambiguity
    if (jacobians[2] != nullptr) {
      Eigen::Map<Eigen::VectorXd> J2(jacobians[2], num_residuals());
      J2 =
          whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
              J_dambiguity_base);
    }
    // Rest satellites ambiguities
    for (size_t i = 0; i < num_residuals(); i++) {
      if (jacobians[3 + i] != nullptr) {
        Eigen::Map<Eigen::VectorXd> J3(jacobians[3 + i],
                                       num_residuals());
        J3 = J_dambiguities.block(0, i, num_residuals(), 1).eval();
      }
    }
  }

  return true;
}
} // namespace gici
