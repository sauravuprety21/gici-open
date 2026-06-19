/**
 * @Function: Multiple pseudorange residuals block for ceres backend
 *
 * @Author  : Saurav Uprety
 * @Email   : 
 *
 * Copyright (C) 2023 by Cheng Chi, All rights reserved.
 **/
#include "gici/gnss/multi_pseudoranges_err.h"

#include "gici/estimate/pose_local_parameterization.h"
#include "gici/gnss/code_phase_maps.h"
#include "gici/gnss/gnss_common.h"
#include "gici/utility/global_variable.h"
#include "gici/utility/transform.h"

namespace gici {

// Construct with measurement and information matrix
template <int... Ns>
MultiPseudorangesError<Ns ...>::MultiPseudorangesError(
    const GnssMeasurement &measurement, const GnssMeasurementsIndexes &indexes,
    const GnssErrorParameter &error_parameter)
    : indexes_(indexes) {

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
    LOG(FATAL) << "MultiPseudorangesError parameter blocks setup invalid!";
  }

  set_num_residuals(indexes_.size());
  *mutable_parameter_block_sizes() = std::vector<int32_t>{Ns...};

  setMeasurement(measurement);

  system_ = measurement_.getSat(indexes_.front()).getSystem();

  if ((system_ == 'G') &&
      checkEqual(measurement_.getObs(indexes_.front()).wavelength,
                 CLIGHT / gnss_common::phaseToFrequency('G', PHASE_L5)))
    correct_ifcb_ = true;

  setInformation(error_parameter);
}

template <int... Ns>
void MultiPseudorangesError<Ns ...>::setInformation(
    const GnssErrorParameter &error_parameter) {
  // General error parameters
  error_parameter_ = error_parameter;
  Eigen::Vector3d factor;
  for (size_t i = 0; i < 3; i++)
    factor(i) = error_parameter_.phase_error_factor[i];
  double ratio = square(error_parameter_.code_to_phase_ratio);
  double system_factor =
      square(error_parameter_.system_error_ratio.at(system_));

  // Covaraince construction
  covariance_.resize(num_residuals(), num_residuals());
  if (!correct_ifcb_)
    covariance_.setConstant(0.0);
  else
    // add IFCB residual error for GPS L5
    covariance_.setConstant(square(error_parameter_.residual_gps_ifcb));

  double timestamp = measurement_.timestamp;

  for (size_t i = 0; i < num_residuals(); i++) {
    const auto &index = indexes_[i];
    Satellite satellite = measurement_.getSat(index);
    Observation observation = measurement_.getObs(index);

    double elevation = gnss_common::satelliteElevation(satellite.sat_position,
                                                       measurement_.position);
    double azimuth = gnss_common::satelliteAzimuth(satellite.sat_position,
                                                   measurement_.position);
    double ephemeris_var, ionosphere_var, troposphere_var;
    // ephemeris error
    if (satellite.sat_type == SatEphType::Broadcast) {
      ephemeris_var = square(error_parameter_.ephemeris_broadcast);
    } else if (satellite.sat_type == SatEphType::Precise) {
      ephemeris_var = square(error_parameter_.ephemeris_precise);
    }
    // ionosphere error
    if (satellite.ionosphere != 0.0 &&
        satellite.ionosphere_type == IonoType::Augmentation) {
      ionosphere_var = square(error_parameter_.ionosphere_augment);
    } else if (satellite.ionosphere != 0.0 &&
               satellite.ionosphere_type == IonoType::DualFrequency) {
      ionosphere_var = square(error_parameter_.ionosphere_dual_frequency);
    } else {
      double ionosphere_delay = gnss_common::ionosphereBroadcast(
          timestamp, measurement_.position, azimuth, elevation,
          observation.wavelength, measurement_.ionosphere_parameters);
      ionosphere_var = square(error_parameter_.ionosphere_broadcast_factor *
                              ionosphere_delay);
    }
    // troposphere error
    double troposphere_delay = gnss_common::troposphereSaastamoinen(
        timestamp, measurement_.position, elevation);
    troposphere_var =
        square(error_parameter_.troposphere_model_factor * troposphere_delay);
    // troposphere wet delay
    if (measurement_.troposphere_wet != 0.0) {
      troposphere_var = square(error_parameter_.troposphere_augment);
    }

    covariance_(i, i) =
        (square(factor(0)) + square(factor(1) / sin(elevation))) * ratio +
        ephemeris_var + ionosphere_var + troposphere_var;
    covariance_(i, i) *= system_factor;
  }
  Eigen::LLT<covariance_t> lltOfCovariance(covariance_);
  CHECK(lltOfCovariance.info() == Eigen::Success)
      << "Multi-pseudoranges covariance - LLT decomp. fail";
  whitening_cholesky_factor_ = lltOfCovariance.matrixL();
}

template <int... Ns>
bool MultiPseudorangesError<Ns ...>::Evaluate(double const *const *parameters,
                                             double *residuals,
                                             double **jacobians) const {
  return EvaluateWithMinimalJacobians(parameters, residuals, jacobians,
                                      nullptr);
}

template <int... Ns>
bool MultiPseudorangesError<Ns ...>::EvaluateWithMinimalJacobians(
    double const *const *parameters, double *residuals, double **jacobians,
    double **jacobians_minimal) const {
  Eigen::Vector3d t_WR_ECEF, t_WS_W, t_SR_S;
  Eigen::Quaterniond q_WS;
  double clock, ifb;

  // Position and clock
  if (!is_estimate_body_) {
    t_WR_ECEF = Eigen::Map<const Eigen::Vector3d>(parameters[0]);
    clock = parameters[1][0];
  } else {
    // pose in ENU frame
    t_WS_W = Eigen::Map<const Eigen::Vector3d>(&parameters[0][0]);
    q_WS = Eigen::Map<const Eigen::Quaterniond>(&parameters[0][3]);

    // relative position
    t_SR_S = Eigen::Map<const Eigen::Vector3d>(parameters[1]);

    // clock
    clock = parameters[2][0];

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

  // Earth tide
  double timestamp = measurement_.timestamp;
  Eigen::Vector3d tide = gnss_common::solidEarthTide(timestamp, t_WR_ECEF);
  t_WR_ECEF += tide;

  // ECEF position
  Eigen::Matrix<double, Eigen::Dynamic, 3> J_t_ECEF =
      Eigen::Matrix<double, Eigen::Dynamic, 3>::Zero(num_residuals(), 3);
  // Pose
  Eigen::Matrix<double, Eigen::Dynamic, 6> J_T_WS =
      Eigen::Matrix<double, Eigen::Dynamic, 6>::Zero(num_residuals(), 6);
  // Intermediate
  Eigen::Matrix<double, Eigen::Dynamic, 6> J0_minimal =
      Eigen::Matrix<double, Eigen::Dynamic, 6>::Zero(num_residuals(), 6);
  // Antenna position in body
  Eigen::Matrix<double, Eigen::Dynamic, 3> J_t_SR_S =
      Eigen::Matrix<double, Eigen::Dynamic, 3>::Zero(num_residuals(), 3);
  // Manifold lift
  Eigen::Matrix<double, 6, 7, Eigen::RowMajor> J_lift =
      Eigen::Matrix<double, 6, 7, Eigen::RowMajor>::Zero();
  // Clock
  Eigen::VectorXd J_clock = -Eigen::VectorXd::Constant(num_residuals(), 1.0);
  // Residuals
  Eigen::VectorXd Residuals = Eigen::VectorXd::Zero(num_residuals());

  for (size_t i = 0; i < num_residuals(); i++) {

    double troposphere_delay, ionosphere_delay;
    double gmf_wet, gmf_hydro;

    const auto &index = indexes_[i];
    const Satellite &satellite = measurement_.getSat(index);
    const Observation &observation = measurement_.getObs(index);

    double rho = gnss_common::satelliteToReceiverDistance(
        satellite.sat_position, t_WR_ECEF);
    double elevation =
        gnss_common::satelliteElevation(satellite.sat_position, t_WR_ECEF);
    double azimuth =
        gnss_common::satelliteAzimuth(satellite.sat_position, t_WR_ECEF);

    // Atmosphere
    if (!is_estimate_atmosphere_) {
      // troposphere hydro-static delay
      troposphere_delay =
          gnss_common::troposphereSaastamoinen(timestamp, t_WR_ECEF, elevation);
      // troposphere wet delay
      if (measurement_.troposphere_wet != 0.0) {
        gnss_common::troposphereGMF(timestamp, t_WR_ECEF, elevation, nullptr,
                                    &gmf_wet);
        troposphere_delay += measurement_.troposphere_wet * gmf_wet;
      }
      // ionosphere
      if (satellite.ionosphere != 0.0 &&
          satellite.ionosphere_type == IonoType::Augmentation) {
        ionosphere_delay = satellite.ionosphere;
        ionosphere_delay = gnss_common::ionosphereConvertFromBase(
            ionosphere_delay, observation.wavelength);
      } else if (satellite.ionosphere != 0.0 &&
                 satellite.ionosphere_type == IonoType::DualFrequency) {
        ionosphere_delay = satellite.ionosphere;
        ionosphere_delay = gnss_common::ionosphereConvertFromBase(
            ionosphere_delay, observation.wavelength);
      } else {
        ionosphere_delay = gnss_common::ionosphereBroadcast(
            timestamp, t_WR_ECEF, azimuth, elevation, observation.wavelength,
            measurement_.ionosphere_parameters);
      }
    }

    // Get estimate derivated measurement
    double pseudorange_estimate = rho + clock - satellite.sat_clock +
                                  troposphere_delay + ionosphere_delay;
    double pseudorange = observation.pseudorange;

    Residuals(i) = pseudorange - pseudorange_estimate;

    if (jacobians != nullptr) {
      // Receiver position in ECEF
      if (jacobians[0] != nullptr) {
        J_t_ECEF.block(i, 0, 1, 3) =
            -((t_WR_ECEF - satellite.sat_position) / rho).transpose();
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

    Eigen::Map<Eigen::VectorXd> weighted_residuals(residuals, num_residuals());
    weighted_residuals =
        whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
            Residuals);

    // Group 1
    if (parameter_block_group_ == 1) {
      // Position
      if (jacobians != nullptr && jacobians[0] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>
            J0(jacobians[0], num_residuals(), 3);
        J0 = whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
            J_t_ECEF);
      }
      if (jacobians_minimal != nullptr && jacobians_minimal[0] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>
            J0_minimal_mapped(jacobians_minimal[0], num_residuals(), 3);
        J0_minimal_mapped =
            whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
                J_t_ECEF);
      }
      // Clock
      if (jacobians != nullptr && jacobians[1] != nullptr) {
        Eigen::Map<Eigen::VectorXd> J1(jacobians[1], num_residuals());
        J1 = whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
            J_clock);
      }
      if (jacobians_minimal != nullptr && jacobians_minimal[1] != nullptr) {
        Eigen::Map<Eigen::VectorXd> J1_minimal_mapped(jacobians_minimal[1],
                                                      num_residuals());
        J1_minimal_mapped =
            whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
                J_clock);
      }
    }
    // Group 2
    if (parameter_block_group_ == 2) {
      J0_minimal =
          whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
              J_T_WS);

      // Pose
      if (jacobians != nullptr && jacobians[0] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 7, Eigen::RowMajor>>
            J0(jacobians[0], num_residuals(), 7);

        // pseudo inverse of the local parametrization Jacobian:
        PoseLocalParameterization::liftJacobian(parameters[0], J_lift.data());

        J0 = J0_minimal * J_lift;
      }
      if (jacobians_minimal != nullptr && jacobians_minimal[0] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 6, Eigen::RowMajor>>
            J0_minimal_mapped(jacobians_minimal[0], num_residuals(), 6);
        J0_minimal_mapped = J0_minimal;
      }

      // Relative position
      if (jacobians != nullptr && jacobians[1] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>
            J1(jacobians[1], num_residuals(), 3);
        J1 = whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
            J_t_SR_S);
      }
      if (jacobians_minimal != nullptr && jacobians_minimal[1] != nullptr) {
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>
            J1_minimal_mapped(jacobians_minimal[1], num_residuals(), 3);
        J1_minimal_mapped =
            whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
                J_t_SR_S);
      }
      // Clock
      if (jacobians != nullptr && jacobians[2] != nullptr) {
        Eigen::Map<Eigen::VectorXd> J2(jacobians[2], num_residuals());
        J2 = whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
            J_clock);
      }
      if (jacobians_minimal != nullptr && jacobians_minimal[2] != nullptr) {
        Eigen::Map<Eigen::VectorXd> J2_minimal_mapped(jacobians_minimal[2],
                                                      num_residuals());
        J2_minimal_mapped =
            whitening_cholesky_factor_.triangularView<Eigen::Lower>().solve(
                J_clock);
      }
    }
  }
  return true;
}
}
