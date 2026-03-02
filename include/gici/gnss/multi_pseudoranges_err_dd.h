/**
 * @Function: Multiple double-differenced pseudorange residuals block.
 *            Differencing operator covariance propagation for ceres backend
 *
 * @Author  : Saurav Uprety
 * @Email   :
 *
 * Copyright (C) 2023 by Cheng Chi, All rights reserved.
 **/
#pragma once

#pragma diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
// Eigen 3.2.7 uses std::binder1st and std::binder2nd which are deprecated since
// c++11 Fix is in 3.3 devel
// (http://eigen.tuxfamily.org/bz/show_bug.cgi?id=872).
#include <Eigen/Core>
#include <ceres/ceres.h>
#pragma diagnostic pop

#include "gici/estimate/error_interface.h"
#include "gici/gnss/geodetic_coordinate.h"
#include "gici/gnss/gnss_types.h"

namespace gici {

// Double-differenced pseudorange errors
// The candidate parameter setups are:
// Group 1: P1. receiver position in ECEF (3)
// Group 2: P1. body pose in ENU (7), P2. relative position from body to
// receiver
//          in body frame (3)
// Group 3: Group 1 + P2. troposphere delay at rov (1), P3. troposphere delay at
// ref (1)
//          P4. ionosphere delay (1), P5. ionosphere delay of base satellite (1)
// Group 4: Group 2 + P3. troposphere delay at rov (1), P4. troposphere delay at
// ref (1),
//          P5. ionosphere delay (1), P6. ionosphere delay of base satellite (1)
template <int... Ns>
class MultiPseudorangesErrorDD : public ceres::CostFunction,
                                 public ErrorInterface {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using ParameterDims = ceres::internal::StaticParameterDims<Ns...>;

  /// \brief The base class type.
  typedef ceres::CostFunction base_t;

  /// \brief The information matrix type (MxM).
  typedef Eigen::MatrixXd information_t;

  /// \brief The covariance matrix type (same as information).
  typedef Eigen::MatrixXd covariance_t;

  /// \brief Default constructor.
  MultiPseudorangesErrorDD();

  /// \brief Construct with measurement and information matrix
  /// @param[in] measurement_rov The measurement at rover.
  /// @param[in] measurement_ref The measurement at reference.
  /// @param[in] error_parameter To compute GNSS information matrix.
  MultiPseudorangesErrorDD(const GnssMeasurement &measurement_rov,
                           const GnssMeasurement &measurement_ref,
                           const GnssMeasurementDDIndexPairs &index_pairs,
                           const GnssMeasurementIndex index_rov_base,
                           const GnssMeasurementIndex index_ref_base,
                           const GnssErrorParameter &error_parameter);

  /// \brief Trivial destructor.
  virtual ~MultiPseudorangesErrorDD() {}

  // setters
  /// \brief Set the measurement.
  /// @param[in] measurement The measurement.
  void setMeasurement(const GnssMeasurement &measurement_rov,
                      const GnssMeasurement &measurement_ref) {
    measurement_rov_ = measurement_rov;
    measurement_ref_ = measurement_ref;
  }

  /// \brief Set the information.
  /// @param[in] information The information (weight) matrix.
  void setInformation(const GnssErrorParameter &error_parameter);

  // Set coordinate for ENU to ECEF convertion
  void setCoordinate(const GeoCoordinatePtr &coordinate) {
    coordinate_ = coordinate;
  }

  // error term and Jacobian implementation
  /**
   * @brief This evaluates the error term and additionally computes the
   * Jacobians.
   * @param parameters Pointer to the parameters (see ceres)
   * @param residuals Pointer to the residual vector (see ceres)
   * @param jacobians Pointer to the Jacobians (see ceres)
   * @return success of th evaluation.
   */
  virtual bool Evaluate(double const *const *parameters, double *residuals,
                        double **jacobians) const;

  /**
   * @brief This evaluates the error term and additionally computes
   *        the Jacobians in the minimal internal representation.
   * @param parameters Pointer to the parameters (see ceres)
   * @param residuals Pointer to the residual vector (see ceres)
   * @param jacobians Pointer to the Jacobians (see ceres)
   * @param jacobians_minimal Pointer to the minimal Jacobians (equivalent to
   * jacobians).
   * @return Success of the evaluation.
   */
  bool EvaluateWithMinimalJacobians(double const *const *parameters,
                                    double *residuals, double **jacobians,
                                    double **jacobians_minimal) const;

  // sizes
  /// \brief Residual dimension.
  size_t residualDim() const { return num_residuals(); }

  /// \brief Number of parameter blocks.
  size_t parameterBlocks() const { return dims_.kNumParameterBlocks; }

  /// \brief Dimension of an individual parameter block.
  size_t parameterBlockDim(size_t parameter_block_idx) const {
    return dims_.GetDim(parameter_block_idx);
  }

  /// @brief Residual block type as string
  virtual ErrorType typeInfo() const {
    return ErrorType::kMultiPseudorangesErrorDD;
  }

  // Convert normalized residual to raw residual
  virtual void deNormalizeResidual(double *residuals) const {
    Eigen::Map<Eigen::Matrix<double, 1, 1>> Residual(residuals);
    Residual = whitening_cholesky_factor_ * Residual;
  }

protected:
  GnssMeasurement measurement_rov_, measurement_ref_;
  GnssMeasurementIndex index_rov_base_, index_ref_base_;
  GnssMeasurementDDIndexPairs index_pairs_;

  Satellite satellite_rov_base_, satellite_ref_base_;
  Observation observation_rov_base_, observation_ref_base_;

  // weighting related
  GnssErrorParameter error_parameter_;
  covariance_t covariance_; ///< The DimxDim covariance matrix.
  information_t whitening_cholesky_factor_;



  // Parameter dimensions
  ceres::internal::StaticParameterDims<Ns...> dims_;

  // Geodetic coordinate
  GeoCoordinatePtr coordinate_;

  // parameter types
  bool is_estimate_body_;
  bool is_estimate_atmosphere_;
  int parameter_block_group_;

  char system_;
};

// Explicitly instantiate template classes
template class MultiPseudorangesErrorDD<3>;    // Group 1
template class MultiPseudorangesErrorDD<7, 3>; // Group 2

} // namespace gici