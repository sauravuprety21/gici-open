/**
* @Function: Align differential measurements
*
* @Author  : Cheng Chi
* @Email   : chichengcn@sjtu.edu.cn
*
* Copyright (C) 2023 by Cheng Chi, All rights reserved.
**/
#pragma once

#include <deque>
#include <iomanip>
#include <sstream>

#include "gici/gnss/gnss_types.h"
#include "gici/estimate/estimator_types.h"
#include "gici/utility/common.h"

namespace gici {

// Estimator
class DifferentialMeasurementsAlign {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // The default constructor
  DifferentialMeasurementsAlign() {}
  ~DifferentialMeasurementsAlign() {}

  // Set measurements
  inline void add(const EstimatorDataCluster& measurement) {
    if (measurement.gnss && measurement.gnss_role == GnssRole::Rover) {
      measurement_rov_.push_back(*measurement.gnss);
      // std::ostringstream stream;
      // stream << std::setprecision(17);
      // stream << "DGNSS align add rover: t=" << measurement.gnss->timestamp
      //        << " tag=" << measurement.gnss->tag
      //        << " queue_rov=" << measurement_rov_.size()
      //        << " queue_ref=" << measurement_ref_.size();
      // DLOG(INFO) << stream.str();
    }
    if (measurement.gnss && measurement.gnss_role == GnssRole::Reference) {
      measurement_ref_.push_back(*measurement.gnss);
      // std::ostringstream stream;
      // stream << std::setprecision(17);
      // stream << "DGNSS align add ref: t=" << measurement.gnss->timestamp
      //        << " tag=" << measurement.gnss->tag
      //        << " queue_rov=" << measurement_rov_.size()
      //        << " queue_ref=" << measurement_ref_.size();
      // DLOG(INFO) << stream.str();
    }
  }

  // Get aligned
  inline bool get(const double max_age, 
    GnssMeasurement& rov, GnssMeasurement& ref) {

    const double offset = 0.0;  // shift a time offset for test

    if (measurement_rov_.size() == 0) return false;
    if (measurement_ref_.size() == 0) {
      measurement_rov_.clear(); return false;
    }

    rov = measurement_rov_.front();
    measurement_rov_.pop_front();
    
    // {
    //   std::ostringstream stream;
    //   stream << std::setprecision(17);
    //   stream << "DGNSS align get: rov_t=" << rov.timestamp
    //          << " remaining_rov=" << measurement_rov_.size()
    //          << " ref_queue=[";
    //   for (size_t i = 0; i < measurement_ref_.size(); ++i) {
    //     if (i != 0) stream << ",";
    //     stream << measurement_ref_.at(i).timestamp;
    //   }
    //   stream << "]";
    //   DLOG(INFO) << stream.str();
    // }

    // get the nearest timestamp
    double min_dt = 1.0e6;
    int index = -1;
    for (int i = measurement_ref_.size() - 1; i >= 0; i--) {
      double dt = fabs(rov.timestamp - measurement_ref_.at(i).timestamp - offset);
      if (min_dt > dt) {
        min_dt = dt; index = i;
      }
    }
    CHECK(index != -1);
    ref = measurement_ref_.at(index);
    // {
    //   std::ostringstream stream;
    //   stream << std::setprecision(17);
    //   stream << "DGNSS align choose ref: rov_t=" << rov.timestamp
    //          << " ref_t=" << ref.timestamp
    //          << " ref_index=" << index
    //          << " min_dt=" << min_dt;
    //   DLOG(INFO) << stream.str();
    // }
    // pop all in front of this
    double cut_timestamp = ref.timestamp;
    while (!checkEqual(cut_timestamp, measurement_ref_.front().timestamp, max_age)) {
      measurement_ref_.pop_front();
    }
    // {
    //   std::ostringstream stream;
    //   stream << std::setprecision(17);
    //   stream << "DGNSS align after pop: chosen_ref_t=" << ref.timestamp
    //          << " remaining_ref_queue=[";
    //   for (size_t i = 0; i < measurement_ref_.size(); ++i) {
    //     if (i != 0) stream << ",";
    //     stream << measurement_ref_.at(i).timestamp;
    //   }
    //   stream << "]";
    //   DLOG(INFO) << stream.str();
    // }

    // Check timestamps
    if (!checkEqual(rov.timestamp, ref.timestamp, max_age)) {
      LOG(WARNING) << "Max age between two measurements exceeded! "
        << "age = " << fabs(rov.timestamp - ref.timestamp)
        << ", max_age = " << max_age;
      return false;
    }

    return true;
  }

protected:
  // Measurement storage for aligning
  std::deque<GnssMeasurement> measurement_rov_;
  std::deque<GnssMeasurement> measurement_ref_;
};

}