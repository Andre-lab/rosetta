// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/StringMetric.hh
///
/// @brief  Base class for core::Real StringMetrics
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_simple_metrics_StringMetric_hh
#define INCLUDED_core_simple_metrics_StringMetric_hh

// Unit headers
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/StringMetric.fwd.hh>

// Core headers

#include <core/pose/Pose.fwd.hh>


// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace simple_metrics {

/// @brief A class that is derived to calculate a string value.
///
///  Apply(pose) method calculates this metric and adds it to the pose score for output.
///  Calculate(pose) method returns a string
///
class StringMetric : public core::simple_metrics::SimpleMetric {

public: // constructors / destructors
	StringMetric();
	~StringMetric() override;

	StringMetric( StringMetric const & other );

	using SimpleMetric::apply; // Add back the base class implementation

	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as out_label
	///
	/// @details Score is added to the SimpleMetricData cache in the pose.
	///           Data is output into the final score file, but can be accessed if needed through the cache.
	void
	apply(
		std::string const & out_label,
		pose::Pose & pose,
		bool override_existing_data=false) const override;

	///@brief Calculate the metric.
	///
	///@details
	/// Calculate and return the value.  This function can not store the value.
	/// Use the apply function to store the value in the pose.
	virtual std::string
	calculate( pose::Pose const & pose ) const = 0;

	///@brief Grab the data from the pose if it exists or calculate the metric
	///
	///@details If use_cache is true, we will attempt to pull the data from the pose.
	/// If fail_on_missing_cache is true, we will fail, otherwise, we will calculate the metric.
	///
	/// This function is meant to support caching metrics, so values do not need to be calculated twice,
	///  for example in SimpleMetricFilter/Features
	///  or code-wise where data takes a while to calculate and can be reused.
	///
	std::string
	cached_calculate(
		pose::Pose const & pose,
		bool use_cache,
		std::string prefix="",
		std::string suffix="",
		bool fail_on_missing_cache=true) const;


public:

	///@brief Name of the class
	std::string
	name() const override = 0;

	///@brief Name of the metric
	std::string
	metric() const override = 0;

	///@brief Get the metric name that this Metric will calculate
	utility::vector1< std::string >
	get_metric_names() const override;

public:
	/// @brief called by parse_my_tag -- should not be used directly
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override = 0;

	SimpleMetricOP
	clone() const override = 0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION


}; //class StringMetrics

} //namespace simple_metrics
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_simple_metrics_StringMetric )
#endif // SERIALIZATION


#endif
