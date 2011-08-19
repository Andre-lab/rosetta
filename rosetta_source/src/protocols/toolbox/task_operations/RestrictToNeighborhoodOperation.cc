// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.cc
/// @brief  TaskOperation class that finds a neighborhood and leaves it mobile in the PackerTask
/// @author Steven Lewis smlewi@unc.edu

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperationCreator.hh>

// Project Headers
#include <core/pose/Pose.hh>

#include <core/pack/task/PackerTask.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
// AUTO-REMOVED #include <basic/MetricValue.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1_bool.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <set>

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.toolbox.TaskOperations.RestrictToNeighborhoodOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

RestrictToNeighborhoodOperation::RestrictToNeighborhoodOperation() {}

///@details this ctor assumes a pregenerated Neighborhood and Neighbors calculators - if you want a particular non-default cutoff distance, assemble those calculators separately then pass them to this operation.
RestrictToNeighborhoodOperation::RestrictToNeighborhoodOperation( std::string const & calculator )
	: parent(), calculator_name_(calculator)
{
	//I suppose you could reasonably create this object BEFORE the calculator was generated/registered
// 	if( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( calculator_name_ ) ){
// 		utility_exit_with_message("In RestrictToNeighborhoodOperation, calculator " + calculator + " does not exist.");
// 	}
}

///@details this ctor generates calculators (easier to use but will rely on defaults, including default distance cutoff)
RestrictToNeighborhoodOperation::RestrictToNeighborhoodOperation( std::set< core::Size > const & central_residues )
	: parent(), calculator_name_("")
{
	make_calculator( central_residues );
}

///@details private helper function to make calculator - runs in the ctor
void RestrictToNeighborhoodOperation::make_calculator( std::set< core::Size > const & central_residues ) {
	make_name( central_residues );

	using namespace core::pose::metrics;
	if( CalculatorFactory::Instance().check_calculator_exists( calculator_name_ ) ){
		Warning() << "In RestrictToNeighborhoodOperation, calculator " << calculator_name_
							<< " already exists, this is hopefully correct for your purposes" << std::endl;
	} else {
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;
	CalculatorFactory::Instance().register_calculator( calculator_name_, new NeighborhoodByDistanceCalculator( central_residues ) );
	}
}

///@details private helper function to name calculator- runs in the ctor
void RestrictToNeighborhoodOperation::make_name( std::set< core::Size > const & central_residues ) {
	calculator_name_ = "RTNhO_calculator";

	for(SizeSet::const_iterator it(central_residues.begin()), end(central_residues.end()) ; it != end; ++it){
		calculator_name_ += '_' + utility::to_string( *it );
	}

}

RestrictToNeighborhoodOperation::~RestrictToNeighborhoodOperation() {}

core::pack::task::operation::TaskOperationOP
RestrictToNeighborhoodOperationCreator::create_task_operation() const
{
	return new RestrictToNeighborhoodOperation;
}

///@details be warned if you use clone that you'll not get a new interface calculator
core::pack::task::operation::TaskOperationOP RestrictToNeighborhoodOperation::clone() const
{
	return new RestrictToNeighborhoodOperation( *this );
}

void
RestrictToNeighborhoodOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{

 	//vector for filling packertask
 	utility::vector1_bool repack(pose.total_residue(), false);

	run_calculator(pose, calculator_name_, "neighbors", repack);

 	task.restrict_to_residues(repack);

}

} //namespace protocols
} //namespace toolbox
} //namespace task_operations
