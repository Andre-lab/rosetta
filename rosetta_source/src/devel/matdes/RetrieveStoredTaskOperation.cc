// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/RetrieveStoredTaskOperation.cc
/// @brief	Restricts specified residue types to only repack, no design.
/// @author Neil King (neilking@uw.edu)

// Unit Headers
#include <devel/matdes/RetrieveStoredTaskOperation.hh>
#include <devel/matdes/RetrieveStoredTaskOperationCreator.hh>

#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/operation/ReplicateTask.hh>
#include <core/pose/symmetry/util.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <devel/matdes/STMStoredTask.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>
#include <boost/functional/hash.hpp>


using basic::Error;
using basic::Warning;
static basic::Tracer TR( "devel.matdes.RetrieveStoredTaskOperation" );

namespace devel {
namespace matdes {

// @brief default constructor	
RetrieveStoredTaskOperation::RetrieveStoredTaskOperation() {}

// @brief destructor
RetrieveStoredTaskOperation::~RetrieveStoredTaskOperation() {}

core::pack::task::operation::TaskOperationOP
RetrieveStoredTaskOperationCreator::create_task_operation() const
{
	return new RetrieveStoredTaskOperation;
}

// @brief copy constructor
core::pack::task::operation::TaskOperationOP RetrieveStoredTaskOperation::clone() const
{
	return new RetrieveStoredTaskOperation( *this );
}

// @brief apply function
void
RetrieveStoredTaskOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		utility_exit_with_message("Your pose does not have CacheableData of type STM_STORED_TASKS");
	} else {
		devel::matdes::STMStoredTask const & stored_tasks = *( static_cast< devel::matdes::STMStoredTask const* >( 
			             pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS )() ) );
		if (!stored_tasks.has_task(task_name_)) {
			utility_exit_with_message("No stored task with the name " + task_name_ + " found");
		} else {
  			task.update_commutative( *( stored_tasks.get_task( task_name_ ) ) );
		}
	}
}

// @brief parse xml
void
RetrieveStoredTaskOperation::parse_tag( TagPtr tag )
{
	task_name_ = tag->getOption< std::string >( "task_name" ) ;
}
	
} //namespace matdes
} //namespace devel
