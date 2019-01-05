// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/TaskFactory.hh
/// @brief  Task class to describe packer's behavior header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#include <utility/exit.hh>

// Unit Headers
#include <core/pack/task/TaskFactory.hh>

// Package Headers
#include <core/pack/task/PackerTask_.hh> // only place in all of mini where this gets #included (except PackerTask_.cc)
#include <core/pack/task/operation/TaskOperation.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#endif // SERIALIZATION



namespace core {
namespace pack {
namespace task {

TaskFactory::TaskFactory() : parent() {}
TaskFactory::TaskFactory( TaskFactory const & src)
:
	parent()
{
	copy_operations( src );
}

TaskFactoryOP TaskFactory::clone() const
{
	return utility::pointer::make_shared< TaskFactory >(*this);
}

TaskFactory::~TaskFactory() = default;

TaskFactory &
TaskFactory::operator = ( TaskFactory const & rhs )
{
	copy_operations( rhs );
	return *this;
}

void
TaskFactory::modify_task( core::pose::Pose const & pose, PackerTaskOP task ) const
{
	runtime_assert( task != nullptr );
	for ( TaskOperationCOP taskop : *this ) {
		taskop->apply( pose, *task );
	}
}

// Non static version.
PackerTaskOP
TaskFactory::create_task_and_apply_taskoperations( pose::Pose const & pose ) const
{
	PackerTaskOP task( new PackerTask_( pose ) );
	modify_task( pose, task );
	return task;
}

// clones the input task, and pushes it back into the list
void
TaskFactory::push_back( TaskOperationCOP taskop )
{
	operations_.push_back( taskop->clone() );
}

TaskFactory::const_iterator
TaskFactory::begin() const
{
	return operations_.begin();
}

TaskFactory::const_iterator
TaskFactory::end() const
{
	return operations_.end();
}

void
TaskFactory::clear()
{
	operations_.clear();
}


PackerTaskOP
TaskFactory::create_packer_task(
	pose::Pose const & pose
)
{
	return utility::pointer::make_shared< PackerTask_ >( pose );
}

void
TaskFactory::copy_operations( TaskFactory const & src )
{
	for ( auto const & taskop_iter : src ) {
		operations_.push_back( taskop_iter->clone() );
	}
}

core::Size
TaskFactory::size() const
{
	core::Size const size( operations_.size() );
	return( size );
}

} //namespace task
} //namespace pack
} //namespace core

#ifdef SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::TaskFactory::save( Archive & arc ) const {
	arc( CEREAL_NVP( operations_ ) ); //OperationList
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::TaskFactory::load( Archive & arc ) {
	arc( operations_ ); //OperationList
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::TaskFactory );
CEREAL_REGISTER_TYPE( core::pack::task::TaskFactory )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_TaskFactory )
#endif // SERIALIZATION
