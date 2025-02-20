// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file InterfaceTaskOperation.cc
/// @brief Base class for TaskOperations that will work with docking to define an interface
/// @author Brian Weitzner (brian.weitzner@jhu.edu)

#include <protocols/simple_task_operations/InterfaceTaskOperation.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


namespace protocols {
namespace simple_task_operations {

//using namespace core;
using namespace core::pack::task::operation;
using namespace utility::tag;

InterfaceTaskOperation::~InterfaceTaskOperation()= default;

/* AMW No creator exists, perhaps because base class.
void InterfaceTaskOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
InterfaceTaskOperation::provide_xml_schema( xsd );
}

std::string InterfaceTaskOperationCreator::keyname() const
{
return InterfaceTaskOperation::keyname();
}
*/

InterfaceTaskOperation::InterfaceTaskOperation() : parent()
{
	movable_jumps_ = utility::vector1_int();
	movable_jumps_.push_back( 1 );
	movable_chains_ = utility::vector1_int();
}
InterfaceTaskOperation::InterfaceTaskOperation( InterfaceTaskOperation const & old_instance ) :
	//utility::VirtualBase(),
	parent( old_instance )
{
	movable_jumps_ = old_instance.movable_jumps_;
	movable_chains_ = old_instance.movable_chains_;
}

void InterfaceTaskOperation::add_movable_jump( int const additional_jump ) {
	movable_jumps_.push_back( additional_jump );
}

void InterfaceTaskOperation::set_movable_jumps( utility::vector1_int const & movable_jumps ) {
	movable_jumps_ = movable_jumps;
}

void InterfaceTaskOperation::add_movable_chain( int const additional_chain ) {
	movable_chains_.push_back( additional_chain );
}

void InterfaceTaskOperation::set_movable_chains( utility::vector1_int const & movable_chains ) {
	movable_chains_ = movable_chains;
}

utility::vector1_int const &
InterfaceTaskOperation::movable_chains() const{
	return movable_chains_;
}

utility::vector1_int const &
InterfaceTaskOperation::movable_jumps() const{
	return movable_jumps_;
}

void InterfaceTaskOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

}
}
