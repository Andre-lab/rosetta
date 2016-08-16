// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#include <--path--/--class--.hh>
#include <--path--/--class--Creator.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>


static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--

--class--::--class--():
	TaskOperation()

{

}

--class--::~--class--() {}

TaskOperationOP
--class--::clone() const {
	return TaskOperationOP( new --class--( *this ) );
}

--class--::--class--( --class-- const & src ):
	TaskOperation(src)

{

}

void
--class--::parse_tag( utility::tag::TagCOP , basic::datacache::DataMap& ){


}



void
--class--::apply( const core::pose::Pose& pose, core::pack::task::PackerTask& ) const {


}

std::string
--class--::keyname() {
	return "--class--";
}

void
--class--::provide_xml_schema( XMLSchemaDefinition & xsd )
{
}

core::pack::task::operation::TaskOperationOP
--class--Creator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new --class-- );
}

std::string
--class--Creator::keyname() const
{
	return --class--::keyname();
}

void
--class--Creator::provide_xml_schema( XMLSchemaDefinition & xsd )
{
	--class--::provide_xml_schema( xsd );
}





--end_namespace--









