// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/AddJobPairData.cc
/// @brief  A really simple mover that takes some data in through xml and appends it to the pose
/// @author Sam DeLuca <Samuel.l.deluca@vanderbilt.edu)

#include <protocols/simple_moves/AddJobPairData.hh>
#include <protocols/simple_moves/AddJobPairDataCreator.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace simple_moves {

std::string AddJobPairDataCreator::keyname() const
{
	return AddJobPairDataCreator::mover_name();
}

moves::MoverOP AddJobPairDataCreator::create_mover() const
{
	return new AddJobPairData;
}

std::string AddJobPairDataCreator::mover_name()
{
	return "AddJobPairData";
}

AddJobPairData::AddJobPairData() :
	string_key_(""),
	string_value_(""),
	real_value_(0.0),
	value_type_(string_value)
{

}

AddJobPairData::AddJobPairData(AddJobPairData const & src) : Mover(src),
	string_key_(src.string_key_),
	string_value_(src.string_value_),
	real_value_(src.real_value_),
	value_type_(src.value_type_)
{

}

AddJobPairData::~AddJobPairData()
{

}

void AddJobPairData::apply( Pose & )
{
	jd2::JobOP current_job = jd2::JobDistributor::get_instance()->current_job();
	if(value_type_ == string_value)
	{
		current_job->add_string_string_pair(string_key_,string_value_);
	}else if(value_type_ == real_value)
	{
		current_job->add_string_real_pair(string_key_,real_value_);
	}else
	{
		//This really shouldn't happen
		assert(false);
	}
}

std::string AddJobPairData::get_name() const
{
	return "AddJobPairData";
}

moves::MoverOP AddJobPairData::clone() const
{
	return new AddJobPairData(*this);
}

moves::MoverOP AddJobPairData::fresh_instance() const
{
	return new AddJobPairData;
}

void AddJobPairData::parse_my_tag(
		TagPtr const tag,
		moves::DataMap &,
		Filters_map const &,
		moves::Movers_map const &,
		Pose const & )
{
	if (!tag->hasOption("value_type"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' mover requires option 'value_type'");
	}
	if (!tag->hasOption("key"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' mover requires option 'key'");
	}
	if (!tag->hasOption("value"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' mover requires option 'value'");
	}

	std::string value_type_string(tag->getOption<std::string>("value_type"));

	if(value_type_string == "string")
	{
		value_type_ = string_value;
	}
	else if(value_type_string == "real")
	{
		value_type_ = real_value;
	}
	else
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'AddJobPairData' option 'value_type' can only be 'string' or 'real'");
	}

	string_key_ = tag->getOption<std::string>("key");
	if(value_type_ == string_value)
	{
		string_value_ = tag->getOption<std::string>("value");
	}else if(value_type_ == real_value)
	{
		real_value_ = tag->getOption<core::Real>("value");
	}else
	{
		// If this happens all the error checking code above is wrong
		assert(false);
	}

}

}
}
