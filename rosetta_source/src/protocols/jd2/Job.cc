// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/Job.cc
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Job classes
/// @author Steven Lewis smlewi@gmail.com

///Unit headers
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///Project headers

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.hh>

#include <core/pose/Pose.hh>

///Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

///C++ headers
#include <string>

#include <utility/vector1.hh>


static basic::Tracer TR("protocols.jd2.Job");

namespace protocols {
namespace jd2 {


////////////////////////////Job/////////////////////////////
Job::Job( InnerJobOP inner_job, core::Size nstruct_index )
	: inner_job_(inner_job),
		nstruct_index_(nstruct_index),
		status_prefix_( "" ),
		completed_(false)
{
	//TR.Trace << "Using Job (base class) for JobDistributor" << std::endl;
}

Job::~Job(){}

///@brief returns a copy of this object whose "output fields" are zeroed out.  Used by the JobDistributor in cases where the job fails and must be retried to prevent accumulation of Job state after a failure.  This implementation was chosen over a clear_all_output function to prevent mover A from deleting mover B's hard work!  You probably should not be trying to call this function.
JobOP Job::copy_without_output() const{
	return new Job(inner_job_, nstruct_index_);
}

InnerJobCOP Job::inner_job() const { return inner_job_; }

std::string const & Job::input_tag() const {
	return inner_job_->input_tag();
}

InnerJobOP Job::inner_job_nonconst() { return inner_job_; }

core::Size Job::nstruct_index() const { return nstruct_index_; }

core::Size Job::nstruct_max() const {
	return inner_job_->nstruct_max();
}

//functions for loading output info into the job
///@brief add an output string
void Job::add_string( std::string const & string_in ){
	long_strings_.push_back(string_in);
}

///@brief adds output strings
void Job::add_strings( Strings const & strings )
{
	long_strings_.insert( long_strings_.end(), strings.begin(), strings.end() );
}

///@brief add a string/string pair
void Job::add_string_string_pair( std::string const & string1, std::string const & string2 ){
	string_string_pairs_.push_back(std::make_pair(string1, string2));
}

///@brief add a string/real pair
void Job::add_string_real_pair( std::string const & string_in, core::Real const real_in ){
	string_real_pairs_[ string_in ] = real_in;
}

///@brief return a COP to the input pose
core::pose::PoseCOP Job::get_pose() const {
	// if pose is loaded into job-object this pose has precedence.
	if ( inner_job_->get_pose() ) {
		return inner_job_->get_pose();
	} else { //if not ask job-inputter
		protocols::jd2::JobDistributor* jd
			= protocols::jd2::JobDistributor::get_instance();
		core::pose::PoseOP aPose = new core::pose::Pose;
		// in the following call we use a copy of the job-output. the idea is that if the inputter has not stored the
		// pose in the Job-Object already he will not do it now... in principle one could also use a const_cast which has
		// basically the same effect -- the copy copies a pointer and an integer...
		if ( jd->job_inputter() ) jd->job_inputter()->pose_from_job( *aPose, copy_without_output() /*or const_cast(this) */ );
		return aPose;
	}
	utility_exit_with_message( "Programming error: you asked for Job::get_pose() but there is neither a job_inputter nor a pose loaded into the Job-Object ");
	return NULL;
}


///@brief in-place copy of input pose
void Job::get_pose( core::pose::Pose& pose ) const {
	// if pose is loaded into job-object this pose has precedence.
	if ( inner_job_->get_pose() ) {
		pose = *( inner_job_->get_pose() );
		return;
	} else { //if not ask Job-inputter
		protocols::jd2::JobDistributor* jd
			= protocols::jd2::JobDistributor::get_instance();
		if ( jd->job_inputter() ) {
			// in the following call we use a copy of the job-output. the idea is that if the inputter has not stored the
			// pose in the Job-Object already he will not do it now... in principle one could also use a const_cast which has
			// basically the same effect -- the copy copies a pointer and an integer...
			jd->job_inputter()->pose_from_job( pose, copy_without_output() );
		}
		return;
	}
	utility_exit_with_message( "Programming error: you asked for Job::get_pose() but there is neither a job_inputter nor a pose loaded into the job-object");
}


//functions for returning output info from the job.  You get iterators so that this interface can stay constant as the underlying implementation changes
Job::Strings::const_iterator Job::output_strings_begin() const
{ return long_strings_.begin(); }

Job::Strings::const_iterator Job::output_strings_end() const
{ return long_strings_.end(); }

Job::StringStringPairs::const_iterator Job::output_string_string_pairs_begin() const
{ return string_string_pairs_.begin(); }

Job::StringStringPairs::const_iterator Job::output_string_string_pairs_end() const
{ return string_string_pairs_.end(); }

Job::StringRealPairs::const_iterator Job::output_string_real_pairs_begin() const
{ return string_real_pairs_.begin(); }

Job::StringRealPairs::const_iterator Job::output_string_real_pairs_end() const
{ return string_real_pairs_.end(); }


bool Job::bad() const {
	return inner_job_->bad();
}

void Job::set_bad(bool value)  {
	inner_job_->set_bad( value );
}


JobOP const JD2_BOGUS_JOB( new Job( (new InnerJob("EMPTY_JOB_use_jd2", 0)), 0) );

} // jd2
} // protocols
