// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/SilentFileJobOutputter.cc
/// @brief
/// @author Oliver Lange

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/pose/Pose.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

///Utility headers
#include <utility/file/FileName.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

#include <string>
#include <algorithm>

static basic::Tracer tr("protocols.jd2.SilentFileJobOutputter");

namespace protocols {
namespace jd2 {

SilentFileJobOutputter::SilentFileJobOutputter() {
	set_defaults();
	read_done_jobs();
}

SilentFileJobOutputter::~SilentFileJobOutputter() {
	//DO NOT PUT THINGS HERE - it is not guarunteed to get called - use flush below instead.
}

void SilentFileJobOutputter::flush() {
	write_all_structs();
}

void SilentFileJobOutputter::write_all_structs() {
	using std::pair;
	using utility::vector1;
	using utility::file::FileName;
	using core::io::silent::SilentStructOP;

	core::io::silent::SilentFileData sfd;
	typedef vector1< std::pair< SilentStructOP, FileName > >::iterator iter;

	// Only write structures if the user hasn't disabled it - otherwise it totally breaks
	// the user's expectation.
	if( !bWriteNoStructures_ ){
		tr.Debug << "writing " << saved_structs_.size() << " structs." << std::endl;
		for ( iter it = saved_structs_.begin(), end = saved_structs_.end();
					it != end; ++it
		) {
			//tr.Debug << "writing struct " << ss->decoy_tag() << std::endl;
			//tr.Debug << "writing struct " << (*it->first)->decoy_tag() << std::endl;
			//SilentStructOP ss = it->first;
			sfd.write_silent_struct( (*it->first), it->second );
		}
	}
	// very important to clear after writing!
	saved_structs_.clear();

	tr.Debug << "currently have " << saved_structs_.size() << " structs."
		<< std::endl;
	tr.flush();
}

void SilentFileJobOutputter::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	silent_file_ = option[ out::file::silent ]();
	if ( silent_file_.relative() ) {
		silent_file_.path( option[ out::path::all ]().path() + "/" + silent_file_.path() );
		//FileName takes care of platform-specific path seperator, i.e.,  "/" or "\" ...
	}
	tr.Debug << "SilentFileJobOutputter setup for file " << silent_file_ << std::endl;

#ifdef USEMPI
	int mpi_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);/* get current process id */
	std::string name = silent_file_;
	// attach mpi rank to out files
	size_t lastslash = name.find_last_of("/\\");
	size_t lastdot   = name.find_last_of('.');
	silent_file_ = name;
#endif

	// dd_parser should use binary silent files as a default. but score-only
	// silent files are also an option
	// this is really stupid. What knowledge do you about the user's intentions
	// that trumps the user's input via the command-line? A better idea would be
	// to exclude certain SilentStruct types, throw a warning or error message,
	// or do something more intelligent like that.
	if ( option[ basic::options::OptionKeys::jd2::dd_parser ]() ) {
		if ( option[ out::file::silent_struct_type ]() != "score")
			option[ out::file::silent_struct_type ].value( "binary" );
	}

	//default is 1
	n_to_buffer_ = option[ basic::options::OptionKeys::jd2::buffer_silent_output ]();

	bWriteIntermediateFiles_ = (
		option[ run::intermediate_scorefiles ]() ||
		option[ run::intermediate_structures ]()
	);

	bWriteIntermediateStructures_ = option[ run::intermediate_structures ]();

	bWriteNoStructures_ = false;

	if ( option[ OptionKeys::out::file::scorefile ].user() ){
		write_separate_scorefile_ = true;
	}	else write_separate_scorefile_ = false;
}

void SilentFileJobOutputter::read_done_jobs() {
	core::io::silent::SilentFileData sfd;
	if ( utility::file::file_exists( silent_file_ ) ) {
		silent_file_tags_ = sfd.read_tags_fast( silent_file_ );
		foreach( std::string & tag, silent_file_tags_ ) {
		/// eliminate the FAILURE_ prefix so that jobs know to start from
		/// the 'next' nstruct on restart. This is important to avoid duplicate
		/// entries
			if( tag.substr( 0, 8 ) == "FAILURE_" )
				tag = tag.substr( 8, 1000 );
		} //foreach
	} //fi
}

void SilentFileJobOutputter::final_pose(
	JobCOP job, core::pose::Pose const & pose
) {
	core::io::silent::SilentStructOP ss =
		dump_pose( silent_file_, job, pose,  bWriteNoStructures_ /*this is always false */ /* bWriteScoresOnly */);


	// only write a scorefile if specified by the user
	using namespace basic::options;

	if ( write_separate_scorefile_ ) {
		//write here to avoid 2x evaluated for same structure
		core::io::silent::SilentFileData sfd;
		sfd.write_silent_struct( *ss, scorefile_name(), true );
	}

	//calling dump_pose twice means calling evaluator twice which can be expensive!!

	//	if ( !bWriteNoStructures_ ) dump_pose( scorefile_name(), job, pose, true );
	//  if (bWriteNoStructures_)  we have already  score-only output ... don't have it redundant
}

/// @brief this function is intended for saving mid-protocol poses; for example
/// the final centroid structure in a combined centroid/fullatom protocol.
/// --->these go to file silent_filename+tag
void SilentFileJobOutputter::other_pose(
	JobCOP job,
	core::pose::Pose const & pose,
	std::string const & tag,
	int copy_count, /*default -1 */
	bool score_only /*default false*/
) {
	utility::file::FileName filename( silent_file_ );
	filename.base( silent_file_.base() +"_"+ tag );

	core::io::silent::SilentStructOP ss;
	if ( bWriteIntermediateFiles_ ) {
		ss=dump_pose( filename, job, pose, !bWriteIntermediateStructures_ || score_only , copy_count );
	}

	if ( write_separate_scorefile_ && ss ) {
		//write here to avoid 2x evaluated for same structure
		core::io::silent::SilentFileData sfd;
		utility::file::FileName filename( scorefile_name() );
		filename.base( scorefile_name().base() +"_"+ tag );
		sfd.write_silent_struct( *ss, filename, true );
	}


}

core::io::silent::SilentStructOP SilentFileJobOutputter::dump_pose(
	utility::file::FileName const & filename,
	JobCOP job,
	core::pose::Pose const & pose_in,
	bool bWriteScoreOnly,
	int copy_count
) {
	PROF_START( basic::JD2_SILENT_OUTPUTTER );
	core::io::silent::SilentFileData sfd;

	using core::io::silent::SilentStructFactory;
	core::io::silent::SilentStructOP ss;
	if ( bWriteScoreOnly ) {
		ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("score");
	} else {
		ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( pose_in );
	}

	std::ostringstream tag;
	tag << output_name( job );
	if ( copy_count>=0 ) tag << '_' << std::setfill('0') << std::setw(8) << copy_count;
	ss->fill_struct( pose_in, tag.str() );
	add_job_data_to_ss( ss, job );

	core::pose::Pose pose( pose_in );
	evaluate( pose, tag.str(), *ss );


	add_silent_struct( ss, filename );
	tr.Debug << "adding struct " << ss->decoy_tag() << std::endl;
	tr.Debug << "have " << saved_structs_.size() << ", buffering " << n_to_buffer_ << std::endl;
	if ( saved_structs_.size() >= n_to_buffer_ ) {
		write_all_structs();
	}

	tr.flush();
	PROF_STOP( basic::JD2_SILENT_OUTPUTTER );
	return ss;
}

void SilentFileJobOutputter::add_silent_struct(
	core::io::silent::SilentStructOP ss,
	utility::file::FileName const & fn
) {
	saved_structs_.push_back( std::make_pair( ss, fn ) );
}

/////////////////////////////////state of output functions/////////////////////////////////
bool SilentFileJobOutputter::job_has_completed( JobCOP job ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// did we complete the job later ?
	if ( job->completed() ) {
		return true;
	}

	// was the job completed beforehand ( in the silent file before the app even
	// started ) ?
	if ( option[ run::multiple_processes_writing_to_one_directory ].value() ) {
		read_done_jobs(); // refresh silent_file_tags_ for parallel processes
	}
	CompareTags predicate( output_name(job) );

	bool const already_written(
		find_if(silent_file_tags_.begin(), silent_file_tags_.end(), predicate) != silent_file_tags_.end()
	);

	using std::pair;
	using std::string;
	using utility::vector1;
	using utility::file::FileName;
	using core::io::silent::SilentStructOP;

	vector1< string > tags;
	typedef vector1< pair< SilentStructOP, FileName > >::const_iterator iter;
	for ( iter it = saved_structs_.begin(), end = saved_structs_.end(); it != end; ++it ) {
		tags.push_back( it->first->decoy_tag() );
	}

	bool const already_buffered(
		find_if( tags.begin(), tags.end(), predicate ) != tags.end()
	);

	return ( already_written || already_buffered );
}

/// @details
/// SilentFile tags should preserve the FULL NAME such that we don't end up with
/// duplicate tags. This will cause problems on BOINC if changed.
std::string SilentFileJobOutputter::output_name( JobCOP job ) {
	return affixed_numbered_name( job );
}


void
SilentFileJobOutputter::set_silent_file_name( utility::file::FileName name ){
	silent_file_ = name;
	read_done_jobs(); //safer to do this again
}

void
SilentFileJobOutputter::set_write_separate_scorefile( bool write_separate_scorefile ){
	write_separate_scorefile_ = write_separate_scorefile;
}

} //jd2
} //protocols
