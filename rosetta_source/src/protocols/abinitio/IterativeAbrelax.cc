// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file IterativeAbrelax
/// @brief iterative protocol starting with abinitio and getting progressively more concerned with full-atom relaxed structures
/// @detailed
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/abinitio/IterativeAbrelax.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <basic/Tracer.hh>

#include <utility/file/FileName.hh>

//// C++ headers
#include <cstdlib>
#include <string>

// Utility headers
#include <basic/options/option_macros.hh>

#include <utility/vector1.hh>


static basic::Tracer tr("protocols.iterative");

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_1GRP_KEY( Boolean, iterative, fullatom )

namespace protocols {
namespace abinitio {
using namespace jd2::archive;

bool IterativeAbrelax::options_registered_( false );


void protocols::abinitio::IterativeAbrelax::register_options() {
	if ( !options_registered_ ) {
		NEW_OPT( iterative::fullatom, "sample fullatom structures in iterative protocol", false );
		options_registered_ = true;
	}
	IterativeBase::register_options();
	IterativeFullatom::register_options();
	ArchiveManager::register_options();
}

IterativeAbrelax::IterativeAbrelax()
	: Parent( "abstract" ),
		centroid_archive_(&fullatom_archive_ ),
		fullatom_( false )
{
	fullatom_ = option[ iterative::fullatom ];
	//not needed anymore: if I initialize base with this finish_stage
	if ( !fullatom_ ) centroid_archive_.set_finish_stage( IterativeBase::LAST_CENTROID_START );
}

void IterativeAbrelax::read_structures( core::io::silent::SilentFileData& sfd, Batch const& batch ) {
	if ( sfd.begin() != sfd.end() ) {
		pose::Pose pose;
		sfd.begin()->fill_pose( pose );
		sfd.begin()->has_energy( "user_tag" );
		std::string centroid_str( sfd.begin()->get_string_value( "user_tag" ) );
		bool isPseudoFullatom ( centroid_str == "centroid" );
		if ( pose.is_fullatom() && !isPseudoFullatom ) {
			runtime_assert( fullatom_ );
			tr.Debug << "reading full-atom structures into fullatom_pool" << std::endl;
			fullatom_archive_.read_structures( sfd, batch );
		} else {
			tr.Debug << "reading " << (isPseudoFullatom ? " pseuod-fullatom " : "centroid" ) <<  " structures into centroid_pool" << std::endl;
			centroid_archive_.read_structures( sfd, batch );
		}
	}
}

bool IterativeAbrelax::finished() const {
	return ( centroid_archive_.finished() && !fullatom_ ) || fullatom_archive_.finished();
}

bool IterativeAbrelax::still_interested( Batch const& batch ) const {
	//	return Parent::still_interested( batch )&&
	return ( centroid_archive_.still_interested( batch ) && ( fullatom_ && fullatom_archive_.still_interested( batch ) ) );
}

///@details ready for new batch .... if queue is empty batch will be generated any way, but otherwise we only generate if this yields true.
///  logic here: new batch at beginning, but only if we are in startup phase ( not a reload of a full archive )
///              otherwise make new batch if sufficiently many structures have been accepted since last batch
// bool IterativeAbrelax::ready_for_batch() const {
// 	return centroid_archive_.ready_for_batch() || ( fullatom_ && fullatom_archive_.ready_for_batch() );
// }

///@details generate new batch...
/// type of batch depends on stage_. we switch to next stage based on some convergence criteria:
/// right now it is how many decoys were accepted from last batch.. if this number drops sufficiently ---> next stage...
///    (maybe need to put a safeguard in here: ratio small but at least XXX decoys proposed since last batch... )
///
void IterativeAbrelax::generate_batch() {
	if ( fullatom_ && fullatom_archive_.ready_for_batch() ) {
		// fixed this already by using the noesy_cst_file from the very first batch that is read in...
		// re-starting assign cycle from cmdline parameter..
// 		if ( fullatom_archive_.first_noesy_fa_cst_file()=="n/a") {
// 			fullatom_archive_.set_first_noesy_fa_cst_file( centroid_archive_.first_noesy_fa_cst_file() );
// 			fullatom_archive_.set_noesy_assign_float_cycle( centroid_archive_.noesy_assign_float_cycle() );
// 		}
		fullatom_archive_.generate_batch();
	} else {
		centroid_archive_.generate_batch();
	}
}

void IterativeAbrelax::set_manager( jd2::archive::ArchiveManagerAP manager ) {
	Parent::set_manager( manager );
	tr.Info << "IterativeAbrelax: set ArchiveManager also for sub-archives... " << std::endl;
	fullatom_archive_.set_manager( manager );
	centroid_archive_.set_manager( manager );
}


void IterativeAbrelax::idle() {
	centroid_archive_.idle();
	if ( fullatom_ ) fullatom_archive_.idle();
}

void IterativeAbrelax::save_to_file( std::string suffix ) {
	centroid_archive_.save_to_file( suffix );
	if ( fullatom_ ) fullatom_archive_.save_to_file( suffix );
}

void IterativeAbrelax::save_status( std::ostream& os ) const {
	centroid_archive_.save_status( os );
	if ( fullatom_ ) fullatom_archive_.save_status( os );
}

bool IterativeAbrelax::restore_from_file() {
	bool	b_have_restored =	centroid_archive_.restore_from_file();
	if ( fullatom_ ) fullatom_archive_.restore_from_file();
	return b_have_restored;
}

void IterativeAbrelax::init_from_decoy_set( core::io::silent::SilentFileData const& sfd ) {
	centroid_archive_.init_from_decoy_set( sfd );
}



} //abinitio
} //protocols
