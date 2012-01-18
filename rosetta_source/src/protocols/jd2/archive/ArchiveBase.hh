// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/jd2/MPIWorkPoolJobDistributor.hh
/// @brief  header for MPIWorkPoolJobDistributor - intended for continuous resamplig jobs  that spawn new jobs based on a pool/archive of
///         structures
/// @author Oliver Lange olange@u.washington.edu

#ifndef INCLUDED_protocols_jd2_archive_ArchiveBase_hh
#define INCLUDED_protocols_jd2_archive_ArchiveBase_hh

// Unit headers
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>
#include <protocols/jd2/archive/ArchiveBase.fwd.hh>


// Package headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>


// Utility headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <list>
#include <deque>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {
namespace archive {
//class ArchiveManager;



///@brief Tags used to tag messeges sent by MPI functions used to decide whether a slave is requesting a new job id or
///flagging as job as being a bad input

///@details This job distributor is meant for running jobs where the machine you are using has a large number of
///processors, the number of jobs is much greater than the number of processors, or the runtimes of the individual jobs
///could vary greatly. It dedicates the head node (whichever processor gets processor rank #0) to handling job requests
///from the slave nodes (all nonzero ranks). Unlike the MPIWorkPartitionJobDistributor, this JD will not work at all
///without MPI and the implementations of all but the interface functions have been put inside of ifdef directives.
///Generally each function has a master and slave version, and the interface functions call one or the other depending
///on processor rank.
class AbstractArchiveBase : public utility::pointer::ReferenceCount {
public:
	AbstractArchiveBase( ArchiveManagerAP ptr ) : manager_( ptr ), name_( "archive" ) {};
	AbstractArchiveBase() : manager_( NULL ), name_( "archive" ) {};

	///@brief is archive converged ?
	virtual bool finished() const = 0;

	//this is probably obsoleted
	//	virtual bool ready_for_batch() const = 0;

	///@brief old-batches might be outdated and should not be computed anymore
	/// return true for this query if this is the case for old_batch
 	virtual bool still_interested( jd2::archive::Batch const& /*old_batch*/ ) const { return true; };

	///@brief create a new batch with manager().start_new_batch() and manager().finalize_batch();
	virtual void generate_batch() = 0;

	///@brief do some computations on archive that can be done while we are waiting
	virtual void idle() = 0;


	///@brief read 'returned_decoys' from 'batch' into archive.
	virtual void read_structures( core::io::silent::SilentFileData& returned_decoys, Batch const& batch ) = 0;

	///@brief save archive to file .. you can put 'suffix' at end of dirname to save other snapshots than the 'current'
	virtual void save_to_file( std::string suffix = "" ) = 0;
	virtual void save_status( std::ostream& ) const = 0;

	///@brief restore archive
	virtual bool restore_from_file() = 0;

	virtual void init_from_decoy_set( core::io::silent::SilentFileData const& sfd ) = 0;

	///@brief set name of archive ( used also for save_to_file and restore_from_file )
	void set_name( std::string const& set )  { name_ = set; };
	std::string const& name() const { return name_; };

	///@brief access to the ArchiveManager (control of batches)
	ArchiveManager& manager() {
		runtime_assert( manager_ );
		return *manager_;
	}

	virtual
	void set_manager( ArchiveManagerAP manager ) {
		manager_=manager;
	}
protected:

	ArchiveManagerAP manager_ptr() {
		return manager_;
	}

private:
	ArchiveManagerAP manager_;
	std::string name_;
};

class ArchiveBase : public AbstractArchiveBase {
	//to make removal of decoys easy, this might better be a map...

protected:
	typedef std::list< core::io::silent::SilentStructOP > SilentStructs;
	typedef SilentStructs::const_iterator const_decoy_iterator;
	typedef SilentStructs::const_iterator decoy_iterator;

	//silent-struct comment identifiers
	static std::string const TAG_IN_FILE;//( "tag_in_file" );
	static std::string const SOURCE_FILE;//( "source_file" );
public:
	ArchiveBase( ArchiveManagerAP ptr=NULL );
	~ArchiveBase();
	static void register_options();
	virtual bool finished() const { return true; };

	//obsolet ?
	//	virtual bool ready_for_batch() const { return false; };


	virtual void generate_batch() = 0;

	///@brief add structure to Archive.. return false if structure is rejected.
	virtual bool add_structure( core::io::silent::SilentStructOP orig_from_batch );

	///@brief how many structures should be in archive .. varies from decoys().size() in startup phase.
	core::Size nstruct() const { return nstruct_; };

	///@brief set target size of pool
	void set_nstruct( core::Size set ) { nstruct_ = set; };

	///@brief save and restore archive to file-system
	virtual void save_to_file( std::string suffix = "" );
	virtual bool restore_from_file();

	///@brief save and restore status of archive to file-system
	virtual void save_status( std::ostream& ) const;
	virtual void restore_status( std::istream& );

	///@brief called when nothing is happening
	virtual void idle() {};

	///@brief read externally provided structures from decoy_file into archive
	virtual void init_from_decoy_set( core::io::silent::SilentFileData const& sfd );

	///@brief SilentFileData contains the new structures belonging to this batch.
	virtual void read_structures( core::io::silent::SilentFileData&, Batch const& batch );

	///
	///---- methods to keep statistics of acceptance
	///
	core::Size& accepts_since_last_batch() { return accepts_since_last_batch_; };
	core::Size accepts_since_last_batch() const { return accepts_since_last_batch_; };

	//	core::Size& proposed_since_last_batch() { return accepts_since_last_batch_; };
	core::Size proposed_since_last_batch() const { return proposed_since_last_batch_; };
	core::Real current_acceptance_ratio() const {
		return floating_acceptance_ratio_; //will always be upper bound of true acceptance rtio
		//	statistics_valid() ? floating_acceptance_ratio_ : 1.0;
		//		return proposed_since_last_batch_ ? 1.0*accepts_since_last_batch_ / proposed_since_last_batch_ : 1.0;
	}

	void reset_accept_counter() {
		total_accepts_+=accepts_since_last_batch_; accepts_since_last_batch_ = 0;
		total_proposed_ += proposed_since_last_batch_; proposed_since_last_batch_ = 0;
		floating_acceptance_ratio_ = 1.0;
		acceptance_history_.clear();
	}

	core::Size total_proposed() { return total_proposed_ + proposed_since_last_batch(); };
	core::Size total_accepts() { return total_accepts_ + accepts_since_last_batch(); };
	bool statistics_valid() { return  acceptance_history_.size() > min_structures_for_acceptance_statistics_; };


protected:
	virtual void count_structure( Batch const& batch, bool accepted  );
	void count_removed_structures( core::Size n_removed );

	SilentStructs& decoys() { return decoys_; };
	SilentStructs const& decoys() const { return decoys_; };

	///@brief call to insert structure at position given by iterator
	void add_structure_at_position( SilentStructs::iterator iss, core::io::silent::SilentStructOP new_decoy );

private:
	core::Size max_nstruct_; //how many structures maximally maintained in archive
	core::Size nstruct_; //how many structures maintained in archive

	SilentStructs decoys_;

  core::Size accepts_since_last_batch_;
  core::Size total_accepts_;

	core::Size proposed_since_last_batch_;
	core::Size total_proposed_;

	typedef std::deque< bool > AcceptHistoryQueue;
	AcceptHistoryQueue acceptance_history_;

	core::Real floating_acceptance_ratio_;
	core::Size min_structures_for_acceptance_statistics_;

	static bool options_registered_;
};

class DebugArchive : public ArchiveBase {
public:
	DebugArchive( ArchiveManagerAP ptr );

	virtual bool add_structure( core::io::silent::SilentStructOP );

	virtual bool finished() const { return ct_batches_ > 4; };
	//	virtual bool ready_for_batch() const { return !finished() && (ct_batches_ < 3 || decoys().size() > 200); };
	virtual void generate_batch();
	virtual void score( core::pose::Pose& pose ) const;

	virtual void save_status(  std::ostream& ) const;
	virtual void restore_status( std::istream& );
private:
	core::Size ct_batches_;
	core::scoring::ScoreFunctionOP cen_score_;
	bool make_mistake_;
};

}//archive
}//jd2
}//protocols


#endif //INCLUDED_protocols_jd2_Archive_HH
