// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_IterativeBase_hh
#define INCLUDED_protocols_abinitio_IterativeBase_hh

// Unit Headers
//#include <protocols/abinitio/IterativeAbrelax.fwd.hh>

// Package Headers
#include <protocols/jd2/archive/EvaluatedArchive.hh>
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>

// Project Headers
#include <protocols/abinitio/PairingStatistics.fwd.hh>
#include <protocols/loops/Loops.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStruct.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <core/fragment/FragSet.fwd.hh>
// AUTO-REMOVED #include <protocols/noesy_assign/NoesyModule.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>


// Third-party Headers
#include <boost/functional/hash.hpp>


//// C++ headers
#include <string>

#include <protocols/noesy_assign/NoesyModule.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class IterativeBase : public jd2::archive::EvaluatedArchive {
	typedef jd2::archive::EvaluatedArchive Parent;
	typedef utility::vector1< core::io::silent::SilentStructOP > SilentStructVector;
public:
	enum IterationStage {
		ENUMERATION = 1,
		TOPO_RESAMPLING,
		PURE_TOPO_RESAMPLING,
		STAGE2_RESAMPLING,
		CEN2FULLATOM,
		//		CEN2FULLATOM_NON_POOL_DECOYS,
		LAST_CENTROID_START = CEN2FULLATOM,
		//		FLEX_CORE_RESAMPLING,
		RIGID_CORE_RESAMPLING,
		FINISHED //keep last
	};

	IterativeBase( jd2::archive::ArchiveManagerAP ptr, std::string name );
	~IterativeBase();

	///@brief archive is finished when at last stage
  virtual bool finished() const { return stage_ >= finish_stage_; };

	///@brief where to stop ?
	void set_finish_stage( IterationStage setting ) {
		finish_stage_ = setting;
	}

	///@brief we are always ready to generate a new batch
	virtual bool ready_for_batch() const { return true; };

	///@brief we are not interested in batches that were generated in old stages
 	virtual bool still_interested( jd2::archive::Batch const& batch ) const;

	///@brief generate a new batch, use different recipe according to current stage

	///@brief generate a new batch, use different recipe according to current stage
  virtual void generate_batch();

	///@brief while waiting for jobs to finish
  virtual void idle();

	virtual void save_status( std::ostream& ) const;
	virtual void restore_status( std::istream& );

	///@brief overloaded to handel special convergence check 'pool_converged_rmsd'
	virtual bool add_structure( core::io::silent::SilentStructOP );

	///@brief setup JumpNrEvaluator
	void setup_default_evaluators();

	///@brief overloaded so we can test for end of IterationStage after reading
	virtual void read_structures( core::io::silent::SilentFileData& returned_decoys, jd2::archive::Batch const& batch );

	///@brief generate flags and stuff for the out-sourced evaluation ---> such that score_final column is returned for each decoy
	/// note needs to be public, since IterativeCentroid calls this from IterativeFullatom to prepare evaluation for soon to be full-atom decoys
	// cen2fullatom-stage ( stage5 )
	virtual void gen_evaluation_output( jd2::archive::Batch& batch, bool fullatom = false );


// ///@brief need to get these from the IterativeCentroid to IterativeFullatom at end of stage5 ;
// 	std::string const& first_noesy_fa_cst_file() const { return first_noesy_fa_cst_file_; }
protected:
 	//void set_first_noesy_fa_cst_file( std::string setting ) { first_noesy_fa_cst_file_ = setting; }

// 	core::Real noesy_assign_float_cycle() const { return noesy_assign_float_cycle_; }
	void set_noesy_assign_float_cycle( core::Real setting ) { noesy_assign_float_cycle_ = setting; }
	bool scored_core_initialized_;
	loops::Loops scored_core_;
	///@brief even in centroid mode the end of abinitio will have a fast relax... enables cs-score and noe-assignment
	bool super_quick_relax_of_centroids_;

	/// ------------- helper functions to be used from generate_batch() --------------------

	void gen_resample_topologies( jd2::archive::Batch& batch );
	void gen_start_structures( jd2::archive::Batch& batch );
	void gen_enumerate_pairings( jd2::archive::Batch& batch );
	void gen_resample_stage2( jd2::archive::Batch& batch );
	void gen_resample_fragments( jd2::archive::Batch& batch );
	void gen_cen2fullatom( jd2::archive::Batch& batch );
	void gen_cen2fullatom_non_pool_decoys( jd2::archive::Batch& batch );
	void add_fullatom_flags( jd2::archive::Batch& batch );

	/// actually run the assignment machinery (only after batch is started to keep archive from hogging the queue... )
	void reassign_noesy_data( jd2::archive::Batch& batch );

	/// generate cst-input from current assigned noesy data
	void gen_noe_assignments( jd2::archive::Batch& batch );

	/// some helpers for the helpers
	PairingStatisticsOP compute_beta_topology();
	void guess_pairings_from_secondary_structure(
	    core::fragment::FragSet const& frags,
			std::string const& out_pairings_file,
			std::string const& out_frag_ss_file
	) const;
	void compute_cores();

	///these are set by the cmd-line options iterative::fa_score and iterative::fa_score_patch
	std::string const& fa_score() const {
		return fa_score_;
	}
	std::string const& fa_score_patch() const {
		return fa_score_patch_;
	}

	///these are set by the cmd-line options iterative::cen_score and iterative::cen_score_patch
	std::string const& cen_score() const {
		return cen_score_;
	}
	std::string const& cen_score_patch() const {
		return cen_score_patch_;
	}

	///OBSOLET cores are computed by compute_cores() in idle()
	loops::Loops const& core( core::Size i ) {
		if ( i == 2 ) { return core2_; };
		if ( i == 3 ) { return core3_; };
		if ( i == 4 ) { return core4_; };
		runtime_assert( false );
		return core2_; //happy compiler
	}

	///@brief current stage?
	IterationStage stage() const {
		return stage_;
	}

	///@brief needed for writing of psi-pred fiels (guess_pairings_from_secondary_structure)
	std::string const& target_sequence() const {
		return target_sequence_;
	}

	void set_stage( IterationStage setting ) {
		stage_ = setting;
	}


	///@brief cluster structures with min_diversity_list_[ stage_ ] as cluster:radius
	void cluster();

	std::string const& chemshift_column() const {
		return chemshift_column_;
	}



private:
	void collect_alternative_decoys( SilentStructs primary_decoys, std::string alternative_decoy_file, SilentStructVector& output_decoys );

	///@brief what is the expected lowest acceptance ratio at the current stage ?
	core::Real target_accept_ratio() const { return target_accept_ratio_[ stage_ ]; }

	///@brief [OBSOLET] add score_coreX and rms_coreX evaluators (and columns) with 0.0 weight
	void add_core_evaluator( loops::Loops const& core, std::string const& core_tag );

	///@brief restrict scoring to core-regions
	void set_scored_core();

	///@brief  calls increment_stage() if appropriate
	void test_for_stage_end();

	///@brief necessary steps to go to next stage... e.g., saving snapshot of archive
	void increment_stage();

	void read_noisy_assing_data_from_last_batch();

private:
	///  ----------------- -- private data members -- --------------------

	/// ------------------------------ stage - control --------------------------
	///@brief current stage
	IterationStage stage_;

	///@brief end-condition
	IterationStage finish_stage_;

	///@brief indices of prominent batches ( STATUS file )
	core::Size first_batch_this_stage_;
	core::Size first_fullatom_batch_;

	/// --------------------------------- other ------------------------------------
	///@brief toggle to keep track of the enumerate-pairings mode .. want to run this only every 2nd batch
	bool bEnumeratedLastTime_;

	///@brief [OBSOLET?] keep track when idle() has been run ...
	core::Size last_accepted_decoys_in_idle_;

	///core-regions --- used in IterativeFullatom for the "rigid-core" sampling step...
	loops::Loops core2_;
	loops::Loops core3_;
	loops::Loops core4_;

	/// ----------- some cmd-line controlled settings -----------------

	///@brief how many structures are maximally produced in stage X
	utility::vector1< int > max_nstruct_list_; //-1 --> skip stage, 0 infinite, N>0 make a maximum of N structures.

	///@brief cluster:radius for minimum diversity in stage X
	utility::vector1< core::Real > min_diversity_list_;

	///@brief minimum acceptance ratio .. current_accept_ratio < target_accept_ratio[ stage_ ] --> increment_stage
	utility::vector1< core::Real > target_accept_ratio_;

	///@brief for RMSD ... e.g., in:file:native
	core::pose::PoseCOP reference_pose_;

	///@brief from cmd-line cen-score and patch names
	std::string cen_score_;
	std::string cen_score_patch_;

	///@brief from cmd-line fa-score and patch names
	std::string fa_score_;
	std::string fa_score_patch_;

	std::string target_sequence_; //read from in:file:fasta in c'stor

	protocols::noesy_assign::NoesyModuleOP noesy_module_;
	core::Real noesy_assign_float_cycle_;

	std::string first_noesy_cst_file_;
	std::string first_noesy_fa_cst_file_;

	std::string current_noesy_sampling_file_;
	bool bCombineNoesyCst_;


	//hash value to see if decoys for noesy assign have changed
	size_t noesy_assign_hash_;
	boost::hash<std::string> hasher;

	///@brief even in centroid mode the end of abinitio will have a fast relax... enables cs-score and noe-assignment
	//bool super_quick_relax_of_centroids_;
	bool recover_centroid_structure_for_pool_;

	std::string chemshift_column_;
	bool bDoBetaJumping_;


	/// ------------------ register cmdline options ---------------------------

private:
	static bool options_registered_;
public:
	static void register_options();

};


}
}

#endif
