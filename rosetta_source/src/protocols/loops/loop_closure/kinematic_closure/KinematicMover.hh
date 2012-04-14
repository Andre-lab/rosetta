// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @brief  kinematic closure move
/// @author Daniel J. Mandell (dmandell@itsa.ucsf.edu)
/// @date   Tues Jan 08 12:08:31 2008
///


#ifndef protocols_loops_loop_closure_kinematic_closure_KinematicMover_HH
#define protocols_loops_loop_closure_kinematic_closure_KinematicMover_HH

// Unit Headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.fwd.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.fwd.hh>

// Package Headers
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <utility/LexicographicalIterator.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace loops {
namespace loop_closure {
namespace kinematic_closure {

///@brief A mover class for performing kinematic loop closure on a peptide segment
///@detail
///////////////////////////////////////////////////////////////////////////////
class KinematicMover : public moves::Mover {

public:

	KinematicMover();
	~KinematicMover();

	virtual std::string get_name() const;

	void set_vary_bondangles( bool vary );
	bool get_vary_bondangles();
	void set_sample_nonpivot_torsions( bool sample );
	bool get_sample_nonpivot_torsions();
	void set_sample_vicinity( bool sample_vicinity );
	void set_degree_vicinity( core::Real degree_vicinity );

	Size start_res() const {
		return start_res_; }

	Size middle_res() const {
		return middle_res_; }

	Size end_res() const {
		return end_res_; }

	Size segment_length() const {
		return seg_len_; }

	core::Real BANGLE_MIN() const {
		return BANGLE_MIN_; }

	core::Real BANGLE_SD() const {
		return BANGLE_SD_; }

	void set_sweep_nonpivot_torsions( bool sweep );
	void set_nonpivot_res_to_sweep( utility::vector1< Size > const & resids );
	void set_nonpivot_bb_torsion_id( utility::vector1< Size > const & bbtorids );
	void set_sweep_start_angle( utility::vector1< core::Real > const & angles_in_degrees );
	void set_sweep_step_size( utility::vector1< core::Real > const & angle_steps_in_degrees );
	void set_sweep_nsteps( utility::vector1< Size > const & nsteps );
	/// @details returns true as long as the Lexicographical iterator has not covered all
	/// angles -- useful in a while-loop.
	bool sweep_incomplete() const;

	void set_filters( utility::vector1< protocols::filters::FilterCOP > const & filters_in ){
		filters_ = filters_in; }

	void add_filter( protocols::filters::FilterCOP filter ){
		filters_.push_back( filter ); }

	void clear_filters() {
		filters_.clear(); }

	virtual void set_rama_check ( bool do_rama_check );
	virtual bool get_rama_check ();
	virtual void set_hardsphere_bump_check( bool do_bump_check );
	virtual void set_do_sfxn_eval_every_iteration( bool do_sfxn_eval );
	virtual bool get_hardsphere_bump_check();
	virtual void set_pivots( Size start_res, Size middle_res, Size end_res );
	virtual void apply( core::pose::Pose & );
	virtual void set_idealize_loop_first( bool idealize );
	virtual bool get_idealize_loop_first();
	void set_temperature(core::Real temp_in);
	void set_sfxn(core::scoring::ScoreFunctionCOP sfxn_in);
	bool check_rama(core::Real old_rama_score, core::Real new_rama_score);
	bool last_move_succeeded();
	void set_perturber( KinematicPerturberOP perturber_in );
	core::Real get_bump_overlap_factor() const
  {
    return bump_overlap_factor_;
  }
  void set_bump_overlap_factor(core::Real bump_overlap_factor)
	{
		bump_overlap_factor_=bump_overlap_factor;
	}

private:

	// pivot residues
	Size start_res_, middle_res_, end_res_;
	Size seg_len_;

	//the perturber that sets/samples the chain angles/torsions
	KinematicPerturberOP perturber_;

	// idealization
	core::Real idl_C_N_CA_; // ideal C_N_CA bond angle from Rosetta for idealized kinematic closure
	core::Real idl_N_CA_C_; // ideal N_CA_C bond angle from Rosetta for idealized kinematic closure
	core::Real idl_CA_C_N_; // ideal CA_C_N bond angle from Rosetta for idealized kinematic closure
	core::Real idl_C_N_;  // ideal C_N  bond length from Rosetta for idealized kinematic closure
	core::Real idl_N_CA_;   // ideal N_CA bond length from Rosetta for idealized kinematic closure
	core::Real idl_CA_C_; // ideal CA_C bond length from Rosetta for idealized kinematic closure

	// bond angle sampling
	core::Real BANGLE_MEAN_; // mean N-CA-C bond angle value from PDB
	core::Real BANGLE_SD_; // N-CA-C bond angle standard deviation from PDB
	core::Real BANGLE_MIN_; // min allowed N-CA-C bond angle sampled (otherwise softly enforced by normal distribution)
	core::Real BANGLE_MAX_; // max allowed N-CA-C bond angle sampled (otherwise softly enforced by normal distribution)

	// omega angle sampling
	core::Real OMEGA_MEAN_; // mean omega torsion from loop set
	core::Real OMEGA_SCALE_FACTOR_; // chosen to reproduce observed double-exponential distribution in pdb loops

	core::Real MAX_SAMPLE_ITS_; // maximum number of iterations in torsion / bond angle sampling loop
	bool vary_bond_angles_; // should we vary bond angles
	bool sample_nonpivot_torsions_; // should we sample non-pivot torsions (restricted to Ramachandran space)
	bool sample_vicinity_; //or should we sample around the input conformation
	core::Real degree_vicinity_; //how far away from the input torsions we will sample

	bool sweep_nonpivot_torsion_; // APL: sweep through non-pivot torsions.
	utility::vector1< core::Size > nonpivot_res_to_sweep_;
	utility::vector1< core::Size > sweep_torsion_ids_;
	utility::vector1< core::Real > sweep_nonpivot_torsion_starts_;
	utility::vector1< core::Real > sweep_step_sizes_;
	utility::LexicographicalIterator sweep_iterator_;

	bool idealize_loop_first_; // should we start with an idealized loop
	bool do_rama_check_; // should we do a rama check before accepting loop closure solutions
	// enables use of soft bumps
	bool do_hardsphere_bump_check_; //hard sphere bump check before accepting
	bool do_sfxn_eval_every_iteration_;
	core::scoring::ScoreFunctionCOP sfxn_;
	bool last_move_succeeded_; // did the last move succeed
	core::Real temperature_; // current tempature of the system
	core::Real bump_overlap_factor_; // reduce sum-squared-distance threshold by this factor in bump check

	//filters that are to be applied to every solution before it gets accepted
	utility::vector1< protocols::filters::FilterCOP > filters_;

	// private functions

	bool pivots_within_vicinity(
							core::pose::Pose const & pose,
							utility::vector1<core::Real> const & t_ang,
							utility::vector1<Size> const & pivots,
							Size const start_res,
							Size const middle_res,
							Size const end_res
							);

	// this version checks rama for all residues in loop segment
	bool perform_rama_check( core::pose::Pose const & pose,
							 utility::vector1<core::Real> const & t_ang,
							 utility::vector1<Size> const & pivots,
							 Size const start_res,
							 Size const seg_len
						   );

	// this version only checks rama for pivot residues
	bool perform_rama_check( core::pose::Pose const & pose,
							 utility::vector1<core::Real> const & t_ang,
							 utility::vector1<Size> const & pivots,
							 Size const start_res,
							 Size const middle_res,
							 Size const end_res
							 );
	// checks for backbone-backbone clashes for loop residues
	bool perform_bump_check ( core::pose::Pose const & pose,
							  Size const start_res,
							  Size const end_res
							  );
	// sets default options
	void set_defaults();

};

} // namespace kinematic_closure
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif
