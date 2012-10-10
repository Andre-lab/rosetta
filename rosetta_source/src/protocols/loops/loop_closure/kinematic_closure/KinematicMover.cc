// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Performs a kinematic closure move on a peptide segment
/// @author Daniel J. Mandell
/// @author Amelie Stein (amelie.stein@ucsf.edu), Oct 2012 -- next-generation KIC

// Unit Headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
//#include <protocols/loops/kinematic_closure/KinematicPerturber.hh>

// Rosetta Headers
#include <core/conformation/Conformation.hh> // DJM: may go away
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh> // DJM: may go away
#include <core/conformation/ResidueFactory.hh> // DJM: may go away
#include <core/conformation/util.hh>
// AUTO-REMOVED #include <core/fragment/BBTorsionAndAnglesSRFD.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/kinematics/DomainMap.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
//#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
//#include <core/scoring/etable/count_pair/CountPairAll.hh>
//#include <core/scoring/etable/count_pair/CountPairFunction.hh>
//#include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/NeighborList.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <numeric/xyzVector.hh>

// option key includes
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>
#include <core/id/TorsionID.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>
#include <utility/vector1.hh>
#include <numeric/conversions.hh>
#include <boost/foreach.hpp>

//Auto Headers
#define foreach BOOST_FOREACH
static numeric::random::RandomGenerator RG(43135);
static basic::Tracer TR("protocols.loops.loop_closure.kinematic_closure.KinematicMover");

using namespace numeric::kinematic_closure;
using namespace core;
using namespace core::kinematics;
using namespace core::scoring;

namespace protocols {
namespace loops {
namespace loop_closure {
namespace kinematic_closure {

// default constructor
KinematicMover::KinematicMover() :
	Mover(),
	start_res_(2),
	middle_res_(3),
	end_res_(4),
	seg_len_(3), // dummy value -- was previously uninitialized (!)
	loop_begin_(2), // AS: start of the full loop -- dummy value
	loop_end_(4),
	idl_C_N_CA_(121.7),
	idl_N_CA_C_(111.2),
	idl_CA_C_N_(116.2),
	idl_C_N_(1.32869),
	idl_N_CA_(1.458),
	idl_CA_C_(1.52326),
	BANGLE_MEAN_(110.86), // from PDB
	//BANGLE_MEAN=107.0;    // from CHARMM (to coincide with mm_bend potential)
	BANGLE_SD_(2.48),
	BANGLE_MIN_(BANGLE_MEAN_ - (0.5 * BANGLE_SD_)),
	//BANGLE_MIN = BANGLE_MEAN - (2.0 * BANGLE_SD);
	//BANGLE_MIN = 102.7284756;
	//BANGLE_MIN = 98.86848; // (to coincide with mm_bend potential)
	BANGLE_MAX_(118.9353622),
	//BANGLE_MAX = 115.07536; // (to coincide with mm_bend potential)
	OMEGA_MEAN_(179.8),
	OMEGA_SCALE_FACTOR_(2.375), // scale factor parameterized to reproduce pdb loop omega angles
	//OMEGA_SD=9.15;
	//OMEGA_MIN = OMEGA_MEAN - (0.5 * OMEGA_SD);
	//MAX_SAMPLE_ITS_(2000), // DJM: deprecated in favor of loops::max_kic_perturber_samples
	vary_bond_angles_(false),
	sample_nonpivot_torsions_(true),
	sweep_nonpivot_torsion_(false),
	idealize_loop_first_(false),
	do_rama_check_(false),
	do_hardsphere_bump_check_(true),
	do_sfxn_eval_every_iteration_(false),
	sfxn_( NULL ),
	last_move_succeeded_ (false),
	temperature_(1.0),
	bump_overlap_factor_(0.49), // 0.6; // 0.8; // 0.8^2, allows some atomic overlap
	taboo_map_max_capacity_(0.95) // how full can the taboo map get before we re-set it? --> the higher, the more likely are almost-endless loops -- but if the value is too low it's not actual taboo sampling

	{
		perturber_ = new kinematic_closure::TorsionSamplingKinematicPerturber( this );
		Mover::type( "KinematicMover" );
		set_defaults(); // doesn't do anything in this implementation
	}

// default destructor
KinematicMover::~KinematicMover() {}

void KinematicMover::set_pivots( Size start_res, Size middle_res, Size end_res )
{
	start_res_	= start_res;
	middle_res_	= middle_res;
	end_res_	= end_res;

	seg_len_ = end_res_ - start_res_ + 1;
	//if these are bad data, we're in trouble!
	runtime_assert(start_res != 0);
	runtime_assert(middle_res != 0);
	runtime_assert(end_res != 0);
	runtime_assert(start_res != middle_res);
	runtime_assert(middle_res != end_res);
	runtime_assert(end_res != start_res);
}

// boundaries of the current loop, needed for torsion bin lookup and offset in torsion-restricted & taboo sampling
void KinematicMover::set_loop_begin_and_end( Size loop_begin, Size loop_end ) 
{
	loop_begin_ = loop_begin;
	loop_end_ = loop_end;
}

void KinematicMover::set_vary_bondangles( bool vary )
{
	KinematicMover::vary_bond_angles_=vary;
}

bool KinematicMover::get_vary_bondangles()
{
	return KinematicMover::vary_bond_angles_;
}

void
KinematicMover::set_perturber( KinematicPerturberOP perturber_in )
{
	perturber_ = perturber_in;
	perturber_->set_kinmover( this );
}

void KinematicMover::set_sample_nonpivot_torsions( bool sample )
{
	KinematicMover::sample_nonpivot_torsions_=sample;
}
	
bool KinematicMover::get_sample_nonpivot_torsions()
{
	return KinematicMover::sample_nonpivot_torsions_;
}

void KinematicMover::set_idealize_loop_first( bool idealize )
{
	KinematicMover::idealize_loop_first_ = idealize;
}

bool KinematicMover::get_idealize_loop_first()
{
	return KinematicMover::idealize_loop_first_;
}

void KinematicMover::set_rama_check( bool do_rama )
{
	KinematicMover::do_rama_check_=do_rama;
}

bool KinematicMover::get_rama_check()
{
	return KinematicMover::do_rama_check_;
}

void KinematicMover::set_hardsphere_bump_check( bool do_bump_check )
{
	do_hardsphere_bump_check_ = do_bump_check;
}

bool KinematicMover::get_hardsphere_bump_check()
{
	return do_hardsphere_bump_check_;
}

void KinematicMover::set_do_sfxn_eval_every_iteration( bool do_sfxn_eval )
{
	do_sfxn_eval_every_iteration_ = do_sfxn_eval;
}

void KinematicMover::set_sfxn(core::scoring::ScoreFunctionCOP sfxn_in)
{
	sfxn_ = sfxn_in;
}



void KinematicMover::set_sweep_nonpivot_torsions( bool sweep ) {
	sweep_nonpivot_torsion_ = sweep;
}

/// @details This must be set before the other properties
/// (torsion angles, starting angles, step sizes, nsteps)
/// are set.
void KinematicMover::set_nonpivot_res_to_sweep( utility::vector1< Size > const & resids )
{
	nonpivot_res_to_sweep_ = resids;
}

void KinematicMover::set_nonpivot_bb_torsion_id( utility::vector1< Size > const & bbtorids )
{
	assert( nonpivot_res_to_sweep_.size() == bbtorids.size() );
	sweep_torsion_ids_ = bbtorids;
}

void KinematicMover::set_sweep_start_angle( utility::vector1< core::Real > const & angles_in_degrees )
{
	assert( nonpivot_res_to_sweep_.size() == angles_in_degrees.size() );
	sweep_nonpivot_torsion_starts_ = angles_in_degrees;
}

void KinematicMover::set_sweep_step_size( utility::vector1< core::Real > const & angle_steps_in_degrees )
{
	assert( nonpivot_res_to_sweep_.size() == angle_steps_in_degrees.size() );
	sweep_step_sizes_ = angle_steps_in_degrees;
}

/// @details Initializes the LexicographicalIterator
void KinematicMover::set_sweep_nsteps( utility::vector1< core::Size > const & nsteps )
{
	assert( nonpivot_res_to_sweep_.size() == nsteps.size() );
	sweep_iterator_.set_dimension_sizes( nsteps );
}

bool KinematicMover::sweep_incomplete() const
{
	return !sweep_iterator_.at_end();
}



void KinematicMover::set_temperature( Real temp_in )
{
	temperature_=temp_in;
}

bool KinematicMover::check_rama(Real old_rama_score, Real new_rama_score) {
	if ( new_rama_score > old_rama_score ) {
		Real const boltz_factor = (old_rama_score - new_rama_score)/temperature_;
		Real const probability = std::exp(std::max(Real(-40.0),boltz_factor) );
		if ( RG.uniform() >= probability ) return( false );
	}
	return( true );
}

bool KinematicMover::last_move_succeeded()
{
	return last_move_succeeded_;
}


void KinematicMover::apply( core::pose::Pose & pose )
{
	//scoring::Ramachandran const & rama( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
	using numeric::conversions::radians;
	using numeric::conversions::degrees;
	using numeric::random::RG;
	using core::id::AtomID;
	last_move_succeeded_=false;

	// DJM: debug
	//TR << "from " << start_res_ << " to " << end_res_ << std::endl;

	// inputs to loop closure
	utility::vector1<utility::vector1<Real> > atoms;
	utility::vector1<Size> pivots (3), order (3);
	// outputs from loop closure
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, save_t_ang, save_b_len, save_b_ang, R0 (3);
	utility::vector1<conformation::ResidueOP> save_residues;

	Size middle_offset = middle_res_ - start_res_; // is used to set central pivot atom
	Size seg_len = 1 + end_res_ - start_res_; // length of the closure chain
	//Size detail_sol = 0; // index of the KIC solution accepted for detailed balance, if > 0
	atoms.resize((seg_len + 2) * 3); // one extra residue on each side to establish the geometric frame

	// KC requires 3 backbone atoms 1 residue N-terminal to start pivot and 1 residue C-terminal to end pivot
	// so check for terminal pivots and compute coordinates from prepended / appended residues if necessary
	if (pose.residue(start_res_).is_lower_terminus()) {
		// make a copy of the lower terminal residue
		conformation::ResidueOP nterm_copy( conformation::ResidueFactory::create_residue(
																						  pose.conformation().residue(start_res_).type(),
																						  pose.conformation().residue(start_res_),
																						  pose.conformation(),
																						  false ) );
		// now create another residue of the same type to prepend
		conformation::ResidueOP pre_nterm( conformation::ResidueFactory::create_residue(
																						pose.residue(start_res_).type() ) );
		// create a new conformation containing the n-term copy
		conformation::Conformation extn;
		extn.append_residue_by_bond( *nterm_copy );
		// attached the pre_nterm residue to the nterm copy
		extn.safely_prepend_polymer_residue_before_seqpos( *pre_nterm, 1, true );
		// set ideal omega angle at junction
		extn.set_torsion( id::TorsionID( 1, id::BB, 3 ),  OMEGA_MEAN_);
		// store the xyz coords for KC from the first 3 atoms of the pre-nterm residue into the first 3 atoms indices
		Size ind=1;
		for (Size j=1; j<=3; j++) {
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (extn.residue(1).xyz(j).x());
			atoms[ind][2] = static_cast<Real> (extn.residue(1).xyz(j).y());
			atoms[ind][3] = static_cast<Real> (extn.residue(1).xyz(j).z());
			ind++;
		}
	}

	// check for upper terminus
	if (pose.residue(end_res_).is_upper_terminus()) {
		// make a copy of the upper terminal residue
		conformation::ResidueOP cterm_copy( conformation::ResidueFactory::create_residue(
																						  pose.conformation().residue(end_res_).type(),
																						  pose.conformation().residue(end_res_),
																						  pose.conformation(),
																						  false ) );

		// now create another residue of the same type to append
		conformation::ResidueOP post_cterm( conformation::ResidueFactory::create_residue(
																						 pose.residue(end_res_).type() ) );
		// create a new conformation containing the c-term copy
		conformation::Conformation extn;
		extn.append_residue_by_bond( *cterm_copy );
		// attached the post_cterm residue to the cterm copy
		extn.safely_append_polymer_residue_after_seqpos( *post_cterm, 1, true );
		// set ideal omega angle at junction
		extn.set_torsion( id::TorsionID( 2, id::BB, 3 ),  OMEGA_MEAN_);
		// store the xyz coords for KC from the first 3 atoms of the post-cterm residue into the last 3 atoms indices
		Size ind=atoms.size()-2;
		for (Size j=1; j<=3; j++) {
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (extn.residue(2).xyz(j).x());
			atoms[ind][2] = static_cast<Real> (extn.residue(2).xyz(j).y());
			atoms[ind][3] = static_cast<Real> (extn.residue(2).xyz(j).z());
			ind++;
		}
	}

	// Get the current coords of the loop closure segment (if a pivot is terminal we took the coords above so skip)
	Size ind = (pose.residue(start_res_).is_lower_terminus() ? 4 : 1);
	for (Size i =  (pose.residue(start_res_).is_lower_terminus() ? start_res_ : start_res_ - 1);
		 i <= (pose.residue(end_res_).is_upper_terminus() ? end_res_ : end_res_ + 1);
		 i++) {
		conformation::Residue res=pose.residue(i);
		for (Size j=1; j<=3; j++) { // DJM: just keeping N, CA, C atoms. We assume these are always the first 3.  BAD -- PROTEIN ONLY ASSUMPTION -- How about metal ions with only 1 atom?
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (res.xyz(j).x());
			atoms[ind][2] = static_cast<Real> (res.xyz(j).y());
			atoms[ind][3] = static_cast<Real> (res.xyz(j).z());
			ind++;
		}
	}

	/// Get alternative conformations from from kinematic loop closure
	// Choose order to solve half-tangents
	order[1]=1;
	order[2]=2;
	order[3]=3;
	// Set the pivot atoms
	Size pvatom1=5; // second C-alpha
	Size pvatom2=5 + (3 * middle_offset); // middle res C-alpha
	Size pvatom3=(3 * (seg_len+1)) - 1; // second-to-last C-alpha
	pivots[1]=pvatom1;
	pivots[2]=pvatom2;
	pivots[3]=pvatom3;


	// chainTORS is used to get dt to remove solutions identical to original torsions
	// DJM: it would be quicker to store the last torsions rather than recomputing them,
	// but we use db_ang and db_len below so it's worth it for now. However, we should
	// be able to get the torsions, bond angles, and bond lengths from the atom_tree now.
	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);

	// save the previous torsions, lengths and angles to restore if no solution is found
	if (! idealize_loop_first_) {
		save_t_ang.resize(dt_ang.size());
		save_b_ang.resize(db_ang.size());
		save_b_len.resize(db_len.size());
		for (Size i=1; i<=dt_ang.size(); i++) {
			save_t_ang[i] = dt_ang[i];
			save_b_ang[i] = db_ang[i];
			save_b_len[i] = db_len[i];
		}
	}

	// idealize the loop if requested
	if (idealize_loop_first_) {
		// save a backbup of the non-ideal residues in case closure fails
		save_residues.resize(seg_len);
		for (Size seqpos=start_res_, i=1; seqpos<=end_res_; seqpos++, i++) {
			conformation::ResidueOP saveres=pose.conformation().residue( seqpos ).clone();
			save_residues[i]=saveres;
		}

		// set all bond lengths, angles, and omegas for closure to ideal values

		// check sizes to prevent memory crashes:
		for (Size i=1; i<=atoms.size(); i+=3) {
			if( db_ang.size() >= (i) )   db_ang[i]=idl_C_N_CA_;
			if( db_ang.size() >= (i+1) ) db_ang[i+1]=idl_N_CA_C_;
			if( db_ang.size() >= (i+2) ) db_ang[i+2]=idl_CA_C_N_;
			if( db_len.size() >= (i) )   db_len[i]=idl_N_CA_;
			if( db_len.size() >= (i+1) ) db_len[i+1]=idl_CA_C_;
			if( db_len.size() >= (i+2) ) db_len[i+2]=idl_C_N_;
			if( dt_ang.size() >= (i+2) ) dt_ang[i+2]=OMEGA_MEAN_;
		}
		// setting pose omegas here, don't set them later!
		for ( Size i= start_res_; i<= end_res_; ++i ) {
			conformation::idealize_position( i, pose.conformation() ); // this is coupled to above values!
			pose.set_omega(i, OMEGA_MEAN_);
		}
	}

	// If preserving detailed balance, compute the Rosenbluth factor for the before state
	//	W_before = rosenbluth_factor( atoms, save_t_ang, save_b_ang, save_b_len, pivots, order, t_ang, b_ang, b_len, nsol );
	// DJM: currently bridgeObjects will be applying save* transformations to atoms, which isn't necessary. Maybe in option in bridge not to?

	for (Size nits=1; nits <= perturber_->max_sample_iterations(); ++nits ) { // try these pivots until a solution passes all filters or nits > perturber_->max_sample_iterations()

		//make sure the perturber can still generate solutions
		if( perturber_->perturber_exhausted() ) break;

		perturber_->perturb_chain( pose, dt_ang, db_ang, db_len );

		// Perform loop closure
		//TR.Debug << "calling bridgeObjects" << std::endl;
		bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);

		// If preserving detailed balance, try performing detailed balance move here
		//if( preserve_detailed_balance_ ) {
		//detail_sol = check_detailed_balance( W_before, t_ang, b_ang, b_len, nsol ); // DJM: should be able to used chainXYZ to recreate the segments
		//if detail_sol != 0 { // accepted a move by detailed balance
				// DJM: next want to pre-screen for detailed balance acceptance, or only consider solns for detailed balance check if pass sterics and Rama
				//last_move_succeeded_ = true;
				//return;
		//}
		//else continue;
		//}

		utility::vector1<Size> pos(nsol);
		for (int i=1; i<=nsol; i++) {
			pos[i]=i;
		}
		
		//std::random__shuffle(pos.begin(), pos.end());
		numeric::random::random_permutation(pos.begin(), pos.end(), RG);
		
		// Look for solutions passing NaN, Rama, bump checks and eventual filters
		for (Size i=nsol; i>=1; i--) {
			//TR << "KINMOVER checking sol " << nsol - i + 1 << " of " << nsol << " .. "; // DJM: debug
			// make sure no torsions are NaN, Inf, or otherwise unreasonable
			if ( ! ( ( -360.0 <= t_ang[pos[i]][pivots[1]-1] ) && (t_ang[pos[i]][pivots[1]-1] <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[1]]   ) && (t_ang[pos[i]][pivots[1]]   <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[2]-1] ) && (t_ang[pos[i]][pivots[2]-1] <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[2]]   ) && (t_ang[pos[i]][pivots[2]]   <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[3]-1] ) && (t_ang[pos[i]][pivots[3]-1] <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[3]]   ) && (t_ang[pos[i]][pivots[3]]   <= 360.0) ) ) {
				TR << "solution " << i << " of " << nsol << " failed NaN / Inf check. Skipping..." << std::endl;
				TR << "failed NaN t_ang matrix: " << std::endl;
				printMatrix(t_ang);
				continue;
			}

			//if (!perform_rama_check(pose, t_ang[i], pivots, start_res, seg_len)) { // checks all loop residues
			if (!perform_rama_check(pose, t_ang[pos[i]], pivots, start_res_, middle_res_, end_res_)) { // checks only pivot residues
				//TR << "Rama continue" << std::endl; // DJM: debug
				continue;
			}
			// place the solution into the pose and bump check+eventual filters
			// set the torsions
			perturber_->set_pose_after_closure( pose, t_ang[pos[i]], b_ang[pos[i]], b_len[pos[i]], true );

			// AS: fetch & set the torsion string for Taboo Sampling 
			insert_sampled_torsion_string_into_taboo_map( torsion_features_string( pose ) );
			
			
			//now check if the pose passes all the filters
			if( do_hardsphere_bump_check_ && !perform_bump_check(pose, start_res_, end_res_) ){
				continue;
			}

			if( do_sfxn_eval_every_iteration_ ) (*sfxn_)( pose );

			if( filters_.size() != 0 ){
				bool all_filters_passed( true );
				foreach(protocols::filters::FilterCOP filter, filters_){
					if( ! filter->apply( pose ) ) {
						all_filters_passed = false;
						break;
					}
				}
				if ( !all_filters_passed ){
					continue;
				}
			}

			last_move_succeeded_ = true;
			//std::cerr << "kinmover success after " << nits << "...   ";
			return;
		}
	} //for (Size nits=1; nits <= MAX_SAMPLE_ITS_; ++nits )

	// prepare to exit here if no solutions
	last_move_succeeded_ = false;
	//TR << "!!!last move not succeeded to pass filters" << std::endl;

	//// no solution found. restore the old pose
	if (! idealize_loop_first_) { // if we didn't idealize, we only changed the torsions and N_CA_C bond angles

		perturber_->set_pose_after_closure( pose, save_t_ang, save_b_ang, save_b_len, false );

	}
	else { // loop residues were idealized so restore from saved non-ideal residues
		for (Size seqpos=start_res_, i=1; seqpos<=end_res_; seqpos++, i++) {
			pose.conformation().replace_residue(seqpos, *save_residues[i], false);
		}
	}

	//TR << "\tnumber of its: " << nits << std::endl;
	return;

} //new_apply

std::string
KinematicMover::get_name() const {
	return "KinematicMover";
}

// checks if pivot solutions are within vicinity degrees of current pivot torsion values
// DJM: currently not used because it significantly reduces the number of accepted moves
	/* AS Oct 03, 2012 -- commenting this out for vicinity refactoring
bool KinematicMover::pivots_within_vicinity(
	core::pose::Pose const & pose,
	utility::vector1<Real> const & t_ang,
	utility::vector1<Size> const & pivots,
	Size const start_res,
	Size const middle_res,
	Size const end_res
	)
{
	Real old_phi, new_phi, old_psi, new_psi; // torsions before and after the move is applied

	// check start_res
	new_phi = t_ang[pivots[1]-1];
	new_psi = t_ang[pivots[1]];
	old_phi = pose.phi(start_res);
	old_psi = pose.psi(start_res);
	if ( ( std::abs( new_phi - old_phi ) > degree_vicinity_ ) || ( std::abs( new_psi - old_psi ) > degree_vicinity_ ) )
		return false;

	// check middle_res
	new_phi = t_ang[pivots[2]-1];
	new_psi = t_ang[pivots[2]];
	old_phi = pose.phi(middle_res);
	old_psi = pose.psi(middle_res);
	if ( ( std::abs( new_phi - old_phi ) > degree_vicinity_ ) || ( std::abs( new_psi - old_psi ) > degree_vicinity_ ) )
		return false;

	// check end_res
	new_phi = t_ang[pivots[3]-1];
	new_psi = t_ang[pivots[3]];
	old_phi = pose.phi(end_res);
	old_psi = pose.psi(end_res);
	if ( ( std::abs( new_phi - old_phi ) > degree_vicinity_ ) || ( std::abs( new_psi - old_psi ) > degree_vicinity_ ) )
		return false;

	// all passed
	return true;
}
*/
	
// this version checks rama for all residues in loop segment
bool KinematicMover::perform_rama_check(
	core::pose::Pose const & pose,
	utility::vector1<Real> const & t_ang,
	utility::vector1<Size> const & pivots,
	Size const start_res,
	Size const seg_len
)
{
	core::scoring::Ramachandran const & rama( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
	Real old_phi, new_phi, old_psi, new_psi; // torsions before and after the move is applied
	Real old_rama_score, new_rama_score; // rama scores before and after the move is applied

	for (Size offset=0; offset<seg_len; offset++) { // offset is natoms into kic segment
		conformation::Residue const & current_rsd( pose.residue(start_res+offset) );
		new_phi = t_ang[pivots[1]+(3*offset)-1];
		new_psi = t_ang[pivots[1]+(3*offset)];
		new_rama_score = rama.eval_rama_score_residue( current_rsd.aa(), new_phi, new_psi );
		if (new_rama_score == 20.0) {
			//std::cout << "Rama fail" << std::endl;
			return false;
		}
		old_psi = pose.psi(start_res+offset);
		old_phi = pose.phi(start_res+offset);
		old_rama_score = rama.eval_rama_score_residue( current_rsd.aa(), old_phi, old_psi );
		if (!check_rama(old_rama_score, new_rama_score)) {
			// DJM: debug
			// TR << "Residue " << start_res + offset << " failed Rama. old_rama_score: " << old_rama_score
			//    << " new_rama_score " << new_rama_score << std::endl;
			///std::cout << "Rama fail" << std::endl;
			return false; // failed rama check
		}
	}
	return true; // passed rama check
}

// this version only checks rama for pivot residues
bool KinematicMover::perform_rama_check(
	core::pose::Pose const & pose,
	utility::vector1<Real> const & t_ang,
	utility::vector1<Size> const & pivots,
	Size const start_res,
	Size const middle_res,
	Size const end_res
)
{
	core::scoring::Ramachandran const & rama( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
	Real old_phi, new_phi, old_psi, new_psi; // torsions before and after the move is applied
	Real old_rama_score, new_rama_score; // rama scores before and after the move is applied

	// check start_res
	conformation::Residue const & start_rsd( pose.residue(start_res) );
	new_phi = t_ang[pivots[1]-1];
	new_psi = t_ang[pivots[1]];
	new_rama_score = rama.eval_rama_score_residue( start_rsd.aa(), new_phi, new_psi );
	if (new_rama_score == 20.0) {
		//TR << "score for pivot 1 was 20. returning false" << std::endl;
		///std::cout << "Rama fail" << std::endl;
		return false;
	}
	old_phi = pose.phi(start_res);
	old_psi = pose.psi(start_res);
	old_rama_score = rama.eval_rama_score_residue( start_rsd.aa(), old_phi, old_psi );
	// DJM: debug
	//TR << "new_phi, new_psi, rama score: " << new_phi << " " << new_psi << " " << new_rama_score << std::endl;
	if (!check_rama(old_rama_score, new_rama_score)) {
		// DJM: debug
		// TR << "Residue " << start_res + offset << " failed Rama. old_rama_score: " << old_rama_score
		//    << " new_rama_score " << new_rama_score << std::endl;
		///std::cout << "Rama fail" << std::endl;
		return false; // failed rama check
	}

	// check middle_res
	conformation::Residue const & middle_rsd( pose.residue(middle_res) );
	new_phi = t_ang[pivots[2]-1];
	new_psi = t_ang[pivots[2]];
	new_rama_score = rama.eval_rama_score_residue( middle_rsd.aa(), new_phi, new_psi );
	if (new_rama_score == 20.0) {
		//TR << "score for pivot 2 was 20. returning false" << std::endl;
		//std::cout << "Rama fail" << std::endl;
		return false;
	}
	old_phi = pose.phi(middle_res);
	old_psi = pose.psi(middle_res);
	old_rama_score = rama.eval_rama_score_residue( middle_rsd.aa(), old_phi, old_psi );
	if (!check_rama(old_rama_score, new_rama_score)) {
		// DJM: debug
		// TR << "Residue " << start_res + offset << " failed Rama. old_rama_score: " << old_rama_score
		//    << " new_rama_score " << new_rama_score << std::endl;
		///std::cout << "Rama fail" << std::endl;
		return false; // failed rama check
	}

	// check end_res
	conformation::Residue const & end_rsd( pose.residue(end_res) );
	new_phi = t_ang[pivots[3]-1];
	new_psi = t_ang[pivots[3]];
	new_rama_score = rama.eval_rama_score_residue( end_rsd.aa(), new_phi, new_psi );
	if (new_rama_score == 20.0) {
		//TR << "score for pivot 3 was 20. returning false" << std::endl;
		//std::cout << "Rama fail" << std::endl;
		return false;
	}
	old_phi = pose.phi(end_res);
	old_psi = pose.psi(end_res);
	old_rama_score = rama.eval_rama_score_residue( end_rsd.aa(), old_phi, old_psi );
	if (!check_rama(old_rama_score, new_rama_score)) {
		// DJM: debug
		// TR << "Residue " << start_res + offset << " failed Rama. old_rama_score: " << old_rama_score
		//    << " new_rama_score " << new_rama_score << std::endl;
		//std::cout << "Rama fail" << std::endl;
		return false; // failed rama check
	}

	// DJM: debug - print out all the rama score
/*
	TR << "All residues passed Rama. Scores follow:" << std::endl;
	for (Size i=start_res; i<=end_res; i++) {
		Real cur_phi = pose.phi(i);
		Real cur_psi = pose.psi(i);
		conformation::Residue const & cur_rsd( pose.residue(i) );
		Real cur_rama_score = rama.eval_rama_score_residue( cur_rsd.aa(), cur_phi, cur_psi);
		TR << "Res " << i << ": " << cur_rama_score << std::endl;
	}
*/
	return true; // passed rama check
}


// perform_bump_check
// 1. loop over all residues in the loop
// 2. for each loop residue, grab the neighbor atom's xyz
//	3. for every other residue res2 in the protein:
//		4. grab the neighbor atom of res2
//		5. grab the neighbor radii of res1.nbratom and res2.nbratom (make sure these are sensible)
//		6. see if res1 and res2 are close enough to be considered neighbors. if so:
//			7. for each main_chain atom in res1,
//				8. for each main_chain atom in res2
//					9. if any of these atoms is within dist_cutoff, clash is true

// returns true if bump check is passed, false if it fails
bool
KinematicMover::perform_bump_check(
	core::pose::Pose const & pose,
	Size const start_res,
	Size const end_res
)
{
	// iterate over loop residues
	//bool passed_but_no_neighbor=true; // DJM: neighbor filter
	//int no_neighbor_res=0; // DJM: neighbor filter
	for (Size i=start_res; i<= end_res; i++) {
		conformation::Residue const & rsd1 = pose.conformation().residue(i);
		// get the xyz vector of the atom of this residue used for neighbor detection
		Vector const & nbr_i( rsd1.xyz( rsd1.nbr_atom() ) );
		// iterate over all other residues in the protein
		//bool had_neighbor = false; // DJM: neighbor filter
		for (Size j=1; j <= pose.total_residue(); j++ ) {
			if ((j == i) || (j == i+1) || (j == i-1)) continue; // don't do same or adjacent residues
			if ((j >= start_res) && (j <= i)) continue; // don't do loop residues multiple times
			conformation::Residue const & rsd2 = pose.conformation().residue(j);
			// get the xyz vector fo the atom of this residue used for neighbor detection
			Vector const & nbr_j( rsd2.xyz( rsd2.nbr_atom() ) );
			// determine the neighbor cutoff from the radii of the neighbor atoms
			Real const nbrcutoff = ( rsd1.nbr_radius() + rsd2.nbr_radius()) ;
			Real const nbrcutoff2 = nbrcutoff * nbrcutoff; // use squared neighbor cutoff
			//if ( (nbr_i - nbr_j).length_squared() < nbrcutoff2 * 2 ) had_neighbor = true; // DJM: neighbor filter
			if ( (nbr_i - nbr_j).length_squared() > nbrcutoff2 ) continue;
			//else had_neighbor = true; // DJM: neighbor filter
			// now check for clashes between the atoms of the two residues
			Size max_res1_heavyatoms = std::min( Size (5), rsd1.nheavyatoms() ); // keep N,CA,C,O (and CB if not Gly)
			for (Size m=1; m<= max_res1_heavyatoms; m++) {
				Size max_res2_heavyatoms;
				if ( rsd2.is_protein() ) max_res2_heavyatoms = std::min( Size (5), rsd2.nheavyatoms() );
				else max_res2_heavyatoms = rsd2.nheavyatoms(); // rsd2 is a ligand, check all atoms
				for (Size n=1; n<= max_res2_heavyatoms; n++) {
					Vector const & atom_i( rsd1.xyz( m ) );
					Vector const & atom_j( rsd2.xyz( n ) );
					// dist cutoff will be squared sum of LJ radii
					Real lj_sum = (rsd1.atom_type(m).lj_radius() + rsd2.atom_type(n).lj_radius());
					//dist_cutoff2 = lj_sum;
					Real dist_cutoff2 = (lj_sum * lj_sum) * bump_overlap_factor_;
					Real dist2 = ( atom_i - atom_j ).length_squared();
					if  (  dist2 <  dist_cutoff2 ) {
						// clash
						//TR << "clash: " << std::endl;
						//TR << "\t Loop residue " << i << ", atom " << m << ", atom_type " << rsd1.atom_type(m).name() << std::endl;
						//TR << "\t Other residue " << j << ", atom " << n << ", atom_type " << rsd2.atom_type(n).name() << std::endl;
						//TR << "dist2 was " << ( atom_i - atom_j ).length_squared() << std::endl;
						//TR << "and dist_cutoff2 was " << dist_cutoff2 << std::endl;
						//pose.dump_pdb("clash_test.pdb");
						//exit(0);
						//break;
						//std::cout << "Bump fail" << std::endl;
						return false;
					}
					//else {
					//	TR << "bump check passed res " << i << ", atom " << m << " against res " << j << ", atom " << n << std::endl;
					//}
				}
			}
		}
		// DJM: neighbor filter
		/*
		if (had_neighbor == false) {
			passed_but_no_neighbor = true;
			no_neighbor_res = i;
			TR << "Pose passed rama and bump check but resid " << no_neighbor_res << " has no neighbors." << std::endl;
			return false;
		}
		 */
	}

	//TR << "\tNO CLASHES!" << std::endl;
	//TR << "t_ang" << std::endl;
	//printVector(t_ang[i]);
	//TR << "b_ang" << std::endl;
	//printVector(b_ang[i]);
	//pose.dump_pdb("no_clash.pdb");
	//exit(0);
	return true; // passed bump check
}

// set default options
// AS -- after vicinity refactoring, only the perturber knows about vicinity sampling -- for now this function doesn't do anything any more... 
void
KinematicMover::set_defaults() {
	using namespace core;
	using namespace basic::options;

	return;
}
	
	/// @brief Taboo Sampling functions
	/// @author Amelie Stein
	/// @date Thu May 17 13:23:52 PDT 2012
	
	void 
	KinematicMover::update_sequence( utility::vector1< core::chemical::AA > const & sequence )
	{
		
		// note: it would be much better to implement this as a vector/map of strings, that then redirect to the corresponding taboo map -- otherwise we lose all information the moment we change loops, which will happen a lot when remodeling multiple loops -- TODO
		// in other words, this implementation is only intended for remodeling single loops
		
		if (sequence != sequence_) { // is there an actual change?
			
			if (sequence_.size() > 0) {
				// store last used taboo_map_ along with the corresponding sequence
				/*
				TR << "storing taboo map for "; // debug/testing info
				for (Size i = 1; i <= sequence_.size(); i++) 
					TR << sequence_[i] << " " ;
				TR << " (" << loop_begin_ << "-" << loop_end_ << ")" << std::endl; // actually I'm not sure when loop_begin_ and loop_end_ are updated... 
				*/
				taboo_master_map_[sequence_] = taboo_map_;
			}
			
			sequence_ = sequence;
			taboo_map_.clear();	
			if (taboo_master_map_.find(sequence) != taboo_master_map_.end()) {
				/*
				TR << "reloading taboo map for "; // debug/testing info
				for (Size i = 1; i <= sequence.size(); i++) 
					TR << sequence[i];
				TR << " (" << loop_begin_ << "-" << loop_end_ << ")" << std::endl;
				*/
				taboo_map_ = taboo_master_map_[sequence]; // reload from previous use, if possible
			}
			//random_torsion_strings_.clear();
			perturber_->clear_torsion_string_stack(); // only affects the TabooSampling perturber for now
		}
	}
	
	// the TabooSamplingKinematicPerturber must access this to derive the appropriate torsion bins
	utility::vector1< core::chemical::AA > 
	KinematicMover::get_loop_sequence() const 
	{
		return sequence_;		
	}
	
	void
	KinematicMover::insert_sampled_torsion_string_into_taboo_map( std::string const & ts )
	{
		if (sequence_.size() > 0 && taboo_map_.size() > pow(4, sequence_.size()) * taboo_map_max_capacity_) {
			TR << "Taboo map filled to " << taboo_map_max_capacity_*100 << "% capacity, clearing the map to avoid almost-endless loops. " << taboo_map_.size() << std::endl;
			taboo_map_.clear();
		}
		
		taboo_map_[ts] = true; 
	}
	
	bool
	KinematicMover::is_in_taboo_map( std::string const & ts ) const
	{
		return taboo_map_.find(ts) != taboo_map_.end();
	}	
	
	
	/// @brief returns the frequency of torsion_bin at position pos in the current taboo map, used for shifting probabilities of sampled torsion bins in later rounds
	core::Real
	KinematicMover::frequency_in_taboo_map( core::Size const & pos, char const & torsion_bin ) const {
		
		if (taboo_map_.size() == 0) 
			return 0;
		
		core::Size counter = 0;
		if (pos >= sequence_.size()) {
			TR << "error -- position " << pos << " exceeds the current taboo map string" << std::endl;
			return counter;
		} // assumption: all strings in the taboo map have the same size, and pos is within these bounds
		for (std::map< std::string, bool >::const_iterator mi = taboo_map_.begin(); mi != taboo_map_.end(); mi++) {
			counter += ( (mi->first)[pos] == torsion_bin ); // this might be slow... 
		}
		return core::Real(counter)/taboo_map_.size();
	}
	
	
	///@brief bin torsion angles as described in http://www.ncbi.nlm.nih.gov/pubmed/19646450
	///@detail generates a string with the torsion angle bins, using uppercase letters as in the publication above for omega ~ 180, and lowercase letters for omega ~ 0; to be used in loop sampling analysis -- same as in the LoopMover implementation... 
	///@author Amelie Stein
	///@date April 26, 2012
	std::string KinematicMover::torsion_features_string( core::pose::Pose const & pose ) const 
	{
		std::string torsion_bins, pos_bin;
		core::Real phi, psi, omega;
		
		// assumption: we have only one loop ... 
		for ( core::Size i = loop_begin_; i <= loop_end_; i++ ) {
			torsion_bins += core::conformation::get_torsion_bin(pose.phi(i), pose.psi(i), pose.omega(i));
		}
		return torsion_bins;
	} // torsion_features_string


} // namespace kinematic_closure
} // namespace loop_closure
} // namespace loops
} // namespace protocols
