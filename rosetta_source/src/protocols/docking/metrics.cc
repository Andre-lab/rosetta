// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file metrics
/// @brief protocols that are specific to docking low resolution
/// @detailed
/// @author Brian Weitzner

// Unit Headers

// Package Headers

// Project Headers

// Utility Headers

// Numeric Headers and ObjexxFCL Headers

// C++ headers

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.hh>

#include <protocols/docking/metrics.hh>

static basic::Tracer TR("protocols.docking.DockingProtocol.metrics");

using namespace core;

namespace protocols {
namespace docking {

core::Real
calc_interaction_energy( const core::pose::Pose & pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps){
	using namespace scoring;

	Real interaction_energy(0);

	core::scoring::ScoreFunctionOP docking_scorefxn;
	core::pose::Pose complex_pose = pose;

	docking_scorefxn = new core::scoring::ScoreFunction( *dock_scorefxn ) ;
    docking_scorefxn->set_weight( core::scoring::atom_pair_constraint, 0.0 );
	/*
	if ( pose.is_fullatom() ){
		docking_scorefxn = new core::scoring::ScoreFunction( *docking_score_high_ ) ;
	} else {
		docking_scorefxn = new core::scoring::ScoreFunction( *docking_score_low_ ) ;
	}
	*/

	// calculate energy of complexed pose
	Real const bound_energy = (*docking_scorefxn)( complex_pose );

	// calculate energy of separated pose over each movable jump
	// ddG is the "right" way to do this, to properly penalize strained rotamers
	// but aroop reports that I_sc yields better results for antibodies
	for( DockJumps::const_iterator it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {
		Size const rb_jump = *it;
		/*
		 Real const threshold = 100000; // dummy threshold
		 Size const repeats = 3;
		 protocols::protein_interface_design::DdgFilter ddg = protocols::protein_interface_design::DdgFilter( threshold, docking_scorefxn, rb_jump, repeats );
		 interaction_energy += ddg.compute( pose );
		 */
		core::pose::Pose unbound_pose = complex_pose;
		Real trans_magnitude = 1000;
		moves::RigidBodyTransMoverOP translate_away ( new moves::RigidBodyTransMover( unbound_pose, rb_jump ) );
		translate_away->step_size( trans_magnitude );
		translate_away->apply( unbound_pose );

		Real const unbound_energy = (*docking_scorefxn)( unbound_pose );

		interaction_energy += (bound_energy - unbound_energy);
	}

	return interaction_energy;
}

core::Real
calc_Lrmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, DockJumps const movable_jumps ){
	using namespace scoring;
	Real Lrmsd(0);

	for( DockJumps::const_iterator it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {
		Size const rb_jump = *it;
		ObjexxFCL::FArray1D_bool temp_part ( pose.total_residue(), false );
		ObjexxFCL::FArray1D_bool superpos_partner ( pose.total_residue(), false );
		/// this gets the wrong partner, therefore it is stored in a temporary
		/// array and then the opposite is put in the actualy array that is used
		/// for superpositioning.  there is probably a better way to do this
		/// need to check TODO
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if (temp_part(i)) superpos_partner(i)=false;
			else superpos_partner(i)=true;
		}

		//Lrmsd += rmsd_no_super_subset( native_pose, pose, superpos_partner, is_protein_CA );
		using namespace core::scoring;
		Lrmsd += core::scoring::rmsd_no_super_subset( native_pose, pose, superpos_partner, is_protein_backbone );
	}
	return Lrmsd;
}

core::Real
calc_Irmsd( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps ){
	using namespace scoring;
	Real Irmsd(0);

	for( DockJumps::const_iterator it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {
		core::Size const rb_jump = *it;
		if (!pose.is_fullatom() || !native_pose.is_fullatom()){
			TR << "Irmsd calc called with non-fullatom pose!!!"<<std::endl;
			return 0.0;
		}

		core::pose::Pose native_docking_pose = native_pose;
		using namespace kinematics;
		FoldTree ft( pose.fold_tree() );
		native_docking_pose.fold_tree(ft);

		//score to set up interface object
		core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction( *dock_scorefxn ) ;
		(*scorefxn)( native_docking_pose );

		protocols::scoring::Interface interface( rb_jump );
		interface.distance( 8.0 );
		interface.calculate( native_docking_pose );
		ObjexxFCL::FArray1D_bool is_interface ( pose.total_residue(), false );

		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if (interface.is_interface(i)) is_interface(i) = true;
		}

		//Irmsd += rmsd_with_super_subset(native_docking_pose, pose, is_interface, is_heavyatom);
		using namespace core::scoring;
		Irmsd += core::scoring::rmsd_with_super_subset(native_docking_pose, pose, is_interface, is_protein_backbone);
	}
	return Irmsd;
}

core::Real
calc_Fnat( const core::pose::Pose & pose, const core::pose::Pose & native_pose, const core::scoring::ScoreFunctionCOP dock_scorefxn, DockJumps const movable_jumps ){
	using namespace scoring;
	using namespace conformation;
	Real Fnat(0);

	for( DockJumps::const_iterator it = movable_jumps.begin(); it != movable_jumps.end(); ++it ) {
		core::Size const rb_jump = *it;

		if (!pose.is_fullatom() || !native_pose.is_fullatom()){
			TR << "Fnat calc called with non-fullatom pose!!!"<<std::endl;
			return 0.0;
		}

		core::pose::Pose native_docking_pose = native_pose;
		using namespace kinematics;
		FoldTree ft( pose.fold_tree() );
		native_docking_pose.fold_tree(ft);
		Real cutoff = 5.0;

		//score to set up interface object
		core::scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction( *dock_scorefxn ) ;
		(*scorefxn)( native_docking_pose );

		ObjexxFCL::FArray1D_bool temp_part ( pose.total_residue(), false );
		pose.fold_tree().partition_by_jump( rb_jump, temp_part );

		utility::vector1< Size > partner1;
		utility::vector1< Size > partner2;

		protocols::scoring::Interface interface( rb_jump );
		interface.distance( 8.0 );
		interface.calculate( native_docking_pose );

		//generate list of interface residues for partner 1 and partner 2
		for ( Size i=1; i <= pose.total_residue(); i++){
			if (interface.is_interface(i)){
				if (!temp_part(i)) partner1.push_back(i);
				if (temp_part(i)) partner2.push_back(i);
			}
		}

		//create contact pair list
		ObjexxFCL::FArray2D_bool contact_list(partner1.size(), partner2.size(), false);

		//identify native contacts across the interface
		//this will probably be changed to use PoseMetrics once I figure out how to use them - Sid
		for ( Size i=1; i<= partner1.size(); i++){
			ResidueOP rsd1 = new Residue(native_docking_pose.residue(partner1[i]));
			for ( Size j=1; j<=partner2.size();j++){
				ResidueOP rsd2 = new Residue(native_docking_pose.residue(partner2[j]));
				contact_list(i,j)=calc_res_contact(rsd1, rsd2, cutoff);
			}
		}

		Real native_ncontact = 0;
		Real decoy_ncontact = 0;

		//identify which native contacts are recovered in the decoy
		for ( Size i=1; i<=partner1.size(); i++){
			for ( Size j=1; j<=partner2.size(); j++){
				if (contact_list(i,j)){
					native_ncontact++;
					ResidueOP rsd1 = new Residue(pose.residue(partner1[i]));
					ResidueOP rsd2 = new Residue(pose.residue(partner2[j]));
					if (calc_res_contact(rsd1, rsd2, cutoff)) decoy_ncontact++;
				}
			}
		}

		Fnat += decoy_ncontact/native_ncontact;
	}
	return Fnat;
}

bool calc_res_contact(
	conformation::ResidueOP rsd1,
	conformation::ResidueOP rsd2,
	Real dist_cutoff
	)
{
	Real dist_cutoff2 = dist_cutoff*dist_cutoff;

	for (Size m=1; m<=rsd1->nheavyatoms(); m++){
		for (Size n=1; n<=rsd2->nheavyatoms(); n++){
			double dist2 = rsd1->xyz(m).distance_squared( rsd2->xyz(n) );  //Is there a reason this is a double?
			if (dist2 <= dist_cutoff2) return true;
		}
	}

	return false;
}

}//docking
}//protocols
