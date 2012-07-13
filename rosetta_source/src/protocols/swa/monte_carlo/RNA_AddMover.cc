// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeleteMover
/// @brief Adds an RNA residue from a chain terminus.
/// @detailed
/// @author Rhiju Das

#include <protocols/swa/monte_carlo/RNA_AddMover.hh>
#include <protocols/swa/monte_carlo/RNA_TorsionMover.hh>
#include <protocols/swa/monte_carlo/SubToFullInfo.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <basic/Tracer.hh>

#include <map>

using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// Removes one residue from a 5' or 3' chain terminus, and appropriately
// updates the pose sub_to_full_info object.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.monte_carlo.rna_add_mover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_AddMover::RNA_AddMover( chemical::ResidueTypeSetCAP rsd_set,
															scoring::ScoreFunctionOP scorefxn ):
		rsd_set_( rsd_set ),
		scorefxn_( scorefxn ),
		presample_added_residue_( true ), // need get/set functions
		start_added_residue_in_aform_( false ), // need get/set functions
		internal_cycles_( 50 ),
		rna_torsion_mover_( new RNA_TorsionMover ),
		sample_range_small_( 5.0 ),
		sample_range_large_( 40.0 ),
		kT_( 0.5 )
	{}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_AddMover::~RNA_AddMover()
  {}

  //////////////////////////////////////////////////////////////////////////
  void
  RNA_AddMover::apply( core::pose::Pose &  )
	{
		std::cout << "not defined yet" << std::endl;
	}

	//////////////////////////////////////////////////////////////////////
  void
  RNA_AddMover::apply( core::pose::Pose & pose, Size const res_to_build_off, MovingResidueCase const moving_residue_case )
	{

		using namespace protocols::moves;

		Size suite_num( 0 ), nucleoside_num( 0 ); // will record which new dofs added.

		//pose.dump_pdb( "before_add.pdb" );
		SubToFullInfo & sub_to_full_info = nonconst_sub_to_full_info_from_pose( pose );
		std::string const  full_sequence  = sub_to_full_info.full_sequence();
		std::map< Size, Size > sub_to_full = sub_to_full_info.sub_to_full();

		if ( moving_residue_case == CHAIN_TERMINUS_3PRIME ){

			runtime_assert( res_to_build_off < pose.total_residue() ); // wait is this necessary?
			runtime_assert( sub_to_full[ res_to_build_off ] < sub_to_full[ res_to_build_off+1 ] -1 );

			Size const res_to_add = res_to_build_off + 1;

			char newrestype = full_sequence[ (sub_to_full[ res_to_build_off ] + 1) - 1 ];
			//std::cout << "I want to add: " << newrestype << std::endl;

			chemical::AA my_aa = chemical::aa_from_oneletter_code( newrestype );
			chemical::ResidueTypeCAPs const & rsd_type_list( rsd_set_->aa_map( my_aa ) );
			// iterate over rsd_types, pick one.
			chemical::ResidueType const & rsd_type = *rsd_type_list[1];
			core::conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

			pose.append_polymer_residue_after_seqpos( *new_rsd, res_to_build_off, true /*build ideal geometry*/ );

			reorder_sub_to_full_info_after_append( pose, res_to_add );

			suite_num = res_to_add-1;
			nucleoside_num = res_to_add;

		} else {

			//Size const res_to_add = moving_res_list[1]; // for now, just build off 5' fragment.
			Size const res_to_add = res_to_build_off;

			runtime_assert( sub_to_full[ res_to_add ] > 1 );

			char newrestype = full_sequence[ (sub_to_full[ res_to_add ] - 1) - 1 ];
			//std::cout << "I want to add: " << newrestype << std::endl;

			chemical::AA my_aa = chemical::aa_from_oneletter_code( newrestype );
			chemical::ResidueTypeCAPs const & rsd_type_list( rsd_set_->aa_map( my_aa ) );
			// iterate over rsd_types, pick one.
			chemical::ResidueType const & rsd_type = *rsd_type_list[1];
			core::conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

			pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_add ); // got to be safe.

			pose.prepend_polymer_residue_before_seqpos( *new_rsd, res_to_add, true /*build ideal geometry*/ );
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_add );

			reorder_sub_to_full_info_after_prepend( pose, res_to_add );

			// initialize with a random torsion... ( how about an A-form + perturbation ... or go to a 'reasonable' rotamer)
			suite_num = res_to_add;
			nucleoside_num = res_to_add;

		}

		if ( start_added_residue_in_aform_ ){
			rna_torsion_mover_->apply_suite_torsion_Aform( pose, suite_num );
			rna_torsion_mover_->apply_nucleoside_torsion_Aform( pose, nucleoside_num );
		} else {
			rna_torsion_mover_->apply_random_nucleoside_torsion( pose, nucleoside_num );
			rna_torsion_mover_->apply_random_suite_torsion( pose, suite_num );
		}

		rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num, sample_range_large_);
		rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_large_);

		///////////////////////////////////
		// Presampling added residue
		///////////////////////////////////
		if ( presample_added_residue_ ){
			TR.Debug << "presampling added residue! " << nucleoside_num << " over " << internal_cycles_ << " cycles " << std::endl;
			MonteCarloOP monte_carlo_internal = new MonteCarlo( pose, *scorefxn_, kT_ );

			std::string move_type( "" );
			for ( Size count_internal = 1; count_internal <= internal_cycles_; count_internal++ ){

				rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num, sample_range_large_);
				rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_large_);
				monte_carlo_internal->boltzmann( pose, move_type );

				rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num, sample_range_small_);
				rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_small_);
				monte_carlo_internal->boltzmann( pose, move_type );
				//std::cout << "During presampling: " << (*scorefxn_)( pose );
			}
		}

	}

	std::string
	RNA_AddMover::get_name() const {
		return "RNA_AddMover";
	}

}
}
}
