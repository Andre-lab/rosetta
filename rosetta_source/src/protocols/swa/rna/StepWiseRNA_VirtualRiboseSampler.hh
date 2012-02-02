// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util.hh
/// @brief
/// @detailed
///
///  @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_SWA_VirtualRiboseSampler_HH
#define INCLUDED_protocols_swa_SWA_VirtualRiboseSampler_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/chemical/AA.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <set>
#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_Bin_Screener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh> //June 02, 2011
#include <core/pose/Pose.hh> //June 02, 2011
#include <core/scoring/rna/RNA_Util.hh> //June 02, 2011


namespace protocols {
namespace swa {
namespace rna {

class FB_Pose_Data{

		public:

		FB_Pose_Data():
		score( 0.0 ),
		tag("tag_blah"),
		base_tag("base_tag_blah"),
		Is_chain_close(false),
		base_rep_score(999999)

	{
	}

	~FB_Pose_Data(){};

	public:

	core::Real score;
	core::pose::PoseOP pose_OP;
	std::string tag;
	std::string base_tag;
	bool Is_chain_close; 
	core::Real base_rep_score; 
	core::kinematics::FoldTree starting_fold_tree;
	core::scoring::constraints::ConstraintSetOP starting_cst_set_OP;

};



class FloatingBaseChainClosureJobParameter{

	public:

	
		FloatingBaseChainClosureJobParameter(core::Size const input_moving_res, core::Size const input_reference_res):
			sample_sugar(true),
			moving_res ( input_moving_res ),
			reference_res( input_reference_res ),
			moving_res_pucker_state( ALL ),
			bulge_res_pucker_state( ALL ),
			moving_res_base_state( BOTH ),
			bulge_res_base_state( BOTH )
		{

			PDL.clear(); //pose_data_list
			Is_prepend= (moving_res>reference_res) ? false: true;
			bulge_res= (Is_prepend) ? reference_res-1: reference_res+1;
			bulge_suite = (Is_prepend) ? bulge_res : bulge_res-1;
			five_prime_chain_break= (Is_prepend) ? moving_res: moving_res-1;

			////////////////////////Consistency check!////////////////////////////////////////
			if(Is_prepend){
				if(moving_res+2!=reference_res){
					std::cout << "moving_res= " << moving_res << " reference_res= " << reference_res;
					utility_exit_with_message("prepend, but moving_res+2!=reference_res!");
				}			
			}else{
				if(moving_res-2!=reference_res){
					std::cout << "moving_res= " << moving_res << " reference_res= " << reference_res;
					utility_exit_with_message("append, but moving_res-2!=reference_res!");
				}			
			}

		}

		FloatingBaseChainClosureJobParameter():
			sample_sugar( false)
		{
			PDL.clear(); //pose_data_list
		};

		~FloatingBaseChainClosureJobParameter(){};


		void
		check_compatibility(core::Size const nres) const{
			using namespace ObjexxFCL;	

			if(moving_res<1 || moving_res> nres) utility_exit_with_message( "moving_res<1 || moving_res> nres. moving_res= " + string_of(moving_res) );
			if(bulge_res <1 || bulge_res> nres) utility_exit_with_message( "bulge_res <1 || bulge_res> nres. bulge_res= " + string_of(bulge_res) );
			if(reference_res<1 || reference_res> nres) utility_exit_with_message( "reference_res<1 || reference_res> nres. reference_res= " + string_of(reference_res) );
	
			//Should check here that moving_res contain virtual ribose and bulge_res is a virtual_rna_residue?
		}

		void
		set_base_and_pucker_state(core::pose::Pose const & pose, StepWiseRNA_JobParametersCOP const & JP){

			////////////////////////June 02, 2011 Add BaseState and PuckerState information///
			moving_res_pucker_state=ALL;
			if( Contain_seq_num(moving_res, JP->working_force_north_ribose_list() ) ) moving_res_pucker_state=NORTH;
			if( Contain_seq_num(moving_res, JP->working_force_south_ribose_list() ) ) moving_res_pucker_state=SOUTH;

			moving_res_base_state = (core::scoring::rna::is_purine( pose.residue( moving_res ) ) ) ? BOTH: ANTI;
			if( Contain_seq_num(moving_res, JP->working_force_syn_chi_res_list()  ) ) moving_res_base_state=SYN;

			bulge_res_pucker_state=ALL;
			if( Contain_seq_num(bulge_res, JP->working_force_north_ribose_list() ) ) bulge_res_pucker_state=NORTH;
			if( Contain_seq_num(bulge_res, JP->working_force_south_ribose_list() ) ) bulge_res_pucker_state=SOUTH;

			bulge_res_base_state = (core::scoring::rna::is_purine( pose.residue( bulge_res ) ) ) ? BOTH: ANTI;
			if( Contain_seq_num(bulge_res, JP->working_force_syn_chi_res_list()  ) ) bulge_res_base_state=SYN;

			////////////////////////Print data!////////////////////////////////////////
			std::cout << "FloatingBaseChainClosureJobParameter: " << std::endl;
			Output_boolean(" Is_prepend= " , Is_prepend);
			std::cout << " reference_res=" << reference_res;

			std::cout <<" moving_res=" << moving_res; 			
			print_base_state("|base_state=", moving_res_base_state );
			print_ribose_pucker_state("|pucker_state=", moving_res_pucker_state);

			std::cout << " bulge_res= " << bulge_res;
			print_base_state("|base_state=", bulge_res_base_state );
			print_ribose_pucker_state("|pucker_state=", bulge_res_pucker_state);

			std::cout << " bulge_suite= " << bulge_suite << " five_prime_chain_break= " << five_prime_chain_break; 
	
		}

	public:
		bool sample_sugar;
		core::Size moving_res;
		core::Size reference_res;
		bool Is_prepend;
		core::Size bulge_res;
		core::Size bulge_suite;
		core::Size five_prime_chain_break;
		PuckerState moving_res_pucker_state;
		PuckerState bulge_res_pucker_state;
		BaseState moving_res_base_state;
		BaseState bulge_res_base_state;
		utility::vector1< pose_data_struct2 > PDL; //pose_data_list of possible ribose conformations.

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
fast_full_atom_VDW_repulsion_screen(core::pose::Pose const & pose, core::Size const res_1, core::Size const res_2, bool const Is_prepend);


//Duplication of Chain_break_screening function from StepWiseRNA_ResidueSampler.cc. NEED TO MERGE THEM BACK TOGETHER AFTER TESTING! Apr 20,2010. Parin S.
bool
floating_base_full_atom_van_der_Waals_screening(core::pose::Pose & current_pose_screen,
																								core::Real const & base_rep_score, 
																								core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn,
																		 						SillyCountStruct & count_data, 
																								bool const verbose);

//Duplication of Chain_break_screening function from StepWiseRNA_ResidueSampler.cc. NEED TO MERGE THEM BACK TOGETHER AFTER TESTING! Apr 20,2010. Parin S.
bool
floating_base_chain_break_screening(core::pose::Pose & chain_break_screening_pose, 
																	core::scoring::ScoreFunctionOP const & chainbreak_scorefxn,
																	SillyCountStruct & count_data, 
																	core::Size const & five_prime_res,
																	std::string const & tag,
                                    bool const verbose);


utility::vector1< FB_Pose_Data >
floating_base_chain_closure_setup(utility::vector1< pose_data_struct2 > const & input_pose_data_list,
																	FloatingBaseChainClosureJobParameter const & FB_job_params,
																	core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn,
																	core::scoring::ScoreFunctionOP const & full_scorefxn,
																	core::pose::Pose & viewer_pose,
																	bool const do_minimize);

void
floating_base_chain_closure_sampling(utility::vector1< FB_Pose_Data > & pose_data_list, 
	                                   core::pose::Pose & viewer_pose,
																		 FloatingBaseChainClosureJobParameter const & FB_job_params,
																	   core::scoring::ScoreFunctionOP const & chainbreak_scorefxn,
													           core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn,
																		 StepWiseRNA_VDW_Bin_ScreenerOP const & VDW_bin_screener,
																	   bool const CCD_grid_index_screen);


utility::vector1< pose_data_struct2 >
floating_base_chain_closure_post_process(utility::vector1< FB_Pose_Data > & pose_data_list, 
	                                       core::pose::Pose & viewer_pose,
																				 core::scoring::ScoreFunctionOP const & sampling_scorefxn,
																		 		 FloatingBaseChainClosureJobParameter const & FB_job_params,
																				 bool const rm_chain_break_jump_point=true);

utility::vector1< pose_data_struct2 > 
sample_virtual_ribose_and_bulge_and_close_chain(core::pose::Pose & viewer_pose, 
																						 FloatingBaseChainClosureJobParameter const & FB_job_params, 
																						 std::string const name,
																						 core::scoring::ScoreFunctionOP const & scorefxn, 
																						 core::scoring::ScoreFunctionOP const & sampling_scorefxn, 
																						 core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn, 
																						 core::scoring::ScoreFunctionOP const & chainbreak_scorefxn,
																						 StepWiseRNA_JobParametersCOP & job_parameters,
																						 bool const virtual_ribose_is_from_prior_step=true);


void
minimize_all_sampled_floating_bases(core::pose::Pose & viewer_pose,
																		utility::vector1<FloatingBaseChainClosureJobParameter> const & FB_JP_list, 
																		utility::vector1< pose_data_struct2 > & pose_data_list,
																		core::scoring::ScoreFunctionOP const & sampling_scorefxn,
																		StepWiseRNA_JobParametersCOP const & job_parameters,
																		bool const virtual_ribose_is_from_prior_step=true);


bool
Is_ribose_virtual( core::pose::Pose const & pose, core::Size const previous_moving_res, core::Size const previous_bulge_res);

void
copy_bulge_res_and_ribose_torsion(FloatingBaseChainClosureJobParameter const & FB_job_params, core::pose::Pose & pose, core::pose::Pose const & template_pose);

void
enumerate_starting_pose_data_list(utility::vector1< pose_data_struct2 > & starting_pose_data_list,
																utility::vector1< FloatingBaseChainClosureJobParameter > const & FB_CC_JP_list, 
																core::pose::Pose const & pose);


utility::vector1< FloatingBaseChainClosureJobParameter > 
setup_FB_CC_JP_list(core::pose::Pose const & pose, utility::vector1< std::string > const & sample_virtual_ribose_string_list, StepWiseRNA_JobParametersCOP & job_parameters);

void
sample_user_specified_virtual_riboses(core::pose::Pose & pose, utility::vector1< std::string > const & sample_virtual_ribose_string_list, 
										 						    StepWiseRNA_JobParametersCOP & job_parameters, core::scoring::ScoreFunctionOP const & scorefxn,
																		std::string const silent_file_out, std::string const input_tag);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


}
}
}

#endif
