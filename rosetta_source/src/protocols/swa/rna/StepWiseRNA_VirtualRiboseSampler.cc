// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_VitualRiboseSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh> //Sept 26, 2011
#include <protocols/swa/rna/StepWiseRNA_VirtualRiboseSampler.hh>
//#include <protocols/swa/rna/StepWiseRNA_FloatingBase_Sampler_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Base_Sugar_Rotamer.hh>
#include <protocols/swa/rna/StepWiseRNA_Base_Sugar_Rotamer.fwd.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_Bin_Screener.hh>
#include <protocols/swa/rna/StepWiseRNA_VDW_Bin_Screener.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/id/TorsionID.hh>
//////////////////////////////////
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh> 
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/util.hh>
//#include <core/id/AtomID_Map.Pose.hh>
#include <set>
#include <numeric/conversions.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>

using namespace core;

static basic::Tracer TR( "protocols.swa.stepwise_rna_virtual_ribose_sampler" );

namespace protocols {
namespace swa {
namespace rna {



//////////////////////////////////////////////////////////////////////////////////////////////////////////////



	bool
	fast_full_atom_VDW_repulsion_screen(core::pose::Pose const & pose, core::Size const res_1, core::Size const res_2, bool const Is_prepend){

		conformation::Residue const & rsd_1=pose.residue(res_1);
		conformation::Residue const & rsd_2=pose.residue(res_2);


		for( Size n_1 = 1; n_1 <= rsd_1.natoms(); n_1++){

			//atom 1-4 are " P  ", " O1P", " O2P" and " O5*"
			Size const act_res_1= (Is_prepend && n_1 <= 4) ? res_1 + 1: res_1;

			if(pose.residue(act_res_1).atom_type(n_1).name()=="VIRT" ) continue;

			for( Size n_2 = 1; n_2 <= rsd_2.natoms(); n_2++){
	
				Size const act_res_2= (Is_prepend && n_2 <= 4) ? res_2 + 1: res_2;

				if(pose.residue(act_res_2).atom_type(n_2).name()=="VIRT" ) continue;

				Real const VDW_radius_1=pose.residue(act_res_1).atom_type(n_1).lj_radius();					
				Real const VDW_radius_2=pose.residue(act_res_2).atom_type(n_2).lj_radius();					

				Real const clash_dist_cutoff=0.8; //Fail van der Waals replusion screen if two atoms radius within 0.5 Angstrom of each other

				Real const clash_radius=VDW_radius_1+VDW_radius_2-clash_dist_cutoff; 

				if( (pose.residue(act_res_1).xyz( n_1) - pose.residue(act_res_2).xyz(n_2) ).length_squared() < clash_radius*clash_radius ){
					return false; //OK consider fail screening if find even one crash...
				}

			}
		}
		return true;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Duplication of Full_atom_van_der_Waals_screening from StepWiseRNA_ResidueSampler.cc. NEED TO MERGE THEM BACK TOGETHER AFTER TESTING! Apr 20,2010. Parin S.
	bool
	floating_base_full_atom_van_der_Waals_screening(core::pose::Pose & current_pose_screen,
																									core::Real const & base_rep_score, 
																									core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn,
																			 						SillyCountStruct & count_data, 
																									bool const verbose){

		using namespace core::scoring;

		(*atr_rep_screening_scorefxn)(current_pose_screen);

		EnergyMap const & energy_map = current_pose_screen.energies().total_energies();
		Real rep_score = atr_rep_screening_scorefxn->get_weight(fa_rep) * energy_map[scoring::fa_rep];

		Real delta_rep_score=rep_score-base_rep_score;

		Real actual_rep_cutoff=10; //defualt

		bool pass_rep_screen=false;

		if( delta_rep_score < actual_rep_cutoff ){
			pass_rep_screen=true;
			count_data.good_rep_rotamer_count++;

			if ( verbose ) {
				std::cout << "rep= " << delta_rep_score;
				std::cout << " rep_n= " << count_data.good_rep_rotamer_count;
				std::cout << " fast_rep_count= " << count_data.fast_full_atom_VDW_replusion_screen;
				std::cout << " bin_rep_count= " << count_data.good_bin_rep_count;
 				std::cout << " tot= " << count_data.tot_rotamer_count << std::endl;
			}
			return true;
		} else {
			return false;
		}
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Duplication of Chain_break_screening function from StepWiseRNA_ResidueSampler.cc. NEED TO MERGE THEM BACK TOGETHER AFTER TESTING! Apr 20,2010. Parin S.
	bool
	floating_base_chain_break_screening( core::pose::Pose & chain_break_screening_pose, 
																		 core::scoring::ScoreFunctionOP const & chainbreak_scorefxn,
																		 SillyCountStruct & count_data, 
																		 core::Size const & five_prime_res,
																		 std::string const & tag,
                                       bool const verbose){

		using namespace core::scoring;
		using namespace core::scoring::rna;

 		static protocols::rna::RNA_LoopCloser rna_loop_closer;

		set_CCD_torsions_to_zero(chain_break_screening_pose, five_prime_res); //Testing Oct 1, 2010
	
		rna_loop_closer.apply( chain_break_screening_pose, five_prime_res);

		(*chainbreak_scorefxn)(chain_break_screening_pose);

		scoring::EMapVector & energy_map= chain_break_screening_pose.energies().total_energies();
		Real angle_score = energy_map[scoring::angle_constraint];
		Real distance_score = energy_map[scoring::atom_pair_constraint];


		if(angle_score<5) count_data.good_angle_count++;
		if(distance_score<5) count_data.good_distance_count++;
		if((angle_score<5) && (distance_score<5)){
			count_data.chain_break_screening_count++;

			if(verbose){
				std::cout << " tag= " << tag;
				std::cout << " chain_closable_count= " << count_data.chain_closable_count;
				std::cout << " angle= " << angle_score << " dist= " << distance_score;
				std::cout << " angle_n= " << count_data.good_angle_count;
				std::cout << " dist_n= " << count_data.good_distance_count;
				std::cout << " chain_break_screening= " << count_data.chain_break_screening_count;
				std::cout << " in_range_CCD_torsion= " << count_data.in_range_CCD_torsion;
				std::cout << " tot= " << count_data.tot_rotamer_count << std::endl;
			}
			return true;
		} else {
			return false;
		}
	}

	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< FB_Pose_Data > 
	floating_base_chain_closure_setup(utility::vector1< pose_data_struct2 > const & input_pose_data_list,
																		FloatingBaseChainClosureJobParameter const & FB_job_params,
																		core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn,
																		core::scoring::ScoreFunctionOP const & full_scorefxn,
																		pose::Pose & viewer_pose,
																		bool const do_minimize){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;
		using namespace core::id;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::scoring::rna;
		using namespace protocols::rna;
		using namespace core::scoring;
		using namespace core::optimization;
		using namespace core::pose;
	
		pose::Pose const viewer_pose_copy=viewer_pose;

		clock_t const time_start( clock() ); 	

		core::scoring::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

		utility::vector1< FB_Pose_Data > pose_data_list; //This is the output pose_data_list;

		if(input_pose_data_list.size()==0) {return pose_data_list;} //return empty list
		/////////////////////////////////////////////////////////////////////////////////////////////////////


		StepWiseRNA_Base_Sugar_RotamerOP base_sugar_rotamer = new StepWiseRNA_Base_Sugar_Rotamer( FB_job_params.moving_res_base_state, FB_job_params.moving_res_pucker_state, rna_fitted_torsion_info);

		//July 28th, 2011 Could set extra_chi here, BUT necessary?

		int num_input_pose_data_pass_screen=0;
		for(Size n=1; n<=input_pose_data_list.size(); n++){

			bool input_pose_data_pass_screen=false;

			pose_data_struct2 input_pose_data=input_pose_data_list[n];
			base_sugar_rotamer->reset();	
			Size count=0;

			while(base_sugar_rotamer->get_next_rotamer()){
				count++;

//				bool verbose = (n==1 && count==1) ? true: false;
				bool verbose = true;


				FB_Pose_Data pose_data;
				pose_data.base_tag=input_pose_data.tag;

				pose_data.tag=input_pose_data.tag + "_" + base_sugar_rotamer->current_tag();


				pose_data.score=0.0;
			
				pose::Pose pose_with_ribose=(*input_pose_data.pose_OP); //Hard copy

				pose_data.starting_fold_tree=pose_with_ribose.fold_tree();

				pose_data.starting_cst_set_OP= pose_with_ribose.constraint_set()->clone();
				assert( pose_data.starting_cst_set_OP ); //if ( !cst_set ) cst_set = new ConstraintSet();



				(*atr_rep_screening_scorefxn)(pose_with_ribose);
				Real const without_ribose_rep_score = atr_rep_screening_scorefxn->get_weight(fa_rep) * pose_with_ribose.energies().total_energies()[ scoring::fa_rep ];

				pose::remove_variant_type_from_pose_residue( pose_with_ribose, "VIRTUAL_RIBOSE", FB_job_params.moving_res );

				Add_harmonic_chainbreak_constraint(pose_with_ribose, FB_job_params.five_prime_chain_break );

				setup_chain_break_jump_point( pose_with_ribose, FB_job_params.moving_res, FB_job_params.reference_res, FB_job_params.five_prime_chain_break, verbose);

				///////////////////////////////////////////////////////////////////////////////////////////////////

									
				if(verbose){	
					std::cout << " 	delta1= " <<  F(8, 3, base_sugar_rotamer->delta()) << " 	chi_1= " <<  F(8, 3, base_sugar_rotamer->chi());
					std::cout << " 	nu2_1= " <<  F(8, 3, base_sugar_rotamer->nu2())    << " 	nu1_1= " <<  F(8, 3, base_sugar_rotamer->nu1());
					std::cout << std::endl;
				}

				pose_with_ribose.set_torsion( TorsionID( FB_job_params.moving_res , id::BB, 4 ) , base_sugar_rotamer->delta());
				pose_with_ribose.set_torsion( TorsionID( FB_job_params.moving_res , id::CHI, 1 ) , base_sugar_rotamer->chi());
				pose_with_ribose.set_torsion( TorsionID( FB_job_params.moving_res , id::CHI, 2 ) , base_sugar_rotamer->nu2());
				pose_with_ribose.set_torsion( TorsionID( FB_job_params.moving_res , id::CHI, 3 ) , base_sugar_rotamer->nu1());

				//testing
				//pose_with_ribose.dump_pdb( "pose_with_ribose_before_minimize_"+ string_of(count) + ".pdb"  );

				//////////////////////////////////////minimize to remove ribose clashes (May 23, 2010)/////////////////////////////////////
				if(do_minimize==true){

					core::scoring::ScoreFunctionOP ribose_scorefxn = full_scorefxn->clone();
					ribose_scorefxn->set_weight( linear_chainbreak, 0.0); 
					ribose_scorefxn->set_weight( angle_constraint, 0.0 );
					ribose_scorefxn->set_weight( atom_pair_constraint, 0.0 );
					ribose_scorefxn->set_weight( coordinate_constraint, 0.1 );
					core::scoring::ScoreFunctionOP ribose_scorefxn_without_ch_bond = ribose_scorefxn->clone();
					ribose_scorefxn_without_ch_bond->set_weight( ch_bond, 0.0 );
					//This makes sure that there are no chain_break score involved.

					//////////////////////Sept 20, 2011 To solve the problem with the floating base res expoding when minimizing/////////////////////////////////////////////////////////////////
					//////////////////////Problem case occur when flaoting base res is the 1st working res///////////////////////////////////////////////////////////////////////////////////////
					//////////////////////Note that this error doesn't seem to occur if virtual_ribose is sampled in same step as SAMPLER/ (no Hbond_tripped!)///////////////////////////////////
					//////////////////////The non-rescale scorefxn does however causes the floating base from moving far away from the starting point even in the same step as SAMPLER case//////
					//////////////////////Also generally the minimizer same and seperate step virtual sampler doesn't give the same results!/////////////////////////////////////////////////////
					std::cout << "--------------START Creating rescaled one_tenth_ribose_score_fxn_without_ch_bond--------------" << std::endl; 
					core::scoring::ScoreFunctionOP rescaled_ribose_score_fxn_without_ch_bond=rescale_scorefxn(ribose_scorefxn_without_ch_bond, 0.1);
					std::cout << "--------------FINISH Creating rescaled one_tenth_ribose_score_fxn_without_ch_bond--------------" << std::endl; 
					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					AtomTreeMinimizer minimizer;
					bool const use_nblist( true );

					//float const dummy_tol( 0.00000025);
					//MinimizerOptions options_standard( "dfpmin", dummy_tol, use_nblist, false, false );

					//float const tolerance= 0.25; //Tested for dfpmin_atol, this lead to different result which is NOT worst than the energy score with 0.00000025 dfpmin_atol in every case!
					//float const tolerance= 0.00000025; //If use this tolerance value, different result for dfp_min and dfpmin_atol.
					float const tolerance= 0.000000000000025; //Sept 21, 2011, same result for  dfp_min and dfpmin_atol| converge to identical energy score with 0.00000025 of dfpmin_atol!

					MinimizerOptions options_standard( "dfpmin_atol", tolerance, use_nblist, false, false );      //Switch to absolute tolerance on Sept 21, 2011 
					MinimizerOptions options_armijo( "dfpmin_armijo_atol", tolerance, use_nblist, false, false ); //Add this on Sept 21, 2011 

					std::cout << "options_standard: min_type= " << options_standard.min_type() << " minimize_tolerance= "  << options_standard.minimize_tolerance() << std::endl;
					std::cout << "options_armijo  : min_type= " << options_armijo.min_type()   << " minimize_tolerance= "  << options_armijo.minimize_tolerance()   << std::endl;

					options_standard.nblist_auto_update( true );
					options_armijo.nblist_auto_update( true );

					core::kinematics::MoveMap mm;

					mm.set_bb(false);	
					mm.set_chi(false);
			

					if(verbose) std::cout << "pose.fold_tree().num_jump()= " << pose_with_ribose.fold_tree().num_jump() << std::endl;

					Size bulge_jump_1, bulge_jump_2, bulge_cutpoint;					
					if(FB_job_params.moving_res>FB_job_params.reference_res){//consistency check (from setup_chain_break_jump_point function)
						bulge_jump_1=FB_job_params.reference_res;
						bulge_jump_2=FB_job_params.moving_res;
						bulge_cutpoint=FB_job_params.moving_res-1;
					}else{
						bulge_jump_1=FB_job_params.moving_res;
						bulge_jump_2=FB_job_params.reference_res;
						bulge_cutpoint=FB_job_params.moving_res;
					}	

					if(verbose) std::cout << "bulge_cutpoint= " << bulge_cutpoint << " bulge_jump_1= " << bulge_jump_1 << " bulge_jump_2= " << bulge_jump_2 << std::endl;


					bool found_desired_jump_ID=false;
					for (Size jump_ID = 1; jump_ID <= pose_with_ribose.fold_tree().num_jump(); jump_ID++ ){
						Size const jump_pos1( pose_with_ribose.fold_tree().upstream_jump_residue( jump_ID ) );
						Size const jump_pos2( pose_with_ribose.fold_tree().downstream_jump_residue( jump_ID ) );
						Size const cutpoint=pose_with_ribose.fold_tree().cutpoint(jump_ID);

						if(verbose) std::cout << "jump at jump_ID= " << jump_ID << " cutpoint= " << cutpoint << " jump_pos1= " << jump_pos1 << " jump_pos2= " << jump_pos2 << std::endl;

						if( (jump_pos1==bulge_jump_1 &&	jump_pos2==bulge_jump_2 ) || (jump_pos1==bulge_jump_2 &&	jump_pos1==bulge_jump_2 ) ) {
							found_desired_jump_ID=true;
							if(verbose) std::cout << "add movemap jump at jump_ID= " << jump_ID << " cutpoint= " << cutpoint << " jump_pos1= " << jump_pos1 << " jump_pos2= " << jump_pos2 << std::endl;
							mm.set_jump( jump_ID, true );	
						}
					}

					if(found_desired_jump_ID==false) utility_exit_with_message( "cannot find desired jump_ID" );
		
					//Rigid body movement...no free torsion!!	
					//	core::scoring::constraints::add_coordinate_constraints( pose_with_ribose );//crap I left this on for May_28_1CSL_SYN_CHI_floating_base_use_May_18_data run!

					viewer_pose= pose_with_ribose;
					std::cout << "removing ribose clashes pose # " << n << " sugar_rotamer # " << count << " " << std::endl;

					/////////////Switch to armijo on Sept 21, 2011///////////////////////////////////////////////////////////////////////////////////////
					/////////////My understanding is that dfpmin_armijo is a "inexact" line search whereas the standard dfpmin is a exact line search///////
					/////////////It seem to indicate the dfpmin should be slower (require more function evaluation) but at the same time more accurate//////
					/////////////See http://www.rosettacommons.org/manuals/archive/rosetta3.3_user_guide/minimization_overview.html for details/////////////
					/////////////However standard dfpmin seem to lead cases where the floating base just "explode" and more far away from starting point////
					/////////////This sometimes lead to the Hbond tripped error/////////////////////////////////////////////////////////////////////////////
					/////////////Side note: switching to dfpmin_atol (atol-> absolute tolerance didn't help!)///////////////////////////////////////////////
					/////////////So switching to dfpmin_armijo which doesn't seem to exhibit this behavior//////////////////////////////////////////////////
					/////////////Note that there is a currently a bug in in dfpmin_armijo:
					/////////////core.optimization.LineMinimizer: Inaccurate G! step= 9.53674e-07 Deriv= -0.0226443 Finite Diff= 0.00628252/////////////////
					/////////////Rhiju mention that this bug is fixed in the latest Rosetta version in trunk////////////////////////////////////////////////
					/////////////So will keep using standard dfp_min except at the first minimiziation step/////////////////////////////////////////////////
					/////////////Also tried two round minimizations with the first using options_armijo. This fix the "explode" bug but led to worst score!/

					minimize_with_constraints(viewer_pose, mm, rescaled_ribose_score_fxn_without_ch_bond, options_armijo ); //Add this round on Sept 20, 2011, Switch to armijo on Sept 21, 2011
					minimizer.run( viewer_pose, mm, *(rescaled_ribose_score_fxn_without_ch_bond), options_armijo );         //Add this round on Sept 20, 2011, Switch to armijo on Sept 21, 2011
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					minimize_with_constraints(viewer_pose, mm, ribose_scorefxn_without_ch_bond, options_standard ); 
					minimizer.run( viewer_pose, mm, *(ribose_scorefxn_without_ch_bond), options_standard );         

					minimize_with_constraints(viewer_pose, mm, ribose_scorefxn, options_standard );
					minimizer.run( viewer_pose, mm, *(ribose_scorefxn), options_standard );

					pose_with_ribose=viewer_pose;
				}

				//testing
				//pose_with_ribose.dump_pdb( "pose_with_ribose_after_minimize_"+ string_of(count) + ".pdb"  );


				////////////////////////////////////////////Screens/////////////////////////////////////////////
				//OK check that with this sugar, the chain can be theoretically closed..
				std::string const moving_atom_name= (FB_job_params.Is_prepend) ? "O3*" : " C5*"; 
				std::string const reference_atom_name= (FB_job_params.Is_prepend) ? " C5*" : "O3*";
				Distance O3i_C5iplus2_distance=(pose_with_ribose.residue(FB_job_params.moving_res).xyz(moving_atom_name) - pose_with_ribose.residue(FB_job_params.reference_res).xyz(reference_atom_name) ).length();

				(*atr_rep_screening_scorefxn)(pose_with_ribose);
				pose_data.base_rep_score = atr_rep_screening_scorefxn->get_weight(fa_rep) * pose_with_ribose.energies().total_energies()[ scoring::fa_rep ];

				std::cout << "tag= " << pose_data.tag << " with_ribose_rep= " << pose_data.base_rep_score << " without_ribose_rep_score= " << without_ribose_rep_score << " O3i_C5iplus2_dist= " << O3i_C5iplus2_distance;
	

				if(O3i_C5iplus2_distance>O3I_C5IPLUS2_MAX_DIST){
					std::cout << " O3i_C5iplus2_dist>O3I_C5IPLUS2_MAX_DIST(" << O3I_C5IPLUS2_MAX_DIST << ") " << std::endl;				
					continue;
				}	

				if((pose_data.base_rep_score-without_ribose_rep_score)>10){
					std::cout << " RIBOSE_rep_score>10! "<< std::endl;
					continue;
				}
				input_pose_data_pass_screen=true;				
				/////////////////////////////////////////////////////////////////////////////////////////////////////

		    // Add to the CCD loop closure phosphate...this is always the 3' chain_break res...

				remove_virtual_rna_residue_variant_type(pose_with_ribose, FB_job_params.bulge_res); 
				pose::add_variant_type_to_pose_residue( pose_with_ribose, "VIRTUAL_PHOSPHATE", FB_job_params.five_prime_chain_break+1 );

				pose_data.pose_OP=new pose::Pose;
				(*pose_data.pose_OP)=pose_with_ribose;
				pose_data.Is_chain_close=false;
				pose_data_list.push_back(pose_data);

			}
			if(input_pose_data_pass_screen) num_input_pose_data_pass_screen++;
		}

		Output_title_text("");

		std::cout << "input_pose_data_list.size()= " << input_pose_data_list.size() << " num_input_pose_data_pass_screen= " << num_input_pose_data_pass_screen; 
		std::cout << " pose_data_list.size()= " << pose_data_list.size() << std::endl;

		std::cout << "Total time in Floating_base_chain_closure SETUP: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

		viewer_pose=viewer_pose_copy;

		return pose_data_list;
	}


	/////////////////////////////////////////////////////////////////////////////////////

	void
	floating_base_chain_closure_sampling(utility::vector1< FB_Pose_Data > & pose_data_list, 
		                                   core::pose::Pose & viewer_pose,
																		 FloatingBaseChainClosureJobParameter const & FB_job_params,
																		 core::scoring::ScoreFunctionOP const & chainbreak_scorefxn,
														          core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn,
																		 StepWiseRNA_VDW_Bin_ScreenerOP const & VDW_bin_screener,
																		 bool const CCD_grid_index_screen){

		Output_title_text("Floating_base_chain_closure SAMPLING");

		if(pose_data_list.size()==0) return; //early return
			
		//std::map<PuckerState, std::map<Base_bin , int , compare_base_bin> > chain_closure_grid_index; //"CCD_grid_index_screen no longer supported!"

		if(CCD_grid_index_screen){
			utility_exit_with_message("CCD_grid_index_screen no longer supported!");
			//chain_closure_grid_index[NORTH]=create_chain_closable_grid_index(FB_job_params, viewer_pose, NORTH);
			//chain_closure_grid_index[SOUTH]=create_chain_closable_grid_index(FB_job_params, viewer_pose, SOUTH);
		}

		clock_t const time_start_sampling( clock() ); 	

		using namespace core::id;
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::scoring::rna;
		using namespace protocols::rna;
		using namespace core::scoring;
		using namespace ObjexxFCL;

		core::scoring::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

		viewer_pose=(*pose_data_list[1].pose_OP);
		pose::Pose screening_pose=viewer_pose;

		///Ok we need to create different rotamer_generator for the 1. backbone/exclude delta, 2. delta. and chi. (depend on delta)
		///This is necessary for optimization since CCD solution doesn't depend on chi torsion value

		std::cout << "setup backbone_rotamer_generator" << std::endl;
		utility::vector1< core::Size > bulge_suite_list;
		bulge_suite_list.clear();
		bulge_suite_list.push_back(FB_job_params.bulge_suite);

		bool sample_sugar_and_base1( false ), sample_sugar_and_base2( false );

		if ( FB_job_params.Is_prepend ) {
			sample_sugar_and_base1 = true;
		} else {
			sample_sugar_and_base2 = true;
		}


		//Will generate both NORTH AND SOUTH PUCKER HERE EVEN IF FB_job_params.bulge_res_pucker_state!=ALL..BUT OK SINCE WILL PERFORM MATCH WITH Base_Sugar_RotamerOP in the actual loop.
		StepWiseRNA_RotamerGenerator_WrapperOP backbone_rotamer_generator = 
						new StepWiseRNA_RotamerGenerator_Wrapper( viewer_pose, bulge_suite_list, sample_sugar_and_base1, sample_sugar_and_base2);

		backbone_rotamer_generator->set_fast( false );
		backbone_rotamer_generator->set_sample_chi_torsion(false);
		backbone_rotamer_generator->set_include_syn_chi(true);
		backbone_rotamer_generator->set_bin_size(20);

		bool more_rotamers=true;
		if(more_rotamers){
			backbone_rotamer_generator->set_extra_epsilon(true);
			//backbone_rotamer_generator->set_extra_beta(true);
		}

		backbone_rotamer_generator->initialize_rotamer_generator_list();

		std::cout << "setup_delta_rotamer_generator" << std::endl;

		StepWiseRNA_Base_Sugar_RotamerOP bulge_base_sugar_rotamer = new StepWiseRNA_Base_Sugar_Rotamer( FB_job_params.bulge_res_base_state, FB_job_params.bulge_res_pucker_state, rna_fitted_torsion_info);

		//July 28th, 2011 Could set extra_chi here, BUT necessary?

		SillyCountStruct count_data;
		RNA_LoopCloser rna_loop_closer;

		Size num_closed_chain_pose=0;
		
		while( backbone_rotamer_generator->has_another_rotamer() ){ //BUT EPSILON DEPENDS on DELTA!
			bulge_base_sugar_rotamer->reset();	

			if(num_closed_chain_pose==pose_data_list.size()){
				break; //early break, all pose had been closed...
			}

			utility::vector1< Torsion_Info > const & BB_rotamer = backbone_rotamer_generator->get_next_rotamer();

			//optimization since CCD solution doesn't depend on chi torsion value, potential 6 fold speed up...
			utility::vector1< bool > CCD_fail_for_BB_rotamer(pose_data_list.size(), false);


			while( bulge_base_sugar_rotamer->get_next_rotamer() ){

				Real BB_delta_value=0.0;

				Size num_BB_delta_ID=0;
				Size num_chi_torsion_ID=0;

				core::id::TorsionID const chi_torsion_ID=TorsionID( FB_job_params.bulge_res , id::CHI, 1 );
				core::id::TorsionID const delta_torsion_ID=TorsionID( FB_job_params.bulge_res , id::BB, 4 );

				for( Size n = 1; n <= BB_rotamer.size(); n++ ){ 

					if(BB_rotamer[n].id==delta_torsion_ID){
						BB_delta_value=BB_rotamer[n].value;
						num_BB_delta_ID++;
					}

					////////Update this on May 1st, 2011, after updating RotamerGeneartor class to not include chi torsion ID when CHI torsion is not sampled.
					if(BB_rotamer[n].id==chi_torsion_ID){
						num_chi_torsion_ID++;
					}

				}
				if(num_BB_delta_ID!=1) utility_exit_with_message( "Error: num_BB_delta_ID in BB_rotamer  !=1" );
				if(num_chi_torsion_ID!=0) utility_exit_with_message( "Error: num_chi_torsion_ID in BB_rotamer !=0" );

				if(BB_rotamer.size()!=8) utility_exit_with_message( "BB_rotamer.size()=" + ObjexxFCL::string_of(BB_rotamer.size()) + "!=8!" ); //consistency_check

				if( (BB_delta_value < (bulge_base_sugar_rotamer->delta()-0.1) ) || ( BB_delta_value > (bulge_base_sugar_rotamer->delta()+0.1) ) ) continue;
									
				utility::vector1< Torsion_Info > current_rotamer = BB_rotamer;

				Torsion_Info chi_torsion_info;
				chi_torsion_info.id=chi_torsion_ID;
				chi_torsion_info.value=bulge_base_sugar_rotamer->chi();
				current_rotamer.push_back(chi_torsion_info);

				if(current_rotamer.size()!=9) utility_exit_with_message( "current_rotamer.size()=" + ObjexxFCL::string_of(current_rotamer.size()) + " !=9!" ); //consistency check

				//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				count_data.tot_rotamer_count++;

				apply_rotamer( screening_pose, current_rotamer); 
			
				apply_rotamer( viewer_pose, current_rotamer);
				viewer_pose.residue(1).xyz("C5*"); //This is just to force the viewer_pose to update...	

								
				//VDW_screen_bin
				if( VDW_bin_screener->VDW_rep_screen(screening_pose, FB_job_params.bulge_res)==false) continue;

				count_data.good_bin_rep_count++;

				for(Size n=1; n<=pose_data_list.size(); n++){

					//if(CCD_fail_for_BB_rotamer[n]==true) continue; //already checked that this BB rotamer doesn't have a CCD solution.  Feb 08, 2012: Commented out until further testing		

					FB_Pose_Data & pose_data=pose_data_list[n];
					if(pose_data.Is_chain_close==true) continue; //chain is already close...

					pose::Pose & current_pose=(*pose_data.pose_OP);
	
					if(FB_job_params.Is_prepend){
						if(Check_chain_closable_floating_base(current_pose, screening_pose, FB_job_params.five_prime_chain_break, 0 )==false) continue;
					}else{
						if(Check_chain_closable_floating_base(screening_pose, current_pose, FB_job_params.five_prime_chain_break, 0 )==false) continue;
					}
					count_data.chain_closable_count++;

					apply_rotamer( current_pose, current_rotamer);

					//NEED THIS:Quick check that there is no VDW clash between bulge_res and moving_res
					if( fast_full_atom_VDW_repulsion_screen(current_pose, FB_job_params.bulge_res, FB_job_params.moving_res, FB_job_params.Is_prepend) == false) continue;
					count_data.fast_full_atom_VDW_replusion_screen++;

					if( floating_base_full_atom_van_der_Waals_screening(current_pose, pose_data.base_rep_score, atr_rep_screening_scorefxn, count_data, true)==false ) continue;

					if( floating_base_chain_break_screening(current_pose, chainbreak_scorefxn, count_data, FB_job_params.five_prime_chain_break, pose_data.tag, true)==false ){

						 //CCD_fail_for_BB_rotamer[n]=true; //Feb 08, 2012. FIX ERROR, used to be == instead of =. Feb 08, 2012: Commented out until further testing

						 continue;
					}

					if( CCD_grid_index_screen){
						utility_exit_with_message("CCD_grid_index_screen no longer supported!");
						//if(floating_base_chain_break_grid_index_screening(current_pose, count_data, FB_job_params, chain_closure_grid_index)==false) continue;
					}

					///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					
//					//Ok now phosphate torsion are defined remove virtual phosphate and check for clashes. Still need to implement this....May 14, 2010
//					pose::remove_variant_type_from_pose_residue( viewer_pose , "VIRTUAL_PHOSPHATE",  FB_job_params.five_prime_chain_break+1 );
//
//					if( floating_base_full_atom_van_der_Waals_screening(current_pose, pose_data.base_rep_score, atr_rep_screening_scorefxn, count_data, true)==false){
//						pose::add_variant_type_to_pose_residue( pose_with_ribose, "VIRTUAL_PHOSPHATE", FB_job_params.five_prime_chain_break+1 ); //Add back the virtual phosphate..
//						continue;
//					}

					/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					pose_data.Is_chain_close=true;
					num_closed_chain_pose++;

					

					viewer_pose=current_pose; //Crap this is hard copy...computationally expensive!!
					viewer_pose.residue(1).xyz("C5*"); //This is just to force the viewer_pose to update...	
				}
			}
		}

		std::cout << " bin_rep_count= " << count_data.good_bin_rep_count;
		std::cout << " fast_rep_count= " << count_data.fast_full_atom_VDW_replusion_screen;
		std::cout << " chain_closable= " << count_data.chain_closable_count;
		std::cout << " angle_n= " << count_data.good_angle_count << " dist_n= " << count_data.good_distance_count;
		std::cout << " rep= " << count_data.good_rep_rotamer_count;
		std::cout << " rmsd= " << count_data.rmsd_count << " tot= " << count_data.tot_rotamer_count << std::endl;
		std::cout << " " << num_closed_chain_pose << " out of " << pose_data_list.size() << " pose were closable" << std::endl;
		std::cout << "Total time in Floating_base_chain_closure SAMPLING: " << static_cast<Real>( clock() - time_start_sampling ) / CLOCKS_PER_SEC << std::endl;
	}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< pose_data_struct2 >
	floating_base_chain_closure_post_process(utility::vector1< FB_Pose_Data > & pose_data_list, 
		                                       core::pose::Pose & viewer_pose,
																				 core::scoring::ScoreFunctionOP const & sampling_scorefxn,
																			 	 FloatingBaseChainClosureJobParameter const & FB_job_params,
																				 bool const	rm_chain_break_jump_point){

			using namespace core::optimization;
			using namespace core::scoring;
			using namespace core::pose;
			using namespace core::io::silent;
			using namespace protocols::rna;
			using namespace core::id;
			using namespace ObjexxFCL;

			Output_title_text("Floating_base_chain_closure POST_PROCESSING");
			clock_t const time_start_post_processing( clock() ); 

			utility::vector1< pose_data_struct2 > output_pose_data_list;

			if( pose_data_list.size()==0) {return output_pose_data_list;} //return empty list

			//Quick minimize to remove error in CCD?///////////////////////////
			AtomTreeMinimizer minimizer;
			float const dummy_tol( 0.00000025);
			bool const use_nblist( true );
			MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
			options.nblist_auto_update( true );

			core::kinematics::MoveMap mm;

			mm.set_bb(false);	
			mm.set_chi(false);

			//The whole richardson's suite ... avoid delta, nu2, nu1, WARNING VIRTUAL_PHOSPHATE is ON at this point...
			mm.set( TorsionID( FB_job_params.bulge_res-1 , id::BB,  5 ), true ); //epsilon
			mm.set( TorsionID( FB_job_params.bulge_res-1 , id::BB,  6 ), true ); //zeta

			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  1 ), true ); //alpha
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  2 ), true ); //beta
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  3 ), true ); //gamma
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  5 ), true ); //epsilon
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  6 ), true ); //zeta
			mm.set( TorsionID( FB_job_params.bulge_res , id::CHI, 1 ), true ); //chi (torsion between base and ribose sugar)

			mm.set( TorsionID( FB_job_params.bulge_res+1 , id::BB,  1 ), true ); //alpha
			mm.set( TorsionID( FB_job_params.bulge_res+1 , id::BB,  2 ), true ); //beta 
			mm.set( TorsionID( FB_job_params.bulge_res+1 , id::BB,  3 ), true ); //gamma

			core::scoring::ScoreFunctionOP  bulge_chain_closure_scorefxn = new ScoreFunction;
			bulge_chain_closure_scorefxn->set_weight( fa_rep  , 0.12 );
			bulge_chain_closure_scorefxn->set_weight( angle_constraint, 1.0 );
			bulge_chain_closure_scorefxn->set_weight( atom_pair_constraint, 1.0 );
			bulge_chain_closure_scorefxn->set_weight( linear_chainbreak, 5.0);

			////Hacky fix: this shouldn't be necessary but is required to fix bug (rna_sugar_close score adn geom_sol doesn't check for virtual_atoms)////////
			////PREVENT RANDOM ENERGETIC PENALTY!////
			//bulge_chain_closure_scorefxn->set_weight( geom_sol, 62); //100X normal weight //MOD OUT July 23, 2011. GEOM_SOL and CI_GEOM_SOL NOW DOES CHECK FOR VIRTUAL_ATOM!
			//bulge_chain_closure_scorefxn->set_weight( CI_geom_sol, 62); //100X normal weight //ADD AND MOD OUT July 23, 2011. GEOM_SOL and CI_GEOM_SOL NOW DOES CHECK FOR VIRTUAL_ATOM!

			bulge_chain_closure_scorefxn->set_weight( rna_sugar_close, 70); //100X normal weight

			mm.set( TorsionID( FB_job_params.bulge_res , id::CHI, 2 ), true ); //nu2
			mm.set( TorsionID( FB_job_params.bulge_res , id::CHI, 3 ), true ); //nu1

			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			for(Size n=1; n<=pose_data_list.size(); n++){
				if(pose_data_list[n].Is_chain_close==false) continue;

				viewer_pose= (*pose_data_list[n].pose_OP);

				std::cout << "POST_PROCESSING pose # " << n << " out of " << pose_data_list.size() << " " << std::endl;;
				minimizer.run( viewer_pose, mm, *(bulge_chain_closure_scorefxn), options );

				viewer_pose.constraint_set( pose_data_list[n].starting_cst_set_OP);

				if(rm_chain_break_jump_point) remove_chain_break_jump_point(viewer_pose, FB_job_params.five_prime_chain_break, pose_data_list[n].starting_fold_tree);
				
				pose::remove_variant_type_from_pose_residue( viewer_pose , "VIRTUAL_PHOSPHATE",  FB_job_params.five_prime_chain_break+1 );
				apply_virtual_rna_residue_variant_type( viewer_pose, FB_job_params.bulge_res);

				pose_data_list[n].score=(*sampling_scorefxn)(viewer_pose); //for output purposes...

				(*pose_data_list[n].pose_OP)=viewer_pose;

				//////////////////////////////////////////				
				//Warning data_data_list and output_pose_data_list doesn't have the same underlying structure!! Apr 20, 2010
				pose_data_struct2 pose_data;
				pose_data.pose_OP=pose_data_list[n].pose_OP;
				pose_data.score=pose_data_list[n].score;
				pose_data.tag=pose_data_list[n].tag;

				output_pose_data_list.push_back(pose_data);
			}
	

			return output_pose_data_list;
			std::cout << "Total time in Floating_base_chain_closure_post_processing: " << static_cast<Real>( clock() - time_start_post_processing ) / CLOCKS_PER_SEC << std::endl;

	}		

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	minimize_all_sampled_floating_bases(core::pose::Pose & viewer_pose,
                                      utility::vector1<FloatingBaseChainClosureJobParameter> const & FB_JP_list, 
																		  utility::vector1< pose_data_struct2 > & pose_data_list,
																			core::scoring::ScoreFunctionOP const & sampling_scorefxn,
																			StepWiseRNA_JobParametersCOP const & job_parameters,
																			bool const virtual_ribose_is_from_prior_step){

		using namespace core::optimization;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::id;
		using namespace ObjexxFCL;

		Output_title_text("Enter minimize_all_sampled_floating_bases");

		pose::Pose const viewer_pose_copy=viewer_pose;


		if(FB_JP_list.size()==0) return;
		if(pose_data_list.size()==0) return;

		AtomTreeMinimizer minimizer;
		float const dummy_tol( 0.00000025);
		bool const use_nblist( true );
		MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
		options.nblist_auto_update( true );


		core::kinematics::MoveMap mm;

		mm.set_bb(false);	
		mm.set_chi(false);

		for(Size n=1; n<=FB_JP_list.size(); n++){

			FloatingBaseChainClosureJobParameter const & FB_job_params=FB_JP_list[n];

			mm.set( TorsionID( FB_job_params.bulge_res-1 , id::BB,  5 ), true ); //epsilon
			mm.set( TorsionID( FB_job_params.bulge_res-1 , id::BB,  6 ), true ); //zeta

			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  1 ), true ); //alpha
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  2 ), true ); //beta
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  3 ), true ); //gamma
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  5 ), true ); //epsilon
			mm.set( TorsionID( FB_job_params.bulge_res , id::BB,  6 ), true ); //zeta
			mm.set( TorsionID( FB_job_params.bulge_res , id::CHI, 1 ), true ); //chi (torsion between base and ribose sugar)

			mm.set( TorsionID( FB_job_params.bulge_res+1 , id::BB,  1 ), true ); //alpha
			mm.set( TorsionID( FB_job_params.bulge_res+1 , id::BB,  2 ), true ); //beta 
			mm.set( TorsionID( FB_job_params.bulge_res+1 , id::BB,  3 ), true ); //gamma

		}

		for(Size n=1; n<=pose_data_list.size(); n++){

			viewer_pose= (*pose_data_list[n].pose_OP);

			utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters->working_moving_partition_pos();
			utility::vector1 < core::Size > already_virtualized_res_list;

			if(virtual_ribose_is_from_prior_step){ //Virtualize the other partition since it doesn't exist in prior step!

				for(Size ii=1; ii<=working_moving_partition_pos.size(); ii++){
					Size const seq_num=working_moving_partition_pos[ii];
					if(viewer_pose.residue(seq_num).has_variant_type("VIRTUAL_RNA_RESIDUE")){
						already_virtualized_res_list.push_back(seq_num);
						continue;
					}

					pose::add_variant_type_to_pose_residue( viewer_pose, "VIRTUAL_RNA_RESIDUE", seq_num );

				}
	
				if (job_parameters->gap_size() == 0) pose::add_variant_type_to_pose_residue( viewer_pose, "VIRTUAL_PHOSPHATE", job_parameters->five_prime_chain_break_res()+1 );
			}

			std::cout << "minimize_all_sampled_floating_bases pose # " << n << " out of " << pose_data_list.size() << " " << std::endl;;
			minimizer.run( viewer_pose, mm, (*sampling_scorefxn), options );
//			o2star_minimize(viwer_pose, sampling_scorefxn)
//			minimizer.run( viewer_pose, mm, (*sampling_scorefxn), options );

			if(virtual_ribose_is_from_prior_step){ //Virtualize the other partition since it doesn't exist in prior step!

				for(Size ii=1; ii<=working_moving_partition_pos.size(); ii++){ 
					Size const seq_num=working_moving_partition_pos[ii];

					if(Contain_seq_num(seq_num, already_virtualized_res_list)) continue;

					pose::remove_variant_type_from_pose_residue( viewer_pose, "VIRTUAL_RNA_RESIDUE", seq_num );
				}
	
				if (job_parameters->gap_size() == 0) pose::remove_variant_type_from_pose_residue( viewer_pose, "VIRTUAL_PHOSPHATE", job_parameters->five_prime_chain_break_res()+1 );
			}


			(*pose_data_list[n].pose_OP)=viewer_pose;
		}

		viewer_pose=viewer_pose_copy;

		Output_title_text("Exit minimize_all_sampled_floating_bases");

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	Is_ribose_virtual( core::pose::Pose const & pose, core::Size const ribose_res, core::Size const bulge_res){
	
		using namespace core::chemical;
		using namespace core::conformation;
		using namespace ObjexxFCL;

		Size const nres= pose.total_residue();

		if( ( ribose_res+1)!= bulge_res && ( ribose_res-1)!= bulge_res ) {
			std::cout << "ribose_res= " << ribose_res << " bulge_res= " << bulge_res << std::endl;
			utility_exit_with_message( "( ribose_res+1)!= bulge_res && ( ribose_res-1)!= bulge_res)" );
		}

		if(ribose_res<1 || ribose_res> nres){
			utility_exit_with_message( "ribose_res<1 || ribose_res> nres("+ string_of(nres) +")!. ribose_res= " + string_of(ribose_res) );
		}

		if(pose.residue(ribose_res).has_variant_type("VIRTUAL_RIBOSE") ) {

			//ok consistency checks:			
			if(bulge_res<1 || bulge_res> nres){
				utility_exit_with_message( "bulge_res<1 || bulge_res> nres("+ string_of(nres) +")!. bulge_res= " + string_of(bulge_res) );
			}

			if(pose.residue(bulge_res).has_variant_type("VIRTUAL_RNA_RESIDUE")==false){
				utility_exit_with_message("pose.residue(bulge_res).has_variant_type(\"VIRTUAL_RNA_RESIDUE\")==false" );
			}

			return true;
		}else{
			return false;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	copy_bulge_res_and_ribose_torsion(FloatingBaseChainClosureJobParameter const & FB_job_params, core::pose::Pose & pose, core::pose::Pose const & template_pose){

		using namespace ObjexxFCL;
		using namespace core::io::silent;
		using namespace core::id;

		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", FB_job_params.moving_res );

		std::map< core::Size, core::Size > res_map;  //This is map from sub numbering to input_res numbering..
		res_map[ FB_job_params.moving_res ] = FB_job_params.moving_res;
		res_map[ FB_job_params.bulge_res ] = FB_job_params.bulge_res;
		res_map[ FB_job_params.reference_res ] = FB_job_params.reference_res;
						
		//copy_dofs( pose, template_pose, res_map, true , false /*not verbose*/ ); //need to make sure this works for every case!

		//Dec 24, 2011 Parin S.:Convert to Rhiju's NEW version
		copy_dofs_match_atom_names( pose, template_pose, res_map, false /*backbone_only*/, false /*ignore_virtual*/);

	}

	//////////////////////////////////////////////////////////////////////////
	bool
	sort_pose_data_by_score(pose_data_struct2  pose_data_1, pose_data_struct2 pose_data_2) {  //Duplicate!
		return (pose_data_1.score < pose_data_2.score);
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////July 21, 2011 Move from StepWiseRNA_ResidueSampler.cc ////////////////////////////////////////////////
	////////Code used to be part of the function StepWiseRNA_ResidueSampler::previous_floating_base_chain_closure();

	//Dec 9, 2010
	//Warning, currently the code doesn't handle the case of i-1, i-3 (and potentially i-5 and so on) bulges very well. Basically will sample the sugar conformation of base i. But even though base i-2, i-4 have virtual riboses, these sugars are not sampled (use the one that currently exist in the pose) which is NOT a complete fail, since the ribose did pass distance check of Check_chain_closable()(correctly implement this on Dec 9,2010) and VDW screening.
	//Possible ways to fix this (HOWEVER lets not fix this until we actual find a REAL existing case where i-1 and i-3 are bulge res):

	//Will have to at least rewrite both floating_base_chain_closure_setup() and floating_base_chain_closure_sampling() to sample both i and i-2 ribose. 
	//Hardest 1. To change this would have to recursively sample the conformation of the bulge and corresponding virtual ribose starting from the bulge that is futhest away. If assuem that each bulge conformation is independent of the other bulges then this is not very computationally expensive. Although this would require sampling two ribose for each bulge except for the bulge that is futhest again.   A trace back is then needed to find combinations of all the bulge conformation which is closable.
	//Medium 2. Another possibility is to simple assume that i-2 base can assume all bulge conformation subjected to check_chain_closable.
	//Easiest 3. Just set chain_closure_sampling to false in this function. This will just keep all possible bulge conformation in pose_data_list that pass O3i_C5iplus2_distance>O3I_C5IPLUS2_MAX_DIST (in setup)
	//However in Easiest 3. case, there is a problem is problem in that O3i_C5iplus2_distance has determined with fixed i-2 ribose.

	utility::vector1< pose_data_struct2 > 
	sample_virtual_ribose_and_bulge_and_close_chain(pose::Pose & viewer_pose, 
																							 FloatingBaseChainClosureJobParameter const & FB_job_params, 
																							 std::string const name,
																							 core::scoring::ScoreFunctionOP const & scorefxn, 
																							 core::scoring::ScoreFunctionOP const & sampling_scorefxn, 
																							 core::scoring::ScoreFunctionOP const & atr_rep_screening_scorefxn, 
																							 core::scoring::ScoreFunctionOP const & chainbreak_scorefxn,
																							 StepWiseRNA_JobParametersCOP & job_parameters,
																							 bool const virtual_ribose_is_from_prior_step){

		using namespace ObjexxFCL;
		using namespace core::io::silent;

		Output_title_text("Enter sample_virtual_ribose_and_bulge_and_close_chain()");	
		clock_t const time_start( clock() ); 

		pose::Pose const viewer_pose_copy=viewer_pose; //BACKUP

		FB_job_params.check_compatibility(viewer_pose.total_residue());  		

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///Switch order so that virtualize before creating VDW_bin_screener July 18, 2010

		pose_data_struct2 pose_data;
		pose_data.pose_OP=new pose::Pose;
		(*pose_data.pose_OP)= viewer_pose;
		pose_data.score = 0.0;
		pose_data.tag="";

		pose::Pose & input_pose= (*pose_data.pose_OP);

		utility::vector1 < core::Size > other_partition_pos; //June 13, 2011..Before used to just use moving_partition_pos here since bulge_res was always in the root_partition!
		utility::vector1 < core::Size > already_virtualized_res_list;
		////////////////////////////////////////////////////////////
		if(virtual_ribose_is_from_prior_step){ //Virtualize the other partition since it doesn't exist in prior step!

			Size const nres = job_parameters->working_sequence().size();
			ObjexxFCL::FArray1D< bool > const & partition_definition = job_parameters->partition_definition();
			bool const bulge_res_partition_value=partition_definition(FB_job_params.bulge_res);

			if( bulge_res_partition_value!=partition_definition(FB_job_params.moving_res) ){ //Check that these three nts are in the same paritition!
				utility_exit_with_message("bulge_res_partition_value!=partition_definition(FB_job_params.moving_res)");
			}

			if( bulge_res_partition_value!=partition_definition(FB_job_params.reference_res) ){ //Check that these three nts are in the same paritition!
				utility_exit_with_message("bulge_res_partition_value!=partition_definition(FB_job_params.reference_res)");
			}

			for(Size seq_num=1; seq_num<=nres; seq_num++){
				if ( partition_definition( seq_num ) != bulge_res_partition_value ) other_partition_pos.push_back( seq_num );
			}

			for(Size ii=1; ii<=other_partition_pos.size(); ii++){
				Size const seq_num=other_partition_pos[ii];
				if(input_pose.residue(seq_num).has_variant_type("VIRTUAL_RNA_RESIDUE")){
					already_virtualized_res_list.push_back(seq_num);
					continue;
				}
				pose::add_variant_type_to_pose_residue( input_pose, "VIRTUAL_RNA_RESIDUE", seq_num );
			}

			if(job_parameters->gap_size() == 0){
				if(input_pose.residue(job_parameters->five_prime_chain_break_res()+1).has_variant_type("VIRTUAL_PHOSPHATE")==true){
					utility_exit_with_message( "input_pose.residue(job_parameters_->five_prime_chain_break_res()+1).has_variant_type(\"VIRTUAL_PHOSPHATE\")==true" );
				}
				pose::add_variant_type_to_pose_residue( input_pose, "VIRTUAL_PHOSPHATE", job_parameters->five_prime_chain_break_res()+1);
			}
		}else{
			if(job_parameters->gap_size() == 0) utility_exit_with_message( "job_parameters_->gap_size() == 0" ); //If virt_ribse is from current step, then it should not be the last step!
		}
		////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////create_VDW_screen_bin//////////////////////////////////////////////////////////////////////////////////////////

		core::kinematics::Stub reference_stub;
		reference_stub.v=core::scoring::rna::get_rna_base_centroid(  input_pose.residue( FB_job_params.reference_res ) , true);
		reference_stub.M=core::scoring::rna::get_rna_base_coordinate_system( input_pose.residue( FB_job_params.reference_res ) , reference_stub.v); 		


		utility::vector1 < core::Size > ignore_res_list;
		ignore_res_list.push_back(FB_job_params.moving_res); 
		ignore_res_list.push_back(FB_job_params.bulge_res); 
		ignore_res_list.push_back(FB_job_params.reference_res);	

 	  StepWiseRNA_VDW_Bin_ScreenerOP prev_floating_base_VDW_bin_screener= new StepWiseRNA_VDW_Bin_Screener();	

		//Feb 21, 2011...one thing is that the 2'-OH hydrogen is not virtualized here...
		//But OK, since consistent with the actual input_pose not having virtualized 2'-OH.
		prev_floating_base_VDW_bin_screener->create_VDW_screen_bin( input_pose, ignore_res_list, FB_job_params.Is_prepend, reference_stub.v, true /*verbose*/ );

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


		utility::vector1< pose_data_struct2 > input_pose_data_list; //Singleton list.
		input_pose_data_list.push_back(pose_data);

		utility::vector1< FB_Pose_Data > pose_data_list= floating_base_chain_closure_setup(input_pose_data_list, FB_job_params, atr_rep_screening_scorefxn, scorefxn, viewer_pose, true /*do_minimize*/);

		bool const chain_closure_sampling=true; 

		if(chain_closure_sampling){
			floating_base_chain_closure_sampling(pose_data_list, viewer_pose, FB_job_params, chainbreak_scorefxn, atr_rep_screening_scorefxn, prev_floating_base_VDW_bin_screener , false);
		}else{
			for(Size n=1; n<=pose_data_list.size(); n++) {pose_data_list[n].Is_chain_close=true; }
		}

		bool const rm_CB_JP_during_post_process=true;
	
		utility::vector1< pose_data_struct2 > final_pose_data_list=floating_base_chain_closure_post_process(pose_data_list, viewer_pose, sampling_scorefxn, FB_job_params, rm_CB_JP_during_post_process);

		SilentFileData silent_file_data;
		
		for(Size n=1; n<=final_pose_data_list.size(); n++){
			pose::Pose & current_pose=(*final_pose_data_list[n].pose_OP);

			if(false) Output_data(silent_file_data, "post_process_" + name + ".out", name + final_pose_data_list[n].tag , false, current_pose, job_parameters->working_native_pose(), job_parameters);		
			////////////////////////////////////////////////////////////
			if(virtual_ribose_is_from_prior_step){ //Virtualize the other partition since it doesn't exist in prior step!

				for(Size ii=1; ii<=other_partition_pos.size(); ii++){
					Size const seq_num=other_partition_pos[ii];
					if(Contain_seq_num(seq_num, already_virtualized_res_list)) continue;
					pose::remove_variant_type_from_pose_residue( current_pose, "VIRTUAL_RNA_RESIDUE", seq_num);
				}
				if (job_parameters->gap_size() == 0) pose::remove_variant_type_from_pose_residue( current_pose, "VIRTUAL_PHOSPHATE", job_parameters->five_prime_chain_break_res()+1 );
			}else{
				if(job_parameters->gap_size() == 0) utility_exit_with_message( "job_parameters->gap_size() == 0" ); //If virt_ribse is from current step, then it should not be the last step!
			}
			////////////////////////////////////////////////////////////

		}

		std::cout << "Time in sample_virtual_ribose_and_bulge_and_close_chain(): " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

		Output_title_text("Exit sample_virtual_ribose_and_bulge_and_close_chain()");

		viewer_pose=viewer_pose_copy;

		return final_pose_data_list;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	enumerate_starting_pose_data_list(utility::vector1< pose_data_struct2 > & starting_pose_data_list,
																	utility::vector1< FloatingBaseChainClosureJobParameter > const & FB_CC_JP_list, 
																	core::pose::Pose const & pose){

		pose::Pose const pose_copy= pose;

		utility::vector1< Size > sugar_ID_counter_list(FB_CC_JP_list.size(), 1);

		while(true){

			pose_data_struct2 start_pose_data; 		

			start_pose_data.pose_OP=new pose::Pose;
			(*start_pose_data.pose_OP)=pose_copy;
			pose::Pose & start_pose=(*start_pose_data.pose_OP);

			start_pose_data.score=0;
			start_pose_data.tag="";

			for(Size n=1; n<=FB_CC_JP_list.size(); n++){

				FloatingBaseChainClosureJobParameter const & curr_FB_JP = FB_CC_JP_list[n];
				Size const sugar_ID=sugar_ID_counter_list[n];

				start_pose_data.tag+= curr_FB_JP.PDL[sugar_ID].tag;
				copy_bulge_res_and_ribose_torsion(curr_FB_JP, start_pose, (*curr_FB_JP.PDL[sugar_ID].pose_OP) );

			}

			starting_pose_data_list.push_back(start_pose_data);

			///////////////////////////Counter/////////////////////////////
			sugar_ID_counter_list[1]++;

			for(Size n=1; n<FB_CC_JP_list.size(); n++){
				if( sugar_ID_counter_list[n]==(FB_CC_JP_list[n].PDL.size()+1) ){
					 sugar_ID_counter_list[n]=1;
					 sugar_ID_counter_list[n+1]++;
				}
			}

			if( sugar_ID_counter_list[FB_CC_JP_list.size()]==(FB_CC_JP_list[FB_CC_JP_list.size()].PDL.size()+1) ) break;
			////////////////////////////////////////////////////////////////

		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< FloatingBaseChainClosureJobParameter > 
	setup_FB_CC_JP_list(pose::Pose const & pose, utility::vector1< std::string > const & sample_virtual_ribose_string_list, StepWiseRNA_JobParametersCOP & job_parameters){

		using namespace ObjexxFCL;

		utility::vector1< FloatingBaseChainClosureJobParameter > FB_CC_JP_list;

		std::string const & working_sequence= job_parameters->working_sequence();

		for(Size n=1; n<=sample_virtual_ribose_string_list.size(); n++){

			utility::vector1< std::string > const tokenize_list=Tokenize(sample_virtual_ribose_string_list[n], "-");
			if(tokenize_list.size()!=2) utility_exit_with_message("tokenize_list!=2");

			if(tokenize_list[2]!="A" && tokenize_list[2]!="P"){
				utility_exit_with_message("tokenize_list[2]!=\"A\" && tokenize_list[2]!=\"P\" (" + tokenize_list[2] + ")" );
			}

			bool const Is_prepend= (tokenize_list[2]=="P") ? true : false;
	
			Size const full_ribose_res=string_to_int( tokenize_list[1] );
			Size const full_bulge_res=(Is_prepend) ? full_ribose_res+1 : full_ribose_res-1;
			Size const full_ref_res=  (Is_prepend) ? full_ribose_res+2 : full_ribose_res-2;

			std::cout << "Case: " << sample_virtual_ribose_string_list[n];
			std::cout << " full_ribose_res= " << full_ribose_res << " full_bulge_res= " << full_bulge_res << " full_ref_res= " << full_ref_res;
			Output_boolean(" Is_prepend= " , Is_prepend); 

			if( pose.total_residue()!=working_sequence.size() ){
				utility_exit_with_message( "pose.total_residue()=("+string_of(pose.total_residue())+")!="+string_of(working_sequence.size())+") working_sequence().size()");
			}

			if(check_is_working_res(full_ribose_res, job_parameters)){

				Size const working_ribose_res=check_validity_and_get_working_res(full_ribose_res, job_parameters);

				bool const ribose_is_virtual=pose.residue(working_ribose_res).has_variant_type("VIRTUAL_RIBOSE"); 

				std::cout << " | working_ribose_res= " << working_ribose_res;
				Output_boolean(" ribose_is_virtual= " , ribose_is_virtual); 

				if(ribose_is_virtual){

					Size const working_bulge_res=check_validity_and_get_working_res(full_bulge_res, job_parameters);
					Size const working_ref_res=check_validity_and_get_working_res(full_ref_res, job_parameters);

					std::cout << " | working_bulge_res= " << working_bulge_res << " working_ref_res= " << working_ref_res;
					
					if(Is_prepend){
						if(working_ribose_res!=(working_bulge_res-1)) utility_exit_with_message("prepend but working_ribose_res!=(working_bulge_res-1)");
						if(working_ribose_res!=(working_ref_res - 2)) utility_exit_with_message("prepend but working_ribose_res!=(working_ref_res - 2)");
					}else{
						if(working_ribose_res!=(working_bulge_res+1)) utility_exit_with_message("append  but working_ribose_res!=(working_bulge_res+1)");
						if(working_ribose_res!=(working_ref_res + 2)) utility_exit_with_message("prepend but working_ribose_res!=(working_ref_res + 2)");
					}

					//ribose_is_virtual=Is_ribose_virtual(pose, working_ribose_res, working_bulge_res);

					if(pose.residue(working_bulge_res).has_variant_type("VIRTUAL_RNA_RESIDUE")==false){
						utility_exit_with_message("pose.residue(working_bulge_res).has_variant_type(\"VIRTUAL_RNA_RESIDUE\")==false" );
					}
				
					FloatingBaseChainClosureJobParameter curr_FB_JP=FloatingBaseChainClosureJobParameter();

					curr_FB_JP=FloatingBaseChainClosureJobParameter(working_ribose_res, working_ref_res);

					FB_CC_JP_list.push_back(curr_FB_JP);

				}


			}else{

				std::cout << " | full_ribose_res is not a working res! ";

			}

			std::cout << std::endl;
		}

		return FB_CC_JP_list;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void
	sample_user_specified_virtual_riboses(pose::Pose & pose, utility::vector1< std::string > const & sample_virtual_ribose_string_list, 
																	    StepWiseRNA_JobParametersCOP & job_parameters , core::scoring::ScoreFunctionOP const & scorefxn, 
																			std::string const silent_file_out, std::string const input_tag){


		using namespace ObjexxFCL;
		using namespace core::io::silent;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::conformation;

		Output_title_text("Enter StepWiseRNA_VirtualRiboseSampler::sample_virtual_ribose");	

		// clock_t const time_start( clock() ); // Unused variable causes warning.
 

		/////////////////////Copy the conformation but nothing else. No energy and no cache data (having cache data can cause problem with column_name order in output silent_file!)//////////////////////////////

		ConformationOP copy_conformation = new Conformation();

		(*copy_conformation)=pose.conformation();

		pose::Pose new_pose;
		new_pose.set_new_conformation( copy_conformation );		

		pose=new_pose;

		///////////////////////////////////////////////////////////////

		pose::Pose const pose_save = pose;
		pose = pose_save; //this recopy is useful for triggering graphics.

		///////////////////////////////////////////////////////////////



		core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn, chainbreak_scorefxn, sampling_scorefxn, o2star_pack_scorefxn;

		initialize_common_scorefxns(scorefxn, sampling_scorefxn, atr_rep_screening_scorefxn, chainbreak_scorefxn, o2star_pack_scorefxn);

		utility::vector1< FloatingBaseChainClosureJobParameter > FB_CC_JP_list=setup_FB_CC_JP_list(pose, sample_virtual_ribose_string_list, job_parameters);

		std::cout << "num_virtual_ribose= " << FB_CC_JP_list.size() << std::endl;

		if(FB_CC_JP_list.size()==0){
			std::cout << "no_virtual_ribose (FB_CC_JP_list.size()==0). EARLY RETURN/NO OUTPUT SILENT_FILE!" << std::endl;			
	
			std::ofstream outfile;
			outfile.open(silent_file_out.c_str()); //Opening the file with this command removes all prior content..
			outfile << "no_virtual_ribose (FB_CC_JP_list.size()==0).\n";
			outfile.flush();
			outfile.close();

			return;
		}


		for(Size n=1; n<=FB_CC_JP_list.size(); n++){

			FloatingBaseChainClosureJobParameter & curr_FB_JP = FB_CC_JP_list[n];

			curr_FB_JP.set_base_and_pucker_state(pose, job_parameters);

			curr_FB_JP.PDL=sample_virtual_ribose_and_bulge_and_close_chain(pose,  curr_FB_JP, "VIRT_RIBOSE_NUM_" + string_of(n), //ACTUAL COMPUTATION OCCUR HERE!!
																																	scorefxn, sampling_scorefxn, atr_rep_screening_scorefxn, chainbreak_scorefxn, job_parameters, 
																																	false /*virtual_ribose_is_from_prior_step*/);

			std::sort(curr_FB_JP.PDL.begin(), curr_FB_JP.PDL.end(), sort_pose_data_by_score);


			if(curr_FB_JP.PDL.size()==0){
				std::cout << "Case n= " << n << " Is_sugar_virt==True but curr_FB_JP.PDL.size()==0. EARLY RETURN!" << std::endl;

				std::ofstream outfile;
				outfile.open(silent_file_out.c_str()); //Opening the file with this command removes all prior content..
				outfile << "num_virtual_ribose != 0 but for one of the sampled virtual_ribose, curr_FB_JP.PDL.size()==0.\n";
				outfile.flush();
				outfile.close();
				return;
			}
		}

		///////////////////////////////////////////////////////////////////////////////////
		utility::vector1< pose_data_struct2 > starting_pose_data_list;

		enumerate_starting_pose_data_list(starting_pose_data_list, FB_CC_JP_list, pose);

		minimize_all_sampled_floating_bases(pose, FB_CC_JP_list, starting_pose_data_list, sampling_scorefxn, job_parameters, false /*virtualize_other_partition*/);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		SilentFileData silent_file_data;
		for(Size n=1; n<=starting_pose_data_list.size(); n++){
			pose=(*starting_pose_data_list[n].pose_OP); //set viewer_pose;

			std::string starting_pose_tag=input_tag + "_sample_ribose" + starting_pose_data_list[n].tag;

			if (job_parameters->gap_size() == 0) utility_exit_with_message( "job_parameters_->gap_size() == 0" );

			(*scorefxn)(pose);

			//std::cout << "starting_pose_tag= " << starting_pose_tag= << std::endl;
			//pose.energies().show();


			Output_data(silent_file_data, silent_file_out , starting_pose_tag , false, pose, job_parameters->working_native_pose(), job_parameters);		

		}

	}


	/////////////////////////////////////////////////////////////////////////////////////


} //rna
} //swa
} // protocols

