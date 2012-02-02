// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_OutputData (Created on Sept 26, 2011)
/// @brief Output silent_file_data functions for Stepwise Assembly RNA.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh> //Oct 22, 2011...Not sure why the code worked without this!
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_FloatingBase_Sampler_Util.hh> //Sept 26, 2011
#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>


#include <protocols/rna/RNA_BasePairClassifier.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/rna/RNA_Util.hh>

#include <core/scoring/ScoreType.hh> //Parin Sept 20, 2011.
//////////////////////////////////

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
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
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <numeric/conversions.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <set>
#include <time.h>
#include <map>

#include <stdio.h> //Sept 26, 2011

//for process_mem_usage:
#include <ios>


using namespace core;

namespace protocols {
namespace swa {
namespace rna {



	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	core::io::silent::BinaryRNASilentStruct
	get_binary_rna_silent_struct_safe(pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::chemical;
		using namespace core::conformation;


		std::string const debug_silent_file=silent_file+"_CONVERSION_DEBUG";
		std::string const debug_tag=tag+"_CONVERSION_DEBUG";

		SilentFileData silent_file_data;

		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

		Size NUM_trails=10;
		Real const local_angle_bin_size=20;
		Real const local_z_bin_size=0.05;

		pose::Pose first_trail_pose_from_silent_file;
		pose::Pose first_trail_pose;



		for(Size trail_num=1; trail_num<=NUM_trails; trail_num++){ //Found that just rigid problem rotation of the pose solves the silent_file conversion problem

			pose::Pose pose=const_pose;

			if(trail_num!=1){

				///////////////////////////get centoid of the structure/////////////////////////////

				numeric::xyzVector<core::Real> centroid=Vector( 0.0, 0.0, 0.0 );
				Size numatoms = 0;

				for(Size seq_num=1; seq_num<=pose.total_residue(); seq_num++){
	
					conformation::Residue const & rsd(pose.residue(seq_num));

					for( Size at = 1; at <= rsd.natoms(); at++){

						if(rsd.is_virtual(at)) continue;
							
			    		centroid += rsd.xyz(at);
 					  	numatoms++;
					}
				}

				if(numatoms==0) utility_exit_with_message("numatoms==0");

				centroid = centroid/numatoms;

				////////////////////////////////////////////////////////////////////////////////


				Euler_angles euler_angles;

				Matrix rotation_matrix;

				euler_angles.alpha=(0.25*trail_num)*local_angle_bin_size*(PI/180); 
				euler_angles.gamma=(0.25*trail_num)*local_angle_bin_size*(PI/180); 

				euler_angles.z=(0.25*trail_num)*local_z_bin_size;	//MAKE SURE THIS DOESN'T GET OUT OF BOUND!
				euler_angles.beta=acos(euler_angles.z);

				convert_euler_to_coordinate_matrix(euler_angles, rotation_matrix);

				//numeric::xyzVector<core::Real> const offset_vector=Vector( 0.0,   0.0,   1.0*(trail_num-1) );

				for(Size seq_num=1; seq_num<=pose.total_residue(); seq_num++){
	
					conformation::Residue const & rsd(pose.residue(seq_num));

					for( Size at = 1; at <= rsd.natoms(); at++){

						id::AtomID const id( at, seq_num);

						pose.set_xyz( id, pose.xyz( id) - centroid ); //This should minimize the error introoduced by the right body rotation!
						pose.set_xyz( id, rotation_matrix * pose.xyz(id));  
					}
				}
			}

			BinaryRNASilentStruct DEBUG_silent_struct( pose, debug_tag );
			BinaryRNASilentStruct const silent_struct( pose, tag );

			if(file_exists(debug_silent_file)) remove_file(debug_silent_file);
	
			silent_file_data.write_silent_struct(DEBUG_silent_struct, debug_silent_file, false);

			///////////////////////////////////////////////////////////////////////////////
	
			core::io::silent::SilentFileData import_silent_file_data;
			import_silent_file_data.read_file( debug_silent_file );
			pose::Pose pose_from_silent_file;

			bool found_tag= false;
			Size num_struct=0;		

			for ( core::io::silent::SilentFileData::iterator iter = import_silent_file_data.begin(), end = import_silent_file_data.end(); iter != end; ++iter ){
				num_struct+=1;
				if ( iter->decoy_tag() != debug_tag ) continue;
				found_tag = true;
				iter->fill_pose( pose_from_silent_file, *rsd_set );
			}

			if(num_struct!=1) utility_exit_with_message( "num_struct=("+string_of(num_struct)+")!=1");
			if(found_tag==false ) utility_exit_with_message( "Could not find specified tag (" + debug_tag + ") in silent file (" + debug_silent_file + ")!" );

			if(file_exists(debug_silent_file)==false){
				utility_exit_with_message("debug_silent_file ("+debug_silent_file+") SHOULD exist!");
			} 

			remove_file(debug_silent_file);

			if(trail_num==1){
				first_trail_pose_from_silent_file=pose_from_silent_file;
				first_trail_pose=pose;
			}

			//pose.dump_pdb( "SILENT_FILE_CONVERSION_TEST_" + tag +"_TRAIL_" + lead_zero_string_of(trail_num, 3) + ".pdb" );
			//pose_from_silent_file.dump_pdb( "IMPORTED_SILENT_FILE_CONVERSION_TEST_" + tag +"_TRAIL_" + lead_zero_string_of(trail_num, 3) + ".pdb" );

			if(check_for_messed_up_structure(pose_from_silent_file, debug_tag)==false){

				return silent_struct;

			}else{
			 	std::cout << "WARNING: Problem with writing pose (" << debug_tag << ") to silent_file [Attempt #" << trail_num << "]" << std::endl;   
			}

		}


		first_trail_pose_from_silent_file.dump_pdb( "SILENT_FILE_CONVERSION_PROBLEM_" + tag + "_pose_from_silent_file.pdb" );
		first_trail_pose.dump_pdb( "SILENT_FILE_CONVERSION_PROBLEM_" + tag + ".pdb" );
		BinaryRNASilentStruct ERROR_silent_struct( first_trail_pose, debug_tag );
		std::string const ERROR_silent_file="SILENT_FILE_CONVERSION_PROBLEM_" + tag + ".out";
		silent_file_data.write_silent_struct(ERROR_silent_struct, ERROR_silent_file, false);

		utility_exit_with_message("Fail to write pose (" + debug_tag + ") to silent_file after "+string_of(NUM_trails)+" trails ");

		////////////This is just to prevent compiler WARNING MESSAGES/////////
		BinaryRNASilentStruct EMPTY_silent_struct;
		return EMPTY_silent_struct;
		//////////////////////////////////////////////////////////////////////

	}


	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	core::io::silent::BinaryRNASilentStruct
	get_binary_rna_silent_struct_safe_wrapper(pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file, bool const write_score_only){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;

		if(write_score_only){

			BinaryRNASilentStruct s( const_pose, tag ); //If write score only, don't have to safe about pose to silent_struct conversion!
			return s;

		}else{

			return (get_binary_rna_silent_struct_safe( const_pose, tag, silent_file));

		}

		////////////This is just to prevent compiler WARNING MESSAGES/////////
		BinaryRNASilentStruct EMPTY_silent_struct;
		return EMPTY_silent_struct;
		//////////////////////////////////////////////////////////////////////

	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//Accept the job_parameter instead.
	void
	Output_data(core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, std::string const & tag, bool const write_score_only, pose::Pose const & pose, core::pose::PoseCOP native_poseCOP, StepWiseRNA_JobParametersCOP job_parameters_){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;

		utility::vector1 < core::Size > const & rmsd_res_list = job_parameters_->rmsd_res_list();
		std::map< core::Size, core::Size > const & full_to_sub = job_parameters_->const_full_to_sub();
		std::map< core::Size, bool > const & Is_prepend_map = job_parameters_->Is_prepend_map();
		bool const Is_prepend(  job_parameters_->Is_prepend() ); // if true, moving_suite+1 is fixed. Otherwise, moving_suite is fixed.
		Size const moving_base_residue( job_parameters_->actually_moving_res() );

		BinaryRNASilentStruct s=get_binary_rna_silent_struct_safe_wrapper( pose, tag, silent_file, write_score_only);

		//s.print_header( std::cout );
		//s.precision(5); REALLY COOL. SET higher precision so that their is no energy rank ambiguity!.


		if ( native_poseCOP ) {
			s.add_energy( "all_rms", rms_at_corresponding_heavy_atoms( pose, *native_poseCOP ) );

			if(write_score_only){ //Basically the optimal alignment, align working_res as well if it is part of the alignment res list.

				// This assumes that pose and native_pose are correctly syperimposed.
				// I added a function in Pose_Setup to make sure this happens. Parin Jan 28, 2010

				s.add_energy( "rmsd", suite_rmsd( pose, *native_poseCOP, moving_base_residue, Is_prepend, false));
				s.add_energy( "loop_rmsd", rmsd_over_residue_list( pose, *native_poseCOP, rmsd_res_list, full_to_sub, Is_prepend_map, false, false) );

				s.add_energy( "V_rms", suite_rmsd( pose, *native_poseCOP, moving_base_residue, Is_prepend, true));
				s.add_energy( "V_loop_rms", rmsd_over_residue_list( pose, *native_poseCOP, rmsd_res_list, full_to_sub, Is_prepend_map, false, true) );

				if(job_parameters_->gap_size()==0){
					s.add_energy( "PBP_rmsd", phosphate_base_phosphate_rmsd( pose, *native_poseCOP, moving_base_residue,  false));
				}else{
					s.add_energy( "PBP_rmsd", 0.0);
				}

			}else{

				pose::Pose current_pose=pose; //hard copy, computationally expensive

				utility::vector1< core::Size > const & working_native_alignment = job_parameters_->working_native_alignment();	
				utility::vector1< core::Size > const & working_best_alignment = job_parameters_->working_best_alignment();

				if(working_native_alignment.size()!=0){ //user specify which residue to align with native.
					align_poses(current_pose, tag, *native_poseCOP, "native", working_native_alignment);
				}else{ //default
					align_poses(current_pose, tag, *native_poseCOP, "native", working_best_alignment);
				}
				s.add_energy( "O_rmsd", suite_rmsd( current_pose, *native_poseCOP, moving_base_residue, Is_prepend, false));
				s.add_energy( "O_loop_rmsd", rmsd_over_residue_list( current_pose, *native_poseCOP, rmsd_res_list, full_to_sub, Is_prepend_map, false, false) );

				s.add_energy( "O_V_rms", suite_rmsd( current_pose, *native_poseCOP, moving_base_residue, Is_prepend, true));
				s.add_energy( "O_V_loop_rms", rmsd_over_residue_list( current_pose, *native_poseCOP, rmsd_res_list, full_to_sub, Is_prepend_map, false, true) );

				if(job_parameters_->gap_size()==0){
					s.add_energy( "O_PBP_rmsd", phosphate_base_phosphate_rmsd( current_pose, *native_poseCOP, moving_base_residue,  false));
				}else{
					s.add_energy( "O_PBP_rmsd", 0.0);
				}
			}

		}


		silent_file_data.write_silent_struct(s, silent_file, write_score_only);
					
	}


}
}
}
