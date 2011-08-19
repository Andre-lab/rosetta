// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
// Protocol Headers
#include <core/init.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <basic/options/option_macros.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <core/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>

#include <basic/options/option_macros.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;


OPT_KEY( String, central_relax_pdb_num )
OPT_KEY( String, input_protein_file )
OPT_KEY( String, input_ligand_file )
OPT_KEY( Integer, num_poses )
OPT_KEY( Integer, num_angles )
OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Boolean, cheat )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, origin_cutoff )
//OPT_KEY( Integer, angle_increment )

int main( int argc, char * argv [] ) {

  NEW_OPT( central_relax_pdb_num, "target residue", "-1" );
  NEW_OPT( input_protein_file, "protein file name", "protein.pdb" );
  NEW_OPT( input_ligand_file, "ligand file name", "ligand.pdb" );
	NEW_OPT( num_poses, "No. of poses to search for initail stochastic search", 1000 );
	NEW_OPT( num_angles, "Number of different pose angles to measure score at", 1);
	NEW_OPT( num_runs, "no. of runs for PSO", 100 );
  NEW_OPT( num_particles, "no. of particles for PSO", 10 );
	NEW_OPT( cheat, "move pocket CoM over Ligand CoM", false );
  NEW_OPT( steric_weight, "steric weight for PSO", 3.0 );
  NEW_OPT( missing_point_weight, "missing point weight", 25.0 );
  NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 5.0 );
	//NEW_OPT( angle_increment, "angle increment", 20.0 );

	core::init(argc, argv);

	std::string const input_protein = option[ input_protein_file ];
	std::string const input_ligand = option[ input_ligand_file ];
	int num_pose_search  = option[ num_poses ];
	std::string const resid = option[ central_relax_pdb_num ];
  int angles = option[ num_angles ];
  int particle_size = option[ num_particles ];
  int run_size = option[ num_runs ];
	core::Real const steric_wt = option[ steric_weight ];
	core::Real const missing_pt_wt = option[ missing_point_weight ];
	core::Real const origin_space = option[ origin_cutoff ];
	//int const ang_inc  = option[ angle_increment ];

	protocols::pockets::NonPlaidFingerprint npf;

	pose::Pose protein_pose;
	core::import_pose::pose_from_pdb( protein_pose, input_protein );

	int  central_relax_residue_number;
  char chain = ' ';
	std::size_t fpos( resid.find(':') );
  if ( fpos != std::string::npos ) {
    central_relax_residue_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
    if (fpos != resid.size()-1 ) {
      chain = resid[ fpos+1 ];
    }
  } else {
    central_relax_residue_number = ObjexxFCL::int_of( resid );
  }
  int seqpos = 0;
  for ( int j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
    if ( protein_pose.pdb_info()->number(j) == central_relax_residue_number ) {
      //seqpos_ = j;
      if (chain != ' '){
        if ( protein_pose.pdb_info()->chain(j) == chain ) {
          seqpos = j;
        }
      }else{
        seqpos = j;
      }
    }
  }
  if ( seqpos == 0 ) {
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
    exit(1);
  }

  utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);
	if(angles <1){
		fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
    return -1;
  }else if (angles > 1)
		{
			core::Real best_vol(0), curr_vol(1);
			for (int i=0; i<angles; ++i){
				core::pose::Pose temp_pose;
				temp_pose = protein_pose;
				core::Real x,y,z;
				x = (int) (numeric::random::uniform() *89 +1);
				y = (int) (numeric::random::uniform() *89 +1);
				z = (int) (numeric::random::uniform() *89 +1);
				numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_degrees(x) );
				numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_degrees(y) );
				numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_degrees(z) );
				numeric::xyzMatrix<core::Real> tot_rot_mat = z_rot_mat * y_rot_mat * x_rot_mat;
				core::Vector v(0,0,0);
				temp_pose.apply_transform_Rx_plus_v(tot_rot_mat, v);
				protocols::pockets::PocketGrid	pg( temp_pose.conformation().residue(seqpos) );
				pg.autoexpanding_pocket_eval( temp_pose.conformation().residue(seqpos), temp_pose ) ;
				curr_vol = pg.netTargetPocketVolume();
				std::cout<<"curr_volume "<<curr_vol<<std::endl;
				if(curr_vol > best_vol){
					best_vol = curr_vol;
					original_pocket_angle_transform[1] = x;
					original_pocket_angle_transform[2] = y;
					original_pocket_angle_transform[3] = z;
				}
			}

			numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_degrees( original_pocket_angle_transform[1] ) );
			numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_degrees( original_pocket_angle_transform[2] ) );
			numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_degrees( original_pocket_angle_transform[3] ) );
			numeric::xyzMatrix<core::Real> bestxyz_rot_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
			core::Vector v(0,0,0);
			protein_pose.apply_transform_Rx_plus_v(bestxyz_rot_mat, v);

			core::pose::Pose best_pose;
			best_pose = protein_pose;

			protocols::pockets::PocketGrid	pg( best_pose.conformation().residue(seqpos) );

			pg.autoexpanding_pocket_eval( best_pose.conformation().residue(seqpos), best_pose ) ;
			std::cout<<"best_volume: "<<pg.netTargetPocketVolume()<<std::endl;
			pg.dumpGridToFile();
			npf.setup_from_PocketGrid( best_pose, pg );

		}else if (angles == 1){
		protocols::pockets::PocketGrid	pg( protein_pose.conformation().residue(seqpos) );
		pg.autoexpanding_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose ) ;
		pg.dumpGridToFile();
		npf.setup_from_PocketGrid( protein_pose, pg );
	}

  int dot_index1 = input_protein.rfind(".", input_protein.size());
  assert(dot_index1 != -1 && "No dot found in filename");
	std::string protein_name = input_protein.substr(0,dot_index1);
  int dot_index2 = input_ligand.rfind(".", input_ligand.size());
  assert(dot_index2 != -1 && "No dot found in filename");
	std::string ligand_name = input_ligand.substr(0,dot_index2);
	std::string tag = ligand_name + "_" + protein_name + "_" + resid;


  // CHEAT
  //Calculate CoM of Ligand
	pose::Pose small_mol_pose;
	core::import_pose::pose_from_pdb( small_mol_pose, input_ligand );
	core::pose::Pose original_pose = small_mol_pose;
	numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_degrees(original_pocket_angle_transform[1] ) );
	numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_degrees(original_pocket_angle_transform[2] ) );
	numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_degrees(original_pocket_angle_transform[3] ) );
	numeric::xyzMatrix<core::Real> bestxyz_rot_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
	core::Vector v(0,0,0);
	small_mol_pose.apply_transform_Rx_plus_v(bestxyz_rot_mat, v);

	core::Size lig_res_num = 0;
  for ( int j = 1, resnum = small_mol_pose.total_residue(); j <= resnum; ++j ) {
    if (!small_mol_pose.residue(j).is_protein()){
      lig_res_num = j;
      break;
    }
  }
  if (lig_res_num == 0){
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
    exit(1);
  }
	numeric::xyzVector<core::Real> ligand_CoM(0.);
	conformation::Residue const & curr_rsd = small_mol_pose.conformation().residue(lig_res_num);
  for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
    ligand_CoM.x() += curr_rsd.atom(i).xyz()(1);
    ligand_CoM.y() += curr_rsd.atom(i).xyz()(2);
    ligand_CoM.z() += curr_rsd.atom(i).xyz()(3);
  }
  ligand_CoM /= curr_rsd.nheavyatoms();
	std::cout<<"Ligand CoM: "<<ligand_CoM.x()<<" "<<ligand_CoM.y()<<" "<<ligand_CoM.z()<<std::endl;
	if (option[ cheat ]()){
		npf.CHEAT_CoM( ligand_CoM );
	}

	std::string np_output_filename = "npf_" + tag + ".txt";
	std::string np_output_pdbname = "npf_" + tag + ".pdb";
  npf.print_to_file(np_output_filename);
  npf.print_to_pdb(np_output_pdbname);
	std::string p_output_filename, p_output_pdbname;
  p_output_filename = "pf_" + tag + ".txt";
  p_output_pdbname = "pf_" + tag + ".pdb";

	protocols::pockets::PlaidFingerprint pf( small_mol_pose, npf );

	numeric::xyzVector<core::Real> pocket_CoM = npf.CoM();
	numeric::xyzVector<core::Real> smallmol_CoM = pf.CoM();
	std::cout<<"Pocket CoM: "<<pocket_CoM.x()<<" "<<pocket_CoM.y()<<" "<<pocket_CoM.z()<<std::endl;
	std::cout<<"smallmol CoM: "<<smallmol_CoM.x()<<" "<<smallmol_CoM.y()<<" "<<smallmol_CoM.z()<<std::endl;
	//numeric::xyzVector<core::Real> pocket_origin = npf.origin();
	//numeric::xyzVector<core::Real> smallmol_origin = pf.origin();
	//std::cout<<"Pocket origin: "<<pocket_origin.x()<<" "<<pocket_origin.y()<<" "<<pocket_origin.z()<<std::endl;
	//std::cout<<"smallmol origin: "<<smallmol_origin.x()<<" "<<smallmol_origin.y()<<" "<<smallmol_origin.z()<<std::endl;


  pf.print_to_file(p_output_filename);
  pf.print_to_pdb(p_output_pdbname);

	std::cout << "SCORE : Unalign1  : " <<  pf.fp_compare( npf, missing_pt_wt, steric_wt ) <<std::endl;
	std::string pose_name = "unal1_pose_" + tag + ".pdb";
	std::string fp_name = "unal1_fp_" + tag + ".pdb";
	pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, 0., 0., 0., original_pocket_angle_transform );

	//std::cout << "SCORE : Unalign2  : " <<  pf.fp_compare( npf, high_rho, low_steric ) << " [high_rho, low_steric] "<<std::endl;
	//pose_name = "unal2_pose_" + tag + ".pdb";
	//fp_name = "unal2_fp_" + tag + ".pdb";
	//pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, 0., 0., 0., original_pocket_angle_transform );

	//	core::Real find_optimal_score, optimal_angle1, optimal_angle2, optimal_angle3, afr_optimal_score;
	//find_optimal_score =  pf.search_random_poses( npf, num_pose_search, optimal_angle1, optimal_angle2, optimal_angle3, low_rho, high_steric );
	//std::cout << "SCORE : OPTIMAL   : " << find_optimal_score << " Optimal_angles: " << optimal_angle1 <<" " << optimal_angle2 <<" " << optimal_angle3 <<std::endl;
	//pose_name = "opt_pose_" + tag + ".pdb";
	//fp_name = "opt_fp_" + tag + ".pdb";
	//pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, optimal_angle1, optimal_angle2, optimal_angle3, original_pocket_angle_transform);

	//	numeric::xyzVector<core::Real> init_origin = pocket_CoM;
	//	utility::vector1<core::Real> vars(6);
	//vars[1] = 0.; // init_origin.x();
	//vars[2] = 0.; // init_origin.y();
	//vars[3] = 0.; // init_origin.z();
	//vars[4]=  optimal_angle1; // * numeric::constants::f::pi_over_180;
	//vars[5] = optimal_angle2; // * numeric::constants::f::pi_over_180;
	//vars[6] = optimal_angle3; // * numeric::constants::f::pi_over_180;

	utility::vector1<core::Real> p_min(6);
	p_min[1] = origin_space * -1; // init_origin.x() - 4.;
	p_min[2] = origin_space * -1; // init_origin.x() - 4.;
	p_min[3] = origin_space * -1; // init_origin.x() - 4.;
	p_min[4] = 0.;
	p_min[5] = 0.;
	p_min[6] = 0.;

	utility::vector1<core::Real> p_max(6);
	p_max[1] = origin_space; // init_origin.x() + 4.;
	p_max[2] = origin_space; // init_origin.y() + 4.;
	p_max[3] = origin_space; // init_origin.z() + 4.;
	p_max[4] = 360.; // * numeric::constants::f::pi_over_180;
	p_max[5] = 360.; // * numeric::constants::f::pi_over_180;
	p_max[6] = 360.; // * numeric::constants::f::pi_over_180;


	ParticleOPs particles;
	protocols::pockets::FingerprintMultifunc fpm(npf, pf, missing_pt_wt, steric_wt);
	core::optimization::ParticleSwarmMinimizer pso(p_min, p_max);
	particles = pso.run(run_size, fpm, particle_size);
	//fpm.dump(vars);
	std::string header = "pso-optimization";
	pso.print_particles(particles, header);

	ParticleOP p = particles[1];
	core::optimization::Particle parti(*p);
	core::Real fit_best = -(parti.fitness_pbest());
	utility::vector1<core::Real> best_vars(6);
	best_vars = parti.pbest();
	std::cout<<"BEST FITNESS ["<<tag<<"]: "<<fit_best<<",   "<<"best dofs:  [  "<<best_vars[1]<<",   "<<best_vars[2]<<",   "<<best_vars[3]<<",   "<<best_vars[4]<<",   "<<best_vars[5]<<",   "<<best_vars[6]<<"   ]"<<std::endl;

	numeric::xyzVector<core::Real> optimized_origin(0.);
	optimized_origin.x() = best_vars[1];
	optimized_origin.y() = best_vars[2];
	optimized_origin.z() = best_vars[3];

	pose_name = "pso_pose_" + tag + ".pdb";
  fp_name = "pso_fp_" + tag + ".pdb";
  pf.dump_oriented_pose_and_fp_to_pdb(pose_name, fp_name, npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin );

	core::pose::Pose oriented_pose = pf.get_oriented_pose(npf, best_vars[4], best_vars[5], best_vars[6], original_pocket_angle_transform, optimized_origin );
	core::Real rmsd_value = pf.rmsd(original_pose, oriented_pose);
	std::cout<<"RMSD ["<<tag<<"]: "<<rmsd_value<<std::endl;
	return 0;

}
