// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

//GPU enabling is not default
//To test how many threads are fastest for your computer,
//use -gpu:threads 1024 (or other number) on the command line

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

// Protocol Headers
#include <devel/init.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <utility/options/StringOption.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( String, central_relax_pdb_num )
OPT_KEY( String, protein )
OPT_KEY( Integer, num_angles )
OPT_KEY( Boolean, cheat )
OPT_KEY( Boolean, trim_pocket )
OPT_KEY( Boolean, adt_grid )
OPT_KEY( Boolean, lig_grid )
OPT_KEY( Boolean, res_grid )
OPT_KEY( Boolean, print_output_complex )
OPT_KEY( Real, gc_x )
OPT_KEY( Real, gc_y )
OPT_KEY( Real, gc_z )
OPT_KEY( Real, gd_x )
OPT_KEY( Real, gd_y )
OPT_KEY( Real, gd_z )
OPT_KEY( Real, gs )
OPT_KEY( String, bound_ligand )
OPT_KEY( Integer, add_grid_size )
OPT_KEY( Integer, add_ext_grid_size )

int main( int argc, char * argv [] ) {

  NEW_OPT( central_relax_pdb_num, "target residue", "-1" );
  NEW_OPT( protein, "protein file name", "protein.pdb" );
  NEW_OPT( num_angles, "no. of angles for rotating the grid", 1 );
  NEW_OPT( cheat, "move pocket CoM over Ligand CoM", false );
  NEW_OPT( trim_pocket, "trim the non-plaid pocket using a known ligand", false );
  NEW_OPT( adt_grid, "resize grid based on user entered AUTODOCK grid values", false );
  NEW_OPT( lig_grid, "resize grid based on bound ligand", false );
  NEW_OPT( res_grid, "resize grid based on a target residue", true );
  NEW_OPT( gc_x, "gid center : X ", 1.0 );
  NEW_OPT( gc_y, "gid center : Y ", 1.0 );
  NEW_OPT( gc_z, "gid center : Z ", 1.0 );
  NEW_OPT( gd_x, "gid dimension : X ", 20.0 );
  NEW_OPT( gd_y, "gid dimension : Y ", 20.0 );
  NEW_OPT( gd_z, "gid dimension : Z ", 20.0 );
  NEW_OPT( bound_ligand, "use bound ligand to set the grid for generating eggshell", "" );
  NEW_OPT( add_grid_size, "add grid dimension along x,y,z axis", 2 );
  NEW_OPT( add_ext_grid_size, "add extra grid dimension along x,y,z axis", 2 );

  devel::init(argc, argv);

  std::string const resid = option[ central_relax_pdb_num ];
  std::string const input_protein = option[ protein ];
  int angles = option[ num_angles ];
  core::Real const grid_cen_x = option[ gc_x ];
  core::Real const grid_cen_y = option[ gc_y ];
  core::Real const grid_cen_z = option[ gc_z ];
  core::Real grid_dim_x = option[ gd_x ];
  core::Real grid_dim_y = option[ gd_y ];
  core::Real grid_dim_z = option[ gd_z ];
  std::string const bound_ligand_file = option[ bound_ligand ];
  int add_grid_dim  = option[ add_grid_size ];
  int add_ext_grid_dim  = option[ add_ext_grid_size ];

  using namespace basic::options;
  core::Real const spacing = option[ OptionKeys::pocket_grid::pocket_grid_spacing ]();

  protocols::pockets::NonPlaidFingerprint npf;
  pose::Pose protein_pose;
  core::import_pose::pose_from_pdb( protein_pose, input_protein );
  utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);

  //find sequence position for the target residue
  int seqpos = 0;
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
  for ( int j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
    if ( protein_pose.pdb_info()->number(j) == central_relax_residue_number ) {
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

  //use autodock grid values for grid size and setup PocketGrid
  if (option[ adt_grid ]()){
    grid_dim_x = (grid_dim_x/2) * spacing;
    grid_dim_y = (grid_dim_y/2) * spacing;
    grid_dim_z = (grid_dim_z/2) * spacing;
    protocols::pockets::PocketGrid	pg( grid_cen_x, grid_cen_y, grid_cen_z, grid_dim_x, grid_dim_y, grid_dim_z );
    numeric::xyzVector<core::Real> grid_center (0.);
    grid_center.x() = grid_cen_x;
    grid_center.y() = grid_cen_y;
    grid_center.z() = grid_cen_z;
    pg.DARC_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose, grid_center ) ;
    npf.setup_from_PocketGrid( protein_pose, pg );
  }

  //use bound ligand to get grid size and setup PocketGrid
  else if (option[ lig_grid ]()){
    if (!bound_ligand_file.empty()){
      std::cout<<"Error, no ligand available for setting the grid" << std::endl;
      exit(1);
    }
    pose::Pose bound_ligand_pose;
    core::import_pose::pose_from_pdb( bound_ligand_pose, bound_ligand_file );
    core::Size lig_res_num = 0;
    for ( int j = 1, resnum = bound_ligand_pose.total_residue(); j <= resnum; ++j ) {
      if (!bound_ligand_pose.residue(j).is_protein()){
       	lig_res_num = j;
	break;
      }
    }
    if (lig_res_num == 0){
      std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
      exit(1);
    }
    numeric::xyzVector<core::Real> input_ligand_CoM(0.);
    conformation::Residue const & curr_rsd = bound_ligand_pose.conformation().residue(lig_res_num);
    core::Real minx(999.), miny(999.), minz(999.), maxx(-999.), maxy(-999.), maxz(-999.);
    for(Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i) {
      if (curr_rsd.atom(i).xyz()(1) > maxx){maxx = curr_rsd.atom(i).xyz()(1);}
      if (curr_rsd.atom(i).xyz()(1) < minx){minx = curr_rsd.atom(i).xyz()(1);}
      if (curr_rsd.atom(i).xyz()(2) > maxy){maxy = curr_rsd.atom(i).xyz()(2);}
      if (curr_rsd.atom(i).xyz()(2) < miny){miny = curr_rsd.atom(i).xyz()(2);}
      if (curr_rsd.atom(i).xyz()(3) > maxz){maxz = curr_rsd.atom(i).xyz()(3);}
      if (curr_rsd.atom(i).xyz()(3) < minz){minz = curr_rsd.atom(i).xyz()(3);}
    }
    core::Real x_from_grd_cen, y_from_grd_cen, z_from_grd_cen, x_halfwidth, y_halfwidth, z_halfwidth;
    core::Real const cen_x = (maxx + minx)/2;
    core::Real const cen_y = (maxy + miny)/2;
    core::Real const cen_z = (maxz + minz)/2;
    x_halfwidth = std::abs(maxx - minx)/2;
    y_halfwidth = std::abs(maxy - miny)/2;
    z_halfwidth = std::abs(maxz - minz)/2;
    numeric::xyzVector<core::Real> const grid_center(cen_x,cen_y,cen_z);
    x_from_grd_cen = x_halfwidth + add_grid_dim;
    y_from_grd_cen = y_halfwidth + add_grid_dim;
    z_from_grd_cen = z_halfwidth + add_grid_dim;

    protocols::pockets::PocketGrid	pg( cen_x, cen_y, cen_z, x_from_grd_cen, y_from_grd_cen, z_from_grd_cen );
    pg.DARC_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose, grid_center ) ;
    x_from_grd_cen = x_halfwidth + add_grid_dim + add_ext_grid_dim;
    y_from_grd_cen = y_halfwidth + add_grid_dim + add_ext_grid_dim;
    z_from_grd_cen = z_halfwidth + add_grid_dim + add_ext_grid_dim;
    protocols::pockets::PocketGrid	ext_grd( cen_x, cen_y, cen_z, x_from_grd_cen, y_from_grd_cen, z_from_grd_cen );
    ext_grd.DARC_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose, grid_center ) ;
    npf.setup_from_PocketGrid( protein_pose, pg, ext_grd );
  }

  //use default grid size centered around target residue to setup pocketgrid
  else if (option[ res_grid ]()){
    if (angles <1){
      //num grid rotation cant be < 1
      fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
      exit(1);
    }
    //rotate grid to choose best pocket volume to choose largest pocket
    else if (angles > 1){
      core::Real best_vol(0), curr_vol(1);
      for (int i=0; i<angles; ++i){
       	core::pose::Pose temp_pose;
	temp_pose = protein_pose;
       	core::Real x = ( numeric::random::uniform() * numeric::constants::r::pi_2 ) + 0.0001;
       	core::Real y = ( numeric::random::uniform() * numeric::constants::r::pi_2 ) + 0.0001;
       	core::Real z = ( numeric::random::uniform() * numeric::constants::r::pi_2 ) + 0.0001;
       	numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_radians(x) );
       	numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_radians(y) );
       	numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_radians(z) );
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
      numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_radians( original_pocket_angle_transform[1] ) );
      numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_radians( original_pocket_angle_transform[2] ) );
      numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_radians( original_pocket_angle_transform[3] ) );
      numeric::xyzMatrix<core::Real> bestxyz_rot_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
      core::Vector v(0,0,0);
      protein_pose.apply_transform_Rx_plus_v(bestxyz_rot_mat, v);
      core::pose::Pose best_pose;
      best_pose = protein_pose;
      protocols::pockets::PocketGrid	pg( best_pose.conformation().residue(seqpos) );
      pg.autoexpanding_pocket_eval( best_pose.conformation().residue(seqpos), best_pose ) ;
      std::cout<<"best_volume: "<<pg.netTargetPocketVolume()<<std::endl;
      npf.setup_from_PocketGrid( best_pose, pg );
    }
    //no grid rotation
    else if (angles == 1){
      protocols::pockets::PocketGrid	pg( protein_pose.conformation().residue(seqpos) );
      pg.autoexpanding_pocket_eval( protein_pose.conformation().residue(seqpos), protein_pose ) ;
      npf.setup_from_PocketGrid( protein_pose, pg );
    }
  }

  //trim pocket based on bound ligand
  if (option[ trim_pocket ]()){
    //calc lig_COM and move pock_COM to lig_com of bound ligand
    pose::Pose bound_ligand_pose;
    core::import_pose::pose_from_pdb( bound_ligand_pose, bound_ligand_file );
    core::Size lig_res_num = 0;
    for ( int j = 1, resnum = bound_ligand_pose.total_residue(); j <= resnum; ++j ) {
      if (!bound_ligand_pose.residue(j).is_protein()){
        lig_res_num = j;
       	break;
      }
    }
    if (lig_res_num == 0){
      std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
      exit(1);
    }
    numeric::xyzVector<core::Real> bound_ligand_CoM(0.);
    conformation::Residue const & curr_rsd = bound_ligand_pose.conformation().residue(lig_res_num);
    for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
      bound_ligand_CoM.x() += curr_rsd.atom(i).xyz()(1);
      bound_ligand_CoM.y() += curr_rsd.atom(i).xyz()(2);
      bound_ligand_CoM.z() += curr_rsd.atom(i).xyz()(3);
    }
    bound_ligand_CoM /= curr_rsd.nheavyatoms();
    npf.CHEAT_CoM( bound_ligand_CoM );
    npf.trim_based_on_known_ligand(bound_ligand_pose);
  }

  //create 'tag' for eggshell output filename
  int dot_index1 = input_protein.rfind(".", input_protein.size());
  assert(dot_index1 != -1 && "No dot found in filename");
  std::string protein_name = input_protein.substr(0,dot_index1);
  std::string eggshell_pdb_tag = "eggshell_" + protein_name + "_" + resid + ".pdb";
  std::string eggshell_triplet_tag = "eggshell_" + protein_name + "_" + resid + ".txt";

  //print the eggshell (ray) files
  npf.write_eggshell_to_pdb_file(eggshell_pdb_tag);
  npf.print_to_file(eggshell_triplet_tag);
  std::cout<< "Written eggshell to pdb file : "<< eggshell_pdb_tag << std::endl;
  std::cout<< "Written eggshell to triplet file: "<< eggshell_triplet_tag << std::endl;
  std::cout<< "DONE!"<< std::endl;

  return 0;
}
