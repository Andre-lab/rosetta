// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file    apps/pilot/bder/SymMetalInterface_TwoZN_setup.cc
/// @brief   Stage 1 of designing a two-zinc metal seeded symmetric interface.
/// @details A two-zinc interface is generated in SymMetalInterface_TwoZN_setup.cc by grafting two 2-residue zinc-binding matches onto the surface, duplicating the pose, rotating 180 degrees about the zinc-zinc axes, then rotating 180 degrees about the line that bisects the two zincs.  At this point, there is one significant degree of freedom to exlore WHILE MAINTAINING SYMMETRY.  This DOF is rotation about the zinc-zinc axis.  The rotation is done in a grid-search manner, and any uneclipsed pose with good metal geometry is dumped/written to disk.  The setup is separate from the design to reduce the computational load, and to provide a debugging checkpoint.
/// @author Bryan Der


// AUTO-REMOVED #include <devel/metal_interface/AddMetalSiteConstraints.hh>
// AUTO-REMOVED #include <devel/metal_interface/FindClosestAtom.hh>
#include <devel/metal_interface/MatchGrafter.hh>
#include <devel/metal_interface/MetalSiteResidue.hh>
#include <devel/metal_interface/ParseMetalSite.hh>
#include <devel/init.hh>

#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/moves/RollMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>
#include <numeric/xyzVector.hh>

#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <sstream>

#include <utility/vector0.hh>



typedef core::pose::Pose Pose;
typedef numeric::xyzVector<core::Real> point;
typedef point axis;

using namespace core;
using basic::Warning;

static basic::Tracer TR("apps.pilot.bder.SymMetalInterface_TwoZN_setup");

basic::options::StringOptionKey const match1( "match1" ); // match refers to 2 residues (Cys/His/Asp/Glu) + zinc, found by the matcher
basic::options::StringOptionKey const match2( "match2" );
basic::options::RealOptionKey const angle_rotation_increment( "angle_rotation_increment" ); // the zinc-zinc axis is rotated 360 degrees (this DOF does maintain dimer symmetry), and this option tells the increment(degrees) in the 360 degree search
basic::options::RealOptionKey const ddG_centroid_cutoff( "ddG_centroid_cutoff" ); // used as a metric to decide if there are irredeemable clashes in the current orientation
basic::options::RealOptionKey const zn_zn_distance_cutoff( "zn_zn_distance_cutoff" ); // will throw out cases where the matches are too close together
basic::options::RealOptionKey const tetrahedral_angle_sumsq_cutoff( "tetrahedral_angle_sumsq_cutoff" ); // how strict you want to be when evaluated metal site geometry.  1800 is for a std dev. of 15 degrees (15 x 15 x 4 angles x 2 zinc sites).  Only 4 angles are considered despite there being 6 tetrahedral angles, because the tetrahedral angle within each match remains fixed


///@brief
class SymMetalInterface_TwoZN_setup : public protocols::moves::Mover {
public:

  SymMetalInterface_TwoZN_setup()
    : msr_1_(5, new devel::metal_interface::MetalSiteResidue), msr_2_(5, new devel::metal_interface::MetalSiteResidue)
  {
    core::import_pose::pose_from_pdb( match1_, basic::options::option[match1].value() );
    core::import_pose::pose_from_pdb( match2_, basic::options::option[match2].value() );

		TR << "//////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;

		TR << "/// @file    apps/pilot/bder/SymMetalInterface_TwoZN_setup.cc" << std::endl;
		TR << "/// @brief   Stage 1 of designing a two-zinc metal seeded symmetric interface." << std::endl;
		TR << "/// @details A two-zinc interface is generated in SymMetalInterface_TwoZN_setup.cc by grafting two matches onto the surface, duplicating the pose, rotating 180 degrees about the zinc-zinc axes, then rotating 180 degrees about the line that bisects the two zincs.  At this point, there one significant degree of freedom to exlore WHILE MAINTAINING SYMMETRY.  This DOF is rotation about the zinc-zinc axes.  The rotation is done in a grid-search manner, and any uneclipsed pose with good metal geometry is dumped/written to disk.  The setup was separate from the design to reduce the computational load, and to provide a debugging checkpoint." << std::endl;

		TR << "Options used in this protocol:" << std::endl;
		TR << "  -s scaffold.pdb      // scaffold used during matching" << std::endl;
		TR << "  -match1 match1.pdb   // match refers to 2 residues (Cys/His/Asp/Glu) + zinc, found by the matcher" << std::endl;
		TR << "  -match2 match2.pdb" << std::endl;
		TR << "  -angle_rotation_increment 5.0  // the zinc-zinc axis is rotated 360 degrees (this DOF does maintain dimer symmetry), and this option tells the increment(degrees) in the 360 degree search" << std::endl;
		TR << "  -ddG_centroid_cutoff  0.0              // used as a metric to decide if there are irredeemable clashes in the current orientation" << std::endl;
		TR << "  -zn_zn_distance_cutoff 10.0            // will throw out cases where the matches are too close together" << std::endl;
		TR << "  -tetrahedral_angle_sumsq_cutoff 1800   // how strict you want to be when evaluated metal site geometry.  1800 is for a std dev. of 15 degrees (15 x 15 x 4 angles x 2 zinc sites).  Only 4 angles are considered despite there being 6 tetrahedral angles, because the tetrahedral angle within each match remains fixed" << std::endl;

		TR << "Steps in the protocol:" << std::endl;
		TR <<   "STEP 1: graft two matches onto monomer scaffold" << std::endl;
		TR <<   "STEP 2: create homodimer pose" << std::endl;
		TR <<   "STEP 3: make homodimer pose an inverse C2 symmetric dimer using rollmoves" << std::endl;
		TR <<   "STEP 4: find the zinc binding sites" << std::endl;
		TR <<   "STEP 5: gridsearch rollmove alignments" << std::endl;
		TR <<   " STEP 6: calculate sum of squares of angles about 2 zincs.  If < 1800 A2 (or other value specified), continue..." << std::endl;
		TR <<   " STEP 7: calculate dG_centroid.  If < 0.0, SUCCESS, output PDB to disk (note: no design!)" << std::endl;

		TR << "//////////////////////////////////////////////////////////////////////////////////////////////" << std::endl << std::endl;
  }

  virtual ~SymMetalInterface_TwoZN_setup(){};



  virtual void
  apply( Pose & pose/*monomer*/ ){

    Pose grafted_monomer = setup_grafted_monomer( pose );

		TR << "Duplicating the monomer (the two chains are overlayed, but fold tree is correct)" << std::endl;
    Pose homodimer = setup_homodimer( grafted_monomer );

		TR << "ROLLMOVE TO INVERSE C2 SYMMETRY:  moves one chain into an inverted symmetric alignment" << std::endl;
    rollmove_to_inverse_C2_symmetry( homodimer );

		TR << "SETUP METALSITES:  locates the two zinc binding sites" << std::endl;
    setup_metalsites( homodimer );

		TR << "SETUP FILTER CLASHES:  makes a scorefunction and an InterfaceAnalyzerMover" << std::endl;
    setup_filter_clashes();


		TR << "GRIDSESARCH ALIGNMENTS:  this is the heart of the protocol, rotates one chain about the zinc-zinc axis and outputs a PDB if there are no serious clashes and the zinc geometry is good." << std::endl;
    gridsearch_symmetric_alignments( homodimer );

    return;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// Setup  ///////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////


  //setup
  virtual Pose
  setup_grafted_monomer( Pose & pose ) {

    zinc1_res_ = pose.total_residue() + 1;
    zinc2_res_ = pose.total_residue() + 2;

    devel::metal_interface::MatchGrafterOP match_grafter = new devel::metal_interface::MatchGrafter;
    Pose first_graft_pose = match_grafter->graft( match1_, pose );
    //first_graft_pose.dump_pdb("first_graft_pose.pdb");

    Pose monomer_grafted_twice = match_grafter->graft( match2_, first_graft_pose );
    //monomer_grafted_twice.dump_pdb("monomer_grafted_twice.pdb");

    utility::vector1<Size> metalsite_seqpos_1;
    metalsite_seqpos_1.push_back(zinc1_res_);
    metalsite_seqpos_1.push_back(match1_.pdb_info()->number(1) );
    metalsite_seqpos_1.push_back(match1_.pdb_info()->number(2) );

    utility::vector1<Size> metalsite_seqpos_2;
    metalsite_seqpos_2.push_back(zinc2_res_);
    metalsite_seqpos_2.push_back(match2_.pdb_info()->number(1) );
    metalsite_seqpos_2.push_back(match2_.pdb_info()->number(2) );

    match_grafter->ensure_proper_his_tautomers( monomer_grafted_twice, metalsite_seqpos_1 );
    match_grafter->ensure_proper_his_tautomers( monomer_grafted_twice, metalsite_seqpos_2 );

    TR << "FOLD TREE OF GRAFTED MONOMER: " << monomer_grafted_twice.fold_tree() << std::endl;

    //monomer_grafted_twice.center(); // important for defining rollmove axes later?  can't remember why I included this

    return monomer_grafted_twice;
  }


  //setup
  virtual Pose
  setup_homodimer( Pose & grafted_monomer ) {

		//the homodimer will be the grafted_monomer duplicated and overlayed on itself, a proper dimer is formed later
    Pose homodimer( grafted_monomer );
    homodimer.append_residue_by_jump(grafted_monomer.residue(1), 1);
    for (core::Size i=2; i<=zinc1_res_ - 1; ++i){
      homodimer.append_residue_by_bond(grafted_monomer.residue(i)); // appending only protein residues, not duplicate zinc atoms
    }
    homodimer.conformation().insert_chain_ending( zinc2_res_ );

    TR << "FOLD TREE HOMODIMER: " << homodimer.fold_tree() << std::endl;

		//homodimer.dump_pdb("homodimer.pdb");
    return homodimer;
  }


  //setup
  virtual void
  rollmove_to_inverse_C2_symmetry( Pose & homodimer ) {

    //First rollmove is along the ZnZn axis
    point zinc1 = homodimer.residue(zinc1_res_).atom(1).xyz();
    point zinc2 = homodimer.residue(zinc2_res_).atom(1).xyz();
    axis const ZnZn_axis = zinc1 - zinc2;
    Size last_residue = homodimer.total_residue(); //((zinc1_res_ - 1) * 2) + 2;
    protocols::moves::RollMoverOP ZnZn_axis_rollmover = new protocols::moves::RollMover( zinc2_res_+1, last_residue, 180, 180, ZnZn_axis, zinc1 );
    ZnZn_axis_rollmover->apply( homodimer );
    //homodimer.dump_pdb("homodimer_ZnZnaxis_rollmove.pdb");


    //Second rollmove is orthogonal to the ZnZn axes, through the ZnZn midpoint
    //note: center of mass of monomer is at the origin
    point ZnZn_midpoint = midpoint(zinc1, zinc2); //#include<numeric/xyzVector.hh>
    axis axis_normal_to_origZn1_origZn2 = cross_product(zinc1, zinc2);
    axis orthogonal_axis = cross_product(zinc1 - zinc2, axis_normal_to_origZn1_origZn2);
    protocols::moves::RollMoverOP orthogonal_axis_rollmover = new protocols::moves::RollMover( zinc2_res_+1, last_residue, 180, 180, orthogonal_axis, ZnZn_midpoint );
    orthogonal_axis_rollmover->apply( homodimer );
    //homodimer.dump_pdb("homodimer_orthaxis_rollmove.pdb");
  }



  //setup
  virtual void
  setup_metalsites( Pose & homodimer ) {
    devel::metal_interface::ParseMetalSiteOP metalsite_parser_1 = new devel::metal_interface::ParseMetalSite( zinc1_res_ );
    devel::metal_interface::ParseMetalSiteOP metalsite_parser_2 = new devel::metal_interface::ParseMetalSite( zinc2_res_ );

    msr_1_ = metalsite_parser_1->parse_metalsite( homodimer );
    msr_2_ = metalsite_parser_2->parse_metalsite( homodimer );

    TR << "zinc1 " << msr_1_[1]->get_seqpos() << std::endl;
    TR << "1_1:  " << msr_1_[2]->get_seqpos() << std::endl;
    TR << "1_2:  " << msr_1_[3]->get_seqpos() << std::endl;
    TR << "1_3:  " << msr_1_[4]->get_seqpos() << std::endl;
    TR << "1_4:  " << msr_1_[5]->get_seqpos() << std::endl;

    TR << "zinc2 " << msr_2_[1]->get_seqpos() << std::endl;
    TR << "2_1:  " << msr_2_[2]->get_seqpos() << std::endl;
    TR << "2_2:  " << msr_2_[3]->get_seqpos() << std::endl;
    TR << "2_3:  " << msr_2_[4]->get_seqpos() << std::endl;
    TR << "2_4:  " << msr_2_[5]->get_seqpos() << std::endl;

  }

  //setup
  virtual void
  setup_filter_clashes() {
    using namespace core::scoring;
    scorefunction_ = core::scoring::ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );
    interface_analyzer_ = new protocols::analysis::InterfaceAnalyzerMover( 3, false, scorefunction_ ); // zinc1 and zinc2 get jumps 1 and 2, chain B gets jump 3
  }



  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// End of Setup  ////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////





  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// Protocol  ////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////




  //protocol
  virtual void
  gridsearch_symmetric_alignments( Pose & homodimer ) {

    //First rollmove is along the ZnZn axis
    point zinc1 = homodimer.residue(zinc1_res_).atom(1).xyz();
    point zinc2 = homodimer.residue(zinc2_res_).atom(1).xyz();

    axis const ZnZn_axis = zinc1 - zinc2;
    Size last_residue = homodimer.total_residue(); //((zinc1_res_ - 1) * 2) + 2;
    Real grid_angle = basic::options::option[angle_rotation_increment].value();

    protocols::moves::RollMoverOP ZnZn_axis_rollmover = new protocols::moves::RollMover( zinc2_res_+1, last_residue, grid_angle, grid_angle, ZnZn_axis, zinc1 );

    for( Size i(1); i <= 360.0/grid_angle; i++ ) {

      ZnZn_axis_rollmover->apply( homodimer );

      Size angle = core::Size(i*grid_angle);
      std::stringstream ss_angle;
      ss_angle << angle;
      std::string name = "homodimer_" + ss_angle.str() + ".pdb";
      //homodimer.dump_pdb( name );

			//FILTER GEOMETRY
      filter_metal_geom( homodimer, ss_angle.str() );
    }

  }


  //protocol
  virtual void
  filter_metal_geom( Pose & homodimer, std::string angle_tag ) {

    utility::vector1< core::Real > angles( 8, 0 );
    utility::vector1< core::Real > anglediffs( 8, 0 );
    core::Real const tetrahedral( 109.5 );
    core::Real sumsq_angle( 0 );

    point zinc1 = msr_1_[1]->get_ligand_atom_xyz();
    point zinc2 = msr_2_[1]->get_ligand_atom_xyz();

    point p1_1 = homodimer.residue( msr_1_[2]->get_seqpos() ).atom( msr_1_[2]->get_ligand_atom_name() ).xyz();
    point p1_2 = homodimer.residue( msr_1_[3]->get_seqpos() ).atom( msr_1_[3]->get_ligand_atom_name() ).xyz();
    point p1_3 = homodimer.residue( msr_1_[4]->get_seqpos() ).atom( msr_1_[4]->get_ligand_atom_name() ).xyz();
    point p1_4 = homodimer.residue( msr_1_[5]->get_seqpos() ).atom( msr_1_[5]->get_ligand_atom_name() ).xyz();

    point p2_1 = homodimer.residue( msr_2_[2]->get_seqpos() ).atom( msr_2_[2]->get_ligand_atom_name() ).xyz();
    point p2_2 = homodimer.residue( msr_2_[3]->get_seqpos() ).atom( msr_2_[3]->get_ligand_atom_name() ).xyz();
    point p2_3 = homodimer.residue( msr_2_[4]->get_seqpos() ).atom( msr_2_[4]->get_ligand_atom_name() ).xyz();
    point p2_4 = homodimer.residue( msr_2_[5]->get_seqpos() ).atom( msr_2_[5]->get_ligand_atom_name() ).xyz();


		TR << "seqpos p1_1 " << msr_1_[2]->get_seqpos() << " " << msr_1_[2]->get_ligand_atom_name() << std::endl;
		TR << "seqpos p1_2 " << msr_1_[3]->get_seqpos() << " " << msr_1_[3]->get_ligand_atom_name() << std::endl;
		TR << "seqpos p1_3 " << msr_1_[4]->get_seqpos() << " " << msr_1_[4]->get_ligand_atom_name() << std::endl;
		TR << "seqpos p1_4 " << msr_1_[5]->get_seqpos() << " " << msr_1_[5]->get_ligand_atom_name() << std::endl;

		TR << "seqpos p2_1 " << msr_2_[2]->get_seqpos() << " " << msr_2_[2]->get_ligand_atom_name() << std::endl;
		TR << "seqpos p2_2 " << msr_2_[3]->get_seqpos() << " " << msr_2_[3]->get_ligand_atom_name() << std::endl;
		TR << "seqpos p2_3 " << msr_2_[4]->get_seqpos() << " " << msr_2_[4]->get_ligand_atom_name() << std::endl;
		TR << "seqpos p2_4 " << msr_2_[5]->get_seqpos() << " " << msr_2_[5]->get_ligand_atom_name() << std::endl;

    // angle_of is a function in numeric::xyzVector
    using numeric::conversions::degrees;
    angles[1] = degrees(angle_of(p1_1, zinc1, p1_3));
    angles[2] = degrees(angle_of(p1_1, zinc1, p1_4));
    angles[3] = degrees(angle_of(p1_2, zinc1, p1_3));
    angles[4] = degrees(angle_of(p1_2, zinc1, p1_4));

    angles[5] = degrees(angle_of(p2_1, zinc2, p2_3));
    angles[6] = degrees(angle_of(p2_1, zinc2, p2_4));
    angles[7] = degrees(angle_of(p2_2, zinc2, p2_3));
    angles[8] = degrees(angle_of(p2_2, zinc2, p2_4));


    for( core::Size j(1); j<= angles.size(); ++j){
      anglediffs[j] = tetrahedral - angles[j];
      sumsq_angle += anglediffs[j] * anglediffs[j];
      TR << "tag " << angle_tag << "  angle " << j << " " << angles[j] << " diff from 109.5 " << anglediffs[j] << std::endl;
    }
    TR << "sumsq_angle: " << sumsq_angle << std::endl;

    if( sumsq_angle < basic::options::option[tetrahedral_angle_sumsq_cutoff].value() ) { // 2*4*(15*15) = 1800
      TR << "Match geometries made the cutoff" << std::endl;

			//FILTER CLASHES
      filter_clashes( homodimer, angle_tag );
    }
  }


  //protocol
  virtual void
  filter_clashes( Pose & homodimer, std::string angle_tag ) {

		TR << "9" << std::endl;
    //Filter Zn-Zn distance
    point zinc1 = msr_1_[1]->get_ligand_atom_xyz();
    point zinc2 = msr_2_[1]->get_ligand_atom_xyz();
    if( zinc1.distance(zinc2) < basic::options::option[zn_zn_distance_cutoff].value() ) {
      TR << "REJECTING, Zincs are too close!" << std::endl;
      return;
    }
    TR << "Zincs are not too close!  Checking centroid ddG..." << std::endl;

    //Filter interchain clashes
    interface_analyzer_->set_use_centroid_dG( true );
    interface_analyzer_->apply( homodimer );
    Real ddG_centroid = interface_analyzer_->get_centroid_dG();
    TR << "ddG_centroid: " << ddG_centroid << std::endl;

    if (ddG_centroid < basic::options::option[ddG_centroid_cutoff].value() ) {
      TR << "SUCCESS, interface clash detection also made the cutoff" << std::endl;
      utility::file::FileName match1_name = match1_.pdb_info()->name();
      utility::file::FileName match2_name = match2_.pdb_info()->name();
      std::string name = match1_name.base() + "." + match2_name.base() + "_" + angle_tag + ".pdb";

			// WE HAVE FOUND A GOOD DIMER STARTING STRUCTURE
      homodimer.dump_pdb( name );
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////// End of Protocol  //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////


	virtual
	std::string
	get_name() const { return "SymMetalInterface_TwoZN_setup"; }


private:
  //options
  Pose match1_;
  Pose match2_;
  Pose homodimer_with_matches_;
  Size symmetric_packrot_iters_;

  Size zinc1_res_;
  Size zinc2_res_;

  utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr_1_;
  utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr_2_;

  protocols::analysis::InterfaceAnalyzerMoverOP interface_analyzer_;
  core::scoring::ScoreFunctionOP scorefunction_;

};

typedef utility::pointer::owning_ptr< SymMetalInterface_TwoZN_setup > SymMetalInterface_TwoZN_setupOP;


int
main( int argc, char* argv[] )
{
  using basic::options::option;
  option.add( match1, "match1" ).def("match1.pdb");
  option.add( match2, "match2" ).def("match2.pdb");
  option.add( angle_rotation_increment, "angle_rotation_increment" ).def(10.0);
  option.add( ddG_centroid_cutoff, "ddG_centroid_cutoff" ).def(10.0);
  option.add( zn_zn_distance_cutoff, "zn_zn_distance_cutoff" ).def(10.0);
  option.add( tetrahedral_angle_sumsq_cutoff, "tetrahedral_angle_sumsq_cutoff" ).def(1800);

  devel::init(argc, argv);
  protocols::jd2::JobDistributor::get_instance()->go( new SymMetalInterface_TwoZN_setup() );

  TR << "************************d**o**n**e**************************************" << std::endl;
  return 0;
}

