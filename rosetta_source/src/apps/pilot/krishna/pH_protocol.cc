// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file src/apps/pilot/krishna/pH_protocol.cc
/// @Calculates the pKa value for a residue

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/pHEnergy.hh>

// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <core/pose/util.hh>

#include <basic/options/option.hh>
#include <core/init.hh>
#include <core/conformation/Residue.hh>

#include <protocols/moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/SidechainMover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

#include <basic/Tracer.hh>
#include <string>

//temp includes krishna
#include <iostream>
//#include <core/io/pdb/pose_io.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/jd2/JobOutputter.hh>
// AUTO-REMOVED #include <protocols/jd2/Job.hh>


// option key includes

#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



static basic::Tracer TR("apps.pilot.krishna.PhProtocol");

class PhProtocol : public protocols::moves::Mover {

public:
  PhProtocol()
  {
    // curr_pH_ = 0.1; //start value for pH
  }

  virtual ~PhProtocol(){};

  virtual
  void
  apply( core::pose::Pose & pose ){

    using namespace core;
    using namespace basic::options;
    using namespace chemical;
    using namespace conformation;

    Real shift_pH = 1.0, pka_value = 1.0, ipka = 1.0;

    //setting up residue details
    std::string pdb_file_name = option[ OptionKeys::in::file::s ]()[1];
    int pdb_res_no = option[ OptionKeys::pH::calc_pka::pka_for_resno]();

    if ( pdb_res_no < 0 ){
      TR << "PLEASE ENTER RESIDUE NO FOR WHICH PKA IS TO BE CALCULATED AND TRY AGAIN" << std::endl;
      return;
    }

    std::string pdb_chain_no = option[ OptionKeys::pH::calc_pka::pka_for_chainno]();
    Size res_no = pose.pdb_info()->pdb2pose( pdb_chain_no[0], pdb_res_no );
    std::string res_name = pose.residue( res_no ).name();
    std::string res_name3 = pose.residue( res_no ).name3();

    //ipKa details
    switch ( pose.residue( res_no ).type().aa() ){
      case aa_asp: ipka = 4.0; break;
      case aa_glu: ipka = 4.4; break;
      case aa_his: ipka = 6.3; break;
      case aa_lys: ipka = 10.4; break;
      case aa_tyr: ipka = 10.0; break;
      default: ipka = 0.0; break;
    }

    TR << " RES NO BEING PACKED IS " << pdb_res_no << " FROM CHAIN " << pdb_chain_no << std::endl;

    pose::Pose old_pose = pose;  //initialize the pose

    for ( Real curr_pH = 0.0; curr_pH <= 12.0; curr_pH += 1.0 ){

      pose::Pose curr_pose = pose;

      TR << " CURRENTLY SIMULATING AT pH " << curr_pH << std::endl;

      pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( curr_pose ));

      for ( Size ii = 1; ii <= curr_pose.total_residue(); ++ii ) {
        if ( ii == res_no )
          task->nonconst_residue_task( ii ).restrict_to_repacking();
        else
          task->nonconst_residue_task( ii ).prevent_repacking();
      }

//      std::cout << *task << std::endl;

      scoring::ScoreFunctionOP score_fxn( scoring::getScoreFunction() );
      scoring::methods::pHEnergy::set_pH ( curr_pH );

      protocols::moves::PackRotamersMoverOP pack_mover( new protocols::moves::PackRotamersMover( score_fxn, task ) );

      pack_mover->apply(curr_pose);

      (*score_fxn)(curr_pose);

      if ( ( old_pose.residue( res_no ).name() != curr_pose.residue( res_no ).name() ) && ( curr_pH != 0.0 ) )
        break;

      old_pose = curr_pose;
      shift_pH = curr_pH;

    }

    for ( Real curr_pH = shift_pH; curr_pH < shift_pH + 1.0; curr_pH += 0.1 ){

      pose::Pose curr_pose = old_pose; //look at smaller intervals of pH taking the previous pose

      pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( curr_pose ));

      for ( Size ii = 1; ii <= curr_pose.total_residue(); ++ii ) {
        if ( ii == res_no )
          task->nonconst_residue_task( ii ).restrict_to_repacking();
        else
          task->nonconst_residue_task( ii ).prevent_repacking();
      }

      scoring::ScoreFunctionOP score_fxn( scoring::getScoreFunction() );
      scoring::methods::pHEnergy::set_pH ( curr_pH );

      protocols::moves::PackRotamersMoverOP pack_mover( new protocols::moves::PackRotamersMover( score_fxn, task ) );

      TR << " CURRENTLY SIMULATING AT pH " << curr_pH << std::endl;

      pack_mover->apply(curr_pose);

      (*score_fxn)(curr_pose);

      if ( old_pose.residue( res_no ).name() != curr_pose.residue( res_no ).name() ){
//        io::pdb::dump_pdb ( curr_pose, "LYS.pdb" );
        break;
      }

      old_pose = curr_pose;
      pka_value = curr_pH;

    }

//    TR << " THE PREDICTED PKA VALUE FOR " << res_name << pdb_res_no << " FROM CHAIN " << pdb_chain_no << " FROM " << pdb_file_name << " IS " << pka_value << std::endl;

    TR << "PKA FOR\t" << res_name3 << "\t" << pdb_res_no << "\t" << pdb_file_name << "\t" << ipka << "\t" << pka_value << std::endl;

  }

	std::string get_name() const { return "PhProtocol"; }

  virtual
  protocols::moves::MoverOP
  fresh_instance() const {
    return new PhProtocol;
  }

  virtual
  bool
  reinitialize_for_each_job() const { return false; }

  virtual
  bool
  reinitialize_for_new_input() const { return false; }

private:

  //core::Real curr_pH_;

};

typedef utility::pointer::owning_ptr< PhProtocol > PhProtocolOP;


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
  core::init(argc, argv);

  PhProtocolOP pH_test(new PhProtocol);

  protocols::jd2::JobDistributor::get_instance()->go(pH_test);

}





