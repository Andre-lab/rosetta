
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief


// libRosetta headers

#include <protocols/viewer/viewers.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/sequence/Sequence.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <protocols/rigid/RigidBodyMover.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.fwd.hh>
// AUTO-REMOVED #include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/loops/loops_main.hh>

#include <basic/options/util.hh>//option.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/scoring/rms_util.hh>

// packing
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
// AUTO-REMOVED #include <core/pack/rotamer_trials.hh>
// AUTO-REMOVED #include <protocols/simple_moves/PackRotamersMover.fwd.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED

// Auto-header: duplicate removed #include <core/sequence/util.hh>
// Auto-header: duplicate removed #include <core/sequence/Sequence.hh>
// Auto-header: duplicate removed #include <core/id/SequenceMapping.hh>

// AUTO-REMOVED #include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/BBTorsionSRFD.hh>
// AUTO-REMOVED #include <core/fragment/util.hh>
// AUTO-REMOVED #include <core/fragment/FragmentIO.hh>

// AUTO-REMOVED #include <basic/prof.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/jobdist/Jobs.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/mpistream.hh>

//#include <apps/benchmark/init_util.hh>  //for core_init_with_additional_options
////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace id;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace optimization;
namespace OK = OptionKeys;
using namespace basic::options::OptionKeys;

using utility::vector1;
using std::string;
using core::import_pose::pose_from_pdb;

basic::Tracer TR( "apps.pilot.lin.symmsilentfile" );


///////////////////////////////////////////////////////////////////////////////
void
SymmSilentFileTest()
{
  using namespace conformation::symmetry;
  using namespace chemical;

  // additional flags needed
  //core_init_with_additional_options( "-in::file::silent_struct_type binary -in::file::fullatom" ); //in
  //core_init_with_additional_options( "-out:file:silent_struct_type binary -out::file::fullatom" ); //out
  //core_init_with_additional_options( "-out:file:output_virtual" ); // debug purpose

  Pose pdb_pose;
  core::import_pose::pose_from_pdb( pdb_pose, start_file() );

  if ( option[ OptionKeys::in::file::psipred_ss2 ].user() ) {
    protocols::loops::set_secstruct_from_psipred_ss2(pdb_pose);
   }

  // set up the symmetric pose
  Pose start_pose = pdb_pose;
  make_symmetric_pose( start_pose );

  ScoreFunctionOP scorefxn3_sym( ScoreFunctionFactory::create_score_function( "score3" ) );
  ScoreFunctionOP scorefxn13_sym( ScoreFunctionFactory::create_score_function( "score13" ) );
  scorefxn3_sym = new core::scoring::symmetry::SymmetricScoreFunction(*scorefxn3_sym);
  scorefxn13_sym = new core::scoring::symmetry::SymmetricScoreFunction(*scorefxn13_sym);
  TR << "start_pose Secondary Structure Set " << start_pose.secstruct() << std::endl;
  TR << "start_pose "<< symmetry_info( start_pose ) << std::endl;
  scorefxn3_sym->show( std::cout, start_pose );
  scorefxn13_sym->show( std::cout, start_pose );
  start_pose.dump_pdb("start.pdb");

  //test: save and restore ProteinSilentStruct
  {
    double rms_threshold = 1e-3;
    pose::Pose restored_pose;
    ResidueTypeSet const & rsd_set( start_pose.residue(1).residue_type_set() );

    // Write out the silent-file
    core::io::silent::SilentFileData sfd;
    std::string silent_outfile = "Silent.out";
    utility::file::file_delete( silent_outfile );
    core::io::silent::BinaryProteinSilentStruct pss,new_pss;
    pss.fill_struct( start_pose, "start_pose" );
    sfd.write_silent_struct( pss, silent_outfile );

    // Read the SilentStruct from the silent-file
    utility::vector1 < std::string > tags;
    tags.push_back( "start_pose" );
    sfd.read_file( silent_outfile, tags );
    core::io::silent::SilentFileData::iterator iter = sfd.begin();
    iter->fill_pose( restored_pose, rsd_set );
 	  TR << "restore_pose Secondary Structure Set " << restored_pose.secstruct() << std::endl;
    TR << "restore_pose "<< symmetry_info( restored_pose ) << std::endl;
    scorefxn3_sym->show( std::cout, restored_pose );
	  scorefxn13_sym->show( std::cout, restored_pose );

    new_pss.fill_struct( restored_pose, "restored_pose" );
    sfd.write_silent_struct( new_pss, silent_outfile );
    restored_pose.dump_pdb("restored.pdb");
    Real rms_to_restored = scoring::CA_rmsd( start_pose, restored_pose );
    TR << "RMS error from save/restore: " << rms_to_restored << std::endl;
    assert( rms_to_restored < rms_threshold );
    //utility::file::file_delete( silent_outfile );
  }

}


///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{
  using namespace basic::options;

  SymmSilentFileTest();

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

  // initialize option and random number system
  devel::init( argc, argv );

  protocols::viewer::viewer_main( my_main );
}
