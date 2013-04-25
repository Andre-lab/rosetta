
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>

#include <boost/cstdint.hpp>
#include <boost/unordered_map.hpp>
#include <core/fragment/picking/VallChunk.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/fragment/picking/VallProvider.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <basic/options/option.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <devel/init.hh>
#include <numeric/HomogeneousTransform.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/match/Hit.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>

// C++ headers
//#include <cstdlib>

#include <iostream>
#include <string>
#include <cstdio>

// option key includes
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/util/SwitchResidueTypeSet.hh>



static basic::Tracer TR("main");

using namespace protocols::moves;
using namespace core::scoring;
using namespace core;
using namespace core::pose;
using namespace conformation;
using namespace kinematics;
using namespace protocols::match;
using namespace core::fragment::picking;


using namespace protocols::loophash;


class LoopHashRelax_Sampler;
typedef utility::pointer::owning_ptr< LoopHashRelax_Sampler > LoopHashRelax_SamplerOP;
typedef utility::pointer::owning_ptr< LoopHashRelax_Sampler const > LoopHashRelax_SamplerCOP;

class LoopHashRelax_Sampler: public protocols::moves::Mover {
public:

  LoopHashRelax_Sampler(
    LoopHashLibraryOP library
  ):
   library_(library)

  {
  }

	virtual void apply( core::pose::Pose& pose );

  virtual protocols::moves::MoverOP clone() const {
		return new LoopHashRelax_Sampler( *this );
	}


	virtual std::string get_name() const {
		return "LoopHashRelax_Sampler";
	}

	virtual	protocols::moves::MoverOP	fresh_instance() const {
		return new LoopHashRelax_Sampler( library_ );
	}

private:
  LoopHashLibraryOP library_;

};

void
LoopHashRelax_Sampler::apply( core::pose::Pose& pose )
{
  if( !library_ ) return;

  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  std::string prefix = option[ out::prefix ]();
  core::Size skim_size = option[ lh::skim_size ]();

  LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );

	core::io::silent::SilentStructOP last_best;

  for(int round = 1; round <= option[ OptionKeys::lh::rounds ]; round ++ ){

		LoopHashSampler  lsampler( library_, simple_inserter );

		core::Size random_start = 1 + rand() % (pose.total_residue() - 10 );
		core::Size random_end = random_start = std::min( int(random_start + 4), int(pose.total_residue() - 5) );
		lsampler.set_start_res( random_start );
		lsampler.set_stop_res ( random_end );
		lsampler.set_max_nstruct( 20 );

		lsampler.set_min_bbrms( 20.0   );
		lsampler.set_max_bbrms( 1400.0 );
		lsampler.set_min_rms( 0.5 );
		lsampler.set_max_rms( 4.0 );

		static int casecount = 0;
    core::pose::Pose opose = pose;
    std::vector< core::io::silent::SilentStructOP > lib_structs;

    TR.Info << "Loophash apply function ! " << std::endl;

    // Set up contraints
    ScoreFunctionOP fascorefxn = core::scoring::getScoreFunction();
    //protocols::relax::FastRelax *qrelax = new protocols::relax::FastRelax( fascorefxn, 1 );
    protocols::relax::FastRelax *relax = new protocols::relax::FastRelax( fascorefxn,  option[ OptionKeys::relax::sequence_file ]() );

    // convert pose to centroid pose:
    core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID);
    core::pose::set_ss_from_phipsi( pose );

    core::Size starttime2 = time(NULL);


    lsampler.build_structures( pose, lib_structs );
    //get_all( pose, lib_structs, 1, 0, 20,1400.0, 0.5, 4.0  );   // old way

    core::Size endtime2 = time(NULL);
    TR.Info << "FOUND " << lib_structs.size() << " alternative states in time: " << endtime2 - starttime2 << std::endl;

    //std::random_shuffle( lib_structs.begin(), lib_structs.end());
    numeric::random::random_permutation(lib_structs.begin(), lib_structs.end(), numeric::random::RG);

    std::vector< core::io::silent::SilentStructOP > select_lib_structs;

    for( core::Size k=0;k< std::min(skim_size, lib_structs.size() ) ;k++){
      select_lib_structs.push_back( lib_structs[k] );
    }

    core::pose::Pose native_pose;
    if( option[ in::file::native ].user() ){
      core::import_pose::pose_from_pdb( native_pose, option[ in::file::native ]() );
    } else {
      utility_exit_with_message("This app requires specifying the -in:file:native flag.");
    }

		core::pose::Pose ref_pose;
    if( option[ lh::refstruct].user() ){
      core::import_pose::pose_from_pdb( ref_pose, option[ lh::refstruct ]() );
    }


    core::Real bestcenscore = MAXIMAL_FLOAT;
    core::Size bestcenindex = 0;

 		if( select_lib_structs.size() == 0 ) continue;

		if( last_best ) {
			select_lib_structs.push_back( last_best );
		}

		for( core::Size h = 0; h < select_lib_structs.size(); h++){
      core::pose::Pose rpose;
      select_lib_structs[h]->fill_pose( rpose );

   		core::Real refrms = 0;
			if( option[ lh::refstruct].user()) {
			 	refrms = scoring::CA_rmsd( ref_pose, rpose );
			}
		  core::Real rms_factor =  5.0;
			core::Real decoy_score = select_lib_structs[h]->get_energy("censcore") + refrms * rms_factor;

			select_lib_structs[h]->add_energy( "refrms",     refrms,      1.0 );
			select_lib_structs[h]->add_energy( "comb_score", decoy_score, 1.0 );
			std::cout << "refrms: " << refrms << "  Energy: " << decoy_score << std::endl;
			if( decoy_score < bestcenscore ){
				bestcenscore = decoy_score;
				bestcenindex = h;
				last_best = select_lib_structs[h];
			}
		}
		std::cout << "Best:" << "  Energy: " << bestcenscore << std::endl;



    if((  option[ OptionKeys::lh::write_centroid_structs ]() ) ||
       (  option[ OptionKeys::lh::centroid_only ]() )){

      core::io::silent::SilentFileData sfd;
      std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
      if( option[ OptionKeys::lh::centroid_only ]() ){
        silent_file_ += ".centroid.out" ;
      }

      for( core::Size h = 0; h < select_lib_structs.size(); h++){
          core::pose::Pose rpose;
          select_lib_structs[h]->fill_pose( rpose );
          core::Real rms = scoring::CA_rmsd( native_pose, rpose );
          select_lib_structs[h]->add_energy( "round", round, 1.0 );
          select_lib_structs[h]->add_energy( "rms", rms, 1.0 );
          select_lib_structs[h]->set_decoy_tag( "S_" + string_of( round ) + "_" + string_of(  h )  );
          if( round >= 2.0 ){
						sfd.write_silent_struct( *(select_lib_structs[h]) , silent_file_ );
      		}
			}

    }

		/// In centroid cleanup mode this is IT
		if( option[ OptionKeys::lh::centroid_only ]() ){
      select_lib_structs[bestcenindex]->fill_pose( pose );
			continue;
		}




		/// For fullatom goodness, continue
    core::Real bestscore = MAXIMAL_FLOAT;
    core::Size bestindex = 0;
    // Batch relax the result:

    core::Size starttime = time(NULL);
    relax->batch_apply( select_lib_structs );
    core::Size endtime = time(NULL);
    TR.Info << "Batchrelax time: " << endtime - starttime << " for " << select_lib_structs.size() << " structures " << std::endl;


    for( core::Size h = 0; h < select_lib_structs.size(); h++){
      TR.Info << "DOING: " << h << " / " << select_lib_structs.size() << std::endl;
      core::pose::Pose rpose;

      select_lib_structs[h]->fill_pose( rpose );

      //core::Real score = scoring::CA_rmsd( native_pose, rpose );
      core::Real score = (*fascorefxn)(rpose);
      TR.Info << "score: " << h << "  " << score << std::endl;

      if( score < bestscore ){
        bestscore = score;
        bestindex = h;
        pose = rpose;
      }
    }
    casecount++;
    //test_loop_sample( pose, pose.total_residue() );

    core::Real bestrms = scoring::CA_rmsd( native_pose, pose );
    TR.Info << "BESTSCORE: " << bestscore << "BESTRMS" << bestrms << std::endl;
    //pose.dump_pdb( "lhb_" + prefix + "_" + utility::to_string( round ) + ".pdb" );


    core::io::silent::SilentFileData sfd;
    std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
    for( core::Size h = 0; h < select_lib_structs.size(); h++){

      if( h == bestindex ) {
        core::pose::Pose rpose;
        select_lib_structs[h]->fill_pose( rpose );
        core::Real rms = scoring::CA_rmsd( native_pose, rpose );
        select_lib_structs[h]->add_energy( "round", round, 1.0 );
        select_lib_structs[h]->add_energy( "rms", rms, 1.0 );
        select_lib_structs[h]->set_decoy_tag( "S_" + string_of( round ) + "_" + string_of(  h )  );
        sfd.write_silent_struct( *(select_lib_structs[h]) , silent_file_ );
      }
    }

  }

}

void run_sandbox( LoopHashLibraryOP /*loop_hash_library*/ ){

  core::pose::Pose tgtpose, srcpose;
	core::import_pose::pose_from_pdb( tgtpose, "input/S_00001_0000001_0_0001.pdb" );
	core::import_pose::pose_from_pdb( srcpose, "input/S_00001_0000001_0_1_0001.pdb" );

  // test silent store class
/*
  protocols::wum::SilentStructStore mystore;

  mystore.add( tgtpose );
  mystore.add( srcpose );


  std::string sdata;
  mystore.serialize( sdata );
  std::cout << "----" << std::endl;
  mystore.print( std::cout );


  protocols::wum::SilentStructStore mystore2;

  mystore2.read_from_cmd_line();
  std::cout << "-------" << std::endl;
  mystore2.print( std::cout );
  std::cout << "-AA-AA-" << std::endl;

  mystore.add( mystore2 );

  std::cout << "-------" << std::endl;

  mystore.print( std::cout );

  std::string serial_form;
  mystore.serialize( serial_form );
  mystore.serialize_to_file( "testoutput.out" );


  protocols::wum::SilentStructStore recovered_store;
  recovered_store.read_from_string( serial_form );

  recovered_store.print( std::cout );

*/
// Loopgraft
  // Test 1
//37 63
  //loop_hash_library->graft_loop( srcpose, tgtpose, protocols::loops::Loop( 37, 63 ) );
  //tgtpose.dump_pdb( "output_graft.pdb" );


  // Test 2

  // results
  //21396 Fullatom protein
  //52499     "    binary
  //13879 Centroid protein
  //23918     "    binary
  // 2424 BackboneSegment float
  // 1212 BackboneSegment word

  {
    core::io::silent::ProteinSilentStruct pss;
    pss.fill_struct(tgtpose, "lala" );

    std::ostringstream ss;
    pss.print_scores( ss );
    pss.print_conformation( ss );
   // std::cout << ss.str() << std::endl;
    std::cout << pss.mem_footprint() << "  " << ss.str().length() << std::endl;
  }

//  {
//    core::io::silent::BinaryProteinSilentStruct pss;
//    pss.fill_struct(tgtpose, "lala" );
//
//    std::ostringstream ss;
//    pss.print_scores( ss );
//    pss.print_conformation( ss );
//   // std::cout << ss.str() << std::endl;
//    std::cout << pss.mem_footprint() << "  " << ss.str().length() << std::endl;
//  }

  core::util::switch_to_residue_type_set( tgtpose, core::chemical::CENTROID);


  {
    core::io::silent::ProteinSilentStruct pss;
    pss.fill_struct(tgtpose, "lala" );

    std::ostringstream ss;
    pss.print_scores( ss );
    pss.print_conformation( ss );
   // std::cout << ss.str() << std::endl;
    std::cout <<pss.mem_footprint() << "  " <<  ss.str().length() << std::endl;
  }

  {
    core::io::silent::BinaryProteinSilentStruct pss;
    pss.fill_struct(tgtpose, "lala" );

    std::ostringstream ss;
    pss.print_scores( ss );
    pss.print_conformation( ss );
   // std::cout << ss.str() << std::endl;
    //std::cout <<pss.mem_footprint() << "  " <<  ss.str().length() << std::endl;
  }
}


int
main( int argc, char * argv [] )
{
    try {
	using namespace protocols;
	using namespace protocols::jd2;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;


	// initialize core
	devel::init(argc, argv);

#ifdef USEMPI
	int mpi_rank_, mpi_npes_;
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank_ ) );
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes_ ) );

	// unless you are rank one - go into infinite sleep loop
	if( mpi_rank_ != 0 ){
		TR << "NOT RANK 0: Sleeping .. " << std::endl;
		while(true){
			sleep( 10 );
		}
	}

#endif


	std::cout << "SIZEOF short: " << sizeof( short ) << std::endl;
	std::cout << "SIZEOF short*: " << sizeof( short * ) << std::endl;
	std::cout << "SIZEOF core::Size: " << sizeof( core::Size ) << std::endl;
	std::cout << "SIZEOF core::Size: " << sizeof( unsigned int ) << std::endl;

	utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
	LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loop_sizes );


  // sandbox mode ?
	if ( option[lh::sandbox]() ){;
		run_sandbox( loop_hash_library );
    return 0;
	}

	// Run simple sampling run test or create the db ?
	if ( option[lh::create_db]() ){;
		loop_hash_library->create_db();
		loop_hash_library->save_db();
		return 0;
	}

  LoopHashRelax_SamplerOP lh_sampler = new LoopHashRelax_Sampler( loop_hash_library );

  // Normal mode with external loophash library
  loop_hash_library->load_db();
  try{
    //protocols::jd2::JobDistributor::get_instance()->go( loop_hash_library );
    protocols::jd2::JobDistributor::get_instance()->go( lh_sampler );
  } catch ( utility::excn::EXCN_Base& excn ) {
    std::cerr << "Exception: " << std::endl;
    excn.show( std::cerr );
    std::cout << "Exception: " << std::endl;
    excn.show( std::cout ); //so its also seen in a >LOG file
  }

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}


