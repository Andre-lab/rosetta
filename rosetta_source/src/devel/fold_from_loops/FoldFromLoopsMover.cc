// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.
// author bcorreia



#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <devel/fold_from_loops/FoldFromLoopsMover.hh>
#include <devel/fold_from_loops/FoldFromLoops_functions.hh>


#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>


// AUTO-REMOVED #include <core/fragment/FragData.hh>

//scoring

#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>


//protocols

#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/loops/SlidingWindowLoopClosure.hh>
#include <protocols/loops/SlidingWindowLoopClosure.fwd.hh>
#include <protocols/jumping/util.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/loops/Exceptions.hh>
// AUTO-REMOVED #include <protocols/relax_protocols.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>


//constraints
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
// Auto-header: duplicate removed #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Func.hh>

//design headers

// AUTO-REMOVED #include <protocols/moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>


//options

#include <basic/options/option.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>



#include <utility/string_util.hh>

// C++ headers
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>
// AUTO-REMOVED #include <ctime>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
// Auto-header: duplicate removed #include <core/pose/util.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <core/util/SwitchResidueTypeSet.hh>


using namespace core;
using namespace kinematics;


namespace devel {
namespace fold_from_loops {

FoldFromLoopsMover::FoldFromLoopsMover()
{}

FoldFromLoopsMover::~FoldFromLoopsMover()
{}


void FoldFromLoopsMover::apply (core::pose::Pose & input_pose )
{

	using namespace core::scoring;
	using namespace kinematics;
	using namespace devel::fold_from_loops;
	using namespace constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::fragment;
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_RETRY;


	basic::Tracer TR("FoldFromLoopsMover");

	protocols::checkpoint::CheckPointer sliding_checkpoint("closing"); // annoying  see if you can take it out of here


	core::pose::Pose archive_pose = input_pose;

	kinematics::FoldTree f;
	f.clear();

	std::vector<Size> cut_points;

	core::util::switch_to_residue_type_set( input_pose, core::chemical::CENTROID );


	core::kinematics::MoveMapOP  movemap = new core::kinematics::MoveMap();


	fold_tree_cutpoints_generator(loops_, cut_points, input_pose , f);



	define_movemap( movemap, input_pose , loops_);


    core::pose::Pose loops_pdb = loops_pdb_;



	new_pose_generator( loops_pdb  , input_pose , loops_, cut_points );



	input_pose = loops_pdb; //update the input pose to the newly generated pose



	input_pose.fold_tree( f );


	refresh_cutpoints(input_pose, cut_points);





	if (option[ OptionKeys::in::file::psipred_ss2 ].user() ) {


			protocols::loops::set_secstruct_from_psipred_ss2( input_pose ); //needs to be done inside the mover



	}




	protocols::abinitio::ClassicAbinitioOP abinitio = new protocols::abinitio::ClassicAbinitio( frag_small_, frag_large_ ,movemap );

		if (option [OptionKeys::fold_from_loops::native_ca_cst].user() ){

				input_pose.constraint_set( ca_csts_ );

				abinitio = new protocols::abinitio::FoldConstraints( frag_small_, frag_large_ ,movemap );

			}

    //input_pose.dump_pdb("extended_pose");


	abinitio->init(input_pose);
	abinitio->apply( input_pose);




	TR << "Pose Secondary Structure After abinitio " << input_pose.secstruct() << std::endl;


	Real rmsd_to_native =  core::scoring::rmsd_with_super(input_pose, archive_pose, is_protein_CA);

	TR << "RMSD TO NATIVE  " <<  rmsd_to_native << std::endl;


	core::scoring::ScoreFunctionOP scorefxn_centr = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std","score4L" );

	pose::setPoseExtraScores( input_pose, "rms",   rmsd_to_native );
	(*scorefxn_centr)(input_pose);





		if ( rmsd_to_native < option[ OptionKeys::fold_from_loops::ca_rmsd_cutoff] )
			{

			set_last_move_status( MS_SUCCESS );

			if ( option[ OptionKeys::loops::ccd_closure ].user()  ){ // to close loops

					Real chain_break_dist = 0; // chain_break eval

					chain_break_dist = input_pose.energies().total_energies()[ scoring::linear_chainbreak ];


					TR << "Chain Break -  " << chain_break_dist << std::endl;




					protocols::loops::SlidingWindowLoopClosureOP closure_protocol = new protocols::loops::SlidingWindowLoopClosure( frag_small_, scorefxn_centr, movemap );


					closure_protocol->scored_frag_cycle_ratio( 0.2 );
					closure_protocol->short_frag_cycle_ratio( 0.1 );

					// TODO: fix the new checkpointer issue
						Real loop_tag = 1;
						try{
							protocols::jumping::close_chainbreaks( closure_protocol, input_pose, sliding_checkpoint ,"sliding", f );
						} catch ( protocols::loops::EXCN_Loop_not_closed& excn ) {
							loop_tag = 0;
							set_last_move_status ( FAIL_RETRY ); // if the loops fail set fail_retry

						}

						TR << "LOOPS_TAG " << loop_tag << std::endl;


		    }



			if ( get_last_move_status() == MS_SUCCESS && option [OptionKeys::fold_from_loops::output_centroid].user() ){

				(*scorefxn_centr)(input_pose);

				TR << "Outputing centroid structure  " << std::endl;

			}


			else if ( get_last_move_status() == MS_SUCCESS  && !option [OptionKeys::fold_from_loops::output_centroid].user() ){

			core::util::switch_to_residue_type_set( input_pose, chemical::FA_STANDARD );


			core::scoring::ScoreFunctionOP scorefxn_fa( ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH ) );
			scorefxn_fa->set_weight( scoring::chainbreak, 1 ); // in order to pipe this scoring function  maybe I have to refresh the cutpoints
			scorefxn_fa->set_weight(scoring::overlap_chainbreak, 1 );




			//if (loops_.size() > 1 ){

				//refresh_cutpoints(input_pose, cut_points);

			//}


			copying_side_chains_swap_loop( loops_pdb_ , input_pose, loops_, movemap );



			(*scorefxn_fa)(input_pose);


			design_excluding_swap_loops (  input_pose, loops_, scorefxn_fa );



			protocols::relax::ClassicRelax relax_protocol( scorefxn_fa, movemap);



			relax_protocol.apply( input_pose );



			std::string decoy_name = option [OptionKeys::out::prefix] + "dr_init";


			protocols::jd2::JobDistributor::get_instance()->current_job()->set_status_prefix( decoy_name ); // setting up struct name



			protocols::jd2::output_intermediate_pose( input_pose, std::string() ); //output it to the same file name



			if ( option [OptionKeys::fold_from_loops::add_relax_cycles ].user() ){

						Size n_relax = option [OptionKeys::fold_from_loops::add_relax_cycles ] ;

						for ( Size i=1; i <= n_relax ; ++i ){


							std::string outfilename_n_rlx =  option [OptionKeys::out::prefix] + "dr_" + right_string_of(i,3,'0')+"_" ;


							design_excluding_swap_loops (  input_pose, loops_, scorefxn_fa );


							relax_protocol.apply( input_pose );


							protocols::jd2::JobDistributor::get_instance()->current_job()->set_status_prefix( outfilename_n_rlx );

							protocols::jd2::output_intermediate_pose( input_pose, std::string() );


						}

				}


			protocols::jd2::JobDistributor::get_instance()->current_job()->set_status_prefix( std::string() );


			}


		} // RMSD cutoff for the relax cycles



		else {

			set_last_move_status( FAIL_RETRY );

		}


	} //mover apply function

std::string
FoldFromLoopsMover::get_name() const {
	return "FoldFromLoopsMover";
}



}
}
