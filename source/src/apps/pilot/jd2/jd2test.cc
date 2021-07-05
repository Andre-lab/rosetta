// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/jd2/jd2test.cc
/// @brief test for new job distributor - demonstrates most of its functionality

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <devel/init.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/Mover.hh>

#include <basic/Tracer.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>


OPT_1GRP_KEY( Boolean, jd2test, set_fail_no_retry )
OPT_1GRP_KEY( Boolean, jd2test, set_fail_bad_input )

static basic::Tracer TR( "jd2test" );

///local mover for testing purposes
class JDtestmover : public protocols::moves::Mover {
public:
	JDtestmover()
	{
		TR << "JDtestmover ctor" << std::endl;
		//setup of packing from fixbb
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		main_task_factory = utility::pointer::make_shared< TaskFactory >();
		main_task_factory->push_back( utility::pointer::make_shared< operation::InitializeFromCommandline >() );
		if ( option[ packing::resfile ].user() ) {
			main_task_factory->push_back( utility::pointer::make_shared< operation::ReadResfile >() );
		}

		score_fxn = core::scoring::get_score_function();

		pack_mover = utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >();
		pack_mover->task_factory( main_task_factory );
		pack_mover->score_function( score_fxn );

		//set up another mover
		sc_mover = utility::pointer::make_shared< protocols::simple_moves::sidechain_moves::SidechainMover >();
		sc_mover->set_task_factory( main_task_factory );

		set_fail_no_retry_ = false;
		set_fail_bad_input_ = false;
	}

	~JDtestmover() override = default;


	void
	apply( core::pose::Pose & pose ) override{
		static int counter(0);
		++counter;

		using protocols::jd2::JobDistributor;
		protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
		std::string me(JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );
		TR << "starting a job, JD says it is: " << me << std::endl;

		//clock_t starttime = clock();
		//TR << starttime << std::endl;

		//trivial changes so that input and output poses are not identical
		pack_mover->apply(pose);

		//this tests that the mid-pose output has the correct scorefunction (meaning that it also includes chainbreak)
		core::scoring::ScoreFunctionOP score_fxn2 = core::scoring::get_score_function();
		score_fxn2->set_weight( core::scoring::chainbreak, 2.0 );
		(*score_fxn2)(pose);
		/// Now handled automatically.  score_fxn2->accumulate_residue_total_energies(pose);

		//this will show up in the intermediate output
		job_me->add_string("somestuff");
		job_me->add_string_string_pair("alhpa", "btea");
		job_me->add_string_real_pair("theanswer", 42.0);

		//test mid-job output
		JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "lookit");
		JobDistributor::get_instance()->job_outputter()->file( job_me, "hey here's some junk from job " + me + "\n" );

		//Now, what happens if we don't want that intermediate output?
		TR << "next line will be odd-looking Job construction" << std::endl;
		protocols::jd2::JobOP job_faux( JobDistributor::get_instance()->current_job()->copy_without_output());
		JobDistributor::get_instance()->job_outputter()->other_pose( job_faux, pose, "no_accessory_data");

		//trivial change to make sure that the "moving" actually occurs
		sc_mover->apply(pose);

		//examples for pose extra scores
		//DEPRECATED
		//   core::pose::setPoseExtraScore(pose, "this_is_a_test" + me, 987654321);
		//   core::pose::setPoseExtraScore(pose, "alpha", 1.000);
		//   core::pose::setPoseExtraScore(pose, "beta", 24000);
		//   core::pose::setPoseExtraScore(pose, "gamma", 3.000);
		//   core::pose::setPoseExtraScore(pose, "delta", 4.000);

		//this will show up in the final output, along with one set of intermediate stuff above
		job_me->add_string("somemorestuff");
		job_me->add_string_string_pair("alpha", "beta");
		job_me->add_string_real_pair("thequestion", 42.0);

		//test mover reporting
		if ( counter == 1 && set_fail_no_retry_ ) set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		if ( counter == 3 ) set_last_move_status(protocols::moves::FAIL_RETRY);
		if ( counter == 1 && set_fail_bad_input_ ) set_last_move_status(protocols::moves::FAIL_BAD_INPUT);
		//MPI testing
		//if(me == "1UBQ_0002") set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		//if(me == "1EM7_0003" && counter <= 40) set_last_move_status(protocols::moves::FAIL_RETRY);
		//if(me.substr(0,5) == "1NLO_") set_last_move_status(protocols::moves::FAIL_BAD_INPUT);
		//while( clock() <= (starttime + 6000000) ) {}
		//TR << clock() << std::endl;

		//necessary for score printout to work right
		(*score_fxn)(pose);
		/// Now handled automatically.  score_fxn->accumulate_residue_total_energies(pose);
		return;
	}


	std::string
	get_name() const override {
		return "JDtestmover";
	}


	protocols::moves::MoverOP
	fresh_instance() const override {
		return utility::pointer::make_shared< JDtestmover >();
	}


	bool
	reinitialize_for_each_job() const override { return false; }


	bool
	reinitialize_for_new_input() const override { return false; }

	bool set_fail_no_retry_;
	bool set_fail_bad_input_;

private:
	//original holdovers
	core::scoring::ScoreFunctionOP score_fxn;
	protocols::minimization_packing::PackRotamersMoverOP pack_mover;
	core::pack::task::TaskFactoryOP main_task_factory;

	//second mover for testing
	protocols::simple_moves::sidechain_moves::SidechainMoverOP sc_mover;

};

using JDtestmoverOP = utility::pointer::shared_ptr<JDtestmover>;


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		NEW_OPT( jd2test::set_fail_no_retry, "Set mover status to fail_no_retry", false );
		NEW_OPT( jd2test::set_fail_bad_input, "Set mover status to fail_bad_input", false );

		devel::init(argc, argv);

		bool set_fail_no_retry = basic::options::option[ basic::options::OptionKeys::jd2test::set_fail_no_retry ];
		bool set_fail_bad_input = basic::options::option[ basic::options::OptionKeys::jd2test::set_fail_bad_input ];

		JDtestmoverOP test_mover(new JDtestmover);


		if ( set_fail_no_retry ) {
			test_mover->set_fail_no_retry_ = true;
			try {
				protocols::jd2::JobDistributor::get_instance()->go(test_mover);
			} catch (protocols::jd2::JD2Failure const & e ) {
				TR << "successfully caught JD2 exception after fail_no_retry" << std::endl;
			}
		} else if ( set_fail_bad_input ) {
			test_mover->set_fail_bad_input_ = true;
			try {
				protocols::jd2::JobDistributor::get_instance()->go(test_mover);
			} catch (protocols::jd2::JD2Failure const & e ) {
				TR << "successfully caught JD2 exception after fail_bad_input" << std::endl;
			}
		} else {
			protocols::jd2::JobDistributor::get_instance()->go(test_mover);
			TR << "*********************successful completion**************************" << std::endl;
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}

