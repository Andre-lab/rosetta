// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

// #include <protocols/abinitio/Templates.hh>
// #include <protocols/abinitio/TemplateJumpSetup.hh>
// #include <protocols/abinitio/PairingStatistics.hh>
// #include <protocols/abinitio/StrandConstraints.hh>
// #include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>
//#include <core/pose/util.hh>
#include <devel/init.hh>
//#include <core/io/pdb/pose_io.hh>

#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
//#include <core/scoring/constraints/BoundConstraint.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/util.hh>
//for derivative check
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
//#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rdc.OptionKeys.gen.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/excn/Exceptions.hh>

//Auto Headers
#include <utility/excn/EXCN_Base.hh>
#include <ObjexxFCL/format.hh>



static basic::Tracer tr("main");

using namespace core;
using namespace protocols;
//using namespace abinitio;
//using namespace jumping;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace scoring;
using namespace ObjexxFCL::fmt;

OPT_1GRP_KEY( Boolean, score_app, linmin )
OPT_KEY( Boolean, dump_all )
OPT_KEY( File, average )
OPT_KEY( File, residue_subset )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//  Templates::register_options();
  OPT( in::file::s ); //should have a Jobdistributor register_options..
	OPT( in::file::silent );
	NEW_OPT( score_app::linmin, "Do a quick round of linmin before reporting the score", false );
	NEW_OPT( dump_all, "show individual RDCs for each structure on screen", false );
	NEW_OPT( average, "compute average RDC and mean_deviation for each residue", "average.dat" );
	NEW_OPT( residue_subset, "loop-file that gives the regions to compute RDC score", "rigid.dat" );
	OPT(  rdc::print_rdc_values );
}

// Forward
class RDCToolMover;

// Types
typedef  utility::pointer::owning_ptr< RDCToolMover >  RDCToolMoverOP;
typedef  utility::pointer::owning_ptr< RDCToolMover const >  RDCToolMoverCOP;

class RDCToolMover : public moves::Mover {
public:
	RDCToolMover() : myRDCs_( NULL ) {};
	virtual void apply( core::pose::Pose& );
	virtual std::string get_name() const {
		return "RDC-Tool";
	}

	void dump_averages();
private:
	ResidualDipolarCouplingOP myRDCs_;
	typedef utility::vector1< core::Real > RDC_Vector;
	RDC_Vector average_RDCs_;
	RDC_Vector var_RDCs_;
	RDC_Vector av_dev_RDCs_;
	Size nstruct_;
};

void RDCToolMover::dump_averages() {
	ResidualDipolarCoupling::RDC_lines data = myRDCs_->get_RDC_data();
	runtime_assert( average_RDCs_.size() == data.size() );
	Real invn = 1.0/nstruct_;
	utility::io::ozstream ofs( option[ average ]() );
	for ( Size d = 1; d<=data.size(); d++ ) {
		ofs << RJ( 6, myRDCs_->get_RDC_data()[ d ].res1() ) << " "
				<< RJ( 6, myRDCs_->get_RDC_data()[ d ].Jdipolar()) << " "
				<< RJ( 6, average_RDCs_[ d ]*invn ) << " "
				<< RJ( 6, sqrt( av_dev_RDCs_[ d ]*invn ) ) << " "
				<< RJ( 6, sqrt(var_RDCs_[ d ]*invn - (average_RDCs_[ d ]*invn)*(average_RDCs_[ d ]*invn))) << std::endl;
	}
}

void RDCToolMover::apply( core::pose::Pose &pose ) {

	if ( !myRDCs_ ) {
		myRDCs_ = new ResidualDipolarCoupling;
		ResidualDipolarCoupling::RDC_lines data = myRDCs_->get_RDC_data();
		if ( option[ average ].user() ) {
			average_RDCs_.resize( data.size(), 0.0 );
			var_RDCs_.resize( data.size(), 0.0 );
			av_dev_RDCs_.resize( data.size(), 0.0 );
			for ( Size d = 1; d<=data.size(); d++ ) {
				average_RDCs_[ d ]=0.0;
				av_dev_RDCs_[d ]=0.0;
				var_RDCs_[d ]=0.0;
			}
			nstruct_=0;
		}
		if ( option[ residue_subset ].user() ) { //filter
			loops::LoopsFileIO loop_file_reader;
			std::ifstream is( option[ residue_subset ]().name().c_str() );
			
			if (!is.good()) {
				utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + option[ residue_subset ]().name() + "'" );
			}
			
			loops::SerializedLoopList loops = loop_file_reader.use_custom_legacy_file_format(is, option[ residue_subset ](), false, "RIGID");
			loops::Loops rigid_core = loops::Loops( loops ); // <==
			ResidualDipolarCoupling::RDC_lines filtered;
			for ( ResidualDipolarCoupling::RDC_lines::const_iterator it = data.begin(); it != data.end(); ++it ) {
				if ( rigid_core.has( it->res1() ) ) {
					filtered.push_back( *it );
				}
			}
			myRDCs_ = new ResidualDipolarCoupling( filtered );
		} // end filter
	} //end init myRDCs_
	store_RDC_in_pose( myRDCs_, pose );

	scoring::ScoreFunction score;
	score.set_weight( scoring::rdc, 1.0 );

	std::string fname( jd2::current_output_name() );
	tr.Info << fname << " " << score( pose ) << std::endl;

	ResidualDipolarCouplingCOP someRDCs = retrieve_RDC_from_pose( pose );
	ResidualDipolarCoupling::RDC_lines data = someRDCs->get_RDC_data();

	if ( option[ dump_all ]() || option[ average ].user() ) {
		Size d = 0;
		for ( ResidualDipolarCoupling::RDC_lines::const_iterator it = data.begin(); it != data.end(); it ++ ) {
			++d;
			if ( option[ dump_all ]() ) std::cout << RJ( 6, it->Jdipolar()) << " " << RJ( 6, it->Jcomputed() ) << std::endl;
			if ( option[ average ].user() ) {
				Real dev = it->Jdipolar() - it->Jcomputed();
				Real dev2 = dev*dev;
				runtime_assert( average_RDCs_.size() >= d );
				average_RDCs_[ d ]+=it->Jcomputed();
				var_RDCs_[ d ]+=it->Jcomputed()*it->Jcomputed();
				av_dev_RDCs_[ d ]+=dev2;
			}
		}
	}
	std::cout << fname << " R: " << someRDCs->R() << " Q: " << someRDCs->Q() << std::endl;
	++nstruct_;
	if ( option[ score_app::linmin ]() ) {
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->set_bb( true ); movemap->set_chi( true );
		protocols::simple_moves::MinMoverOP minmover = new protocols::simple_moves::MinMover(
					 movemap, score.clone(), "linmin", 1e-4,
					 false /*use_nblist*/, true /*deriv_check*/, true /*verbose driv check*/ );
		minmover->apply( pose );
	}


}

void run() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !basic::options::option[ basic::options::OptionKeys::in::path::database ].user() ) {
		basic::options::option[ basic::options::OptionKeys::in::path::database ].def( "/work/olange/minirosetta_database");
	}

	RDCToolMoverOP cst_tool =  new RDCToolMover;
	protocols::jd2::JobDistributor::get_instance()->go( cst_tool, new jd2::NoOutputJobOutputter );


	if ( option[ average ].user() ) {
		cst_tool->dump_averages();
	}

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	register_options();
	devel::init( argc, argv );

	//	try{
		run();
		//	} catch ( utility::excn::EXCN_Base& anExcn ) {
		//		anExcn.show( std::cerr );
		//	}

	return 0;
}


