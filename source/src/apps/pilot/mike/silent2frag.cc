// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

// C++ headers
#include <iostream>
#include <string>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <protocols/idealize/IdealizeMover.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


using basic::Error;
using basic::Warning;


using namespace core;
using namespace core::fragment;
using namespace protocols;
using namespace utility;

using utility::vector1;

namespace protocols {
namespace moves {

class FragCaptureMover : public moves::Mover
{
public:
	/// @brief
	///  empty constructor fills values with the values
	///  read in from the commandline
	FragCaptureMover( core::Size length) :
		Mover(),
		fragset_( length )
	{
		Mover::type( "FragCaptureMover" );
	}


	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return "FragCaptureMover"; }
	void test_move( core::pose::Pose & pose ) override{apply(pose);}

	void write( const std::string &filename );


private:  //  Data

	ConstantLengthFragSet fragset_;


};

using FragCaptureMoverOP = utility::pointer::shared_ptr<FragCaptureMover>;
using FragCaptureMoverCOP = utility::pointer::shared_ptr<const FragCaptureMover>;


void FragCaptureMover::apply( core::pose::Pose & pose )
{
	std::cout << "Stealing .. " << std::endl;
	steal_constant_length_frag_set_from_pose (  pose, fragset_ );
}

void FragCaptureMover::write( const std::string &filename )
{
	FragmentIO fragio;
	fragio.write_data( filename, fragset_ );
}

} // moves
} // protocols


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::moves;

		devel::init(argc, argv);
		jd2::register_options();

		FragCaptureMoverOP capture3mers( new FragCaptureMover(3) );
		FragCaptureMoverOP capture9mers( new FragCaptureMover(9) );

		SequenceMoverOP seqmov( new SequenceMover );
		if ( option[ run::idealize_before_protocol ]() ) {
			seqmov->add_mover( utility::pointer::make_shared< protocols::idealize::IdealizeMover >() );
		}

		seqmov->add_mover( capture3mers );
		seqmov->add_mover( capture9mers );
		MoverOP mover = seqmov;
		protocols::jd2::JobDistributor::get_instance()->go( seqmov );

		std::string prefix = "";
		if ( option[ out::prefix ].user() ) prefix = option[ out::prefix ]() + "_" ;
		capture3mers->write( prefix + "aa3mer.1_3" );
		capture9mers->write( prefix + "aa9mer.1_3" );

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
