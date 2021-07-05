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

#include <core/types.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

//#include <boost/spirit/core.hpp>


#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>

#include <boost/spirit/home/classic/actor/push_back_actor.hpp> // AUTO IWYU For push_back_a, push_back_action
#include <boost/spirit/home/classic/core/composite/actions.hpp> // AUTO IWYU For action
#include <boost/spirit/home/classic/core/composite/kleene_star.hpp> // AUTO IWYU For kleene_star, operator*
#include <boost/spirit/home/classic/core/composite/sequence.hpp> // AUTO IWYU For sequence, operator>>
#include <boost/spirit/home/classic/core/primitives/numerics.hpp> // AUTO IWYU For real_parser, real_p

using namespace boost::spirit::classic;

bool parse_numbers( char const* str, utility::vector1< core::Real > & v ) {
	return parse(str,
		// Begin grammar
		( real_p[push_back_a(v)] >> *( ',' >> real_p[push_back_a(v)] ) ),
		// End grammar
		space_p
		).full;
}

//bool parse_range(
// char const* str,
// utility::vector1< std::pair< core::Real, core::Real > > & v
//) {
// return parse(str,
//  // Begin grammar
//  ( real_p[ push_back_a(v) ] >> *( '-' >> real_p[push_back_a(v)] ) ),
//  // End grammar
//  space_p
// ).full;
//}

int
main( int argc, char * argv [] ) {
	try {

		using std::pair;
		using core::Real;
		using utility::vector1;
		devel::init( argc, argv );
		// parsing a list of comma-separated numbers
		std::string s( "5.0,3,2.1,0" );
		utility::vector1< core::Real > numbers;
		parse_numbers( s.c_str(), numbers );
		for ( vector1< Real >::const_iterator it = numbers.begin(),
				end = numbers.end(); it != end; ++it
				) {
			std::cout << *it << std::endl;
		}

		// parse a list of ranges
		//s = "1-5,3-9,11-27";
		//vector1< pair< Real, Real > > ranges;
		//for ( vector1< pair< Real, Real > >::const_iterator it = ranges.begin(),
		//   end = ranges.end(); it != end; ++it
		//) {
		// std::cout << it->first() << " => " << it->secong() << std::endl;
		//}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
} // main
