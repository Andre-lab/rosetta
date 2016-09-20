// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/string_util.cxxtest.hh
/// @brief  string_util.cxxtest: test suite for utility::string_util
/// @author James Thompson
/// @author Christopher Miles (cmiles@uw.edu)
/// @author Rhiju Das

// Testing headers
#include <cxxtest/TestSuite.h>


// Project headers
#include <utility/string_util.hh>
#include <utility/tools/make_vector1.hh>
#include <core/types.hh>

// C/C++ headers
#include <iostream>
#include <string>
#include <sstream>

//#include <basic/Tracer.hh>

using std::endl;
using std::string;
using std::stringstream;

//static THREAD_LOCAL basic::Tracer TR("StringUtil");

class StringUtilTests : public CxxTest::TestSuite {
 public:

	// duplicated implementation in FileName... now it's just an alias for utility::basename
	void test_file_basename() {
		TS_ASSERT_EQUALS(utility::file_basename("core/scoring/ScoreFunction.cc"),
						 "ScoreFunction.cc");
	}

	void test_filename() {
		TS_ASSERT_EQUALS(utility::filename("/foo/bar/baz"), "baz");
	}

	void test_pathname() {
		TS_ASSERT_EQUALS(utility::pathname("/foo/bar/baz"), "/foo/bar/");
	}

	void test_same_ignoring_spaces() {
		TS_ASSERT(utility::same_ignoring_spaces("CA", "CA"));
		TS_ASSERT(utility::same_ignoring_spaces(" CA", "CA"));
		TS_ASSERT(utility::same_ignoring_spaces("CA ", "CA"));
		TS_ASSERT(utility::same_ignoring_spaces("   CA     ", "  CA  "));
	}

	void test_string2uint() {
		string s = "42";
		unsigned int i;
		utility::string2uint(s, &i);
		TS_ASSERT_EQUALS(i, 42);
	}

	void test_string_to_sha1() {
		std::string test_string = "hello world";
		std::string real_hash = "2aae6c35c94fcfb415dbe95f408b9ce91ee846ed";

		std::string test_hash = utility::string_to_sha1(test_string);
		TS_ASSERT_EQUALS(real_hash,test_hash);
	}


	void test_reschain_to_string() {
		utility::vector1< int >  res_vector = utility::tools::make_vector1( -5, -4, -3, 1, 2, 3, 1, 2);
		utility::vector1< char > chain_vector = utility::tools::make_vector1( 'A','A','A','A','A','A',' ','B' );
		std::string tag = make_tag_with_dashes( res_vector, chain_vector );
		TS_ASSERT_EQUALS( tag, "A:-5--3 A:1-3 1 B:2" );
	}

	void run_test_of_get_resnum_and_chain( std::string const & tag ) {
		bool ok;
		std::pair< std::vector<int>, std::vector<char> > resnum_chain = utility::get_resnum_and_chain( tag, ok );

		TS_ASSERT( ok );
		utility::vector1<int>  resnum( resnum_chain.first );
		utility::vector1<char> chains( resnum_chain.second );
		TS_ASSERT_EQUALS( resnum.size(), 8 );
		TS_ASSERT_EQUALS( resnum[1], -5 );
		TS_ASSERT_EQUALS( resnum[2], -4 );
		TS_ASSERT_EQUALS( resnum[3], -3 );
		TS_ASSERT_EQUALS( resnum[4],  1 );
		TS_ASSERT_EQUALS( resnum[5],  2 );
		TS_ASSERT_EQUALS( resnum[6],  3 );
		TS_ASSERT_EQUALS( resnum[7],  1 );
		TS_ASSERT_EQUALS( resnum[8],  2 );

		TS_ASSERT_EQUALS( chains.size(), 8 );
		TS_ASSERT_EQUALS( chains[1], 'A' );
		TS_ASSERT_EQUALS( chains[2], 'A' );
		TS_ASSERT_EQUALS( chains[3], 'A' );
		TS_ASSERT_EQUALS( chains[4], 'A' );
		TS_ASSERT_EQUALS( chains[5], 'A' );
		TS_ASSERT_EQUALS( chains[6], 'A' );
		TS_ASSERT_EQUALS( chains[7], ' ' );
		TS_ASSERT_EQUALS( chains[8], 'B' );
	}

	void test_string_to_reschain() {
		std::string tag( "hello world" );
		bool ok;
		utility::get_resnum_and_chain( tag, ok );
		TS_ASSERT( !ok );

		tag = "A:-5--3 A:1-3 1 B:2";
		run_test_of_get_resnum_and_chain( tag );

		// try an edge case
		tag = "  A-5--3,A1-3 :1 B2-2";
		run_test_of_get_resnum_and_chain( tag );
	}

	void test_make_segtag() {
		utility::vector1< int >  res_vector = utility::tools::make_vector1( -5, -4, 2, 3, 0, 0, 1, 2);
		utility::vector1< std::string > segid_vector = utility::tools::make_vector1( "    ","    ",
																																								 "   A","   A",
																																								 "X   ","Y   ",
																																								 "BLAH","BLAH" );
		std::string tag = make_segtag_with_dashes( res_vector, segid_vector );
		// comma delimiters would be easier to see...
		TS_ASSERT_EQUALS( tag, "    :-5--4    A:2-3 X   :0 Y   :0 BLAH:1-2" );
	}

	void run_test_of_get_resnum_and_segid( std::string const & tag )
	{
		bool ok;
		std::pair< std::vector<int>, std::vector<std::string> > resnum_segid = utility::get_resnum_and_segid( tag, ok );
		TS_ASSERT( ok );
		utility::vector1<int>  resnum( resnum_segid.first );
		utility::vector1<string> segids( resnum_segid.second );
		TS_ASSERT_EQUALS( resnum.size(), 8 );
		TS_ASSERT_EQUALS( resnum[1], -5 );
		TS_ASSERT_EQUALS( resnum[2], -4 );
		TS_ASSERT_EQUALS( resnum[3],  2 );
		TS_ASSERT_EQUALS( resnum[4],  3 );
		TS_ASSERT_EQUALS( resnum[5],  0 );
		TS_ASSERT_EQUALS( resnum[6],  0 );
		TS_ASSERT_EQUALS( resnum[7],  1 );
		TS_ASSERT_EQUALS( resnum[8],  2 );

		TS_ASSERT_EQUALS( segids.size(), 8 );
		TS_ASSERT_EQUALS( segids[1], "    " );
		TS_ASSERT_EQUALS( segids[2], "    " );
		TS_ASSERT_EQUALS( segids[3], "   A" );
		TS_ASSERT_EQUALS( segids[4], "   A" );
		TS_ASSERT_EQUALS( segids[5], "X   " );
		TS_ASSERT_EQUALS( segids[6], "Y   " );
		TS_ASSERT_EQUALS( segids[7], "BLAH" );
		TS_ASSERT_EQUALS( segids[8], "BLAH" );
	}

	void test_string_to_res_and_segid() {
		std::string tag( "hello world" );
		bool ok;
		utility::get_resnum_and_segid( tag, ok );
		TS_ASSERT( !ok );

		tag = "    :-5--4    A:2-3 X   :0 Y   :0 BLAH:1-2";
		run_test_of_get_resnum_and_segid( tag );

		tag = "    :-5--4,   A:2-3,X   :0,Y   :0,BLAH:1-2";
		run_test_of_get_resnum_and_segid( tag );

		// how the tag might look from, e.g., a silent file -- weird spaces.
		tag = ":-5--4   A:2-3 X   :0    Y   :0    BLAH:1-2";
		run_test_of_get_resnum_and_segid( tag );

	}

	void test_padding() {
		std::string sL = "padL";
		std::string new_strL = utility::pad_left(sL, 6, ' ');
		//std::cout << ":" << sL << ":" << std::endl;
		TS_ASSERT_EQUALS( new_strL, "  padL");

		std::string sR = "padR";
		std::string new_strR = utility::pad_right(sR, 6, ' ');
		//std::cout << ":" << sR << ":" << std::endl;
		TS_ASSERT_EQUALS( new_strR, "padR  ");

	}

	void test_real_format() {
		core::Real n = 1.23456;
		core::Real o = 11;
		//core::Real n2 = 22.23456;


		//std::cout << "2Decimal:"<<utility::Real2string(n, 2) << std::endl;
		//std::cout << "3Decimal:"<<utility::Real2string(n2, 3) << std::endl;
		//std::cout << "padded2:"<<  utility::fmt_real(n, 2, 3) << std::endl;
		TS_ASSERT_EQUALS( utility::Real2string(n, 2), "1.23" );
		TS_ASSERT_EQUALS( utility::fmt_real(n, 2, 3), " 1.235" );
		TS_ASSERT_EQUALS( utility::fmt_real(o, 2, 3), "11.000" );

	}

	void test_quoted_split() {
		utility::vector1< std::string > v;

		v = utility::quoted_split("one");
		TS_ASSERT_EQUALS( v.size(), 1 );
		TS_ASSERT_EQUALS( v[1], "one" );

		v = utility::quoted_split("one two");
		TS_ASSERT_EQUALS( v.size(), 2 );
		TS_ASSERT_EQUALS( v[1], "one" );
		TS_ASSERT_EQUALS( v[2], "two" );

		v = utility::quoted_split("    multiple   spaces  ");
		TS_ASSERT_EQUALS( v.size(), 2 );
		TS_ASSERT_EQUALS( v[1], "multiple" );
		TS_ASSERT_EQUALS( v[2], "spaces" );

		v = utility::quoted_split("'single quotes' are 'working well '");
		TS_ASSERT_EQUALS( v.size(), 3 );
		TS_ASSERT_EQUALS( v[1], "'single quotes'" );
		TS_ASSERT_EQUALS( v[2], "are" );
		TS_ASSERT_EQUALS( v[3], "'working well '" );

		v = utility::quoted_split("\"double quotes\" are \" awesome \"  ");
		TS_ASSERT_EQUALS( v.size(), 3 );
		TS_ASSERT_EQUALS( v[1], "\"double quotes\"" );
		TS_ASSERT_EQUALS( v[2], "are" );
		TS_ASSERT_EQUALS( v[3], "\" awesome \"" );

		v = utility::quoted_split(" \"it's nested\"and con-'\"jugated'");
		TS_ASSERT_EQUALS( v.size(), 2 );
		TS_ASSERT_EQUALS( v[1], "\"it's nested\"and" );
		TS_ASSERT_EQUALS( v[2], "con-'\"jugated'" );

		v = utility::quoted_split("'es\\\"ca\\ \\'pe\\\\d' ");
		TS_ASSERT_EQUALS( v.size(), 1 );
		TS_ASSERT_EQUALS( v[1], "'es\\\"ca\\ \\'pe\\\\d'" );

		v = utility::quoted_split("'quotes will'keep \"things together"); // deliberately missing end quote
		TS_ASSERT_EQUALS( v.size(), 2 );
		TS_ASSERT_EQUALS( v[1], "'quotes will'keep" );
		TS_ASSERT_EQUALS( v[2], "\"things together" );


	}
};
