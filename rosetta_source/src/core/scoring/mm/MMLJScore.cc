// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMLJScore.cc
/// @brief  Molecular mechanics lj score class
/// @author P. Douglas Renfrew (renfrew@unc.edu)

// Unit headers
#include <core/scoring/mm/MMLJScore.hh>
#include <core/scoring/mm/MMLJLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/scoring/ScoringManager.hh>

// AUTO-REMOVED #include <basic/prof.hh>

// Utility header
// AUTO-REMOVED #include <utility/keys/Key4Tuple.hh>
// AUTO-REMOVED #include <utility/keys/Key3Tuple.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>
#include <string>
#include <map>
#include <math.h>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace mm {

MMLJScore::MMLJScore() :
  mm_lj_library_( scoring::ScoringManager::get_instance()->get_MMLJLibrary() )
{ }

MMLJScore::MMLJScore( MMLJLibrary const & mmljl ) :
  mm_lj_library_( mmljl )
{ }

MMLJScore::~MMLJScore() {}

/// @details blah
Energy
MMLJScore::score( Size atom1, Size atom2, Size path_distance, Real distance ) const
{
  // lookup params
  mm_lj_param_set atom1_params, atom2_params;
  if( path_distance == 3 )
    {
      atom1_params = mm_lj_library_.lookup_three_bond( atom1 );
      atom2_params = mm_lj_library_.lookup_three_bond( atom2 );
    }
  else
    {
      atom1_params = mm_lj_library_.lookup( atom1 );
      atom2_params = mm_lj_library_.lookup( atom2 );
    }

  // calc score
  Real epsilon( sqrt( atom1_params.key2() * atom2_params.key2() ) );
  Real rmin_over_dist( ( atom1_params.key1() + atom2_params.key1() ) / distance );
  Real score( epsilon * ( pow( rmin_over_dist, 12 ) - 2 * pow( rmin_over_dist, 6 ) ) );

	return score;
}

/// @details blah
Energy
MMLJScore::deriv_score( Size atom1, Size atom2, Size path_distance, Real distance ) const
{
  // lookup params
  mm_lj_param_set atom1_params, atom2_params;
  if( path_distance == 3 )
    {
      atom1_params = mm_lj_library_.lookup_three_bond( atom1 );
      atom2_params = mm_lj_library_.lookup_three_bond( atom2 );
    }
  else
    {
      atom1_params = mm_lj_library_.lookup( atom1 );
      atom2_params = mm_lj_library_.lookup( atom2 );
    }

	// calc deriv
	Real epsilon( sqrt( atom1_params.key2() * atom2_params.key2() ) );
	Real rmin ( atom1_params.key1() + atom2_params.key1() );
	Real deriv( 12 * epsilon * ( ( pow( rmin, 6 ) / pow( distance, 7 ) ) - ( pow( rmin, 12 ) / pow( distance, 13 ) ) ) );

 	return deriv;
}

/// @derails blah
Real
MMLJScore::min_dist( Size atom1, Size atom2, Size path_distance ) const
{
  // lookup params
  mm_lj_param_set atom1_params, atom2_params;
  if( path_distance == 3 )
    {
      atom1_params = mm_lj_library_.lookup_three_bond( atom1 );
      atom2_params = mm_lj_library_.lookup_three_bond( atom2 );
    }
  else
    {
      atom1_params = mm_lj_library_.lookup( atom1 );
      atom2_params = mm_lj_library_.lookup( atom2 );
    }

	// calc min
	Real rmin ( atom1_params.key1() + atom2_params.key1() );

	return rmin;
}

} // namespace mm
} // namespace scoring
} // namespace core
