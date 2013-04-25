// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/conformation/ResidueFactory.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_ScoreFunctionFactory_hh
#define INCLUDED_core_scoring_ScoreFunctionFactory_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <string>

#ifdef WIN32 //VC++ needs full class declaration
#include <core/scoring/ScoreFunction.hh> // WIN32 INCLUDE
#endif
// #include <core/chemical/ResidueType.fwd.hh>
// #include <core/conformation/Conformation.fwd.hh>
// #include <core/conformation/Residue.fwd.hh>
// #include <core/chemical/AtomTypeSet.fwd.hh>
// #include <core/chemical/MMAtomTypeSet.fwd.hh>

namespace core {
namespace scoring {

/// @brief a collection of functions making a single score_function
class ScoreFunctionFactory
{
public:

	typedef core::scoring::symmetry::SymmetricScoreFunction SymmetricScoreFunction;

	/// @brief Returns a ScoreFunction from the database weights file  <weights_tag>
	///
	/// example(s):
	///     scorefxn = create_score_function('standard')
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.show
	///     ScoreFunction.weights
	///     ScoreType
	static
	ScoreFunctionOP
	create_score_function( std::string weights_tag );

	/// @brief Returns a ScoreFunction from the database weights file  <weights_tag>
	/// with the patch <patch_tag>
	///
	/// example(s):
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.show
	///     ScoreFunction.weights
	///     ScoreType
	static
	ScoreFunctionOP
	create_score_function( std::string weights_tag, std::string const & patch_tag );

	/// @brief Returns a ScoreFunction from the database weights file  <weights_tag>  with patches in <patch_tags>
	static
	ScoreFunctionOP
	create_score_function( std::string weights_tag, utility::vector1< std::string > patch_tags );

private:

	/// @brief Applies user defined re-weighting from the options system. Reweights are applied as a
	/// factor of the original, so -rg_reweight 0.5 would result in half of the previously defined
	/// rg weight.
	static void apply_user_defined_reweighting_( core::scoring::ScoreFunctionOP scorefxn );

	static void load_weights_file( std::string weights_tag, ScoreFunctionOP scorefxn );

};

extern std::string const STANDARD_WTS;
extern std::string const SCORE13;
extern std::string const CENTROID_WTS;
extern std::string const SOFT_REP_WTS;
extern std::string const SOFT_REP_DESIGN_WTS;
extern std::string const DNA_INT_WTS;
extern std::string const DNA_INT_WTS_GB;
extern std::string const MM_STD_WTS;
extern std::string const RNA_LORES_WTS;
extern std::string const RNA_HIRES_WTS;
extern std::string const RNA_LORES_PLUS_HIRES_WTS;
extern std::string const MEMB_HIGHRES_WTS; //pba

extern std::string const SCORE12_PATCH;
extern std::string const DOCK_PATCH;
extern std::string const DOCK_LOW_PATCH;



/// @brief A helper function which returns a scoring function owning pointer according to the
/// user's command line parameters -score:weights and -score:patch
/// By default it returns weights=standard and patch=score12 for fullatom,
/// and weights=cen_std and patch="" for centroid
core::scoring::ScoreFunctionOP getScoreFunction( bool const is_fullatom = true );

/// @brief use the logic of getScoreFunction to get the name.
/// The 	name format is <weights_tag>[_<patch_tag> ... ]
std::string
getScoreFunctionName(
	bool const is_fullatom = true );

} // namespace scoring
} // namespace core

#endif
