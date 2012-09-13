// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  ScoreFunction class definition. A easy way to deal with the low-res score function problem in replica docking
/// @author Zhe Zhang


// Unit headers
#include <core/scoring/MinScoreScoreFunction.hh>

// Package headers

// // Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <basic/Tracer.hh>

static basic::Tracer tr("core.scoring.MinScoreScoreFunction");

namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////
MinScoreScoreFunction::MinScoreScoreFunction():
	ScoreFunction()
	{}

///////////////////////////////////////////////////////////////////////////////
ScoreFunctionOP
MinScoreScoreFunction::clone() const
{
	return new MinScoreScoreFunction( *this );
}

///////////////////////////////////////////////////////////////////////////////
MinScoreScoreFunction &
MinScoreScoreFunction::operator=( MinScoreScoreFunction const & src )
{
	ScoreFunction::operator=( src );
	return *this;
}

///////////////////////////////////////////////////////////////////////////////
MinScoreScoreFunction::MinScoreScoreFunction( MinScoreScoreFunction const & src ):
	ScoreFunction( src ),
	min_score_( src.min_score_ )
{}

MinScoreScoreFunction::MinScoreScoreFunction( ScoreFunction const & src, core::Real min_score ):
  ScoreFunction( src )
{
	min_score_ = min_score;
}

MinScoreScoreFunction::MinScoreScoreFunction( ScoreFunctionOP src, core::Real const min_score ):
  ScoreFunction( *src )
{
	min_score_ = min_score;
}

MinScoreScoreFunction::MinScoreScoreFunction( ScoreFunctionCOP src, core::Real const min_score ):
  ScoreFunction( *src )
{
	min_score_ = min_score;
}

///////////////////////////////////////////////////////////////////////////////

// to start out, just thinking fullatom energies
//
// NOTE: no freakin rotamer trials inside scoring!
Real
MinScoreScoreFunction::operator()( pose::Pose & pose ) const
{
 	ScoreFunction::operator()( pose ); //score -- but without atom_pair_constraints..
	// is probably cheaper to not apply a completely new scorefunction...
	EnergyMap cst_free_weights( weights() );
	Real cst_weight = cst_free_weights.get( atom_pair_constraint );
	cst_free_weights[ atom_pair_constraint ] = 0;
	Real uncst_energy = pose.energies().total_energies().dot( cst_free_weights );
	Real min_energy = uncst_energy < min_score_ ? min_score_ : uncst_energy;
	tr.Debug << "uncst_energy: " << uncst_energy << " min_energy: " << min_energy << std::endl;
	pose.energies().total_energies()[ total_score ] = min_energy + pose.energies().total_energies()[ atom_pair_constraint ]*cst_weight;
	pose::setPoseExtraScores( pose, "min_score", uncst_energy );
	return pose.energies().total_energies()[ total_score ];
}


///////////////////////////////////////////////////////////////////////////////
} // namespace scoring
} // namespace core
