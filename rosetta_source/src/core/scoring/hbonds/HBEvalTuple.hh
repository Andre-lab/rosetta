// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/hbonds/HBEvalTuple.hh
/// @brief  Tuple describing data about the donor and acceptor in a single hbond
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_hbonds_HBEvalTuple_HH
#define INCLUDED_core_scoring_hbonds_HBEvalTuple_HH

// Unit headers
#include <core/scoring/hbonds/HBEvalTuple.fwd.hh>

// Package headers
#include <core/scoring/hbonds/types.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>

namespace core {
namespace scoring {
namespace hbonds {

class HBEvalTuple 
{
private:
	HBDonChemType don_type_;
	HBAccChemType acc_type_;
	HBSeqSep      seq_sep_;
	HBEvalType    eval_type_;

public:
	HBEvalTuple() :
		don_type_( hbdon_NONE ),
		acc_type_( hbacc_NONE ),
		seq_sep_( seq_sep_other ),
		eval_type_( hbe_NONE )
	{}
	
	HBEvalTuple(
		int const datm,
		core::conformation::Residue const & don_rsd,
		int const aatm,
		core::conformation::Residue const & acc_rsd
	);
	
	HBEvalTuple(
		HBDonChemType don,
		HBAccChemType acc,
		HBSeqSep sequence_sep
	);

	HBEvalTuple( HBEvalTuple const & src ) :
		don_type_( src.don_type_ ),
		acc_type_( src.acc_type_ ),
		seq_sep_( src.seq_sep_ ),
		eval_type_( src.eval_type_ )
	{}
	
	HBEvalTuple const & operator = ( HBEvalTuple const & rhs ) {
		if ( this != & rhs ) {
			don_type_ = rhs.don_type_;
			acc_type_ = rhs.acc_type_;
			seq_sep_ = rhs.seq_sep_;
			eval_type_ = rhs.eval_type_;
		}
		return *this;
	}

	inline ~HBEvalTuple() {}

	void don_type( HBDonChemType don );
	void acc_type( HBAccChemType acc );
	void sequence_sep( HBSeqSep seqsep );

	inline HBDonChemType don_type() const { return don_type_; }
	inline HBAccChemType acc_type() const { return acc_type_; }
	inline HBSeqSep sequence_sep() const { return seq_sep_; }
	inline HBEvalType eval_type() const { return eval_type_; }	


private:
	void update_hbevaltype();

};

} // namespace hbonds
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_hbonds_HBEvalTuple_HH
