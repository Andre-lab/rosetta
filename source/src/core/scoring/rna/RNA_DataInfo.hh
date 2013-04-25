// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_DataPotential.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_rna_RNA_DataInfo_hh
#define INCLUDED_core_scoring_rna_RNA_DataInfo_hh

#include <core/types.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>

// C++

namespace core {
namespace scoring {
namespace rna {

class RNA_Datum {
public:

	RNA_Datum( Size const position, Size const edge, Real const weight){
		position_ = position;
		edge_ = edge;
		weight_ = weight;
	}

	Size position() const { return position_; };
	Size edge() const { return edge_; };
	Real weight() const { return weight_; };

private:
	Size position_;
	Size edge_;
	Real weight_;
};

typedef utility::vector1< RNA_Datum > RNA_Data;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
class RNA_DataInfo : public basic::datacache::CacheableData {

public:

	RNA_DataInfo(){}

  RNA_DataInfo( RNA_DataInfo const & src );

  basic::datacache::CacheableDataOP
  clone() const
  {
    return new RNA_DataInfo( *this );
  }

	RNA_DataInfo &
	operator=( RNA_DataInfo const & src );

  Size
  size() const {
    return rna_data_.size();
  }

  void
  zero();

	RNA_Data const & rna_data() const { return rna_data_; }

	void
	add_datum( RNA_Datum const & rna_datum ){ rna_data_.push_back( rna_datum ); }

	ObjexxFCL::FArray1D< bool > const & backbone_burial() const { return backbone_burial_; }

	void set_backbone_burial( ObjexxFCL::FArray1D< bool > const & backbone_burial ){ backbone_burial_ = backbone_burial; }

	ObjexxFCL::FArray1D< bool > const & backbone_exposed() const { return backbone_exposed_; }

	void set_backbone_exposed( ObjexxFCL::FArray1D< bool > const & backbone_exposed ){ backbone_exposed_ = backbone_exposed; }

private:

	RNA_Data rna_data_;
	ObjexxFCL::FArray1D< bool > backbone_burial_;
	ObjexxFCL::FArray1D< bool > backbone_exposed_;

};

}
}
}

#endif
