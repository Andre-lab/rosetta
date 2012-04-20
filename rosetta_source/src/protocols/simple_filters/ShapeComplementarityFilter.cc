// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/filters/ShapeComplementarityFilter.cc
/// @brief  Filter structures by shape complementarity and/or interface area
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

// Unit Headers
#include <protocols/simple_filters/ShapeComplementarityFilter.hh>
#include <protocols/simple_filters/ShapeComplementarityFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static basic::Tracer tr("protocols.filters.ShapeComplementarityFilter");

namespace protocols {
namespace simple_filters {

// @brief default constructor
ShapeComplementarityFilter::ShapeComplementarityFilter():
	Filter( "ShapeComplementarity" ),
	filtered_sc_( 0.50 ),
	filtered_area_( 250 ),
	jump_id_( 1 ),
	quick_( false ),
	verbose_( false ),
	residues1_( ),
	residues2_( )
{}


// @brief constructor with arguments
ShapeComplementarityFilter::ShapeComplementarityFilter( Real const & filtered_sc, Real const & filtered_area,
	Size const & jump_id, Size const & quick, Size const & verbose):
	Filter( "ShapeComplementarity" ),
	filtered_sc_( filtered_sc ),
	filtered_area_( filtered_area ),
	jump_id_( jump_id ),
	quick_( quick ),
	verbose_( verbose ),
	residues1_( ),
	residues2_( )
{}

// @brief copy constructor
ShapeComplementarityFilter::ShapeComplementarityFilter( ShapeComplementarityFilter const & rval ):
	Super( rval ),
	filtered_sc_( rval.filtered_sc_ ),
	filtered_area_( rval.filtered_area_ ),
	jump_id_( rval.jump_id_ ),
	quick_( rval.quick_ ),
	verbose_( rval.verbose_ ),
	residues1_( rval.residues1_ ),
	residues2_( rval.residues2_ )
{}

void ShapeComplementarityFilter::filtered_sc( Real const & filtered_sc ) { filtered_sc_ = filtered_sc; }
void ShapeComplementarityFilter::filtered_area( Real const & filtered_area ) { filtered_area_ = filtered_area; }
void ShapeComplementarityFilter::jump_id( Size const & jump_id ) { jump_id_ = jump_id; }
void ShapeComplementarityFilter::quick( Size const & quick ) { quick_ = quick; }
void ShapeComplementarityFilter::verbose( Size const & verbose ) { verbose_ = verbose; }

/// @brief
core::Size ShapeComplementarityFilter::compute( Pose const & pose ) const
{
	if(scc_.GetResults().valid)
		return 1;

	if(!scc_.Init())
		return 0;
	if(quick_)
		scc_.settings.density = 5.0;
	scc_.Reset();

	if(!residues1_.empty() && !residues2_.empty()) {
		for(utility::vector1<Size>::const_iterator r = residues1_.begin();
			r != residues1_.end(); ++r)
				scc_.AddResidue(0, pose.residue(*r));

		for(utility::vector1<Size>::const_iterator r = residues2_.begin();
			r != residues2_.end(); ++r)
				scc_.AddResidue(1, pose.residue(*r));

		if(!scc_.Calc())
			return 0;

	} else {
		int sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num(pose, jump_id_);
		if(!scc_.Calc( pose, sym_aware_jump_id ))
			return 0;
	}

	core::scoring::sc::RESULTS const &r = scc_.GetResults();
	if(verbose_) {

		// Verbose view
		tr << "==================================================" << std::endl;
		tr << std::endl;
		for(int i = 0; i <= 2; i++) {
			if(i < 2)
				tr << "Molecule " << (i+1) << ":" << std::endl;
			else
				tr << "Total/Average for both molecules:" << std::endl;

			tr << "          Total Atoms: " << r.surface[i].nAtoms << std::endl;
			tr << "         Buried Atoms: " << r.surface[i].nBuriedAtoms << std::endl;
			tr << "        Blocked Atoms: " << r.surface[i].nBlockedAtoms << std::endl;
			tr << "           Total Dots: " << r.surface[i].nAllDots << std::endl;
			tr << " Trimmed Surface Dots: " << r.surface[i].nTrimmedDots << std::endl;
			tr << "         Trimmed Area: " << r.surface[i].trimmedArea << " (avg) " << std::endl;
			tr << std::endl;
    }
		tr << std::endl;

		for(int i = 0; i <= 2; i++) {
			if(i < 2)
				tr << "Molecule " << (i+1) << "->" << ((i+1)%2+1) << ": " << std::endl;
			else
				tr << "Average for both molecules:" << std::endl;
			tr << "      Mean Separation: " << r.surface[i].d_mean << std::endl;
			tr << "    Median Separation: " << r.surface[i].d_median << std::endl;
			tr << "    Mean Shape Compl.: " << r.surface[i].s_mean << std::endl;
			tr << "  Median Shape Compl.: " << r.surface[i].s_median << std::endl;
			tr << std::endl;
		}

	}

	tr << "Shape complementarity: " << r.sc << std::endl;
	tr << "Interface area: " << r.area << std::endl;
	tr << "Interface seperation: " << r.distance << std::endl;

	return 1;
}

/// @brief
core::Real ShapeComplementarityFilter::report_sm( Pose const & pose ) const
{
	if(compute( pose ))
		return scc_.GetResults().sc;
	return -1;
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose has high enough shape
// complementarity.
bool ShapeComplementarityFilter::apply( Pose const & pose ) const
{
	scc_.Reset();

	if(!compute( pose ))
		return false;

	Real sc = scc_.GetResults().sc;
	Real area = scc_.GetResults().area;

	if( sc < filtered_sc_ ) {
		tr << "Filter failed current < threshold sc: " << sc << " < " << filtered_sc_ << std::endl;
		return false;
	}

	if( area < filtered_area_ ) {
		tr << "Filter failed current < threshold interface area: " << area << " < " << filtered_area_ << std::endl;
		return false;
	}

	tr << "Successfully filtered: " << sc << std::endl;
	return true;
} // apply_filter

/// @brief parse xml
void
ShapeComplementarityFilter::parse_my_tag(
	TagPtr const tag,
	DataMap &,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & pose )
{
	filtered_sc_ = tag->getOption<Real>( "min_sc", 0.50 );
	filtered_area_ = tag->getOption<Real>( "min_interface", 0 );
	verbose_ = tag->getOption<Size>( "verbose", false );
	quick_ = tag->getOption<Size>( "quick", false );
	jump_id_ = tag->getOption<Size>( "jump", 1 );

	if(tag->hasOption("residues1")) {
		residues1_ = protocols::rosetta_scripts::get_resnum_list(tag, "residues1", pose);
		if(residues1_.empty())
			tr.Warning << "Failed to parse residue range: " << tag->getOption<std::string> ("residues1") << ". Using default." << std::endl;
	}
	if(tag->hasOption("residues2")) {
		residues2_ = protocols::rosetta_scripts::get_resnum_list(tag, "residues2", pose);
		if(residues2_.empty())
			tr.Warning << "Failed to parse residue range: " << tag->getOption<std::string> ("residues2") << ". Using default." << std::endl;
	}

	tr.Info << "Structures with shape complementarity < " << filtered_sc_ << ", interface area < " <<
		filtered_area_ << " A^2 will be filtered." << std::endl;

	if(quick_)
		tr.Info << "Calculating shape complementarity in quick mode with less accuracy." << std::endl;
	if(!residues1_.empty() && !residues2_.empty()) {
		tr.Info << "Using residues for molecule surface (rosetta numbering):" << std::endl;
		tr.Info << "  Surface 1: ";
		for(utility::vector1<Size>::const_iterator r = residues1_.begin(); r != residues1_.end(); ++r)
	                tr.Info << (r == residues1_.begin() ? "" : ", ") << *r;
		tr.Info << std::endl;
		tr.Info << "  Surface 2: ";
		for(utility::vector1<Size>::const_iterator r = residues2_.begin(); r != residues2_.end(); ++r)
	                tr.Info << (r == residues2_.begin() ? "" : ", ") << *r;
		tr.Info << std::endl;
	} else {
		if(!residues1_.empty() || !residues2_.empty())
			tr.Warning << "Ignoring residue range selection since residues" << (residues1_.empty() ? 1 : 2) << " is empty." << std::endl;
		if(jump_id_ != 1)
			tr.Info << "Using Jump ID " << jump_id_ << " to define surfaces." << std::endl;
	}
}

filters::FilterOP
ShapeComplementarityFilterCreator::create_filter() const { return new ShapeComplementarityFilter; }

std::string
ShapeComplementarityFilterCreator::keyname() const { return "ShapeComplementarity"; }


} // filters
} // protocols
