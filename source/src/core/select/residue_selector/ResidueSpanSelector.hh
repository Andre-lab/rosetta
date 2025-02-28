// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueSpanSelector.hh
/// @brief  The ResidueSpanSelector selects residues using a span designation
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_select_residue_selector_ResidueSpanSelector_HH
#define INCLUDED_core_select_residue_selector_ResidueSpanSelector_HH

// Unit headers
#include <core/select/residue_selector/ResidueSpanSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/ResidueIndexDescription.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

/// @brief The ResidueSpanSelector returns a ResidueSubset, i.e. a utility::vector1< bool > containing
/// 'true' for residue positions which are within a given span.
/// A span is defined by start and end point, each as
/// either Rosetta indices (e.g. 10) or PDB numbers (e.g. 10A, residue 10 of chain A).
/// Detection and mapping from PDB to Rosetta residue numbers is done internally.
/// (The span is defined based on Rosetta indicies, even if residues are given by PDB numbers.)
class ResidueSpanSelector : public ResidueSelector {
public:
	// derived from base class
	ResidueSpanSelector();

	/// @brief Copy constructor
	///
	ResidueSpanSelector( ResidueSpanSelector const & src);

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the new object.
	ResidueSelectorOP clone() const override;

	ResidueSpanSelector(
		core::pose::ResidueIndexDescriptionCOP start,
		core::pose::ResidueIndexDescriptionCOP end);

	ResidueSpanSelector( std::string const & start_str, std::string const & end_str );
	ResidueSpanSelector( core::Size start, core::Size end );
	~ResidueSpanSelector() override;

	ResidueSubset apply( core::pose::Pose const & pose ) const override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	get_name() const override;

	static std::string class_name();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void set_span( std::string const & start, std::string const & end);

private: // data members
	core::pose::ResidueIndexDescriptionCOP start_;
	core::pose::ResidueIndexDescriptionCOP end_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_ResidueSpanSelector )
#endif // SERIALIZATION


#endif
