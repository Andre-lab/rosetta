// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/IdealAbegoGenerator.hh
/// @brief Logic for selection of abego values
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_components_IdealAbegoGenerator_hh
#define INCLUDED_protocols_denovo_design_components_IdealAbegoGenerator_hh

// Unit headers
#include <protocols/denovo_design/components/IdealAbegoGenerator.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>
#include <map>
#include <string>

namespace protocols {
namespace denovo_design {
namespace components {

///@brief Logic for selection of abego values
class IdealAbegoGenerator : public utility::VirtualBase {
public:
	typedef Segment Motif;
	typedef SegmentOP MotifOP;
	typedef utility::vector1< MotifOP > MotifOPs;
	typedef std::set< core::Size > LengthSet;

	typedef std::string Abego;
	typedef utility::vector1< std::string > Abegos;
	typedef std::map< std::string, Abegos > SecstructAbegoMap;

public:
	IdealAbegoGenerator( std::string const & id_val );
	~IdealAbegoGenerator() override;

	IdealAbegoGeneratorOP
	clone() const;

	/// @brief Given desired lengths, compute a set of idealized loop motifs via Nobu/Rie/YuRu rules
	/// @param[in] abego1 Abego of residue immediately before the loop
	/// @param[in] abego2 Abego of residue immediately after the loop
	/// @param[in] lenset Set of allowed loop lengths
	/// @param[in] cutpoint_set Set of allowed cutpoint residues.  A value of N indicates the Nth residue
	///                         of the loop is a LOWER_CUTPOINT and the (N+1)th residue is an UPPER_CUTPOINT
	/// @returns MotifOPs of all ideal loops with given lengths/cutpoints connecting two residues of the given
	///          abegos
	/// @details If extend_ss is true, "extension" of SS elements by adding residues of the same abego type
	///          is allowed. For example, connecting A-->A, you might get a 4-residue loop
	///          of type A-ABGA-A. If false, only different abegos are allowed, and you might
	///          get a 4-residue loop of type A-BGAB-A
	MotifOPs
	generate(
		char abego1,
		char abego2,
		LengthSet const & lenset,
		LengthSet const & cutpoint_set ) const;

	/// @brief Returns true if the given SS type should be extended or not
	bool
	extend_ss( char const secstruct ) const;

	/// @brief if true, secondary structure may be extended to close the loop.  Overall size of the insert
	///        will not change.  Therefore, a 4 residue loop might extend a helix by 2 residues and have
	///        a 2-residue ideal loop if extend_ss is true.  If extend_ss is false, a 4 residue loop must
	///        have four loop residues.
	void
	set_extend_ss( std::string const & extend_ss );

	/// @brief if true, the connection rules from https://doi.org/10.1101/2021.03.10.434454 will be used.. If false, all possible abego combinations will be sampled
	void
	set_hh_rules_2021( bool const use_hh_rules ) { use_hh_rules_2021_ = use_hh_rules; }

	// public member functions
	Abegos
	retrieve_loop_abegos( char const abego1, char const abego2 ) const;

private:
	// private functions

	/// @brief takes a list of motifs without cutpoints and generates all permutations with cutpoints from cutpoint_set
	/// @param[in] orig         List of motifs without cutpoints
	/// @param[in] cutpoint_set Set of loop-relative indices of cutpoint residues
	/// @returns MotifOPs of all motifs in orig with all possible cutpoints in cutpoint_set
	MotifOPs
	add_cutpoints( MotifOPs const & orig, LengthSet const & cutpoint_set ) const;

	MotifOPs
	extract_ideal_motifs(
		Abegos const & abegolist,
		LengthSet const & lenset ) const;

	void
	note_fwd_element_extensions(
		Motif & motif,
		char const sschar,
		char const element_abego ) const;

	void
	note_rev_element_extensions(
		Motif & motif,
		char const sschar,
		char const element_abego ) const;

private:
	// private data
	std::string segment_name_;
	std::string extend_ss_ = "1";
	bool match_ss_to_abego_ = false;
	bool use_hh_rules_2021_ = true;
	SecstructAbegoMap extended_abegos_;

private:
	// private constants
	static Abegos const
		ab_no_extension;

	static Abegos const
		ba_no_extension;

	static Abegos const
		bb_no_extension;

	static Abegos const
		aa_no_extension;

	static Abegos const
		aa_no_extension_all;

	static std::map< char, core::Size > const
		extension_lengths;

	static Abegos
	generate_extended_abegos( Abegos const & no_extension_loops, char const abego_n, char const abego_c );

	// hide default constructor
private:
	IdealAbegoGenerator() {};

};

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_IdealAbegoGenerator_hh
