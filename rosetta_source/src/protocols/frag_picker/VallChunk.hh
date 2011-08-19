// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/VallChunk.hh
/// @brief  a contiguous chunk of residues taken from a vall.
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_protocols_frag_picker_VallChunk_hh
#define INCLUDED_protocols_frag_picker_VallChunk_hh

// unit headers
#include <protocols/frag_picker/VallChunk.fwd.hh>

// package headers
#include <protocols/frag_picker/VallProvider.fwd.hh>
#include <protocols/frag_picker/VallResidue.hh>

// mini
#include <core/sequence/SequenceProfile.fwd.hh>

//Auto Headers
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace frag_picker {

/// @brief  represents a chunk of residues extracted from a vall.
/// @detailed VallChunk contains a vector of VallResidue objects and provides a basic ways to access them
class VallChunk: public utility::pointer::ReferenceCount {
public:

	VallChunk(VallProviderAP provider);

	/// @brief  returns a PDB id (a string of four letters, e.g. "4mba")
	inline std::string get_pdb_id() const {
		return at(1)->id().substr(0, 4);
	}

	/// @brief  returns protein chain ID
	inline char get_chain_id() const {
		return at(1)->id()[4];
	}

	/// @brief  returns integer key of this chunk, which is the key of this chunk's first residue
	inline Size key() const {
		return at(1)->key();
	}

	/// @brief  returns the size of this chunk i.e. the number of residues stored in there
	inline Size size() const {
		return residues_.size();
	}

	/// @brief  returns i-th residue form this chunk. The first residue has index 1
	inline VallResidueOP at(Size index) const {
		return residues_.at(index);
	}

	/// @brief  appends a residue to this chunk
	inline void push_back(VallResidueOP what) {
		residues_.push_back(what);
	}

	/// @brief  returns amino acid sequence of this chunk
	std::string& get_sequence();

	/// @brief  returns amino acid profile of this chunk
	/// @detailed the profile object is created when this function is called for the first time
	/// and then cached within a VallProvider object.
	/// Every time this method is called for a new chunk, VallProvider caches new data
	core::sequence::SequenceProfileOP get_profile();

	/// @brief  returns a pose created for this chunk
	/// @detailed the pose object is created when this function is called for the first time
	/// and then cached within a VallProvider object
	/// Every time this method is called for a new chunk, VallProvider caches new data
	core::pose::PoseOP get_pose();

	/// @brief  returns a string that is unique for each chunk in vall
	std::string & chunk_key() { if(!has_key_) create_key(); return chunk_key_; }
private:
	utility::vector1<VallResidueOP> residues_;
	std::string sequence_;
	VallProviderAP my_provider_;
	std::string chunk_key_;
	bool has_key_;

	void create_key();
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_VallChunk_HH */
