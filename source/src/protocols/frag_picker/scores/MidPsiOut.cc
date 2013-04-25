// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite && is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions && developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.cc
/// @brief  Object that scores a fragment by root mean square deviation of Psi && Psi dihedrals
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/MidPsiOut.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
//#include <protocols/frag_picker/PsiPsiTalosIO.hh>
// AUTO-REMOVED #include <protocols/frag_picker/scores/FragmentCrmsd.hh>
// AUTO-REMOVED #include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
// AUTO-REMOVED #include <protocols/frag_picker/FragmentPicker.hh>

// option key includes
// AUTO-REMOVED #include <core/init.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>

// AUTO-REMOVED #include <protocols/Protocol.hh>

// AUTO-REMOVED #include <numeric/model_quality/rms.hh>

// utils
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <basic/prof.hh>

// C++
// AUTO-REMOVED #include <math.h>
// AUTO-REMOVED #include <stdio.h>
// AUTO-REMOVED #include <string.h>

// Boost
// AUTO-REMOVED #include <boost/algorithm/string.hpp>
// AUTO-REMOVED #include <boost/tuple/tuple.hpp>

#include <utility/vector1.hh>

#ifdef WIN32
	#include <protocols/frag_picker/FragmentPicker.hh>
#endif

namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

MidPsiOut::MidPsiOut(Size priority, Real lowest_acceptable_value, bool use_lowest) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "MidPsiOut") {
}

void MidPsiOut::do_caching(VallChunkOP current_chunk) {

	chunk_psi_.redimension(current_chunk->size());
	for (Size i = 1; i <= current_chunk->size(); ++i) {
		VallResidueOP r = current_chunk->at(i);
		chunk_psi_(i) = r->psi();
		//chunk_psi_(i) = r->psi();
	}

}

bool MidPsiOut::score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool MidPsiOut::cached_score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {

	std::string tmp = fragment->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0) {
		do_caching(fragment->get_chunk());
		cached_scores_id_ = tmp;
	}

	Size offset_v = fragment->get_first_index_in_vall() - 1;

	Real r_psi = static_cast< Real > ( fragment->get_length() );
	r_psi = (r_psi/2)+0.5+offset_v;
	Size i_psi = static_cast< Size > ( r_psi );

	Real psi = chunk_psi_(i_psi);

	scores->set_score_component( psi, id_);

	return true;
}

void MidPsiOut::clean_up() {
}

/// @brief Creates a MidPsiOut scoring method
/// @param priority - priority of the scoring method. The higher value the earlier the score
///		will be evaluated
/// @param lowest_acceptable_value - if a calculated score is higher than this value,
///		fragment will be neglected
/// @param FragmentPickerOP object - not used
/// @param line - the relevant line extracted from the scoring configuration file that defines this scoring method
/// 		It could look like: "MidPsiOut                140     -5.0     100.0 additional_string"
///		where 140, -5.0 && 100.0 are priority, weight && treshold, respectively.
///		The additional string may be:
///		- empty: then the maker tries to create a scoring object from a TALOS file
///			trying in::file::talos_psi_psi flag. If fails, will try to use a pose from in::file::s
///		- a pdb file, pdb extension is necessary. This will create a pose && steal Psi && Psi
///		- a TALOS file with Psi/Psi prediction (tab extension is necessary)
FragmentScoringMethodOP MakeMidPsiOut::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
		, std::string) {

	return (FragmentScoringMethodOP) new MidPsiOut(priority,
																								 lowest_acceptable_value, use_lowest);
}

} // scores
} // frag_picker
} // protocols
