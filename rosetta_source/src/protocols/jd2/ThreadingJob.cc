// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/ThreadingJob.cc
/// @author James Thompson

#include <protocols/jd2/ThreadingJob.hh>
#include <protocols/jd2/ThreadingJob.fwd.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/id/SequenceMapping.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/loops/Loops.hh>

#include <basic/Tracer.hh>

#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

ThreadingJob::ThreadingJob(
	core::pose::PoseCOP template_pdb,
	core::sequence::SequenceAlignmentCOP alignment,
	std::string const & input_tag,
	core::Size nstruct_max
) :	InnerJob( template_pdb, input_tag, nstruct_max ) {
	alignment_ = alignment;  // the alignment from the input file

	//fpd  Make a transitve map from the template alignment to the template pdb seq
	//fpd  Accounts for missing density in the PDB
	static basic::Tracer tr("protocols.jd2.ThreadingJob");
	using namespace core::sequence;

	SequenceOP query_sequence(
		new Sequence(
			alignment->sequence( 1 )->ungapped_sequence(),
			alignment->sequence( 1 )->id(),
			alignment->sequence( 1 )->start()
		)
	);

	SequenceOP aligned_template(
		alignment->sequence( 2 )->clone()
	);

	SequenceOP t_align_seq(
		new Sequence(
			aligned_template->ungapped_sequence(),
			aligned_template->id(),
			aligned_template->start()
		)
	);


	SequenceOP t_pdb_seq(
		new Sequence (
			template_pdb->sequence(),
			alignment->sequence( 2 )->id(),
			1
		)
	);

	// construct an intermediate alignment of the sequence from the alignment
	// to the sequence in the PDB file.
	SWAligner sw_align;
	ScoringSchemeOP ss( new SimpleScoringScheme( 120, 0, -100, 0 ) );

	fasta2template_ = new core::sequence::SequenceAlignment( sw_align.align( t_align_seq, t_pdb_seq, ss ) );
	nres_template_  = aligned_template->ungapped_length();

	// std::cerr << "OLD" << std::endl;
	// std::cerr << *alignment << std::endl;
	// std::cerr << "NEW" << std::endl;
	// std::cerr << *alignment_ << std::endl;
}

///@brief returns the "standard" loop definition (as conservative as possible)
protocols::loops::Loops ThreadingJob::loops( core::Size nres ) const {
	using namespace protocols::loops;
	using namespace protocols::comparative_modeling;
	Loops loops = loops_from_transitive_alignments( nres, *alignment_, nres_template_, *fasta2template_, 0 );

	// remove loops that overlap stolen residues
	Loops valid_loops;
	for ( Size jj = 1; jj <= loops.size(); ++jj ) {
		bool valid(true);
		for ( Size ii = 1; ii <= extra_residues_to_steal_.size(); ++ii ) {
			Size const extra_res( extra_residues_to_steal_[ii] );
			if ( loops[jj].start() <= extra_res && extra_res <= loops[jj].stop() ) {
				valid = false;
			}
		}
		if ( valid ) valid_loops.add_loop( loops[jj] );
	}

	return valid_loops;
}

utility::vector1< core::Size > const & ThreadingJob::extra_residues_to_steal() const {
	return extra_residues_to_steal_;
}

void ThreadingJob::extra_residues_to_steal( utility::vector1< core::Size > const & res ) {
	extra_residues_to_steal_ = res;
}

} // namespace jd2
} // namespace protocols
