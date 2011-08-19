// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <devel/init.hh>
#include <core/types.hh>

#include <protocols/jobdist/standard_mains.hh>
// AUTO-REMOVED #include <protocols/rbsegment_Moves/RBSegmentMover.hh>
#include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/ProteinSilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>

// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <basic/Tracer.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/jobdist/Jobs.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <numeric/random/random_permutation.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <core/scoring/rms_util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <protocols/evaluation/RmsdEvaluator.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/chemical/ChemicalManager.hh>


// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>

#include <list>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/rbsegment_moves/RBSegment.hh>



int main(int argc, char **argv) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;

	using core::Size;
	using std::string;

	devel::init( argc,argv );

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose native_pose, current_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, *rsd_set, option[ in::file::native ]() );
	}

	// rbsegfile
	utility::vector1< protocols::RBSegment::RBSegment > rbsegs;
	protocols::loops::Loops loops;
	std::string filename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );
	protocols::RBSegment::read_RBSegment_file( rbsegs, loops, filename );
 	std::list< core::Size > core_reses;
	for ( core::Size i=1; i <= rbsegs.size(); ++i )
		for ( core::Size j=1, j_end=rbsegs[i].nContinuousSegments(); j<=j_end; ++j)
			for ( core::Size k=rbsegs[i][j].start(), k_end=rbsegs[i][j].end() ; k<=k_end; ++k)
				core_reses.push_back( k );



	MetaPoseInputStream input = streams_from_cmd_line();

	SilentFileData sfd_out;
	while( input.has_another_pose() ) {
		input.fill_pose( current_pose, *rsd_set );

		//
		if ( option[ in::file::native ].user() ) {
			core::Real CA_rmsd = core::scoring::CA_rmsd( native_pose, current_pose );
			core::Real CA_core_rmsd = core::scoring::CA_rmsd( native_pose, current_pose , core_reses );

			std::cout << protocols::jobdist::extract_tag_from_pose( current_pose ) << "  " << CA_core_rmsd << "  " << CA_rmsd << std::endl;
		}

	} // while( input.has_another_pose() )

	return 0;
}
