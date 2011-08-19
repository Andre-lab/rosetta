// :noTabs=false:tabSize=4:indentSize=4:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>




#include <core/types.hh>
#include <devel/init.hh>

#include <basic/ComparingTracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


int main( int argc, char * argv [] )
{
	devel::init(argc, argv);

	/// if command line option -out:dry_run set to 'false' (default) output will be compare to file "core.tr"
	/// if -out:dry_run set to 'true' then no asserts generated and output will be save to "core.tr"
	basic::otstreamOP tr = new basic::ComparingTracer("core.tr");

	// Set to monitor space separated list of chanels (hierarchy of chanels works here)
	//basic::Tracer::set_ios_hook(tr, "core protocols");
	basic::Tracer::set_ios_hook(tr, basic::Tracer::AllChannels);

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose, "test_in.pdb");

	basic::Tracer::set_ios_hook(0, "");  // Turn off channels monitoring

	return 0;
}



