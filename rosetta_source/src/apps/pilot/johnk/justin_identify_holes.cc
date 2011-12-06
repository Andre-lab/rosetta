// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/database/open.hh>
//#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//OPT_KEY( String, ref_decoy )

static basic::Tracer TR( "apps.pilot.justin_identify_holes.main" );

/// General testing code
int
main( int argc, char * argv [] )
{

	//	NEW_OPT( ref_decoy, "the structure to compute RMSD and relative score to", "" );

	core::init(argc, argv);

	TR << "Starting to identify holes" << std::endl;

	//	std::string const ref_decoy_fname = option[ ref_decoy ];

	// Read holes params
	TR << "Reading input params" << std::endl;
	core::scoring::packing::HolesParams hp_resl,hp_dec,hp_dec15;
	hp_resl.read_data_file(basic::database::full_name("scoring/rosettaholes/resl.params"));
	hp_dec.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy25.params"));
	hp_dec15.read_data_file(basic::database::full_name("scoring/rosettaholes/decoy15.params"));

	// scoring function
	scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(STANDARD_WTS, SCORE12_PATCH) );

	// Open output file, generate the header line (save it for printing in the log later), print to file
	std::string outfname = "list_of_holes.out";
  utility::io::ozstream outstream;
	outstream.open(outfname, std::ios::out);

	outstream << "fname max_void_volume" << std::endl;

	for (core::Size f=1; f <= basic::options::start_files().size(); f++) {

		std::string const curr_decoy_fname = basic::options::start_files().at(f);
		TR << "Processing decoy " << curr_decoy_fname << std::endl;

		pose::Pose curr_pose;
		core::import_pose::pose_from_pdb( curr_pose, curr_decoy_fname );
		(*scorefxn)(curr_pose);

		core::scoring::packing::HolesResult holes_result = core::scoring::packing::compute_rosettaholes_score(curr_pose,hp_resl,hp_dec,hp_dec15);

		TR << "RosettaHoles: " << holes_result.score << " " << holes_result.resl_score << " " << holes_result.decoy_score << std::endl;
		TR << "Rosetta scoring: " << (*scorefxn)(curr_pose) << std::endl;

		std::string out_cav_fname = "cav_"+curr_decoy_fname;
		utility::io::ozstream out_cavpdb_stream;
		out_cavpdb_stream.open( out_cav_fname, std::ios::out);
		curr_pose.dump_pdb(out_cavpdb_stream);
		core::scoring::packstat::output_packstat_pdb( curr_pose, out_cavpdb_stream );

		//		core::Real max_cluster_volume = core::scoring::packstat::output_packstat_pdb( curr_pose, out_cavpdb_stream );
		//		outstream << curr_decoy_fname << ' ' << max_cluster_volume << std::endl;
		//		TR << "Found max void volume " << max_cluster_volume << std::endl;

	}

	TR << "Done identifying holes" << std::endl;

	outstream.close();
	outstream.clear();

	return 0;

}



