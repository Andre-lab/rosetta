// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief Rescore membrane protein test
/// @author Bjorn Wallner

// libRosetta headers


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <ObjexxFCL/formatted.o.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

//Auto Headers
#include <core/import_pose/import_pose.hh>


int
main( int argc, char* argv [] )
{
	try {
	// options, random initialization
	devel::init( argc, argv );

	//	using namespace core;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;
	//using namespace core::chemical;

	// setup residue types
	//		core::chemical::ResidueTypeSetCAP rsd_set;
	//if ( option[ in::file::fullatom ]() ) {
	//	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	//} else {
	//	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	//	rsd_set=core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	//	}

	// configure score function
	//core::scoring::ScoreFunctionOP scorefxn;
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("score_membrane");
	//scorefxn->set_weight(core::scoring::Mlipo,1.0);
	// configure silent-file data object
	core::io::silent::SilentFileData sfd;
	core::pose::Pose pose;
	std::string infile  = *(option[ in::file::silent ]().begin());
	std::string const spanfile = option[ in::file::spanfile ]();
	std::string outfile = option[ out::file::silent ]();
	utility::io::ozstream output;
	output.open( outfile.c_str() );
	std::cout << "spanfile: " << spanfile << "\n";
	//	core::scoring::MembraneTopologyOP topology=new core:scoring:MembraneTopology;

	core::scoring::MembraneTopologyOP topologyOP = new core::scoring::MembraneTopology;
	pose.data().set( MEMBRANE_TOPOLOGY, topologyOP );
	//	core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( pose.data().get_ptr( basic::MEMBRANE_TOPOLOGY )() ));
	core::scoring::MembraneTopology & topology=*( static_cast< core::scoring::MembraneTopology * >( pose.data().get_ptr( MEMBRANE_TOPOLOGY )() ));
	topology.initialize(spanfile);


	//	core::pose::metrics::PoseMetricCalculatorOP center_normal= new protocols::toolbox::pose_metric_calculators::MembraneCenterNormal;
	//core::pose::metrics::CalculatorFactory::Instance().register_calculator("MCN",center_normal);
	//basic::MetricValue<core::Vector> center;

	//std::cout << pose.print_metric("MembraneCenterNormal","center");

	//	topology.attach_to_pose(pose);

	if ( option[ in::file::silent ].user() ) {
		sfd.read_file( infile );
	}
	utility::vector1< std::string > tag_list;
	tag_list=sfd.tags();

	//core::Size ntimes = option[ in::file::fragA_size ]();

	//	std::exit(1);
	//	for ( core::Size i = 1; i <= ntimes; ++i ) {
		clock_t start_time = clock();
		core::Size nscores = 0;
		//core::Real weight = 1 / ( (double) sfd.size() * ntimes );

		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end();
					iter != end;
					++iter
		) {
			std::string tag(iter->decoy_tag());
			//	core::Real memb_env=iter->get_energy("MEMB_ENV");
			//core::Real memb_pair=iter->get_energy("MEMB_PAI");
			//core::Real memb_cb=iter->get_energy("MEMB_CB");
			//iter->fill_pose( pose,*rsd_set) //, *rsd_set );
			//core::import_pose::pose_from_pdb(pose,*rsd_set,"S_49.pdb.mini");
				iter->fill_pose( pose,
							 option[ in::file::fullatom ] ?
												 *(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )) :
												 *(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ) ));


			pose.data().set( MEMBRANE_TOPOLOGY, topologyOP );
			//			pose.metric("MCN","center",center);

			(*scorefxn)(pose);
			//			scorefxn->show(std::cout,pose);
			//			std::cout << A( 15, "S CORE: TOTAL");
			//scorefxn->show_line_headers(std::cout);
	/*		if(nscores==0) {
				std::cout << "SCORE : TOTAL ";
				scorefxn->show_line_headers(std::cout);
				std::cout << "\n";
			}
			std::cout << "SCORE: ";
			scorefxn->show_line(std::cout,pose);
			std::cout << " " << tag << "\n";
	 */
			//std::cout << "center: " << center.value();
			//	std::cout << "\n";
			//std::cout << "SCORE OLD: MEM_ENV " << iter->get_energy("MEMB_ENV") << " MEMB_CB " << iter->get_energy("MEMB_CB") << " MEMB_PAIR " << iter->get_energy("MEMB_PAI") << " " << tag << "\n";


			//std::cout << "SCORE NEW:" << (*scorefxn)(pose) << " " << tag << "\n" ;
			output << "SCORE: TOTAL\t" ; //A( 15, "SCORE: TOTAL"); //"SCORE: TOTAL ";
			scorefxn->show_line_headers(output);
			output << "tag" << "\n"; //A( 15, "tag") << "\n";
			output << "SCORE: "; //A( 15, "SCORE NEW: ");
			scorefxn->show_line(output,pose);
			output << " " << tag << "\n"; //A(15,tag) << "\n";
			//output << "SCORE OLD: MEM_ENV " << iter->get_energy("MEMB_ENV") << " MEMB_CB " << iter->get_energy("MEMB_CB") << " MEMB_PAIR " << iter->get_energy("MEMB_PAI") << " rms " << iter->get_energy("rms") << " MAXSUB " << iter->get_energy("MAXSUB") << " " << tag << "\n";
			++nscores;
			if(nscores % 100 == 0)
			{
				std::cout << nscores << "\n";
			}
			//std::cout << scorefxn->get_energy("Menv") << "\n";
			if( option[out::pdb].user() )
			{
				std::cout << "Outputting "<< tag << "\n";
				pose.dump_pdb(tag+".pdb");
			}

			pose.energies().clear();
		} // for sfd
		clock_t end_time = clock();
		//clock_t avg_time = (start_time - endtime) / nscores;
		std::cout << "Total time: "
							<< (double(end_time) - double(start_time)) / ( CLOCKS_PER_SEC )
							<< " seconds for " << nscores << " score operations." << std::endl;

		//	} // for ntimes
		output.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
