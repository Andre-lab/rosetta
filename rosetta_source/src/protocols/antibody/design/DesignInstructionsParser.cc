// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/DesignInstructionsParser.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/DesignInstructionsParser.hh>

#include <protocols/antibody/design/AntibodyGraftDesigner.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/CDRClusterEnumManager.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <utility/string_util.hh>
#include <utility/PyAssert.hh>
#include <basic/Tracer.hh>
//#include <utility/io/izstream.hh>

#include <basic/database/open.hh>

#include <iostream>
#include <fstream>
#include <utility/io/izstream.hh>

#include <boost/algorithm/string.hpp>

static basic::Tracer TR("protocols.antibody.design.DesignInstructionsParser");

namespace protocols {
namespace antibody {
namespace design {

	using namespace boost;
	using namespace protocols::antibody;
	using std::string;
	
	typedef std::map< CDRNameEnum, CDRGraftInstructions > GraftInstructions;

DesignInstructionsParser::DesignInstructionsParser(AntibodyInfoOP const & ab_info, string const path){
	ab_info_=ab_info;
	instructions_path_ = path;
	//This will be refactored for singletons.
	ab_manager_ = ab_info_->get_antibody_enum_manager();
	cluster_manager_ = ab_info_->get_cdr_cluster_enum_manager();
}

DesignInstructionsParser::~DesignInstructionsParser(){}

void
DesignInstructionsParser::check_path(){
	using namespace std;
	ifstream check( instructions_path_.c_str(), ifstream::in);
	if (check.good()){return;}
	else{
		ifstream check2((basic::database::full_name(instructions_path_, false)).c_str(), ifstream::in);
		if (check2.good()){
			instructions_path_ = basic::database::full_name(instructions_path_);
			return;
		}
	}
}

void
DesignInstructionsParser::read_cdr_graft_instructions(GraftInstructions & instructions){
	using namespace utility;
	using namespace std;
	
	check_path();
	//This is straight from C++ tutorials.  
	string line;
	ifstream instruction_file(instructions_path_.c_str());//This may need to change to izstream.  Not sure what the difference is really besides compression.
	PyAssert((instruction_file.is_open()), "Unable to open grafting instruction file.");
	while (getline(instruction_file, line)){
		
		//Skip any comments + empty lines
		if (startswith(line, "#")){
			continue;
		}
		if (startswith(line, "\n")){
			continue;
		}
		
		utility::trim(line, "\n"); //Remove trailing line break
		vector1< string > lineSP = string_split_multi_delim(line); //Split on space or tab
		
		//Everything besides comments needs to have a CDR associated with it.
		PyAssert((ab_manager_->cdr_name_is_present(lineSP[1])), "Unrecognized CDR:"+lineSP[1]);
		
		
		//Here we match.  This is rather ugly, as I don't have much C++ expereince in this.  Python however...

		if (lineSP.size() == 2){
			parse_cdr_graft_general_options(instructions, lineSP);
		}
		else if (lineSP.size() == 3){
			parse_cdr_graft_mintype(instructions, lineSP);
		}
		else if (lineSP.size() >= 4){
			parse_cdr_graft_type_options(instructions, lineSP);
		}

	}
	instruction_file.close();
	TR << "Graft Instructions read successfully" <<std::endl;
}

void
DesignInstructionsParser::parse_cdr_graft_type_options(GraftInstructions& instructions, vector1<string> & lineSP){
	
	CDRNameEnum cdr = ab_manager_->cdr_name_string_to_enum(lineSP[1]);
	
	string type  = lineSP[2];
	string includes = lineSP[3];
	boost::to_upper(type);
	boost::to_upper(includes);
	
	if (type == "CLUSTERS") {
		if (includes == "EXCLUDE"){
			for (core::Size i=4; i<=lineSP.size(); ++i){
				//utility::PyAssert((cluster_manager_->cdr_cluster_is_present(lineSP[i])), "Cluster not found:"+lineSP[i]);
				instructions[cdr].leave_out_clusters.push_back(cluster_manager_->cdr_cluster_string_to_enum(lineSP[i]));
			}
		}
		else if (includes == "INCLUDEONLY"){
			for (core::Size i=4; i<=lineSP.size(); ++i){
				instructions[cdr].include_only_clusters.push_back(cluster_manager_->cdr_cluster_string_to_enum(lineSP[i]));	
			}
		}
		else{utility_exit_with_message("Unrecognized cluster option");}	
	}
	else if (type == "PDBIDS" || type == "PDBID" || type == "PDB") {
		
		if (includes == "EXCLUDE"){
			for (core::Size i=4; i<=lineSP.size(); ++i){
				instructions[cdr].leave_out_pdb_ids.push_back(lineSP[i]);
			}
		}
		else if (includes == "INCLUDEONLY"){
			for (core::Size i=4; i<=lineSP.size(); ++i){
				instructions[cdr].include_only_pdb_ids.push_back(lineSP[i]);	
			}
		}
		else{utility_exit_with_message("Unrecognized pdbid option "+lineSP[3]);}
	}
	else if (type == "TYPES") {
		if (includes == "EXCLUDE"){
			for (core::Size i=4; i<=lineSP.size(); ++i){
				instructions[cdr].cluster_types[utility::string2int(lineSP[i])]=false;
			}
		}
		else if (includes == "INCLUDEONLY"){
			for (core::Size i=4; i<=lineSP.size(); ++i){
				instructions[cdr].cluster_types[utility::string2int(lineSP[i])]=true;	
			}
		}
		else{utility_exit_with_message("Unrecognized type option: "+lineSP[3]);}
	}
	else if (type == "LENGTH"){
		if (includes == "MIN"){
			instructions[cdr].min_length = utility::string2int(lineSP[4]);
		}
		else if (includes == "MAX"){
			instructions[cdr].max_length = utility::string2int(lineSP[4]);
		}
		else{utility_exit_with_message("Unrecognized length option: "+lineSP[3]);}
	}
	else{utility_exit_with_message("Unrecognized option: "+lineSP[2]);}
}
	
void
DesignInstructionsParser::parse_cdr_graft_mintype(GraftInstructions & instructions, vector1<string> & lineSP){
	CDRNameEnum cdr = ab_manager_->cdr_name_string_to_enum(lineSP[1]);
	
	string setting  = lineSP[2];
	string mintype = lineSP[3];
	boost::to_upper(setting);
	boost::to_upper(mintype);
	
	//This may be refactored to enums.  
	
	if (setting == "MINTYPE"){
		if (mintype == "RELAX"){
			instructions[cdr].mintype = relax;
		}
		else if (mintype == "MIN" || mintype == "MINIMIZE" || mintype == "MINIMIZER"){
			instructions[cdr].mintype = minimize;
		}
		else if (mintype == "REPACK"){
			instructions[cdr].mintype = repack;
		}
		else if (mintype == "CENTROIDRELAX" || mintype == "CENTROID_RELAX" || mintype == "CENRELAX"){
			instructions[cdr].mintype = centroid_relax;
		}
		else if (mintype == "NONE"){
			instructions[cdr].mintype = no_min;
		}
		else{utility_exit_with_message("Unrecognized mintype option: "+lineSP[3]);}
	}
	else{utility_exit_with_message("Unrecognized setting: "+lineSP[2]);}
}


void
DesignInstructionsParser::parse_cdr_graft_general_options(GraftInstructions & instructions, vector1<string> & lineSP){
	CDRNameEnum cdr = ab_manager_->cdr_name_string_to_enum(lineSP[1]);
	string option  = lineSP[2];
	boost::to_upper(option);

	if (option == "FIX"){
		instructions[cdr].graft=false;
	}
	else if (option == "STAYNATIVECLUSTER"){
		instructions[cdr].stay_native_cluster=true;
	}
	else if (option == "CENTERSONLY"){
		instructions[cdr].cluster_centers_only=true;
	}
	else{utility_exit_with_message("Unrecognized option: "+option);}
}

} //design
} //antibody
} //protocols
