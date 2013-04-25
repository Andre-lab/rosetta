// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Jacob Bale ( balej@uw.edu )

// Unit headers
#include <devel/matdes/ExtractSubpose.hh>
#include <devel/matdes/ExtractSubposeCreator.hh>


// Project Headers
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <devel/matdes/util.hh>

using basic::Warning;
static basic::Tracer TR("devel.matdes.ExtractSubpose");

namespace devel {
namespace matdes {

using namespace core;
using namespace utility;

// -------------  Mover Creator -------------
std::string
ExtractSubposeCreator::keyname() const
{
	return ExtractSubposeCreator::mover_name();
}

protocols::moves::MoverOP
ExtractSubposeCreator::create_mover() const {
	return new ExtractSubpose;
}

std::string
ExtractSubposeCreator::mover_name()
{
	return "ExtractSubpose";
}
// -------------  Mover Creator -------------

ExtractSubpose::ExtractSubpose() :
  sym_dof_names_(""),
	prefix_(""),
	suffix_(""),
  contact_dist_(10.0),
  extras_( 0 )
{ }

ExtractSubpose::ExtractSubpose(const ExtractSubpose& rval) :
	protocols::moves::Mover(),
  sym_dof_names_( rval.sym_dof_names_ ),
	prefix_( rval.prefix_ ),
	suffix_( rval.suffix_ ),
  contact_dist_( rval.contact_dist_ ),
  extras_( rval.extras_ )
{ }

protocols::moves::MoverOP 
ExtractSubpose::clone() const {
	return new ExtractSubpose( *this );
}

protocols::moves::MoverOP 
ExtractSubpose::fresh_instance() const {
	return new ExtractSubpose();
}

void
ExtractSubpose::apply(Pose & pose) {
	using namespace basic;
	using namespace pose;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace scoring;
	typedef vector1<Size> Sizes;

	core::pose::Pose pose_out;
	utility::vector1<core::Size> resis;
	utility::vector1<std::string> sym_dof_name_list = utility::string_split( sym_dof_names_ , ',' );

	if ( extras_ ) {

		// Find out which positions are near the inter-subunit interfaces
		// These will be further screened below, then passed to design()
		SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
		Real const contact_dist_sq = contact_dist_ * contact_dist_;
		utility::vector1<std::string> primary_sub_strings, secondary_sub_strings, tertiary_sub_strings, all_sub_strings;
		utility::vector1<std::string> all_names = sym_dof_names(pose);
		std::map<char,Sizes> comp_subs;
		utility::vector1<char> secondary_comps;
		Size n_sec_tert_subs = 0;

		for (Size i = 1; i <= all_names.size(); i++) {
			comp_subs.insert(std::make_pair(get_jump_name_to_components(pose,all_names[i])[1],get_jump_name_to_subunits(pose, all_names[i])));
		}

		// Use first sym_dof by default.
		if ( sym_dof_name_list.size() == 0) sym_dof_name_list.push_back(all_names[1]);

		// Get all of the indices for the subunits associated with each sym_dof
		for (Size i = 1; i <= sym_dof_name_list.size(); i++) {
			Sizes intra_subs = get_jump_name_to_subunits(pose,sym_dof_name_list[i]);
			for (Size j = 1; j <= intra_subs.size(); j++) {
				std::string s = ObjexxFCL::string_of(intra_subs[j]) + get_jump_name_to_components(pose,sym_dof_name_list[i])[1];
				primary_sub_strings.push_back(s);
				TR << "Primary subunit " << s << " extracted from pose" << std::endl;
			}
		}

		//Get primary subs and all contacting subunits
		for(Size i=1; i<=sym_info->num_total_residues_without_pseudo(); i++) {
			if(find(primary_sub_strings.begin(),primary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(i)) + ObjexxFCL::string_of(get_component_of_residue(pose,i)) )==primary_sub_strings.end()) continue;
			std::string atom_i = (pose.residue(i).name3() == "GLY") ? "CA" : "CB";
			for(Size j=1; j<=sym_info->num_total_residues_without_pseudo(); j++) {
				if(find(primary_sub_strings.begin(),primary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(j)) + ObjexxFCL::string_of(get_component_of_residue(pose,j)) )!=primary_sub_strings.end()) continue;
				if(find(secondary_sub_strings.begin(),secondary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(j)) + ObjexxFCL::string_of(get_component_of_residue(pose,j)) )!=secondary_sub_strings.end()) continue;
				std::string atom_j = (pose.residue(j).name3() == "GLY") ? "CA" : "CB";
				if(pose.residue(i).xyz(atom_i).distance_squared(pose.residue(j).xyz(atom_j)) <= contact_dist_sq) {
					std::string s = ObjexxFCL::string_of(sym_info->subunit_index(j)) + ObjexxFCL::string_of(get_component_of_residue(pose,j)) ;
					secondary_sub_strings.push_back(s);
					TR << "Secondary subunit " << s << " extracted from pose" << std::endl;
					if(find(secondary_comps.begin(),secondary_comps.end(),get_component_of_residue(pose,j))!=secondary_comps.end()) continue;
					secondary_comps.push_back(get_component_of_residue(pose,j));
					n_sec_tert_subs += primary_sub_strings.size() * comp_subs.find(get_component_of_residue(pose,j))->second.size();
					}
				}
			}
		
		//Get other subunits of secondary components
		for(Size i=1; i<=sym_info->num_total_residues_without_pseudo(); i++) {
			if(find(secondary_sub_strings.begin(),secondary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(i)) + ObjexxFCL::string_of(get_component_of_residue(pose,i)) )==secondary_sub_strings.end()) continue;
			if(find(primary_sub_strings.begin(),primary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(i)) + ObjexxFCL::string_of(get_component_of_residue(pose,i)) )!=primary_sub_strings.end()) continue;
			std::string atom_i = (pose.residue(i).name3() == "GLY") ? "CA" : "CB";
			//Stop extracting subunits once all COMPONENTS in contact with the primary subunits have been extracted.
			if( (secondary_sub_strings.size() + tertiary_sub_strings.size())==n_sec_tert_subs) break;
			for(Size j=1; j<=sym_info->num_total_residues_without_pseudo(); j++) {
				//Stop extracting subunits once all COMPONENTS in contact with the primary subunits have been extracted.
				if( (secondary_sub_strings.size() + tertiary_sub_strings.size())==n_sec_tert_subs) break;
				if(get_component_of_residue(pose,i)!=get_component_of_residue(pose,j)) continue;
				if(find(secondary_sub_strings.begin(),secondary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(j)) + ObjexxFCL::string_of(get_component_of_residue(pose,j)) )!=secondary_sub_strings.end()) continue;
				if(find(primary_sub_strings.begin(),primary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(j)) + ObjexxFCL::string_of(get_component_of_residue(pose,j)) )!=primary_sub_strings.end()) continue;
				if(find(tertiary_sub_strings.begin(),tertiary_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(j)) + ObjexxFCL::string_of(get_component_of_residue(pose,j)) )!=tertiary_sub_strings.end()) continue;
				std::string atom_j = (pose.residue(j).name3() == "GLY") ? "CA" : "CB";
				if(pose.residue(i).xyz(atom_i).distance_squared(pose.residue(j).xyz(atom_j)) <= contact_dist_sq) {
					std::string s = ObjexxFCL::string_of(sym_info->subunit_index(j)) + ObjexxFCL::string_of(get_component_of_residue(pose,j));
					tertiary_sub_strings.push_back(s);
					TR << "Tertiary subunit " << s << " extracted from pose" << std::endl;
				}
			}
		}

		for(Size i=1; i<=primary_sub_strings.size() ; i++) {
			all_sub_strings.push_back(primary_sub_strings[i]);
		}	
		for(Size i=1; i<=secondary_sub_strings.size() ; i++) {
			all_sub_strings.push_back(secondary_sub_strings[i]);
		}	
		for(Size i=1; i<=tertiary_sub_strings.size() ; i++) {
			all_sub_strings.push_back(tertiary_sub_strings[i]);
		}	

		//Generate a new pose with the residues from the primary, secondary, and tertiary subunits

		for(Size i=1; i<=sym_info->num_total_residues_without_pseudo(); i++) {
			if(find(all_sub_strings.begin(),all_sub_strings.end(), ObjexxFCL::string_of(sym_info->subunit_index(i)) + ObjexxFCL::string_of(get_component_of_residue(pose,i)) )==all_sub_strings.end()) continue;
			resis.push_back(i);
			}
	} else {
		resis = devel::matdes::get_neighbor_sub_resis(pose, contact_dist_, sym_dof_name_list[1]);
	}
	core::io::pdb::pose_from_pose(pose_out, pose, resis);
	pose_out.dump_pdb(prefix_ + protocols::jd2::JobDistributor::get_instance()->current_output_name() + suffix_ + ".pdb");
}

void 
ExtractSubpose::parse_my_tag( utility::tag::TagPtr const tag,
										 protocols::moves::DataMap & /*data*/,
										 protocols::filters::Filters_map const &,
										 protocols::moves::Movers_map const &,
										 core::pose::Pose const & ) {

	sym_dof_names_ = tag->getOption< std::string >( "sym_dof_names","" );
	prefix_ = tag->getOption< std::string >( "prefix", "" );
	suffix_ = tag->getOption< std::string >( "suffix", "" );
	contact_dist_ = tag->getOption<core::Real>("contact_dist", 10.0);
	extras_ = tag->getOption<bool>("extras", 0 );

}


} // matdes
} // devel
