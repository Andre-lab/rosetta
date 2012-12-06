// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Jacob Bale (balej@uw.edu)
#include <devel/matdes/InterfacePackingFilter.hh>
#include <devel/matdes/InterfacePackingFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <basic/database/open.hh>
#include <core/id/AtomID_Map.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/io/pdb/file_data.hh>
#include <devel/matdes/util.hh>
#include <protocols/jd2/JobDistributor.hh>

namespace devel {
namespace matdes {

static basic::Tracer TR( "devel.matdes.InterfacePackingFilter" );

///@brief default ctor
InterfacePackingFilter::InterfacePackingFilter() :
	parent( "InterfacePacking" ),
	distance_cutoff_( 9.0 ),
	contact_dist_( 10.0 ),
	lower_threshold_( -5 ),
	upper_threshold_( 5 ),
	sym_dof_names_( "" )
{}

core::Real
InterfacePackingFilter::distance_cutoff() const{
	return distance_cutoff_;
}

core::Real
InterfacePackingFilter::contact_dist() const{
	return contact_dist_;
}

core::Real
InterfacePackingFilter::lower_threshold() const{
	return lower_threshold_;
}

core::Real
InterfacePackingFilter::upper_threshold() const{
	return upper_threshold_;
}

std::string
InterfacePackingFilter::sym_dof_names() const{
	return sym_dof_names_;
}

void
InterfacePackingFilter::distance_cutoff( core::Real const d ){
	distance_cutoff_ = d;
}

void
InterfacePackingFilter::contact_dist( core::Real const c ){
	contact_dist_ = c;
}

void
InterfacePackingFilter::lower_threshold( core::Real const l ){
	lower_threshold_ = l;
}

void
InterfacePackingFilter::upper_threshold( core::Real const u ){
	upper_threshold_ = u;
}

void
InterfacePackingFilter::sym_dof_names( std::string const s ){
	sym_dof_names_ = s;
}

bool
InterfacePackingFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real packing_score(compute( pose ));
	if( (packing_score >= lower_threshold_) && (packing_score <= upper_threshold_) ){
		TR<<"passing."<<std::endl;
		return true;
	} 
	else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

core::Real
InterfacePackingFilter::compute( core::pose::Pose const & pose ) const{

	using namespace core::pose::symmetry;

	core::pose::Pose sub_pose; 
  core::scoring::packing::HolesParams hp(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	core::Real cutoff2 = distance_cutoff_*distance_cutoff_;
  core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
  core::Size monomer_lower_bound, monomer_upper_bound, ir, jr, base;
  utility::vector1<core::Size> sub_pose_resis, neighbor_resis;
  core::Size count = 0; Real if_score = 0;
	utility::vector1<std::string> sym_dof_name_list;
	std::string sym_dof_name;

	sym_dof_name_list = utility::string_split( sym_dof_names_ , ',' );

	for (Size i = 1; i <= sym_dof_name_list.size(); i++) {

		sym_dof_name = sym_dof_name_list[i];

		if (sym_dof_name == "") {
			utility_exit_with_message("The required argument, sym_dof_names, for the InterfacePackingFilter was now set properly.");
		}

		int	sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name );
		ObjexxFCL::FArray1D_bool is_upstream ( pose.total_residue(), false );
		pose.fold_tree().partition_by_jump( sym_aware_jump_id, is_upstream );
		sub_pose_resis = devel::matdes::get_neighbor_sub_resis(pose, contact_dist_, sym_dof_name);
		core::io::pdb::pose_from_pose(sub_pose, pose, sub_pose_resis);
		//pose.dump_pdb("pose_" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb");
		//sub_pose.dump_pdb("sub_pose_" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb");
	 	core::scoring::packing::HolesResult hr(core::scoring::packing::compute_holes_score(sub_pose, hp));
		TR << "computed_holes" << std::endl;

		Size nres_monomer = 0;
		bool start = true;
		for (Size i=1; i<=symm_info->num_independent_residues(); ++i) {
			if ( is_upstream(i) ) continue;
			if ( start ) monomer_lower_bound = i;
			start = false;
			monomer_upper_bound = i;
			nres_monomer++;
		}
		TR << "nres_monomer: " << nres_monomer << " for sym_dof_name: " << sym_dof_name << std::endl;

		for(Size r=1; r<=sub_pose_resis.size(); r++) {
			if (monomer_lower_bound == sub_pose_resis[r]) { 
				base = r;
				break;
			}
			TR.Debug << "r: " << r << " sub_pose_resi: " << sub_pose_resis[r] << std::endl;
		}
		TR.Debug << "base residue of monomer: " << base << std::endl; 
 
		for (core::Size sir=base; sir<=base+nres_monomer-1; sir++) {
			ir = sub_pose_resis[sir];
			for (core::Size ia = 1; ia<=pose.residue(ir).nheavyatoms(); ia++) {
				bool contact = false;
				for (core::Size sjr=1; sjr<=sub_pose_resis.size(); sjr++) {
					jr = sub_pose_resis[sjr];
					if ( !is_upstream(jr) ) continue; 
					for (core::Size ja = 1; ja<=pose.residue(jr).nheavyatoms(); ja++) {
						if (pose.residue(ir).xyz(ia).distance_squared(pose.residue(jr).xyz(ja)) <= cutoff2)  {
							contact = true;
							break; // ja
						}
					} // ja
					if (contact == true) {
						TR.Debug << "ir: " << ir << " jr: " << jr << std::endl;
						break;
					}
				} // jr
				if (contact == true) {
					count++;
					if_score += hr.atom_scores[core::id::AtomID(ia, sir)];
					TR.Debug << "count: " << count << " atom_score: " << hr.atom_scores[core::id::AtomID(ia, sir)] << " if_score: " << if_score << std::endl;
				}
			} // ia
		} // ir
	}
		TR << "final if_score / count = " << if_score << " / " << count << " = " << (if_score / (Real)count) << std::endl;
		return if_score / (Real)count;
}

core::Real
InterfacePackingFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real packing_score(compute( pose ));
	return( packing_score );
}

void
InterfacePackingFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"InterfacePackingFilter returns "<<compute( pose )<<std::endl;
}

void
InterfacePackingFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	TR << "InterfacePackingFilter"<<std::endl;
	distance_cutoff( tag->getOption< core::Real >( "distance_cutoff", 9.0 ) );
	contact_dist( tag->getOption<core::Real>("contact_dist", 10.0));
	lower_threshold( tag->getOption< core::Real >( "lower_cutoff", -5 ) );
	upper_threshold( tag->getOption< core::Real >( "upper_cutoff", 5 ) );
	sym_dof_names( tag->getOption< std::string >( "sym_dof_names", "" ) );
	TR<<"with options lower_threshold: "<<lower_threshold()<<", upper_threshold: "<<upper_threshold()<<", distance_cutoff: "<<distance_cutoff()<<"contact_dist: "<<contact_dist()<<"sym_dof_names: "<<sym_dof_names()<<std::endl;
}

void InterfacePackingFilter::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & ,
				utility::lua::LuaObject const & ) {
	TR << "InterfacePackingFilter"<<std::endl;
	distance_cutoff( def["distance_cutoff"] ? def["distance_cutoff"].to<core::Real>() : 9.0 );
	lower_threshold( def["lower_cutoff"] ? def["lower_cutoff"].to<core::Real>() : -5 );
	upper_threshold( def["upper_cutoff"] ? def["upper_cutoff"].to<core::Real>() : 5 );
	TR<<"with options lower_threshold: "<<lower_threshold()<<", upper_threshold: "<<upper_threshold()<<", and distance_cutoff: "<<distance_cutoff()<<std::endl;
}
protocols::filters::FilterOP
InterfacePackingFilter::fresh_instance() const{
	return new InterfacePackingFilter();
}

InterfacePackingFilter::~InterfacePackingFilter(){}

protocols::filters::FilterOP
InterfacePackingFilter::clone() const{
	return new InterfacePackingFilter( *this );
}

protocols::filters::FilterOP
InterfacePackingFilterCreator::create_filter() const { return new InterfacePackingFilter; }

std::string
InterfacePackingFilterCreator::keyname() const { return "InterfacePacking"; }

} // matdes
} // devel
