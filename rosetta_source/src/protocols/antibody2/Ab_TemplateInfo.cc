// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocls/antibody2/Ab_TemplateInfo.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/Ab_TemplateInfo.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/PDBInfo.hh>

static basic::Tracer TR("antibody2.Ab_TemplateInfo");

namespace protocols {
namespace antibody2 {

/// @details Auto-generated virtual destructor
Ab_TemplateInfo::~Ab_TemplateInfo() {}

using namespace core;


/// default constructor
Ab_TemplateInfo::Ab_TemplateInfo() {
    set_default( false/*camelid*/ );
}


/// constructor with arguments  1
Ab_TemplateInfo::Ab_TemplateInfo(bool load_L1, bool load_L2, bool load_L3,
                                 bool load_H1, bool load_H2, bool load_H3)
{
    load_templates_from_pdbs(load_L1, load_L2, load_L3,
                             load_H1, load_H2, load_H3, false/*camelid*/ );
}


/// constructor with arguments  1
Ab_TemplateInfo::Ab_TemplateInfo(bool load_L1, bool load_L2, bool load_L3,
                                 bool load_H1, bool load_H2, bool load_H3, bool camelid)
{
    load_templates_from_pdbs(load_L1, load_L2, load_L3,
                             load_H1, load_H2, load_H3, camelid );
}








void
Ab_TemplateInfo::load_templates_from_pdbs(bool load_L1, bool load_L2, bool load_L3,
                                          bool load_H1, bool load_H2, bool load_H3, bool camelid)
{
    set_default(camelid);

    load_L1_ = load_L1;
    load_L2_ = load_L2;
    load_L3_ = load_L3;
    load_H1_ = load_H1;
    load_H2_ = load_H2;
    load_H3_ = load_H3;


    if(  camelid && (load_L1_||load_L2_||load_L3_)   ){
        utility_exit_with_message("This is Camelid antibody, No Light Chain !!!");
    }


    if(!camelid_){
        if(load_L1_){
            import_pose::pose_from_pdb( L1_t_pose_, "input/l1.pdb" );
            templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("L1", L1_t_pose_) );
        }
        if(load_L2_){
            import_pose::pose_from_pdb( L2_t_pose_, "input/l2.pdb" );
            templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("L2", L2_t_pose_) );
        }
        if(load_L3_){
            import_pose::pose_from_pdb( L3_t_pose_, "input/l3.pdb" );
            templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("L3", L3_t_pose_) );
        }
    }


    if(load_H1_){
        import_pose::pose_from_pdb( H1_t_pose_, "input/h1.pdb" );
        templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("H1", H1_t_pose_) );
    }
    if(load_H2_){
        import_pose::pose_from_pdb( H2_t_pose_, "input/h2.pdb" );
        templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("H2", H2_t_pose_) );
    }
    if(load_H3_){
        import_pose::pose_from_pdb( H3_t_pose_, "input/h3.pdb" );
        templates_poses_.insert(  std::pair<std::string, core::pose::Pose> ("H3", H3_t_pose_) );
    }


}



void
Ab_TemplateInfo::set_default( bool camelid )
{
    camelid_ = camelid;

    load_L1_ = true;
    load_L2_ = true;
    load_L3_ = true;

    load_H1_ = true;
    load_H2_ = true;
    load_H3_ = true;
}




pose::Pose Ab_TemplateInfo::get_one_template_pose(std::string cdr_name){

    TemplatePoseMap::iterator iter = templates_poses_.begin();
	iter = templates_poses_.find(cdr_name);
	if ( iter != templates_poses_.end() ) {return iter->second;}

}



void Ab_TemplateInfo::obtain_templates_names(){
    std::ifstream inf;
    inf.open("query.matches");
    if(!inf.is_open()) {utility_exit_with_message("Cannot open 'query.matches' file!!");}

    std::string temp,tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    LightHeavy_t_name_ = tttt;

    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;

    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    L1_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    L2_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    L3_t_name_ = tttt;

    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    H1_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    H2_t_name_ = tttt;
    inf>>temp>>temp>>temp>>temp>>tttt>>temp>>temp;
    H3_t_name_ = tttt;
}






/// @details  Show the complete setup of the Ab_TemplateInfo
void
Ab_TemplateInfo::show( std::ostream & out ) {
    //	if ( !flags_and_objects_are_in_sync_ ){
    //		sync_objects_with_flags();
    //	}
    out << *this;
}

std::ostream & operator<<(std::ostream& out, const Ab_TemplateInfo & ab_t_info )
{
    using namespace ObjexxFCL::fmt;
    // All output will be 80 characters - 80 is a nice number, don't you think?
    std::string line_marker = "///";
    out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
        out << line_marker << A( 47, "Rosetta Antibody Template Info" ) << space( 27 ) << line_marker << std::endl;
        out << line_marker << space( 74 ) << line_marker << std::endl;

        if(ab_t_info.load_L1_){
            out << line_marker << " L1 template: "<<std::endl;
            out << line_marker << "   template_name:  "<<ab_t_info.L1_t_name_<<std::endl;
            out << line_marker << ab_t_info.L1_t_pose_ << std::endl;
        }

        if(ab_t_info.load_L2_){
            out << line_marker << " L2 template: "<<std::endl;
            out << line_marker << "   template_name:  "<<ab_t_info.L2_t_name_<<std::endl;
            out << line_marker << ab_t_info.L2_t_pose_<<std::endl;
        }

        if(ab_t_info.load_L3_){
            out << line_marker << " L3 template: "<<std::endl;
            out << line_marker << "   template_name:  "<<ab_t_info.L3_t_name_<<std::endl;
            out << line_marker << ab_t_info.L3_t_pose_<<std::endl;
        }

        if(ab_t_info.load_H1_){
            out << line_marker << " H1 template: "<<std::endl;
            out << line_marker << "   template_name:  "<<ab_t_info.H1_t_name_<<std::endl;
            out << line_marker << ab_t_info.H1_t_pose_<<std::endl;
        }

        if(ab_t_info.load_H2_){
            out << line_marker << " H2 template: "<<std::endl;
            out << line_marker << "   template_name:  "<<ab_t_info.H2_t_name_<<std::endl;
            out << line_marker << ab_t_info.H2_t_pose_<<std::endl;
        }

        if(ab_t_info.load_H3_){
            out << line_marker << " H3 template: "<<std::endl;
            out << line_marker << "   template_name:  "<<ab_t_info.H3_t_name_<<std::endl;
            out << line_marker << ab_t_info.H3_t_pose_<<std::endl;
        }

        // Close the box I have drawn
        out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
        return out;
    }






} //namespace antibody2
} //namespace protocols
