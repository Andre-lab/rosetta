// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef INCLUDED_devel_domain_assembly_domain_assembly_setup_hh
#define INCLUDED_devel_domain_assembly_domain_assembly_setup_hh

#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


// AUTO-REMOVED #include <string>

// a class to hold information about the domains that will be assembled together
// in addition to the coordinates of the domains it holds the linker sequences
// that will be attached to each domain and instructions for trimming the N or C terminus
class DomainInfo {
  core::pose::Pose input_pose_;  //not intended to include linkers
  core::pose::Pose processed_pose_;  // includes linkers and truncations
  std::string Nterm_linker_;
  std::string Cterm_linker_;
  core::Size trim_nterm_;
  core::Size trim_cterm_;
  core::Size domain_begin_; //beginning of core domain (not linkers) in the context of the full chain
  core::Size domain_end_; // end ''

public:
  DomainInfo() {
    Nterm_linker_ = Cterm_linker_ = "";
    trim_nterm_ = trim_cterm_ = 0;
    domain_begin_ = domain_end_ = 0;
  }

  void set_input_pose (core::pose::Pose & pose) { input_pose_ = pose; }
  core::pose::Pose get_input_pose () { return input_pose_; }

  void set_processed_pose (core::pose::Pose & pose) { processed_pose_ = pose; }
  core::pose::Pose get_processed_pose () { return processed_pose_; }

  void set_Nterm_linker ( std::string & linker ) { Nterm_linker_ = linker; }
  std::string get_Nterm_linker () { return Nterm_linker_; }

  void set_Cterm_linker ( std::string & linker ) { Cterm_linker_ = linker; }
  std::string get_Cterm_linker () { return Cterm_linker_; }

  void set_trim_nterm ( core::Size & i ) { trim_nterm_ = i; }
  core::Size get_trim_nterm () { return trim_nterm_; }

  void set_trim_cterm ( core::Size & i ) { trim_cterm_ = i; }
  core::Size get_trim_cterm () { return trim_cterm_; }

  void set_domain_begin ( core::Size & i ) { domain_begin_ = i; }
  core::Size get_domain_begin () { return domain_begin_; }

  void set_domain_end ( core::Size & i ) { domain_end_ = i; }
  core::Size get_domain_end () { return domain_end_; }

  ///@brief truncates and adds linker to a domain
  ///  fills processed_pose
  void process_domain( );
};

///@brief builds a full length pose from a set of input pdbs
void assemble_domains_setup();

///@brief optimizes linkers in a multidomain protein
void assemble_domains_optimize();


#endif
