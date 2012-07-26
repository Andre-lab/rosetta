// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data_options.hh
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_io_pdb_file_data_options_HH
#define INCLUDED_core_io_pdb_file_data_options_HH


// Unit headers
#include <core/io/pdb/file_data_options.fwd.hh>

// Basic headers
#include <basic/resource_manager/ResourceOptions.hh>

// C++ headers
#include <string>


namespace core {
namespace io {
namespace pdb {

class FileDataOptions : public basic::resource_manager::ResourceOptions
{
public:
	FileDataOptions();

	virtual ~FileDataOptions();

	virtual
	void parse_my_tag( utility::tag::TagPtr tag );

	virtual
	std::string type() const;

	// accessors
	bool check_if_residues_are_termini() const;
	bool exit_if_missing_heavy_atoms() const;
	bool ignore_unrecognized_res() const;
	bool ignore_waters() const;
	bool ignore_zero_occupancy() const;
	bool keep_input_protonation_state() const;
	bool preserve_header() const;
	bool randomize_missing_coords() const;
	bool remember_unrecognized_res() const;
	bool remember_unrecognized_water() const;
	std::string const & chains_whose_residues_are_separate_chemical_entities() const;

	// mutators
	void set_check_if_residues_are_termini( bool check_if_residues_are_termini );
	void set_exit_if_missing_heavy_atoms( bool exit_if_missing_heavy_atoms );
	void set_ignore_unrecognized_res( bool ignore_unrecognized_res );
	void set_ignore_waters( bool ignore_waters );
	void set_ignore_zero_occupancy( bool ignore_zero_occupancy );
	void set_keep_input_protonation_state( bool keep_input_protonation_state );
	void set_preserve_header( bool preserve_header );
	void set_randomize_missing_coords( bool randomize_missing_coords );
	void set_remember_unrecognized_res( bool remember_unrecognized_res );
	void set_remember_unrecognized_water( bool remember_unrecognized_water );
	void set_chains_whose_residues_are_separate_chemical_entities( std::string const & setting );

private:
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

private:
	bool check_if_residues_are_termini_; //termini
	bool exit_if_missing_heavy_atoms_;
	bool ignore_unrecognized_res_;
	bool ignore_waters_;
	bool ignore_zero_occupancy_;
	bool keep_input_protonation_state_;
	bool preserve_header_;
	bool randomize_missing_coords_;
	bool remember_unrecognized_res_;
	bool remember_unrecognized_water_;

	std::string chains_whose_residues_are_separate_chemical_entities_; //treat_residues_in_these_chains_as_separate_chemical_entities
};

} // namespace pdb
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_file_data_options_HH
