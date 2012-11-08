// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley

/// NOTE-- this file includes both string and map, use .fwd.hh if
/// you can!


#ifndef INCLUDED_core_scoring_methods_EnergyMethodOptions_hh
#define INCLUDED_core_scoring_methods_EnergyMethodOptions_hh

// Unit headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

#include <core/types.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.fwd.hh> // FA_STANDARD_DEFAULT
#include <core/scoring/SecondaryStructureWeights.hh> /// REPLACE THIS WITH .fwd.hh
#include <core/scoring/hbonds/HBondOptions.fwd.hh>

// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh> // CENTROID

#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>

/// Utility headers
// AUTO-REMOVED #include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {

/// add more options here
/// NOTE: If you add an option, make sure you also update the constructor,
/// the assignment operator, the == comparison operator, and the show method in the .cc file!
/// right now this class should be pretty light-weight since a copy is held inside ScoreFunctionInfo
///


class EnergyMethodOptions : public utility::pointer::ReferenceCount {
public:
	///
	EnergyMethodOptions();

	/// @brief Initialize a new EnergyMethodOptions with defaults from the command line.
	void
	initialize_from_options();

	/// copy constructor
	EnergyMethodOptions( EnergyMethodOptions const & src );

	virtual
	~EnergyMethodOptions();

	/// copy operator
	EnergyMethodOptions const &
	operator=( EnergyMethodOptions const & src );

	///
	std::string const &
	etable_type() const;

	///
	void
	etable_type( std::string const & type );

	bool analytic_etable_evaluation() const { return analytic_etable_evaluation_; }
	void analytic_etable_evaluation( bool setting ) { analytic_etable_evaluation_ = setting; }

	///
	std::string const &
	unfolded_energies_type() const;

	///
	void
	unfolded_energies_type( std::string const & type );

	///
	bool
	exclude_protein_protein_hack_elec() const;

	///
	void
	exclude_protein_protein_hack_elec( bool const setting );

	///
	bool
	exclude_monomer_hack_elec() const;

	///
	void
	exclude_monomer_hack_elec( bool const setting );

	/// @brief The maximum (all atom) distance at which hack_elec is non-zero
	core::Real
	hackelec_max_dis() const;

	void
	hackelec_max_dis( core::Real setting );

	/// @brief The minimium (all atom) distance for which hack_elec changes with distances
	core::Real
	hackelec_min_dis() const;

	void
	hackelec_min_dis( core::Real setting );

	/// @brief The dielectric used for the hack_elec term
	core::Real
	hackelec_die() const;

	void
	hackelec_die( core::Real setting );

	/// @brief Should hack_elec use a constant (non-distance dependant) dielectric?
	bool
	hackelec_no_dis_dep_die() const;

	void
	hackelec_no_dis_dep_die( bool setting );

	///
	bool
	exclude_DNA_DNA() const;

	///
	void
	exclude_DNA_DNA( bool const setting );

 	/// @brief Read access to the hbond options object
 	hbonds::HBondOptions const &
 	hbond_options() const;

 	/// @brief non-const access to the hbond options object
 	hbonds::HBondOptions &
 	hbond_options();

 	/// @breif Set the hbond options object -- makes a deep copy
 	void
 	hbond_options( hbonds::HBondOptions const & opts );

	/// @brief  This is used in the construction of the VDW_Energy's AtomVDW object
	std::string const &
	atom_vdw_atom_type_set_name() const;

	///
	void
	atom_vdw_atom_type_set_name( std::string const & setting );

	///
	core::Size
	cst_max_seq_sep() const;

	///
	void
	cst_max_seq_sep( Size const setting );

	/// deprecated
	utility::vector1<std::string> const &
	bond_angle_central_atoms_to_score() const;

	/// depricated
	void
	bond_angle_central_atoms_to_score(
		utility::vector1<std::string> const & atom_names);

	scoring::mm::MMBondAngleResidueTypeParamSetOP
	bond_angle_residue_type_param_set();

	scoring::mm::MMBondAngleResidueTypeParamSetCOP
	bond_angle_residue_type_param_set() const;

	void
	bond_angle_residue_type_param_set(
    scoring::mm::MMBondAngleResidueTypeParamSetOP param_set);

	void set_strand_strand_weights(
		int ss_lowstrand,
		int ss_cutoff);

	///
	SecondaryStructureWeights const &
	secondary_structure_weights() const;

	///
	SecondaryStructureWeights &
	secondary_structure_weights();

	///
	bool
	has_method_weights( ScoreType const & type ) const;

	///
	utility::vector1< Real > const &
	method_weights( ScoreType const & type ) const;

	///
	void
	set_method_weights(
		ScoreType const & type,
		utility::vector1< Real > const & wts);

	///@brief get the harmonic bond angle and bond-length spring constants
	void
	get_cartesian_bonded_parameters( Real &len, Real &ang, Real &tors, Real &proton , Real &imp ) const {
		len=cartbonded_len_;
		ang=cartbonded_ang_;
		tors=cartbonded_tors_;
		proton=cartbonded_proton_;
		imp=cartbonded_improper_;
	}

	///@brief set the harmonic bond angle and bond-length spring constants
	void
	set_cartesian_bonded_parameters( Real len, Real ang, Real tors, Real proton, Real imp ) {
		cartbonded_len_=len;
		cartbonded_ang_=ang;
		cartbonded_tors_=tors;
		cartbonded_proton_=proton;
		cartbonded_improper_=imp;
	}

	///@brief get the harmonic bond angle and bond-length spring constants
	bool get_cartesian_bonded_linear() const {
		return cartbonded_linear_;
	}

	///@brief set the harmonic bond angle and bond-length spring constants
	void set_cartesian_bonded_linear( bool lin_in ) {
		cartbonded_linear_ = lin_in;
	}

	/// used inside ScoreFunctionInfo::operator==
	friend
	bool
	operator==( EnergyMethodOptions const & a, EnergyMethodOptions const & b );

	/// used inside ScoreFunctionInfo::operator==
	friend
	bool
	operator!=( EnergyMethodOptions const & a, EnergyMethodOptions const & b );

	///
	void
	show( std::ostream & out ) const;

	static
	void
	write_score_function_method_options_table_schema(
		utility::sql_database::sessionOP db_session
	);

	void
	insert_score_function_method_options_rows(
		Size batch_id,
		std::string const & score_function_name,
		utility::sql_database::sessionOP db_session) const;

private:
	/// expand this to a class and include ss weights inside
	typedef 	std::map< ScoreType, utility::vector1< Real > > MethodWeights;

private:

	/////////////////////////////////////////////////
	// IMPORTANT NOTE!  If you add an option, make sure you also update the constructor,
	// the assignment operator, the == comparison operator, and the show method in the .cc file!
	/////////////////////////////////////////////////
	std::string etable_type_;
	bool analytic_etable_evaluation_;
	std::string atom_vdw_atom_type_set_name_;
	std::string unfolded_energies_type_;
	MethodWeights method_weights_;
	SecondaryStructureWeights ss_weights_;
	bool exclude_protein_protein_hack_elec_;
	bool exclude_monomer_hack_elec_;
	core::Real hackelec_max_dis_;
	core::Real hackelec_min_dis_;
	core::Real hackelec_die_;
	bool hackelec_no_dis_dep_die_;
	bool exclude_DNA_DNA_;
	hbonds::HBondOptionsOP hbond_options_;

	core::Size cst_max_seq_sep_;
	core::Real cartbonded_len_, cartbonded_ang_, cartbonded_tors_, cartbonded_proton_, cartbonded_improper_;
	bool cartbonded_linear_;

	/// deprecated
	utility::vector1<std::string> bond_angle_central_atoms_to_score_;
	core::scoring::mm::MMBondAngleResidueTypeParamSetOP bond_angle_residue_type_param_set_;
};


std::ostream &
operator<< ( std::ostream & out, EnergyMethodOptions const & options );

}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
