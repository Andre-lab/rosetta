// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Patrick Barth

#ifndef INCLUDED_core_scoring_memb_etable_MembEtable_hh
#define INCLUDED_core_scoring_memb_etable_MembEtable_hh

// Unit Headers
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/memb_etable/MembEtable.fwd.hh>

// Package Headers
#include <core/scoring/etable/EtableOptions.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <utility/pointer/access_ptr.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
//#include <ObjexxFCL/FArray5D.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>

// C++ Headers
//#include <iosfwd>
//#include <list>
//#include <string>


namespace core {
namespace scoring {
namespace etable {


/// jk Class definition for Etable
class MembEtable : public Etable {

public:
  typedef Etable Etable;

	///  constructor
	MembEtable(
		chemical::AtomTypeSetCAP atom_set_in, // like etable namespace
		EtableOptions const & options,
		std::string const alternate_parameter_set = ""
  );

  //pba needed ?
  void copy_from( Etable const * source );

	/// const access to the arrays
/*  ObjexxFCL::FArray3D< Real > const &
	ljatr() const
	{
		return ljatr_;
	}

	ObjexxFCL::FArray3D< Real > const &
	ljrep() const
	{
		return ljrep_;
	}
*/
	ObjexxFCL::FArray3D< Real > const &
	solv1() const
	{
		return solv1_;
	}

	ObjexxFCL::FArray3D< Real > const &
	solv2() const
	{
		return solv2_;
	}

  //pba
  ObjexxFCL::FArray3D< Real > const &
  memb_solv1() const
  {
    return memb_solv1_;
  }

  //pba
  ObjexxFCL::FArray3D< Real > const &
  memb_solv2() const
  {
    return memb_solv2_;
  }

	/// const access to the deriv arrays
/*	ObjexxFCL::FArray3D< Real > const &
	dljatr() const
	{
		return dljatr_;
	}

	ObjexxFCL::FArray3D< Real > const &
	dljrep() const
	{
		return dljrep_;
	}
*/
	/// @brief return the solvation derivative table for the desolvation of atom1 by atom2
	ObjexxFCL::FArray3D< Real > const &
	dsolv1() const
	{
		return dsolv1_;
	}

	/// @brief return the solvation derivative table that combines atom1 and atom2's desolvations
	ObjexxFCL::FArray3D< Real > const &
	dsolv2() const
	{
		return dsolv2_;
	}

  /// pba @brief return the solvation derivative table for the desolvation of atom1 by atom2
  ObjexxFCL::FArray3D< Real > const &
  memb_dsolv1() const
  {
    return memb_dsolv1_;
  }

  /// pba @brief return the solvation derivative table that combines atom1 and atom2's desolvations
  ObjexxFCL::FArray3D< Real > const &
  memb_dsolv2() const
  {
    return memb_dsolv2_;
  }

	Real
	max_dis() const
	{
		return max_dis_;
	}

	Real
	get_safe_max_dis2() const
	{
		return safe_max_dis2;
	}

	int
	get_bins_per_A2() const
	{
		return bins_per_A2;
	}

	chemical::AtomTypeSetCAP
	atom_set() const
	{
		return atom_set_;
	}

	Real
	hydrogen_interaction_cutoff2() const
	{
		return hydrogen_interaction_cutoff2_;
	}

	///
	Real
	nblist_dis2_cutoff_XX() const
	{
		return nblist_dis2_cutoff_XX_;
	}

	///
	Real
	nblist_dis2_cutoff_XH() const
	{
		return nblist_dis2_cutoff_XH_;
	}

	///
	Real
	nblist_dis2_cutoff_HH() const
	{
		return nblist_dis2_cutoff_HH_;
	}

	/// @brief Returns the maximum lj radius for any non-hydrogen
	/// atom as defined by the atom-type-set used to create this Etable.
	Real
	max_non_hydrogen_lj_radius() const;

	/// @brief Returns the maximum lj radius for any hydrogen atom as
	/// defined by the input atom-type-set used to create this Etable.
	Real
	max_hydrogen_lj_radius() const;

	/// set these up in the ctor
	inline
	Real
	lj_radius( int const i ) const
	{
		return lj_radius_[i];
	}
/*
	///
	Real
	lj_wdepth( int const i ) const
	{
		return lj_wdepth_[i];
	}
*/
	///
	Real
	lk_dgfree( int const i ) const
	{
		return lk_dgfree_[i];
	}

	///
	Real
	lk_volume( int const i ) const
	{
		return lk_volume_[i];
	}

	///
	Real
	lk_lambda( int const i ) const
	{
		return lk_lambda_[i];
	}

  //pba
  Real
  memb_lk_dgfree( int const i ) const
  {
    return memb_lk_dgfree_[i];
  }

  //pba
  /*Real
  lk_dgrefce( int const i ) const
  {
    return lk_dgrefce_[i];
  }*/

  ObjexxFCL::FArray1D< Real > const &
  lk_dgrefce() const
  {
    return lk_dgrefce_;
  }

  //pba
  /*Real
  memb_lk_dgrefce( int const i ) const
  {
    return memb_lk_dgrefce_[i];
  }*/

  ObjexxFCL::FArray1D< Real > const &
  memb_lk_dgrefce() const
  {
    return memb_lk_dgrefce_;
  }

private:

	void
	output_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		std::string label,
		std::ostream & out
	);

	void
	input_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		const std::string label,
		std::istream & in
	);

private:

	chemical::AtomTypeSetCAP atom_set_;

	// parameters:
	int const n_atomtypes;

	// from options
	Real const max_dis_;
	int const bins_per_A2;
	Real const Wradius; // global mod to radii
	Real const lj_switch_dis2sigma; // actual value used for switch
	Real const max_dis2;
	int const etable_disbins;

	// hard-coded for now
	bool const lj_use_lj_deriv_slope;
	Real const lj_slope_intercept;
	bool const lj_use_hbond_radii;
	Real const lj_hbond_dis;
	Real const lj_hbond_hdis;
	Real const lj_hbond_accOch_dis;
	Real const lj_hbond_accOch_hdis;
	bool const lj_use_water_radii;
	Real const lj_water_dis;
	Real const lj_water_hdis;
	Real const lk_min_dis2sigma;
	Real const min_dis;
	Real const min_dis2; // was double
	bool const add_long_range_damping;
	Real const long_range_damping_length;
	Real const epsilon;
	Real const safe_max_dis2;
	Real hydrogen_interaction_cutoff2_;
	Real const nblist_dis2_cutoff_XX_;
	Real const nblist_dis2_cutoff_XH_;
	Real const nblist_dis2_cutoff_HH_;
	Real max_non_hydrogen_lj_radius_;
	Real max_hydrogen_lj_radius_;

	// these three derived from other data
/*	Real lj_switch_sigma2dis;
	Real lj_switch_value2wdepth;
	Real lj_switch_slope_sigma2wdepth;*/

	//
	utility::vector1< Real > lj_radius_;
	/*utility::vector1< Real > lj_wdepth_;*/
	utility::vector1< Real > lk_dgfree_;
	utility::vector1< Real > lk_volume_;
	utility::vector1< Real > lk_lambda_;
  utility::vector1< Real > memb_lk_dgfree_; //pba
  //utility::vector1< Real > lk_dgrefce_;     //pba
  //utility::vector1< Real > memb_lk_dgrefce_;   //pba
  ObjexxFCL::FArray1D< Real > lk_dgrefce_; //was FArray1D
  ObjexxFCL::FArray1D< Real > memb_lk_dgrefce_; //was FArray1D

	// the etables themselves
/*	ObjexxFCL::FArray3D< Real > ljatr_;
	ObjexxFCL::FArray3D< Real > ljrep_;*/
	ObjexxFCL::FArray3D< Real > solv1_;
	ObjexxFCL::FArray3D< Real > solv2_;
/*	ObjexxFCL::FArray3D< Real > dljatr_;
	ObjexxFCL::FArray3D< Real > dljrep_;*/
	ObjexxFCL::FArray3D< Real > dsolv1_;
	ObjexxFCL::FArray3D< Real > dsolv2_;
  ObjexxFCL::FArray3D< Real > memb_solv1_; //pba
  ObjexxFCL::FArray3D< Real > memb_solv2_; //pba
  ObjexxFCL::FArray3D< Real > memb_dsolv2_; //pba
  ObjexxFCL::FArray3D< Real > memb_dsolv1_; //pba

	// private methods

	// get the AtomType object corresponding to a give type index
	chemical::AtomType const &
	atom_type( int const type )
	{
		return (*atom_set_)[ type ];
	}

	void smooth_etables();
	void modify_pot();
	void make_pairenergy_table();

	// helper functions
	void
	precalc_etable_coefficients(
		ObjexxFCL::FArray2< Real > & lj_sigma,
		/*ObjexxFCL::FArray2< Real > & lj_r6_coeff,
		ObjexxFCL::FArray2< Real > & lj_r12_coeff,
		ObjexxFCL::FArray2< Real > & lj_switch_intercept,
		ObjexxFCL::FArray2< Real > & lj_switch_slope,*/
		ObjexxFCL::FArray1< Real > & lk_inv_lambda2,
		ObjexxFCL::FArray2< Real > & lk_coeff,
    ObjexxFCL::FArray2< Real > & memb_lk_coeff, //pba
		ObjexxFCL::FArray2< Real > & lk_min_dis2sigma_value,
    ObjexxFCL::FArray2< Real > & memb_lk_min_dis2sigma_value //pba
	);

	void
	calc_etable_value(
		Real & dis2,
		int & atype1,
		int & atype2,
/*		Real & atrE,
		Real & d_atrE,
		Real & repE,
		Real & d_repE,*/
		Real & solvE1,
		Real & solvE2,
		Real & dsolvE1,
		Real & dsolvE2,
		ObjexxFCL::FArray2< Real > & lj_sigma,
		/*ObjexxFCL::FArray2< Real > & lj_r6_coeff,
		ObjexxFCL::FArray2< Real > & lj_r12_coeff,
		ObjexxFCL::FArray2< Real > & lj_switch_intercept,
		ObjexxFCL::FArray2< Real > & lj_switch_slope,*/
		ObjexxFCL::FArray1< Real > & lk_inv_lambda2,
		ObjexxFCL::FArray2< Real > & lk_coeff,
		ObjexxFCL::FArray2< Real > & lk_min_dis2sigma_value,
    Real & memb_solvE1, //pba
    Real & memb_solvE2,
    ObjexxFCL::FArray2< Real > & memb_lk_coeff,
    ObjexxFCL::FArray2< Real > & memb_lk_min_dis2sigma_value,
    Real & memb_dsolvE1,
    Real & memb_dsolvE2
	);

/*	void
	zero_hydrogen_and_water_ljatr();*/

};

} // etable
} // scoring
} // core

#endif
