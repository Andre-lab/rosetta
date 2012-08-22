// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran.hh
/// @brief  Ramachandran potential class delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_Ramachandran_hh
#define INCLUDED_core_scoring_Ramachandran_hh

// Unit Headers
#include <core/scoring/Ramachandran.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/interpolation/spline/Bicubic_spline.hh>


// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>


namespace core {
namespace scoring {



class Ramachandran : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	Ramachandran();

	Ramachandran(
		std::string const & rama_map_filename,
		bool use_bicubic_interpolation);

	~Ramachandran() {}

	Real
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi
	) const;

	void
	eval_rama_score_residue(
		conformation::Residue const & res,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	eval_rama_score_residue(
		bool use_bicubic_interpolation,
		bool rama_not_squared,
		AA const res_aa,
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	random_phipsi_from_rama(
		AA const res_aa,
		Real & phi,
		Real & psi
	) const;


	///////////////////////////////
	// unused??
	void
	eval_rama_score_all(
		Pose & pose,
		ScoreFunction const & scorefxn
	) const;

	void
	write_rama_score_all(
		Pose const & pose
	) const;

	//Real get_rama_score_residue_deriv( int res, Pose const & a_pose, ProteinTorsion torsion ) const;
	void eval_procheck_rama( Pose const & a_pose,
		Real & favorable, Real & allowed, Real & generous ) const;


private:

	void read_rama(
		std::string const & rama_map_filename,
		bool use_bicubic_interpolation);

	void init_rama_sampling_table();
/*
	static
	void
	procheck_map_initializer( FArray1D_string & map );
*/
	static bool rama_initialized_;
	//static utility::vector1< utility::vector1 < FArray2D< Real > > > > ram_probabil_;
	static ObjexxFCL::FArray4D< Real > ram_probabil_;
	//static utility::vector1< utility::vector1 < FArray2D_int > > > ram_counts_;
	static ObjexxFCL::FArray4D_int ram_counts_;
	//static utility::vector1< utility::vector1 < FArray2D< Real > > > > ram_energ_;
	static ObjexxFCL::FArray4D< Real > ram_energ_;
	utility::vector1< numeric::interpolation::spline::BicubicSpline > rama_energy_splines_; // loops only

	static int const n_phi_ = 36;
	static int const n_psi_ = 36;
	static Real const binw_; // 360 / n_phi_ = 10;
	static Real const rama_sampling_thold_;
	static Real const rama_sampling_factor_;
	static int const n_aa_ = 20; // Ramachandran score defined for the cananical AAs only.
	static ObjexxFCL::FArray2D< Real > ram_entropy_;
	utility::vector1< utility::vector1< utility::vector1< Real > > > rama_sampling_table_; // vector of allowed phi/psi pairs for each residue

/*
	static bool procheck_map_initialized_
	static FArray1D_string procheck_map_;
*/

};

}
}

#endif
