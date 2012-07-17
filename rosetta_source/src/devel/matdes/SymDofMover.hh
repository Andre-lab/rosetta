
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

#ifndef INCLUDED_devel_matdes_SymDofMover_HH
#define INCLUDED_devel_matdes_SymDofMover_HH

// Unit headers
#include <devel/matdes/SymDofMover.fwd.hh>

// Package headers

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

namespace devel {
namespace matdes {

class SymDofMover : public protocols::moves::Mover {
	typedef core::pose::Pose Pose;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;
	typedef utility::tag::TagPtr TagPtr;

public:
  SymDofMover();
	//SymDofMover(const SymDofMover& rval);

  // --- virtual functions from mover ---
  virtual std::string get_name() const { return "SymDofMover"; }
  virtual void apply(Pose& pose);
	virtual void trans_pose( Pose & pose, numeric::xyzVector<Real> const & trans, Size start, Size end );
	virtual void rot_pose( Pose & pose, numeric::xyzMatrix<Real> const & rot, Size start, Size end );
	virtual void rot_pose( Pose & pose, numeric::xyzMatrix<Real> const & rot, numeric::xyzVector<Real> const & cen, Size start, Size end );
	virtual void rot_pose( Pose & pose, numeric::xyzVector<Real> const & axis, double const & ang, Size start, Size end );
	virtual void rot_pose( Pose & pose, numeric::xyzVector<Real> const & axis, double const & ang, numeric::xyzVector<Real> const & cen, Size start, Size end );
	virtual void alignaxis(core::pose::Pose & pose, numeric::xyzVector<Real> newaxis, numeric::xyzVector<Real> oldaxis, numeric::xyzVector<Real> cen , Size start, Size end );

	// --- virtual copy constructors 
	virtual MoverOP clone() const;


	/// @brief create this type of object
	virtual MoverOP fresh_instance() const;


	virtual void parse_my_tag( TagPtr const tag,
														 DataMap & data,
														 Filters_map const &,
														 Movers_map const &,
														 Pose const & );

private:
	utility::vector1<std::string> get_sym_dof_names();
	utility::vector1<Real> get_radial_disps();
	utility::vector1<Real> get_angles();

private:
	std::string symm_file_;
	utility::vector1<std::string> sym_dof_names_;
	utility::vector1<Real> radial_disps_;
	utility::vector1<Real> angles_;
	utility::vector1<Real> radial_disps_range_min_;
	utility::vector1<Real> radial_disps_range_max_;
	utility::vector1<Real> angles_range_min_;
	utility::vector1<Real> angles_range_max_;
	utility::vector1<Real> radial_disp_steps_;
	utility::vector1<Real> angle_steps_;
	utility::vector1<Real> radial_disp_deltas_;
	utility::vector1<Real> angle_deltas_;
	std::string sampling_mode_;
	
};

} // matdes
} // devel
#endif
