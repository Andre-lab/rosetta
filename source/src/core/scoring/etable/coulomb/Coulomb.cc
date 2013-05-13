// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/coulomb/Coulomb.hh
/// @brief  Evaluate Coulombic potential
/// @author Phil Bradley, modifed by James Gleixner
/// @author Matthew O'Meara

// unit headers
#include <core/scoring/etable/coulomb/Coulomb.hh>

// project headers
#include <core/scoring/methods/EnergyMethodOptions.hh>

// numeric headers
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>

namespace core {
namespace scoring {
namespace etable {
namespace coulomb {


////////////////////////////////////////////////////////////////////////////
Coulomb::Coulomb( methods::EnergyMethodOptions const & options ):
	max_dis_( options.hackelec_max_dis() ),
	min_dis_( options.hackelec_min_dis() ),
	smooth_hack_elec_( options.smooth_hack_elec() ),
	die_( options.hackelec_die() ),
	no_dis_dep_die_( options.hackelec_no_dis_dep_die() )
{
	initialize();
}


////////////////////////////////////////////////////////////////////////////
Coulomb::Coulomb( Coulomb const & src ): ReferenceCount(),
	max_dis_( src.max_dis_ ),
	min_dis_( src.min_dis_ ),
	smooth_hack_elec_( src.smooth_hack_elec_ ),
	die_( src.die_ ),
	no_dis_dep_die_( src.no_dis_dep_die_ )
{
	initialize();
}

void
Coulomb::initialize() {
	// Must have already initialized max_dis_, min_dis_, die_, and no_dis_dep_die_

	//max_dis_ = 5.5;
	max_dis2_ = max_dis_ * max_dis_;
	//min_dis_ = 1.5;
	min_dis2_ = min_dis_ * min_dis_ ;

	// default dielectric is 10r
	//die_ = 10.0;
	//no_dis_dep_die_ = false;

	C0_ = 322.0637 ;
	C1_ = C0_ / die_ ;
	if( no_dis_dep_die_ ) {
		C2_ = C1_ / max_dis_ ;
		min_dis_score_ = C1_ / min_dis_ - C2_ ;
		dEfac_ = -1.0 * C0_ / die_ ;
	} else {
		C2_ = C1_ / max_dis2_ ;
		min_dis_score_ = C1_ / min_dis2_ - C2_ ;
		dEfac_ = -2.0 * C0_ / die_ ;
	}

	if ( smooth_hack_elec_ ) {


		low_poly_start_ = min_dis_ - 0.25;
		low_poly_end_   = min_dis_ + 0.25;
		low_poly_start2_ = low_poly_start_ * low_poly_start_;
		low_poly_end2_   = low_poly_end_ * low_poly_end_;
		low_poly_width_ = low_poly_end_ - low_poly_start_;
		low_poly_invwidth_ = 1.0 / low_poly_width_;

		// scope low polynomial
		{
			using namespace numeric::interpolation::spline;
			Real low_poly_end_score(0.0), low_poly_end_deriv(0.0);
			if ( no_dis_dep_die_ ) {
				low_poly_end_score = C1_ / low_poly_end_ - C2_;
				low_poly_end_deriv = -1 * C1_ / low_poly_end2_ ;

			} else {
				low_poly_end_score = C1_ / low_poly_end2_ - C2_;
				low_poly_end_deriv = -2 * C1_ / ( low_poly_end2_ * low_poly_end_ );
			}
			SplineGenerator gen_low_poly(
				low_poly_start_, min_dis_score_, 0,
				low_poly_end_, low_poly_end_score, low_poly_end_deriv );
			InterpolatorOP interp_low( gen_low_poly.get_interpolator() );
			SimpleInterpolatorOP sinterp_low = dynamic_cast< SimpleInterpolator * > (interp_low() );
			if ( ! sinterp_low ) {
				utility_exit_with_message( "Hack Elec created non-simple-interpolator in initialize()" );
			}
			low_poly_.ylo  = sinterp_low->y()[ 1 ];
			low_poly_.yhi  = sinterp_low->y()[ 2 ];
			low_poly_.y2lo = sinterp_low->ddy()[ 1 ];
			low_poly_.y2hi = sinterp_low->ddy()[ 2 ];
		}

		hi_poly_start_    = max_dis_ - 1.0;
		hi_poly_end_      = max_dis_;
		hi_poly_start2_   = hi_poly_start_ * hi_poly_start_;
		hi_poly_end2_     = hi_poly_end_ * hi_poly_end_;
		hi_poly_width_    = hi_poly_end_ - hi_poly_start_;
		hi_poly_invwidth_ = 1.0 / hi_poly_width_;

		// scope hi polynomial
		{
			using namespace numeric::interpolation::spline;
			Real hi_poly_start_score(0.0), hi_poly_start_deriv(0.0);
			if ( no_dis_dep_die_ ) {
				hi_poly_start_score = C1_ / hi_poly_start_ - C2_;
				hi_poly_start_deriv = -1 * C1_ / hi_poly_start2_ ;

			} else {
				hi_poly_start_score = C1_ / hi_poly_start2_ - C2_;
				hi_poly_start_deriv = -2 * C1_ / ( hi_poly_start2_ * hi_poly_start_ );
			}

			SplineGenerator gen_hi_poly(
				hi_poly_start_, hi_poly_start_score, hi_poly_start_deriv,
				hi_poly_end_, 0, 0 );
			InterpolatorOP interp_hi( gen_hi_poly.get_interpolator() );
			SimpleInterpolatorOP sinterp_hi = dynamic_cast< SimpleInterpolator * > (interp_hi() );
			if ( ! sinterp_hi ) {
				utility_exit_with_message( "Hack Elec created non-simple-interpolator in initialize()" );
			}
			hi_poly_.ylo  = sinterp_hi->y()[ 1 ];
			hi_poly_.yhi  = sinterp_hi->y()[ 2 ];
			hi_poly_.y2lo = sinterp_hi->ddy()[ 1 ];
			hi_poly_.y2hi = sinterp_hi->ddy()[ 2 ];
		}
	} else {
		low_poly_start_ = min_dis_;     low_poly_start2_ = std::pow( low_poly_start_, 2 );
		low_poly_end_   = min_dis_ / 2; low_poly_end2_   = std::pow( low_poly_end_, 2 );
		hi_poly_start_  = max_dis_;     hi_poly_start2_  = std::pow( hi_poly_start_, 2 );
		low_poly_width_ = 0; low_poly_invwidth_ = 0;
		hi_poly_width_ = 0; hi_poly_invwidth_ = 0;
	}


	//low_fade_start_ = min_dis_;
	//low_fade_start2_ = low_fade_start_ * low_fade_start_;
	//low_fade_end_ = min_dis_ + 0.75;
	//low_fade_end2_ = low_fade_end_ * low_fade_end_;
	//low_fade_d0_ = min_dis_ + 0.25;
	//low_fade_K_ = 16;
	//high_fade_start_ = max_dis_ - 1.0;
	//high_fade_start2_ = high_fade_start_ * high_fade_start_;
	//high_fade_end_ = max_dis_;
	//high_fade_end2_ = high_fade_end_ * high_fade_end_;
	//high_fade_K_ = -12;
	//high_fade_d0_ = max_dis_ - 0.25;

}

}
}
}
}
