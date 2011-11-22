// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


// Rosetta Headers
#include <basic/basic.hh>
#include <basic/interpolate.hh>

#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iosfwd>
#include <limits>
#include <string>

//Auto Headers
#include <execinfo.h>

//#include "wobble.h"

// apl -- Can this be removed?  We have numeric/interpolation/periodic_range etc


namespace basic {

////////////////////////////////////////////////////////////////////////////////
/// @begin interpolate_bilinear_by_value
///
/// @brief get bilinear interpolate, given four points of a 2d periodic function
///
/// @detailed
///
///     Value and derivatives of bilinear interpolation between four
///     specified input values, which are the function values
///     corresponding to four points which bracket the interpolated point
///     in 2-dimensions.  (see diagram below) (see Num. Recipes v2, sec
///     3.6, "Interpolation in two or more dimensions")
///
///     Note that it is automatically assumed that the arguments are
///     periodic, that is f(a,b),a,b=periodic value; the *value* of the
///     function can be treated as periodic ( fixed to a periodicity of
///     360.) as well when the input parameter angle is set to true
///
///     The derivatives are of the bilinear interpolation only. This is
///     not the same value as a numerical derivative of function
///     values. The derivatives here are discontinuous every time you
///     cross over a bin boundary in either dimension
///
///
/// -  bilinear interpolation:
///
///    x0y1 +----------+ x1y1
///         |          |
///   1-yd  |    (x,y) |
///         - - +      |
///     yd  |   .      |
///         |   .      |
///    x0y0 +---|------+ x1y0
///          xd   1-xd
///
///
///
/// @param[in]   x0y0 - in - input bracketing function value (see diagram)
/// @param[in]   x1y0 - in - " "
/// @param[in]   x0y1 - in - " "
/// @param[in]   x1y1 - in - " "
/// @param[in]   xd - in - error value in the 1st dim between lower lst dim
///                   bracket value 'x0' and point to be interpolated, 'x'
/// @param[in]   yd - in - " 2nd dim "
/// @param[in]   binrange - in - range of bin in angles ( for both dimensions )
/// @param[in]   angles - in - if true, treat the the *values* of xy_func as
///                       as having a periodicity of 360. ( note that
///                       the bin ranges are already assumed periodic )
/// @param[out]   val - out - biliinear interpolated value
/// @param[out]   dval_dx - out - derivative of 'val' in the 1st ('x') dim
/// @param[out]   dval_dy - out - " 2nd dim "
///
/// @remarks
///
/// @references
///
/// @authors ctsa 8-19-03
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void
interpolate_bilinear_by_value(
	double const x0y0,
	double const x1y0,
	double const x0y1,
	double const x1y1,
	double const xd,
	double const yd,
	double const binrange,
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy
)
{
	double const w1 = ( 1.0f - xd ) * ( 1.0f - yd );
	double const w2 = xd * ( 1.0f - yd );
	double const w3 = ( 1.0f - xd ) * yd;
	double const w4 = xd * yd;

	if ( angles ) {
		//// ctsa - hack method for weighted angle averaging, if we could burn
		////  cpu time, then the right thing to do would probably be to
		////  break func into x and y components, add up, interpolate,
		////  and find the direction

		double const w12 = w1 + w2;
		double a12;
		if ( w12 != 0.0 ) {
			a12 = ( w1 * x0y0 + w2 * ( subtract_degree_angles(x1y0,x0y0) + x0y0 ) ) / w12;
		} else {
			a12 = 0.0;
		}

		double const w34 = w3 + w4;
		double a34;
		if ( w34 != 0.0 ) {
			a34 = ( w3 * x0y1 + w4 * ( subtract_degree_angles(x1y1,x0y1) + x0y1 ) ) / w34;
		} else {
			a34 = 0.0;
		}

		val = w12 * a12 + w34 * ( subtract_degree_angles(a34,a12) + a12 );
		angle_in_range(val);

		dval_dx = ( ( 1.0f - yd ) * subtract_degree_angles(x1y0,x0y0) +
		 yd * subtract_degree_angles(x1y1,x0y1) ) / binrange;
		dval_dy = ( ( 1.0f - xd ) * subtract_degree_angles(x0y1,x0y0) +
		 xd * subtract_degree_angles(x1y1,x1y0) ) / binrange;
	} else {
		val = x0y0 * w1 + x1y0 * w2 + x0y1 * w3 + x1y1 * w4;

		dval_dx = ( ( 1.0f - yd ) * ( x1y0 - x0y0 ) + yd * ( x1y1 - x0y1 ) ) / binrange;
		dval_dy = ( ( 1.0f - xd ) * ( x0y1 - x0y0 ) + xd * ( x1y1 - x1y0 ) ) / binrange;
	}
}


} // namespace basic
