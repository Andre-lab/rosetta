// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/epr_deer/metrics/DEERJSMethod.fwd.hh
/// @brief   Container for scoring method for DEER distributions
/// @author   Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_metrics_DEERJSMethod_fwd_hh
#define INCLUDED_core_scoring_epr_deer_metrics_DEERJSMethod_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace epr_deer {
namespace metrics {

class DEERJSMethod;

typedef utility::pointer::shared_ptr< DEERJSMethod > DEERJSMethodOP;
typedef utility::pointer::shared_ptr< DEERJSMethod const > DEERJSMethodCOP;

} // namespace metrics
} // namespace epr_deer
} // namespace scoring
} // namespace core

#endif
