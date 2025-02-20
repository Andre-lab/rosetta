// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/JobSummary.hh
/// @brief  The definition for class protocols::jd3::JobSummary
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_JobSummary_HH
#define INCLUDED_protocols_jd3_JobSummary_HH

// Unit headers
#include <protocols/jd3/JobSummary.fwd.hh>

// Package headers

//utility headers
#include <utility/VirtualBase.hh>

//C++ headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

/// @brief %JobSummary class holds the output that's generated by a Job over the course
/// of its execution.  The %JobSummary is handed by the JobQueen to the JobOutputWriter
/// objects, each of which have the opportunity to pull data out of the JobSummary
/// class.
class JobSummary : public utility::VirtualBase
{
public:

	JobSummary();
	~JobSummary() override;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // JobSummary

} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_JobSummary )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_JobSummary_HH
