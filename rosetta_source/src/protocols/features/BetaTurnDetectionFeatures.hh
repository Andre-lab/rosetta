// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/LoopAnchorFeatures.hh
/// @brief  report comments stored with each pose
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_features_BetaTurnDetectionFeatures_hh
#define INCLUDED_protocols_features_BetaTurnDetectionFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/LoopAnchorFeatures.fwd.hh>

//External
#include <boost/uuid/uuid.hpp>

// Project Headers
#include <core/types.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <numeric/HomogeneousTransform.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols{
namespace features{

class BetaTurnDetectionFeatures : public FeaturesReporter {
public:
	BetaTurnDetectionFeatures();

	BetaTurnDetectionFeatures( BetaTurnDetectionFeatures const & );

	virtual ~BetaTurnDetectionFeatures();

	///@brief return string with class name
	virtual std::string
	type_name() const;

	///@brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	///@brief generate the beta_turns table schema
	void
	write_beta_turns_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	virtual utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief collect all the feature data for the pose
	virtual core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & /*relevant_residues*/,
		boost::uuids::uuid struct_id,
		utility::sql_database::sessionOP db_session);
private:
	static void setup_conformation_to_turn_type_map();

		bool all_turn_residues_are_on_the_same_chain( core::pose::Pose const & pose, Size first_residue ) const;

		bool residue_range_is_relevant( utility::vector1< bool > const & relevant_residues, Size range_begin, Size range_end ) const;

		bool residue_range_is_protein( core::pose::Pose const & pose, Size range_begin, Size range_end ) const;

		bool beta_turn_present( core::pose::Pose const & pose, Size first_residue ) const;

		std::string const & beta_turn_type( core::pose::Pose const & pose, Size first_residue ) const;
	std::string determine_ramachandran_hash( core::Real phi, core::Real psi, core::Real omega ) const;

private:
		static bool initialized_;
	static Size const beta_turn_length;
		static core::Real const beta_turn_distance_cutoff;
		static std::map< std::string, std::string > conformation_to_turn_type_;

};

} // features namespace
} // protocols namespace

#endif //INCLUDED_protocols_features_BetaTurnDetectionFeatures_hh
