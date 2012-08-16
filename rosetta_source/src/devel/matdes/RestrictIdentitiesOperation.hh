// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/RestrictIdentitiesOperation.hh
/// @brief  TaskOperation class that restricts a vector of Size defined residues to repacking
///			when parsed, it takes in a string and splits by ","
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_devel_matdes_RestrictIdentitiesOperation_hh
#define INCLUDED_devel_matdes_RestrictIdentitiesOperation_hh

// Unit Headers
#include <devel/matdes/RestrictIdentitiesOperation.fwd.hh>
#include <devel/matdes/RestrictIdentitiesOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>

namespace devel {
namespace matdes {

///@details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictIdentitiesOperation : public core::pack::task::operation::TaskOperation {
public:

	RestrictIdentitiesOperation();
	RestrictIdentitiesOperation( utility::vector1 < std::string > identities, bool prevent_repacking );
	virtual core::pack::task::operation::TaskOperationOP clone() const;
	virtual ~RestrictIdentitiesOperation();

	// @brief getters
	utility::vector1< std::string > identities() const;
	bool prevent_repacking() const;

	// @brief setters
	void identities( utility::vector1 < std::string > residues_vec );
	void prevent_repacking( bool const prevent_repacking);

	virtual void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	virtual void parse_tag( TagPtr );

private:
	std::string unparsed_identities_;
	utility::vector1 < std::string > identities_;
	bool prevent_repacking_;
	
};

} //namespace matdes
} //namespace devel

#endif // INCLUDED_devel_matdes_RestrictIdentitiesOperation_HH
