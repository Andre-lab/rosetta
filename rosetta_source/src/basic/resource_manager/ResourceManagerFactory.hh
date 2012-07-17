// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceManagerFactory.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceManagerFactory_HH
#define INCLUDED_basic_resource_manager_ResourceManagerFactory_HH

//unit headers
#include <basic/resource_manager/ResourceManagerFactory.fwd.hh>

// package headers
#include <basic/resource_manager/ResourceManagerCreator.fwd.hh>

// Utility headers
#include <utility/factory/WidgetRegistrator.hh>

//C++ headers
#include <map>
#include <string>

namespace basic {
namespace resource_manager {

class ResourceManagerFactory {
private:
	ResourceManagerFactory(); // singleton, private constructor

public:

	/// Should only be called by the ResourceManager in its singleton construction!
	ResourceManager *
	create_resource_manager_from_options_system() const;

	static ResourceManagerFactory * get_instance();

	void
	factory_register( ResourceManagerCreatorOP creator );

private:
	static ResourceManagerFactory * instance_;

	typedef std::map< std::string, ResourceManagerCreatorOP > ResourceManagerCreatorsMap;
	ResourceManagerCreatorsMap creators_map_;

};

template < class T >
class ResourceManagerRegistrator : public utility::factory::WidgetRegistrator< ResourceManagerFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ResourceManagerFactory, T > parent;
public:
	ResourceManagerRegistrator() : parent() {}
};


} // namespace resource_manager
} // namespace basic

#endif //INCLUDED_basic_resource_manager_ResourceManager_HH
