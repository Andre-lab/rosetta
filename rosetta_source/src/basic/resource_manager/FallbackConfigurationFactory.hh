// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FallbackConfigurationFactory.hh
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configuration_factory_HH
#define INCLUDED_basic_resource_manager_fallback_configuration_factory_HH

//unit headers
#include <basic/resource_manager/FallbackConfigurationFactory.fwd.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfigurationCreator.fwd.hh>

// Utility headers
#include <utility/factory/WidgetRegistrator.hh>

//C++ headers
#include <map>
#include <string>

namespace basic {
namespace resource_manager {

/// @brief The FallbackConfigurationFactory is a singleton factory with which
/// FallbackConfigurationCreators may (should) be registered.  This factory
/// serves as a mechanism for preserving command-line functionality even
/// while switching more protocols from requesting Resources from the
/// ResourceManager rather than reading directly from the command line.
///
/// When a resource with a particular resource-description string is requested
/// from the ResourceManager, and the ResourceManager ddoes not have any instructions
/// on how to load a resource matching that description, then the ResourceManager
/// will then ask the FallbackConfigurationFactory for help.  The
/// FallbackConfigurationFactory can first answer "does this resource description
/// match any resource descriptions for which you have a registered 
/// FallbackConfigurationCreator?" and if the answer is "yes", then the
/// ResourceManager can request that FallbackConfigurationCreator and ask,
/// first, if it is able to construct a resource (FallbackConfigurations
/// often read from the command line, and it's possible that the required
/// flags have not been put on the command line), and if so, then request data
/// from the FallbackConfiguration needed to contruct that resource.
class FallbackConfigurationFactory {
private:
	FallbackConfigurationFactory(); // singleton, private constructor

public:

	FallbackConfigurationOP
	create_fallback_configuration( std::string const & resource_description	) const;

	static FallbackConfigurationFactory * get_instance();

	void
	factory_register( FallbackConfigurationCreatorOP creator );
	
	bool
	has_fallback_for_resource( std::string const & desc ) const;

	/// @brief Only useful for unit testing.  Since factory registration happens (sometimes) at
	/// load time, there may be no one to catch a thrown exception in the event of a name collision
	/// two FallbackConfigurationCreators that register for the same 
	void
	set_throw_on_double_registration();

private:
	static FallbackConfigurationFactory * instance_;

	typedef std::map< std::string, FallbackConfigurationCreatorOP > FallbackConfigurationCreatorsMap;
	bool throw_on_double_registration_;
	FallbackConfigurationCreatorsMap creators_map_;

};

template < class T >
class FallbackConfigurationRegistrator : public utility::factory::WidgetRegistrator< FallbackConfigurationFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< FallbackConfigurationFactory, T > parent;
public:
	FallbackConfigurationRegistrator() : parent() {}
};


} // namespace resource_manager
} // namespace basic

#endif // INCLUDED_basic_resource_manager_fallback_configuration_factory_HH
