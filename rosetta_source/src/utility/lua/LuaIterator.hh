// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=1 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/lua/LuaIterator.hh
/// @brief A wrapper around a luabind::iterator
/// has skey() and ikey() functions that return a string or int conversion of the key
// 
// this is really just a convenience class
/// @author Ken Jung

#ifndef INCLUDED_utility_lua_LuaIterator_hh
#define INCLUDED_utility_lua_LuaIterator_hh

#include <utility/lua/LuaIterator.fwd.hh>
#include <utility/lua/LuaObject.fwd.hh>
#ifdef USELUA
#include <lua.hpp>
#include <luabind/luabind.hpp>

#include <boost/mpl/assert.hpp>

namespace utility {
namespace lua {

class LuaIterator {

		public:
    LuaIterator( luabind::iterator iterator) : 
						iterator_ (iterator) {}

    LuaIterator(){}

				~LuaIterator(){}

				std::string skey();

				int ikey();

				luabind::iterator raw();

				bool operator==(LuaIterator & other);

				bool operator!=(LuaIterator & other);

				LuaObject operator * ();
				LuaObject *  operator -> ();

				LuaIterator operator++(int);

				LuaIterator & operator++();

		private:
				luabind::iterator iterator_;

};

} //lua
} //utility
#endif
#endif
