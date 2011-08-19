// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/string_util.hh
///
/// @brief  Some std::string helper functions.
/// @author Sergey Lyskov
#ifndef INCLUDED_utility_string_util_hh
#define INCLUDED_utility_string_util_hh

// Utility headers
#include <utility/vector1.hh>
#include <utility/stream_util.hh>

// Boost headers
#include <boost/algorithm/string/erase.hpp>

// C++ headers
#include <list>
#include <sstream>
#include <string>
#include <vector>

namespace utility {

/// @brief Reads the contents of <filename> into <contents>, preserving newline
/// characters. Aborts if an error is encoutered.
void ReadFromFileOrDie(const std::string& filename, std::string* contents);
	
/// @brief split given std::string using ' ' symbol.
std::vector< std::string > split(const std::string &s);

///@brief combine strings with anything
std::string join(utility::vector1<std::string> const & s, std::string const & connector);

///@brief combine strings with anything
std::string join(std::vector<std::string> const & s, std::string const & connector);

/// @brief join space separations in a string with a connector such as '_'
std::string join(std::string const & string_w_spaces, std::string const & connector);

/// @brief split given std::string using ' ' symbol.
std::list< std::string > split_to_list(const std::string &s);

/// @details split to vector< std::string > using arbitrary split character
std::vector< std::string >
string_split( std::string const & in, char splitchar = ' ' );

/// @brief convert a string to a float
float string2float( std::string st );

/// @brief convert a string to an int
int string2int( std::string st );

// @brief Reads an unsigned int from string <x>, writing the result
// to output parameter <y>, which must be non-NULL. If the read was not
// successful, this function call has no effect on the value of <y> that
// was present prior to invokation.
void string2uint(const std::string& x, unsigned int *y);

/// @brief True iff haystack starts with needle
bool startswith(std::string const & haystack, std::string const & needle);

/// @brief True iff haystack ends with needle
bool endswith(std::string const & haystack, std::string const & needle);

void slurp(std::istream & in, std::string & out);

void trim( std::string & s, const std::string & drop = " " );

/// @brief create a new string that drops all the unwanted substrings of
/// the original string.
std::string
trim( std::string const & s, std::string const & drop = " " );

/// @brief compares two strings ignoring leading and trailing spaces
bool trimmed_compare( std::string const & s1, std::string const & s2 );

/// @brief adds spaces to a left aligned string until a given length is reached
void add_spaces_left_align( std::string & st, std::size_t const newlen );

/// @brief adds spaces to a right aligned string until a given length is reached
void add_spaces_right_align( std::string & st, std::size_t const newlen );

// @brief return true of the string has only [0-9], ,'+','-','.' or '[Ee]'
bool is_string_numeric(std::string const & input);

/// @brief Read the entire contents of a file into a string.  All end-of-line characters are replaced
/// by "\n".  Throws a utility::excn::EXCN_msg_exception if the file cannot be opened.
std::string file_contents( std::string const & file_name );

std::string file_basename( std::string const & full_path );

// "/foo/bar/baz" => "baz"
// "/foo/bar/baz.cc" => "baz.cc"
std::string filename(const std::string& path);

// "/foo/bar/baz" => "/foo/bar/"
std::string pathname(const std::string& path);


///@brief find all environment variables with the form ${VARIABLE}
/// and replace with the contents of that environment variable.
/// if the environment variable does not exist, return string::npos
std::string replace_environment_variables(std::string input);

/// @brief Compares two strings, ignoring spaces.  Useful for comparing atom
/// name strings which have pdb-alignment built into them.  Slightly dangerous
/// if you consider the fact that atom names in the PDB are different for
/// different indentation rules: ' CA ' is c-alpha.  'CA  ' is calcium.
inline
bool same_ignoring_spaces( std::string const & s1, std::string const & s2 ) {
	std::string t1 = boost::algorithm::erase_all_copy(s1, " ");
	std::string t2 = boost::algorithm::erase_all_copy(s2, " ");
	return t1 == t2;
}

inline
void replace_in( std::string & s, const char from, const char *to )
{
	// fix string
	for ( unsigned int c = 0; c < s.length(); ++c ) {
		if( s[c] == from ) s.replace(c,1,to);
	}
}

template <class T>
inline std::string to_string (const T & t)
{
	std::ostringstream ss;
	ss << t;
	return ss.str();
}

template <class T>
inline T const from_string (std::string const & s, T )
{
	T t;
	std::istringstream ss(s);
	ss >> t;
	return t;
}

}  // namespace utility

#endif  // INCLUDED_utility_string_util_HH
