// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.hh
/// @brief  header of classes for resfile options
/// @author Steven Lewis
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_task_ResfileReader_hh
#define INCLUDED_core_pack_task_ResfileReader_hh

// Unit Headers
#include <core/pack/task/ResfileReader.fwd.hh>

// Package Headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pose/Pose.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>

// Utility Headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// STL Headers
#include <iosfwd>
#include <list>
#include <map>
#include <string>

//Auto Headers
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {



class ResfileReaderException
{

public:

	ResfileReaderException( std::string message){
		message_ = message;
	}

	std::string get_message(){
		return message_;
	}

private:
	std::string message_;
};

class ResfileContents : public utility::pointer::ReferenceCount
{
public:

	ResfileContents( pose::Pose const & pose, std::istream & istream );
	virtual ~ResfileContents();

	std::list< ResfileCommandCOP > const &
	default_commands() const;

	bool specialized_commands_exist_for_residue( Size resid ) const;

	std::list< ResfileCommandCOP > const &
	commands_for_residue( Size resid ) const;


private:
	std::list< ResfileCommandCOP > default_commands_;
	utility::vector1< std::list< ResfileCommandCOP > > commands_;
};

///@brief abstract/interface class for Resfile reader command objects
class ResfileCommand : public utility::pointer::ReferenceCount
{
public:
	virtual ResfileCommandOP clone() const = 0;

	// @brief Read the contents of the Resfile and store the state
	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	) = 0;

	///@brief Modify the packer task with the command that was read in
	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const = 0;

};

///////////////////////////////////////////////////////////////////////
//In this section, list mode-style options (NATAA, etc)

///@brief NATRO disables packing and designing at a position, the residue
///will be totally unchanged
class NATRO : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new NATRO; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "NATRO";}
};

///@brief NATAA allows repacking but no sequence changes (all rotamers are of the original residue)
class NATAA : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new NATAA; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "NATAA";}
};

///@brief ALLAA is deprecated; allows repacking and designing to any canonical residue (default state of PackerTask)
class ALLAA : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new ALLAA; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "ALLAA";}
};

///@brief ALLAAxc allows repacking and designing to any canonical noncysteine residue
class ALLAAxc : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new ALLAAxc; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "ALLAAXC";}
};

///@brief allows repacking and designing to any canonical residue (default state of PackerTask)
class ALLAAwc : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new ALLAAwc; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "ALLAAWC";}
};

///@brief PIKAA allows residues specifed in a following string and packing
class PIKAA : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new PIKAA; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "PIKAA";}
private:
	utility::vector1< bool > keep_canonical_aas_;
	std::list< chemical::AA > na_allowed_; // nucleic acids which are allowed.
};

///@brief PIKNA allows nucleic acid residues specifed in a following string
class PIKNA : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new PIKNA; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "PIKNA";}
private:
	utility::vector1< chemical::AA > keep_nas_;
};

///@brief PIKNA allows nucleic acid residues specifed in a following string
class PIKRNA : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new PIKRNA; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "PIKRNA";}
private:
	utility::vector1< chemical::AA > keep_rnas_;
};

///@brief NOTAA disallows residues specified in a following string, and allows packing
class NOTAA : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new NOTAA; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "NOTAA";}
private:
	utility::vector1< bool > keep_aas_;
};

///@brief EMPTY disallows all canonical residues but leaves packing and designing unchanged (for noncanonicals)
class EMPTY : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new EMPTY; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "EMPTY";}
};

///@brief POLAR allows polar residues and packing
class POLAR : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new POLAR; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "POLAR";}

};

///@brief APOLAR allows nonpolar residues and packing
class APOLAR : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new APOLAR; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "APOLAR";}
};

///@brief APOLA is deprecated, it calls APOLAR to allow nonpolar residues and packing
class APOLA : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new APOLA; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "APOLA";}
};

////////end of mode-style options////////////////////////

////////in this section list other options///////////////

///@brief EX handles the various extrachi options
class EX : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new EX; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "EX";}
private:
	bool           aro_specified_;
	Size           which_chi_;
	ExtraRotSample chi_sample_level_;
};

///@brief NC handles explicit allowance of noncanonical residue types
class NC : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new NC; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "NC";}
private:
	std::string nc_to_include_;
};

///@brief EX_CUTOFF allows setting of the extrachi_cutoff (for determining burial for extra rotamers)
class EX_CUTOFF : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new EX_CUTOFF; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "EX_CUTOFF";}

private:
	Size ex_cutoff_;
};

///@brief USE_INPUT_SC turns on inclusion of the current rotamer for the packer
class USE_INPUT_SC : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new USE_INPUT_SC; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "USE_INPUT_SC";}
};

///@brief AUTO suggests that a packer can/should reconsider the design setting at a/each residue
///@details This is a protocol-level flag to be used in non-vanilla packers. For example, one may want an ALLAA tag to be effective only if the residue is in an automatically-determined region of interest, without knowing which residues qualify a-priori
///@author ashworth
class AUTO : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new AUTO; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "AUTO";}
};

///@brief SCAN suggests to some packing routines that if there are multiple type choices for this residue, then each of them should be considered explicitly in one way or another
///@details This is a protocol-level flag to be used in non-vanilla packers
///@author ashworth
class SCAN : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new SCAN; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "SCAN";}
};

///@brief TARGET flags the position as "targeted", and can optionally specify a "targeted" type
///@details This is a protocol-level flag to be used in non-vanilla packers--positions flagged as "targeted" may be treated in a special fashion
///
///  The optional specification of a target type is be useful for multistate considerations:
///  multistate protocols need 1) rotamers and energies for all possible states, and 2) a target state
///  The target type must be a member of PackerTask's allowed_types_
///
///@author ashworth
class TARGET : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new TARGET; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "TARGET";}
private:
	std::string argstring_;
};

///@brief NO_ADDUCTS will disable adducts, assuming they exist
///@detailed This command exists because if adducts exist, then they are enabled by default for all residues.
///@author ashworth
class NO_ADDUCTS : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new NO_ADDUCTS; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "NO_ADDUCTS";}
};

///@brief FIX_HIS_TAUTOMER: when a histidine is present when the PackerTask is initialized, this flag will fix its tautomer (whether its hydrogen is on ND1 or NE2.  Does nothing if not histidine at initialization (meaning if it mutates to histidine later this flag will have no effect).
class FIX_HIS_TAUTOMER : public ResfileCommand
{
public:
	virtual ResfileCommandOP clone() const { return new FIX_HIS_TAUTOMER; }

	virtual
	void initialize_from_tokens(
		utility::vector1< std::string > const & tokens,
		Size & which_token,
		Size resid
	);

	virtual
	void residue_action(
		PackerTask &,
		Size resid
	) const;

	static std::string name() {return "FIX_HIS_TAUTOMER";}
};

///////////end of other options//////////////////////////
///////////utility functions for resfile reader//////////

///@brief utility function to increment next token to be parsed
	std::string
get_token( const Size which_token, const utility::vector1<std::string> & tokens);

///@brief utility function for resfile reader
utility::vector1< std::string >
tokenize_line( std::istream & inputstream );

///@brief utility for resfile reader, commands MUST be entered into this hard-coded map
std::map< std::string, ResfileCommandOP >
create_command_map();

///@brief utility function for resfile reader (checks for a leading # signaling a comment)
bool
comment_begin( utility::vector1< std::string > const & tokens, Size which_token );

///@brief changes the state of the given PackerTask according to the commands in the resfile at read in from the -pack:resfile option system.
void
parse_resfile(
	pose::Pose const & pose,
	PackerTask & the_task);

///@brief changes the state of the given PackerTask according to the commands in the resfile at filename
void
parse_resfile(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string filename );

///@brief changes the state of the given PackerTask according to the commands in the resfile.
void
parse_resfile_string(
	pose::Pose const & pose,
	PackerTask & the_task,
	std::string const & resfile_string ) throw(ResfileReaderException);

void
onError( std::string message );


}//namespace task
}//namespace pack
}//namespace core

#endif
