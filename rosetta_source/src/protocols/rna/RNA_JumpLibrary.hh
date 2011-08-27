// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 1.1.2.1 $
//  $Date: 2005/11/07 21:05:35 $
//  $Author: rhiju $
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_rna_RNA_JumpLibrary_HH
#define INCLUDED_protocols_rna_RNA_JumpLibrary_HH

// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>
// C++ Headers
#include <string>
#include <map>


namespace protocols {
namespace rna {

	///////////////////////////////////////////////////////////////////////////////////
	class BasePairType{
	public:
		char aa1;
		char aa2;
		char edge1;
		char edge2;
		char orientation;

		BasePairType( char const aa1_in, char const aa2_in,
									char const edge1_in, char const edge2_in,
									char const orientation_in){
			aa1 = aa1_in;		aa2 = aa2_in;
			edge1 = edge1_in;		edge2 = edge2_in;
			orientation = orientation_in;
		}

		friend
		bool operator < (BasePairType const & lhs, BasePairType const & rhs )
		{
			//There must be a more elegant way to do this...
			if( lhs.aa1 < rhs.aa1 )
				return true;
			else if ( lhs.aa1 == rhs.aa1 )
				if ( lhs.aa2 < rhs.aa2 )
					return true;
				else if ( lhs.aa2 == rhs.aa2 )
					if ( lhs.edge1 < rhs.edge1 )
						return true;
					else if ( lhs.edge1 == rhs.edge1 )
						if ( lhs.edge2 < rhs.edge2 )
							return true;
						else
							if ( lhs.edge2 == rhs.edge2)
								return ( lhs.orientation < rhs.orientation);
			return false;
		}

	};

	///////////////////////////////////////////////////////////////////////////////////
	class RNA_PairingTemplate : public utility::pointer::ReferenceCount {

	public:

		RNA_PairingTemplate( core::kinematics::Jump const j, std::string const atom_name1, std::string const atom_name2 );

		RNA_PairingTemplate( core::kinematics::Jump const j1, core::kinematics::Jump const j2, std::string const atom_name1, std::string const atom_name2 );

		core::kinematics::Jump const &
		jump() const { return jump_forward_; }

		core::kinematics::Jump const &
		jump_forward() const { return jump_forward_; }

		core::kinematics::Jump const &
		jump_backward() const { return jump_backward_; }

		std::string const &
		atom_name1() const {return atom_name1_; }

		std::string const &
		atom_name2() const {return atom_name2_; }

	private:
		core::kinematics::Jump const jump_forward_;
		core::kinematics::Jump const jump_backward_;
		std::string const atom_name1_;
		std::string const atom_name2_;
	};

	typedef utility::pointer::owning_ptr< RNA_PairingTemplate > RNA_PairingTemplateOP;

	///////////////////////////////////////////////////////////////////////////////////

	typedef utility::vector1< RNA_PairingTemplateOP > RNA_PairingTemplateList;
	typedef std::map< BasePairType, RNA_PairingTemplateList > RNA_PairingTemplateMap;

	///////////////////////////////////////////////////////////////////////////////////
	class RNA_JumpLibrary : public utility::pointer::ReferenceCount {
	public:
		RNA_JumpLibrary( std::string const filename ){ read_jumps_from_file( filename ); };

		void
		read_jumps_from_file( std::string const & jump_library_filename );

		void
		check_forward_backward(
													 std::string & atom_name,
													 bool const forward,
													 core::kinematics::Jump & j,
													 RNA_PairingTemplateOP const & t) const;

		core::kinematics::Jump
		get_random_base_pair_jump(
															char const aa1,
															char const aa2,
															char const edge1,
															char const edge2,
															char const orientation,
															std::string & atom_name1,
															std::string & atom_name2,
															bool & success,
															bool const forward1 = true,
															bool const forward2 = true ) const;

	private:

		void
		save_in_jump_library( core::Size const reschar1, core::Size const reschar2,
													char const edgechar1, char const edgechar2,
													char const orientation,
													std::string const & atom_name1,
													std::string const & atom_name2,
													core::kinematics::Jump const & jump1,
													core::kinematics::Jump const & jump2 );


		RNA_PairingTemplateMap rna_pairing_template_map_;

	};

	typedef utility::pointer::owning_ptr< RNA_JumpLibrary > RNA_JumpLibraryOP;

} //rna
} //protocols

#endif
