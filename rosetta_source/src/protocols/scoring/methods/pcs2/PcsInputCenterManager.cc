// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @begin
 ///
 /// @file protocols/scoring/methods/pcs2/PcsInputCenterManager.cc
 ///
 /// @brief Singleton that hold everything about the input PCS
 /// This avoid multiple reading of the input file.
 ///
 /// @detailed
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references
 ///
 /// @authorsv Christophe Schmitz
 ///
 /// @last_modified February 2010
 ////////////////////////////////////////////////


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsInputCenterManager.hh>
// AUTO-REMOVED #include <protocols/scoring/methods/pcs2/PcsInputCenter.fwd.hh>
//#include <protocols/scoring/methods/pcs2/PcsGridSearchParameterManager.hh>
//#include <protocols/scoring/methods/pcs2/PcsEnergyParameterManager.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers
#include <iostream>

#include <utility/vector1.hh>


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

basic::Tracer TR_PcsInputCenterManager("protocols.scoring.methods.pcs.PcsInputCenterManager");

PcsInputCenterManager::PcsInputCenterManager(){
	TR_PcsInputCenterManager << "Empty constructor called" << std::endl;
}

PcsInputCenterManager *
PcsInputCenterManager::get_instance(){
	if ( instance_ == 0 ){
		 instance_ = new PcsInputCenterManager();
	}
	return instance_;
}

std::ostream &
operator<<(std::ostream& out, const PcsInputCenterManager &me){

	core::Size i, n;
	std::map<std::string, PcsInputCenter> junk = me.PcsInputCenter_all_;
	std::map< std::string, PcsInputCenter >::iterator it;
	n = me.PcsInputCenter_all_.size();

	out << "Number of paramagnetic center: "<<  n << std::endl;

	i = 1;
	for ( it = junk.begin(); it != junk.end(); ++it) {
		out << "Paramagnetic center " << i << " / " << n << std::endl;
		out << "Filename(s) " << it->first <<std::endl;
		out << it->second ;
		i++;
	}
	//	out << std::endl;

	return out;
}


void
PcsInputCenterManager::re_init(){

	//	core::Size n(PcsInputCenter_all_.size());
	//	core::Size i(1);

	PcsInputCenter_all_.clear();
	//	std::cerr <<"CHECKING aa 0 = " << PcsInputCenter_all_.size();

}

	//TODO Why don't I give back a reference?
PcsInputCenter
PcsInputCenterManager::get_PcsInputCenter_for(utility::vector1<std::string> const & filenames, utility::vector1<core::Real> const & weight){
	std::string id;
	core::Size i;

	for(i = 1; i <= filenames.size(); ++i){
		id += filenames[i];
	}

	std::map< std::string, PcsInputCenter >::iterator it;

	for ( it = PcsInputCenter_all_.begin(); it != PcsInputCenter_all_.end(); ++it ) {
		if(it->first == id){
			return(it->second);
		}
	}

	//	PcsInputCenter pcs_i_c;
	//  pcs_i_c = new PcsInputCenter(filenames, weight);

	PcsInputCenter pcs_i_c(filenames, weight);

	//	PcsInputCenter pcs_i_c = new PcsInputCenter(filenames, weight);

	it = PcsInputCenter_all_.begin();
	PcsInputCenter_all_.insert(it, std::pair< std::string , PcsInputCenter >( id ,pcs_i_c));

	//	TR_PcsInputCenterManager << pcs_i_c << std::endl;

	return(pcs_i_c);
}

PcsInputCenterManager * PcsInputCenterManager::instance_( 0 );

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
