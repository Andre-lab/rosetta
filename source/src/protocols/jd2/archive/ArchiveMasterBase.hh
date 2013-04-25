// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/jd2/MpiFileBuffer.hh
/// @brief  header file for MPISilentFileJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @detail this outputter will send silentstructs via MPI to dedicated node that will collect all structures
/// @author Oliver Lange olange@u.washington.edu


#ifndef INCLUDED_protocols_jd2_archive_ArchiveMasterBase_hh
#define INCLUDED_protocols_jd2_archive_ArchiveMasterBase_hh

#ifdef USEMPI
#include <mpi.h>
#endif


//unit headers
#include <protocols/mpi/MpiFileBuffer.fwd.hh>

//project headers
#include <core/types.hh>

//utility headers
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>

//C++ headers

#include <string>
#include <map>

namespace protocols {
namespace mpi {

class ArchiveMasterBase {
  //singleton?
public:
  //ArchiveMasterBase();
  virtual void get_new_decoys( SilentStructOPs& new_structures );
private:
}

}
}

#endif
