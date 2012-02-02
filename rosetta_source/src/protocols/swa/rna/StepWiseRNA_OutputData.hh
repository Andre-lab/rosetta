// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_SWA_RNA_OutputData_HH
#define INCLUDED_protocols_swa_SWA_RNA_OutputData_HH


#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh> /*For PuckerState and Torsion_Info*/
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh> //Sept 26, 2011 Parin S.

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/chemical/AA.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <set>


typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {


core::io::silent::BinaryRNASilentStruct
get_binary_rna_silent_struct_safe(core::pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file);

core::io::silent::BinaryRNASilentStruct
get_binary_rna_silent_struct_safe_wrapper(core::pose::Pose const & const_pose, std::string const & tag, std::string const & silent_file, bool const write_score_only);

void
Output_data(core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, std::string const & tag, bool const write_score_only, core::pose::Pose const & pose, core::pose::PoseCOP native_poseCOP, StepWiseRNA_JobParametersCOP job_parameters_);


}
}
}

#endif
