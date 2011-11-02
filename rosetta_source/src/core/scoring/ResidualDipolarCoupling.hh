// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ResidualDipolarCoupling.hh
/// @brief  Uses NMR RDC for scoring
/// @author Srivatsan Raman

#ifndef INCLUDED_core_scoring_ResidualDipolarCoupling_hh
#define INCLUDED_core_scoring_ResidualDipolarCoupling_hh

#include <core/scoring/ResidualDipolarCoupling.fwd.hh>
// AUTO-REMOVED #include <core/scoring/methods/ResidualDipolarCouplingEnergy.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <basic/datacache/CacheableData.hh>
#include <numeric/numeric.functions.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {

void store_RDC_in_pose(ResidualDipolarCouplingOP, core::pose::Pose&);
ResidualDipolarCouplingOP retrieve_RDC_from_pose(core::pose::Pose&);
ResidualDipolarCouplingCOP retrieve_RDC_from_pose(core::pose::Pose const&);

///@brief ResidualDipolarCouplings are mainly handled by this class
///@detail related classed: RDC --- a single line in an RDC file - representing a single dipolar coupling
///                         ResidualDipolarCouplingEnergy -- an energy method which triggers computations handled by this class.
///
///
class ResidualDipolarCoupling: public basic::datacache::CacheableData {
	//	friend class ResidualDipolarCouplingEnergy;
public:
	// typedefs
	typedef core::Real Real;
	typedef core::Size Size;
	typedef utility::vector1<core::scoring::RDC> RDC_lines;
	typedef core::Real Tensor[3][3];
  typedef core::Real rvec[3];

public:
	/// @brief standard c'stor -- will access option -in:file:rdc to read RDC data
// 	ResidualDipolarCoupling() :
// 		nex_(0), nrows_(0) {
// 		reserve_buffers();
// 		read_RDC_file();
// 	}

	ResidualDipolarCoupling( std::string const& filename = "" ) :
		nex_(0), nrows_(0) {
		reserve_buffers();
		if ( filename != "" ) {
			read_RDC_file( 1, filename );
			release_buffers();
			preprocess_data();
			reserve_buffers();
		} else {
			read_RDC_file();
		}
	}

	/// @brief alternative c'stor if you have a list of RDC lines
	ResidualDipolarCoupling(RDC_lines data_in) :
		All_RDC_lines_(data_in) {
		preprocess_data();
		reserve_buffers();
	}

	//explicit copy c'stor to initialize buffers
	ResidualDipolarCoupling(ResidualDipolarCoupling const& other);
	//explicit assignment operator to initialize buffers
	ResidualDipolarCoupling& operator=(ResidualDipolarCoupling const & other);

	//explicit destructor because we use raw pointers for buffers
	virtual ~ResidualDipolarCoupling();

	//this class lives in the PoseCache.... need to provide clone()
	basic::datacache::CacheableDataOP clone() const {
		return new ResidualDipolarCoupling(*this);
	}

	///@brief compute dipolar score for given pose
	/// will set alignment tensor and force-fields in RDC
	core::Real compute_dipscore(core::pose::Pose const& pose);

	/// fit rdc using RDC data
	core::Real compute_dipscore_nls(core::pose::Pose const& pose);
	core::Real compute_dipscore_nlsDa(core::pose::Pose const& pose, utility::vector1<Real> const tensorDa);
	core::Real compute_dipscore_nlsR(core::pose::Pose const& pose, utility::vector1<Real> const tensorR);
	core::Real compute_dipscore_nlsDaR(core::pose::Pose const& pose, utility::vector1<Real> const tensorDa, utility::vector1<Real> const tensorR);

	//wRDC (like wRMSD .. iter i + 1 tensor weights are ~exp(  - dev^2/sigma ))
	Real iterate_tensor_weights(core::pose::Pose const& pose,
			core::Real sigma2, core::Real tolerance, bool reset);

	void show(std::ostream&) const;

	///@brief read RDC data from file
	void read_RDC_file();

	///@brief fill internal buffers... call always when RDC lines change.
	void preprocess_data();

	///@brief free memory of buffers
	void release_buffers();

	///@brief get memory for buffers
	void reserve_buffers();

	///@brief get the raw RDC data
	inline RDC_lines const& get_RDC_data() const {
		return All_RDC_lines_;
	}

	///@brief return the Q value ( cornilescu )  --- only valid after compute_dipscore
	Real Q() const {
		return R() * sqrt(2.0f);
	}

	///@brief return the R value ( M Clore )  --- only valid after compute_dipscore
	Real R() const {
		return R_;
	}

  //@brief return the number of experiments
	core::Size get_n_alignments() const {return nex_;}

	///@brief return tensor of certain experiment... exp_id starts at 1
	Tensor* tensor( ) {
		return S_;
	}


	/*	core::Real get_fractional_anisotropy(core::Size ex) const {
		//		std::cout << "FFA " << FA_  <<std::endl;
		return FA_[ex] != NULL ? FA_[ex] : -1;
		}*/

	core::Real get_al_tensor_trace(core::Size ex) const {
		//		std::cout << "FFA " << FA_  <<std::endl;
		return trace_[ex];
	}


 core::Real get_al_tensor_max_z(core::Size ex) const {
		//		std::cout << "FFA " << FA_  <<std::endl;
		return maxz_[ex];
	}

	void compute_tensor_stats();
	void show_tensor_stats( std::ostream&, core::Size ex ) const;
	void show_tensor_matrix( std::ostream&, core::Size ex ) const;
	void show_rdc_values( std::ostream&, core::Size ex ) const;

	void show_tensor_stats_nls( std::ostream&, core::Size ex, const double *par) const;
	//void show_tensor_stats_nlsDa( std::ostream&, core::Size ex, const double tensorDa, const double *par) const;
	//void show_tensor_stats_nlsR( std::ostream&, core::Size ex, const double tensorR, const double *par) const;
	//void show_tensor_stats_nlsDaR( std::ostream&, core::Size ex, const double tensorDa, const double tensorR, const double *par) const;
private:
	///@brief read RDC data from file
	void read_RDC_file( Size nex, std::string const& filename );

	///@brief non-const reference to RDC data     private use only.
	RDC_lines& get_RDC_data_ref();

private:
	/// some internal buffers in
	typedef core::Real rvec5[5];
	typedef core::Real Tensor5[5][5];
	RDC_lines All_RDC_lines_;
	rvec* EV_; // rvec3
	rvec5* D_; //rvec5 x nrows
	rvec5* rhs_; //rvec5 x nrows
	Tensor* S_; // 3 x 3 x nex
	Tensor5* T_; // 5 x 5 x nex
	core::Size nex_; //number of alignment media, i.e., how many tensors
	core::Size nrows_;
	core::Real R_; //clore R factor only valid after compute_dipscore...
	core::Real rmsd_;

	//stuff only computed for information purposes -- call compute_tensor_stats()
	Tensor* SD_;//3 x 3 x nex
	Tensor* EIG_;
	core::Real* FA_;//Fractional Anisotropy of the diffusion tensor - NGS
	core::Real* trace_;
	core::Real* maxz_;

  //stuff for nls
	core::Real* r0_;
	core::Real* r1_;
	core::Real* r2_;
	core::Real* exprdc_;
	core::Real* rdcconst_;
	core::Size* lenex_;
};

/////////////////////////////////////////////////
//@brief short class that stores the RDC data lines
/////////////////////////////////////////////////
class RDC {

public:
	enum RDC_TYPE {
		RDC_TYPE_NH = 1, RDC_TYPE_NC, RDC_TYPE_CH, RDC_TYPE_CC
	};

	RDC() {
	}

	RDC(Size res1, std::string const& atom1, Size res2,
			std::string const& atom2, Real Jdipolar, Real weight = 1.0,//for alignment calculation
			Size expid = 1
	//		core::Real Reduced_Jdipolar
	) :
		type_(get_RDC_data_type(atom1, atom2)), res1_(res1), res2_(res2),
				Jdipolar_(Jdipolar), weight_(weight), Jdipolar_computed_(-999),
				expid_(expid), atom1_(atom1), atom2_(atom2)
	//		Reduced_Jdipolar_( Reduced_Jdipolar )
	{
	}

	RDC_TYPE type() const {
		return type_;
	}

	///@brief which type of RDC pairing are we dealing with ?
	RDC_TYPE get_RDC_data_type(std::string const & atom1,
			std::string const & atom2);

	inline Size expid() const {
		return expid_;
	}

	inline Size res1() const {
		return res1_;
	}

	inline Size res2() const {
		return res2_;
	}

	inline Real Jdipolar() const {
		return Jdipolar_;
	}

	Real const& Jcomputed() const {
		return Jdipolar_computed_;
	}

	Real& Jcomputed() {
		return Jdipolar_computed_;;
	}

	Vector fij() const {
		return fij_;
	}

	inline Real fixed_dist() const {
		runtime_assert( 0 ); ///don't use this
		if (type() == RDC_TYPE_NH)
			return 1.02;
		else if (type() == RDC_TYPE_NC)
			return 1.341;
		else if (type() == RDC_TYPE_CH)
			return 1.08;
		else if (type() == RDC_TYPE_CC)
			return 1.525;
		//should never get here...
		runtime_assert( 0 );
		return 0;
	}

	inline Real Reduced_Jdipolar() const {
		using namespace numeric;
		return Jdipolar_ * Dconst() * numeric::cube(fixed_dist());

	}

	Real weight() const {
		return weight_;
	}

	void weight(Real w_in) {
		weight_ = w_in;
	}

	inline Real Dconst() const {
		core::Real Dcnst(0.0);
		if (type() == RDC_TYPE_NH)
    //Concepts in Magnetic Resonance Part A, 21, 10-21
			Dcnst = 36.5089;
		if (type() == RDC_TYPE_NC)
			Dcnst = 9.179532;
		if (type() == RDC_TYPE_CH)
			Dcnst = 90.55324;
			//Dcnst = -90.55324;
		if (type() == RDC_TYPE_CC)
			Dcnst = 22.76805;
			//Dcnst = -22.76805;
		//		std::cout << Dcnst <<std::endl;
				runtime_assert( Dcnst != 0 ); // at this position the type() should already
		return Dcnst;
	}
	friend class ResidualDipolarCoupling;

	std::string const& atom1() const {
		return atom1_;
	}

	std::string const& atom2() const {
		return atom2_;
	}

	void show(std::ostream&) const;


private:
	RDC_TYPE type_;
	Size res1_, res2_;
	Real Jdipolar_, Reduced_Jdipolar_;
	Real weight_;

public :
	Real Jdipolar_computed_;
	core::Vector fij_;

private:
	Size expid_;
	std::string atom1_;
	std::string atom2_;
};

extern std::ostream& operator<<(std::ostream&, ResidualDipolarCoupling const&);
extern std::ostream& operator<<(std::ostream&, RDC const&);

} //scoring
} //core

#endif
