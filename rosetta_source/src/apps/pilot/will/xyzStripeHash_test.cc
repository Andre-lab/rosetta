// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <numeric/geometry/hashing/xyzStripeHash.hh>
#include <protocols/sic_dock/types.hh>
#include <protocols/sic_dock/xyzStripeHashPose.hh>
//#include <apps/pilot/will/gpu/gpu_refold.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <devel/init.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.io.hh>
#include <fstream>
#include <protocols/scoring/ImplicitFastClashCheck.hh>

#include <core/kinematics/Stub.hh>

using core::Size;
using core::id::AtomID;
typedef numeric::xyzVector<double> Vec;

#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#ifdef __MACH__
#include <mach/mach_time.h>
#endif
double time_highres() {
#ifdef __MACH__
  mach_timebase_info_data_t info;
  mach_timebase_info(&info);
  return mach_absolute_time() / 1000000000.0;
  //uint64_t duration = mach_absolute_time();
  //duration *= info.numer;
  //duration /= info.denom;
#else
  timespec tp;
  clock_gettime(CLOCK_REALTIME, &tp);
  return tp.tv_sec + tp.tv_nsec/1000000000.0;
#endif
  return 0;
}

OPT_KEY( Boolean, dump_hash )
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( dump_hash,"dump hash data and quit", false );
}



void dump_points_pdb(utility::vector1<Vec> & p, std::string fn) {
	using namespace ObjexxFCL::fmt;
	std::ofstream o(fn.c_str());
	for(Size i = 1; i <= p.size(); ++i) {
		std::string rn = "VIZ";
		o<<"HETATM"<<I(5,i)<<' '<<" CA "<<' '<<rn<<' '<<"A"<<I(4,i)<<"    "<<F(8,3,p[i].x())<<F(8,3,p[i].y())<<F(8,3,p[i].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	}
	o.close();
}


int main(int argc, char *argv[]) {
	using numeric::geometry::hashing::xyzStripeHash;
	register_options();
	devel::init(argc,argv);

	std::cout << "sizeof(utility::vector1< xyzVector<double> >) " 
	          << sizeof(utility::vector1< numeric::xyzVector<double> >) << " " 
	          << sizeof(xyzStripeHash<double>::ushort2)
	          << std::endl;

	double const DIST(3.5);
	
	core::pose::Pose p;
	core::import_pose::pose_from_pdb(p,basic::options::option[basic::options::OptionKeys::in::file::s]()[1]);
	if(true) {
		for(Size ir = 1; ir <= p.n_residue(); ++ir) {
			if( p.residue(ir).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(p,ir);
			if( p.residue(ir).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(p,ir);		
		}
	}
	
	protocols::scoring::ImplicitFastClashCheck ifc(p,DIST);

	protocols::sic_dock::xyzStripeHashPose xyzhash(DIST,p,protocols::sic_dock::HVY);
	xyzhash.sanity_check();

	// xyzStripeHash<double>::float4 const * j = xyzhash.grid_atoms_;
	// for(xyzStripeHash<double>::xyz_iterator i = xyzhash.xyz_begin(); i != xyzhash.xyz_end(); ++i){
	// 	numeric::xyzVector<double> const & a = *i;
	// 	xyzStripeHash<double>::float4 const & b = *(j++);
	// 	std::cout << a.x() <<" "<< b.x <<" "<< a.y() <<" "<< b.y <<" "<< a.z() <<" "<< b.z << std::endl;
	// 	if( a.x() != b.x && a.y() != b.y && a.z() != b.z ){			
	// 		utility_exit_with_message("BAD!!!!!");
	// 	}
	// }

	// utility_exit_with_message("FOO");

	if(basic::options::option[basic::options::OptionKeys::dump_hash]()) {
		using namespace ObjexxFCL::fmt;
		std::cout << xyzhash.grid_size() << std::endl;
		std::cout << xyzhash.natom() << std::endl;
		std::cout << xyzhash.xdim() << " " << xyzhash.ydim() << " " << xyzhash.zdim() << std::endl;
		for(int j = 0; j < (int)xyzhash.natom(); ++j) {
			std::cout << F(12,7,xyzhash.grid_atoms()[j].x) << " ";
			std::cout << F(12,7,xyzhash.grid_atoms()[j].y) << " ";
			std::cout << F(12,7,xyzhash.grid_atoms()[j].z) << std::endl;
		}
		for(int j = 0; j < (int)xyzhash.xdim()*xyzhash.ydim()*xyzhash.zdim(); ++j) {
			std::cout << xyzhash.grid_stripe()[j].x << " ";
			std::cout << xyzhash.grid_stripe()[j].y << std::endl;
		}
		return 0;
	}

	for(Size ir = 1; ir <= p.n_residue(); ++ir) {
		for(Size ia = 1; ia <= p.residue_type(ir).natoms(); ++ia) {
			p.set_xyz(AtomID(ia,ir),p.xyz(AtomID(ia,ir))+xyzhash.translation());
		}
	}
	// p.dump_pdb("input_trans.pdb");
	

	Vec mn(9e9),mx(-9e9);
	int real_natom = 0.0;
	for(Size ir = 1; ir <= p.n_residue(); ++ir) {
		for(Size ia = 1; ia <= p.residue(ir).natoms(); ++ia) {
			real_natom++;
			// mn.x() = numeric::min(mn.x(),p.xyz(AtomID(ia,ir)).x());
			// mn.y() = numeric::min(mn.y(),p.xyz(AtomID(ia,ir)).y());
			// mn.z() = numeric::min(mn.z(),p.xyz(AtomID(ia,ir)).z());
			// mx.x() = numeric::max(mx.x(),p.xyz(AtomID(ia,ir)).x());
			// mx.y() = numeric::max(mx.y(),p.xyz(AtomID(ia,ir)).y());
			// mx.z() = numeric::max(mx.z(),p.xyz(AtomID(ia,ir)).z());
		}
	}
	// mn -= Vec(15,15,15);
	// mx += Vec(15,15,15);
	mn = Vec(-10,-10,-10);
	mx = Vec( 50, 50, 50);

	std::cout << "ncells: " << xyzhash.xdim()*xyzhash.ydim()*xyzhash.zdim() << std::endl;
	std::cout << "Bounds: " << mn << " " << mx << std::endl;
	std::cout << "Natom: " << real_natom << " " << xyzhash.natom() << std::endl;
	{
		utility::vector1<Vec> hashpts;
		for(int i = 0; i < (int)xyzhash.natom(); ++i) {
			hashpts.push_back( Vec( xyzhash.grid_atoms()[i].x, xyzhash.grid_atoms()[i].y, xyzhash.grid_atoms()[i].z ) );
		}
		// dump_points_pdb(hashpts,"input_xyzhash.pdb");
	}
	
	
	double tifc=0.0,th=0.0,ts=0.0,t=0.0;
	int tot = 0;
	for(int i = 0; i < (int)basic::options::option[basic::options::OptionKeys::out::nstruct](); ++i) {
		Vec rv( numeric::random::uniform(),numeric::random::uniform(),numeric::random::uniform() );
		rv.x() = (mx.x()-mn.x()) * rv.x() + mn.x();
		rv.y() = (mx.y()-mn.y()) * rv.y() + mn.y();
		rv.z() = (mx.z()-mn.z()) * rv.z() + mn.z();		
		// float const rx = rv.x();
		// float const ry = rv.y();
		// float const rz = rv.z();

		numeric::geometry::hashing::Counter<double> counter;

		int hash_nbcount;
		t = time_highres();
		hash_nbcount = xyzhash.nbcount(rv);
		th += time_highres()-t;

		Vec tmp = rv+xyzhash.translation();
		t = time_highres();
		int safe_nbcount = 0.0;
		for(int j = 0; j < (int)xyzhash.natom(); ++j) {
			double const & hx = xyzhash.grid_atoms()[j].x;
			double const & hy = xyzhash.grid_atoms()[j].y;
			double const & hz = xyzhash.grid_atoms()[j].z;
			if( (hx-tmp.x())*(hx-tmp.x()) + (hy-tmp.y())*(hy-tmp.y()) + (hz-tmp.z())*(hz-tmp.z()) <= DIST*DIST ) safe_nbcount++;
		}
		ts += time_highres()-t;

		t = time_highres();
		int ifc_nbcount = ifc.clash_count(rv);
		tifc += time_highres()-t;

		tot += hash_nbcount;
		if( safe_nbcount != hash_nbcount /*|| safe_nbcount != ifc_nbcount*/ ) {
//			if(rv.x() > 36.0) continue; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! BUG!!!!!!!!!!!!!!!!!!!!!!!!!
			std::cout << rv << std::endl;
			std::cout << safe_nbcount << " " << hash_nbcount << " " << ifc_nbcount << std::endl;			
			utility_exit_with_message("FAIL after "+ObjexxFCL::fmt::I(10,i)+"!!!");
		}
		
	}
	std::cout << ObjexxFCL::fmt::I(10,basic::options::option[basic::options::OptionKeys::out::nstruct]())+" nb counts of random points match. Woot! " << tot << std::endl;
	std::cout << "hash speedup over N^2 count: " << ts  /th << std::endl;
	std::cout << "hash speedup over IFC count: " << tifc/th << std::endl;

	return 0;
}
