// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:

#include <protocols/sic_dock/scores/MotifHashRigidScore.hh>
	#include <protocols/sic_dock/util.hh>
	#include <protocols/motif_hash/motif_hash_stuff.hh>

	#include <basic/options/keys/sicdock.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/option.hh>
	#include <basic/options/option_macros.hh>
	#include <numeric/constants.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/string_util.hh>
	#include <core/id/AtomID_Map.hh>
	#include <core/pose/util.hh>
	#include <basic/Tracer.hh>
	#include <boost/foreach.hpp>
	#include <protocols/fldsgn/topology/SS_Info2.hh>
	#include <core/scoring/dssp/Dssp.hh>

namespace protocols {
namespace sic_dock {
namespace scores {

static basic::Tracer TR( "protocols.sic_dock.scores.MotifHashRigidScore" );

using core::Size;
	using core::Real;
	using numeric::min;
	using core::id::AtomID;
	using std::cout;
	using std::endl;
	using utility::vector1;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef Xform Xform;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	typedef utility::vector1<Vec> Vecs;
	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef numeric::Xforms Xforms;
	typedef utility::vector1<RigidScoreCOP> Scores;


MotifHashRigidScore::MotifHashRigidScore(
	Pose const & _pose1,
	Pose const & _pose2
):
	pose1_(_pose1),pose2_(_pose2),mh_(NULL),xs_(NULL),xsee_(NULL),xseh_(NULL),xshe_(NULL),xshh_(NULL),xspp_(NULL)
{
	core::scoring::dssp::Dssp(pose1_).insert_ss_into_pose_no_IG_helix(pose1_);
	core::scoring::dssp::Dssp(pose2_).insert_ss_into_pose_no_IG_helix(pose2_);
	xs_   = protocols::motif_hash::MotifHashManager::get_instance()->xform_score_from_cli();
	xshh_ = protocols::motif_hash::MotifHashManager::get_instance()->xform_score_hh_from_cli();
	xseh_ = protocols::motif_hash::MotifHashManager::get_instance()->xform_score_eh_from_cli();
	xshe_ = protocols::motif_hash::MotifHashManager::get_instance()->xform_score_he_from_cli();
	xsee_ = protocols::motif_hash::MotifHashManager::get_instance()->xform_score_ee_from_cli();
	xspp_ = protocols::motif_hash::MotifHashManager::get_instance()->xform_score_sspair_from_cli();
	mh_ = protocols::motif_hash::MotifHashManager::get_instance()->motif_hash_from_cli();
	for(Size ir = 1; ir <= pose1_.total_residue(); ++ir){
		// if( pose1_.secstruct(ir)=='L' ) continue;
		Vec  N = pose1_.residue(ir).xyz(1);
		Vec CA = pose1_.residue(ir).xyz(2);
		Vec  C = pose1_.residue(ir).xyz(3);
		bbx1_.push_back(Xform(CA,N,CA,C));
	}
	for(Size ir = 1; ir <= pose2_.total_residue(); ++ir){
		// if( pose2_.secstruct(ir)=='L' ) continue;
		Vec  N = pose2_.residue(ir).xyz(1);
		Vec CA = pose2_.residue(ir).xyz(2);
		Vec  C = pose2_.residue(ir).xyz(3);
		bbx2_.push_back(Xform(CA,N,CA,C));
	}
	// cout << pose1_.secstruct() << endl;
	// cout << pose2_.secstruct() << endl;
	// ssinfo1_ = new protocols::fldsgn::topology::SS_Info2(pose1_,pose1_.secstruct());
	// ssinfo2_ = new protocols::fldsgn::topology::SS_Info2(pose2_,pose2_.secstruct());
	// nss1_=0;
	// for(Size ir = 1; ir <= pose1_.n_residue(); ++ir){
	// 	// cout << ir << " " << ssinfo1_->ss_element_id(ir) << endl;
	// 	nss1_ = std::max(nss1_,ssinfo1_->ss_element_id(ir));
	// }
	// // utility_exit_with_message("foo");
	// nss2_=0;
	// for(Size ir = 1; ir <= pose2_.n_residue(); ++ir){
	// 	// cout << ir << " " << ssinfo2_->ss_element_id(ir) << endl;
	// 	nss2_ = std::max(nss2_,ssinfo2_->ss_element_id(ir));
	// }
}

MotifHashRigidScore::~MotifHashRigidScore(){
	if(ssinfo1_) delete ssinfo1_;
	if(ssinfo2_) delete ssinfo2_;
}


core::Real
MotifHashRigidScore::score_meta( Xforms const & x1s, Xforms const & x2s, int & nsheet, Real & rawscore, Real & sselem_score, Real & coverage, Real & res_score, Real & sheetsc, int &nres, int &Nhh, int &Nee, int &Neh, int &Nh, int &Ne, int &Nl ) const {
	using namespace protocols::motif_hash;
	utility::vector1<Real> sselemsc1(nss1_,0.0),sselemsc2(nss2_,0.0);
	Real tot_weighted=0.0, totscore=0.0;
	std::map<Size,Real> ssp1,ssp2,mres1,mres2;
	std::set<Size> tres1,tres2;
	Nhh=0; Neh=0; Nee=0; nres=0;
	BOOST_FOREACH(Xform const & x1,x1s){
		BOOST_FOREACH(Xform const & x2,x2s){
			for(Size ir = 1; ir <= bbx1_.size(); ++ir){
				char ss1 = pose1_.secstruct(ir);
				Xform const xb1 = x1 * bbx1_[ir];
				for(Size jr = 1; jr <= bbx2_.size(); ++jr){
					Xform const xb2 = x2 * bbx2_[jr];
					Xform x = ~xb1 * xb2;
					if(x.t.length_squared() > 100.0 ) continue;
					if( x.t.length_squared() < 64.0 ){ tres1.insert(ir); tres2.insert(jr); }
					char ss2 = pose2_.secstruct(jr);
					protocols::motif_hash::XformScoreCAP xscore_for_this_ss = xs_;
					if(option[mh::score::noloops]() && ss1=='L') continue;
					if(option[mh::score::noloops]() && ss2=='L') continue;
					if(ss1=='H'&&ss2=='H') xscore_for_this_ss = xshh_;
					if(ss1=='E'&&ss2=='H') xscore_for_this_ss = xseh_;
					if(ss1=='H'&&ss2=='E') xscore_for_this_ss = xshe_;
					if(ss1=='E'&&ss2=='E') xscore_for_this_ss = xsee_;
					if(x.t.length_squared() <  49.0 ) tot_weighted -= 1.0;
					Vec e = x.euler_angles_deg();
					Real6 rt6;
					rt6[1] = x.t.x(); rt6[2] = x.t.y(); rt6[3] = x.t.z();
					rt6[4] =   e.x(); rt6[5] =   e.y(); rt6[6] =   e.z();
					Real raw = min(500.0,xscore_for_this_ss->score_of_bin(rt6));
					if(raw <= 0) continue;
					totscore += raw;
					tot_weighted += sqrt( raw );
					// sselemsc1[ssinfo1_->ss_element_id(ir)] += sqrt(raw)+raw/30.0;
					// sselemsc2[ssinfo2_->ss_element_id(jr)] += sqrt(raw)+raw/30.0;//!!!!!!!!!!!!!!!!!!!!!!!!
					if( x.t.length_squared() <  64.0 ){
						if(mres1.find(ir)==mres1.end()) mres1[ir]=0; mres1[ir] += raw;
						if(mres2.find(jr)==mres2.end()) mres2[jr]=0; mres2[jr] += raw;
					}
					if( xspp_ && pose1_.secstruct(ir)=='E' && pose2_.secstruct(jr)=='E' ){
						Real ss = xspp_->score_of_bin(rt6);
						if(ss>100.0){
							if(ssp1.find(ir)==ssp1.end()) ssp1[ir]=0; ssp1[ir] += ss;
							if(ssp2.find(jr)==ssp2.end()) ssp2[jr]=0; ssp2[jr] += ss;
						}
					}
					if(ss1=='H'&&ss2=='H') ++Nhh;
					if(ss1=='E'&&ss2=='H') ++Neh;
					if(ss1=='H'&&ss2=='E') ++Neh;
					if(ss1=='E'&&ss2=='E') ++Nee;
				}
			}
		}
	}

	Real ssscore = 0;
	// BOOST_FOREACH(Real s,sselemsc1) ssscore += sqrt(s) + s/15.0;
	// BOOST_FOREACH(Real s,sselemsc2) ssscore += sqrt(s) + s/15.0;

	res_score = tres1.size()+tres2.size();
	res_score = -res_score * 0.666;
	typedef std::map<Size,Real>::value_type  MapVal;
	BOOST_FOREACH(MapVal v,mres1) res_score += sqrt(v.second);
	BOOST_FOREACH(MapVal v,mres2) res_score += sqrt(v.second);

	Nh=0; Ne=0; Nl=0;
	BOOST_FOREACH(Size n,tres1){ Nh += 'H'==pose1_.secstruct(n); Ne += 'E'==pose1_.secstruct(n); Nl += 'L'==pose1_.secstruct(n); }
	BOOST_FOREACH(Size n,tres2){ Nh += 'H'==pose2_.secstruct(n); Ne += 'E'==pose2_.secstruct(n); Nl += 'L'==pose2_.secstruct(n); }

	coverage = Real(mres1.size()+mres2.size())/Real(tres1.size()+tres2.size());

	Size nsheetres = ssp1.size() + ssp2.size();
	Real nsheetsc = nsheetres / 2.0;

	sselem_score = ssscore;
	rawscore = tot_weighted;
	nsheet = nsheetres;

	sheetsc = option[basic::options::OptionKeys::mh::score::strand_pair_weight]() * ( 30*nsheetsc );

	return res_score + sheetsc;
}

core::Real
MotifHashRigidScore::score( Xforms const & x1, Xforms const & x2 ) const {
	int a,f,g,h,i,j,k,l;
	Real b,c,d,e,m;
	return score_meta(x1,x2,a,b,c,d,e,m,f,g,h,i,j,k,l);
}

int
MotifHashRigidScore::dump_matching_motifs( Pose const & pose1, Pose const & pose2, std::ostream & out, int & count, xyzStripeHashPoseCAP clash_check, bool print ) const {
	return mh_->dump_matching_motifs(pose1,pose2,option[mh::motif_match_radius](),out,count,clash_check,print) +
	       mh_->dump_matching_motifs(pose2,pose1,option[mh::motif_match_radius](),out,count,clash_check,print) ;
}
int
MotifHashRigidScore::dump_matching_motifs( Xforms const & x1s, Xforms const & x2s, std::ostream & out, xyzStripeHashPoseCAP clash_check, bool print ) const {
	if(!mh_) mh_ = protocols::motif_hash::MotifHashManager::get_instance()->motif_hash_from_cli();
	int nhit = 0;
	BOOST_FOREACH(Xform const & x1,x1s){
		Pose pose1(pose1_);
		xform_pose(pose1,x1);
		BOOST_FOREACH(Xform const & x2,x2s){
			Pose pose2(pose2_);
			xform_pose(pose2,x2);
			nhit += dump_matching_motifs(pose1,pose2,out,nhit,clash_check,print);
		}
	}
	return nhit;
}

void
MotifHashRigidScore::show(
	std::ostream & out,
    int width
) const {
	if(!mh_) mh_ = protocols::motif_hash::MotifHashManager::get_instance()->motif_hash_from_cli();
	Pose pose1,pose2;
	Stats stats;
	mh_->stat_matching_motifs(pose1,pose2,stats);
	out << " " << ObjexxFCL::fmt::RJ(width,"MH_RAW");
	out << " " << ObjexxFCL::fmt::RJ(width,"MH_SPREAD");
	out << " " << ObjexxFCL::fmt::RJ(width,"MH_RES");
	// out << " " << ObjexxFCL::fmt::RJ(width,"MH_SSPAIR");
	out << " " << ObjexxFCL::fmt::RJ(2,"EE");
	out << " " << ObjexxFCL::fmt::RJ(width,"MH_COVER");
	out << " " << ObjexxFCL::fmt::RJ(4,"NRES");
	out << " " << ObjexxFCL::fmt::RJ(3,"NH"  );
	out << " " << ObjexxFCL::fmt::RJ(3,"NE"  );
	out << " " << ObjexxFCL::fmt::RJ(3,"NL"  );
	out << " " << ObjexxFCL::fmt::RJ(5,"P_EE"  );
	out << " " << ObjexxFCL::fmt::RJ(5,"P_EH"  );
	out << " " << ObjexxFCL::fmt::RJ(5,"P_HH"  );
	BOOST_FOREACH(Stats::value_type val,stats){
		out << " " << ObjexxFCL::fmt::RJ(width,val.first);
	}

}
void
MotifHashRigidScore::show(
	std::ostream & out,
	Xforms const & x1s,
	Xforms const & x2s,
	xyzStripeHashPoseCAP clash_check,
	int width
) const {
	if(!mh_) mh_ = protocols::motif_hash::MotifHashManager::get_instance()->motif_hash_from_cli();

	Real rawscore,ssscore,coverage,res_score,sheetsc;
	int nsheetres,nres,Nhh,Nee,Neh,Nh,Ne,Nl;
	score_meta(x1s,x2s,nsheetres,rawscore,ssscore,coverage,res_score,sheetsc,nres,Nhh,Nee,Neh,Nh,Ne,Nl);
	out << " " << ObjexxFCL::fmt::F(width,4,rawscore);
	out << " " << ObjexxFCL::fmt::F(width,4,ssscore);
	out << " " << ObjexxFCL::fmt::F(width,4,res_score);
	// out << " " << ObjexxFCL::fmt::F(width,4,sheetsc);
	out << " " << ObjexxFCL::fmt::I(2,nsheetres);
	out << " " << ObjexxFCL::fmt::F(width,4,coverage);
	out << " " << ObjexxFCL::fmt::I(4,nres);
	out << " " << ObjexxFCL::fmt::I(3,Nh);
	out << " " << ObjexxFCL::fmt::I(3,Ne);
	out << " " << ObjexxFCL::fmt::I(3,Nl);
	out << " " << ObjexxFCL::fmt::F(5,3,Real(Nee)/Real(Nee+Neh+Nhh));
	out << " " << ObjexxFCL::fmt::F(5,3,Real(Neh)/Real(Nee+Neh+Nhh));
	out << " " << ObjexxFCL::fmt::F(5,3,Real(Nhh)/Real(Nee+Neh+Nhh));

	Stats stats;
	BOOST_FOREACH(Xform const & x1,x1s){
		BOOST_FOREACH(Xform const & x2,x2s){
			Pose pose1(pose1_),pose2(pose2_);
			xform_pose(pose1,x1);
			xform_pose(pose2,x2);
			mh_->stat_matching_motifs(pose1,pose2,stats,clash_check,option[mh::motif_match_radius]());
			mh_->stat_matching_motifs(pose2,pose1,stats,clash_check,option[mh::motif_match_radius]());
		}
	}


	stats["M_EE"     ] = stats["M_EE"     ]/stats["M_NUM"];
	stats["M_HE"     ] = stats["M_HE"     ]/stats["M_NUM"];
	// stats["M_EL"     ] = stats["M_EL"     ]/stats["M_NUM"];
	stats["M_HH"     ] = stats["M_HH"     ]/stats["M_NUM"];
	// stats["M_HL"     ] = stats["M_HL"     ]/stats["M_NUM"];
	// stats["M_LL"     ] = stats["M_LL"     ]/stats["M_NUM"];
	// stats["M_SSPAIR"] = stats["M_SSPAIR"]/stats["M_NUM"];

	BOOST_FOREACH(Stats::value_type val,stats){
		out << " " << ObjexxFCL::fmt::F(width,5,val.second);
	}
}
void
MotifHashRigidScore::show(
	std::ostream & out,
	Xforms const & x1,
	Xforms const & x2,
	int width
) const {
	show(out,x1,x2,NULL,width);
}



} // namespace scores
} // namespace sic_dock
} // namespace protocols
