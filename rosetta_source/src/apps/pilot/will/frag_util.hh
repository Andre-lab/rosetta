#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>

void read_fragdata( vector1< core::fragment::FragDataOP > & fds, std::istream & in ) {
	using namespace core::fragment;
	Size n,count=0;
	while( in >> n ) {
	 	string pdb;
		char buf[999];
		FragDataOP fd = new FragData;
		for( Size i = 1; i <= n; ++i ) {
			utility::pointer::owning_ptr<SingleResidueFragData> srfd;
			srfd = new BBTorsionSRFD;
			in >> pdb;
			in.getline(buf,999);
			std::istringstream iss(buf);
			iss >> *srfd;
			fd->add_residue(srfd);
		}
		fd->set_valid(true);
		fds.push_back(fd);
		count++;
		// if( count >= nread ) break;
	}
}

std::map<string, vector1<core::fragment::FragDataOP> >
get_frags_map( ) {
	using namespace core::fragment;
	//cout << "reading frags" << std::endl;
	utility::io::izstream in;
	std::map<string,vector1<FragDataOP> > fds;
	basic::database::open(in,"ss_fragfiles/EEE.fragfile"); read_fragdata(fds["EEE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/EEH.fragfile"); read_fragdata(fds["EEH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/EEL.fragfile"); read_fragdata(fds["EEL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/EHH.fragfile"); read_fragdata(fds["EHH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/ELE.fragfile"); read_fragdata(fds["ELE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/ELH.fragfile"); read_fragdata(fds["ELH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/ELL.fragfile"); read_fragdata(fds["ELL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HEE.fragfile"); read_fragdata(fds["HEE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HEH.fragfile"); read_fragdata(fds["HEH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HEL.fragfile"); read_fragdata(fds["HEL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HHE.fragfile"); read_fragdata(fds["HHE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HHH.fragfile"); read_fragdata(fds["HHH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HHL.fragfile"); read_fragdata(fds["HHL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HLE.fragfile"); read_fragdata(fds["HLE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HLH.fragfile"); read_fragdata(fds["HLH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/HLL.fragfile"); read_fragdata(fds["HLL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LEE.fragfile"); read_fragdata(fds["LEE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LEH.fragfile"); read_fragdata(fds["LEH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LEL.fragfile"); read_fragdata(fds["LEL"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LHH.fragfile"); read_fragdata(fds["LHH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LLE.fragfile"); read_fragdata(fds["LLE"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LLH.fragfile"); read_fragdata(fds["LLH"],in); in.close();
	basic::database::open(in,"ss_fragfiles/LLL.fragfile"); read_fragdata(fds["LLL"],in); in.close();
	return fds;
}


core::fragment::FragSetOP make_frag_set(std::string ss, std::map<string, vector1<core::fragment::FragDataOP> > & fds) {
	using namespace core::fragment;
	FragSetOP frags = new ConstantLengthFragSet();
	int const stop = ss.size() + 1 - 3;
	if((int)1 >= stop) return NULL;
	for( Size i = 1; i <= (Size)stop; ++i ) {
		string ss3 = ss.substr(i-1,3);
		bool mkframe = true;
		for(Size j = 0; j < ss3.size(); ++j) if(ss3[j]!='H'&&ss3[j]!='E'&&ss3[j]!='L'&&ss3[j]!='_') mkframe = false;
		if(!mkframe) utility_exit_with_message("oaisnrt");
		FrameOP frame = new Frame(i,3);
		vector1<char> ss0,ss1,ss2;
		if('_'==ss3[0]) { ss0.push_back('H'); ss0.push_back('E'); ss0.push_back('L'); } else ss0.push_back(ss3[0]);
		if('_'==ss3[1]) { ss1.push_back('H'); ss1.push_back('E'); ss1.push_back('L'); } else ss1.push_back(ss3[1]);
		if('_'==ss3[2]) { ss2.push_back('H'); ss2.push_back('E'); ss2.push_back('L'); } else ss2.push_back(ss3[2]);				
		int nfrag = 0;
		for( Size j = 1; j <= ss0.size(); ++j ) {
		for( Size k = 1; k <= ss1.size(); ++k ) {
		for( Size l = 1; l <= ss2.size(); ++l ) {
			string ss=""; ss+=ss0[j]; ss+=ss1[k]; ss+=ss2[l];
			//cout << "adding ss " << ss << " '" << ss0[j] << "' '" << ss1[k] << "' '" << ss2[l] << "'" << std::endl;
			vector1<FragDataOP>::iterator beg = fds[ss].begin();
			vector1<FragDataOP>::iterator end = fds[ss].end();			
			for( vector1<FragDataOP>::iterator fi = beg; fi != end; ++fi ) {
				frame->add_fragment(*fi);
				nfrag++;
			}
		}}}
		frags->add(frame);
		//coutcout << "make frag " << i << ": " << ss3 << std::endl;
		//cout << ss3 << " " << nfrag << " fragments" << endl;
	}
	if(frags->size() == 0) utility_exit_with_message("no fragments!!!");
	return frags;
}

core::fragment::FragSetOP make_frag_set(vector1<char> ss, std::map<string, vector1<core::fragment::FragDataOP> > & fds) {
	string s = "";
	for(Size i = 1; i <= ss.size(); ++i) s += ss[i];
	return make_frag_set(s,fds);
}

core::fragment::FragSetOP make_frag_set_9mers(Size nres) {
	using namespace core::fragment;
	vector1<core::fragment::FragDataOP> fds9;
	std::ifstream in("input/loop_helix.9mers");
	read_fragdata(fds9,in);

	FragSetOP frags = new ConstantLengthFragSet();
	for( Size i = 1; i <= nres-8; ++i ) {
		FrameOP frame = new Frame(i,9);
		for( vector1<FragDataOP>::iterator fi = fds9.begin(); fi != fds9.end(); ++fi ) {
			frame->add_fragment(*fi);
		}
		frags->add(frame);
	}
	if(frags->size() == 0) return NULL;
	return frags;
}
