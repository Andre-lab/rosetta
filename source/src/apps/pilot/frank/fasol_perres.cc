// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>

#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/residue_datacache.hh>


#include <basic/options/option.hh>

#include <devel/init.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>


#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


#include <utility/excn/Exceptions.hh>


#include <basic/options/option_macros.hh>


// C++ headers
//#include <cstdlib>
#include <iostream>
#include <string>

//Auto Headers
#include <basic/Tracer.hh>

#include <protocols/moves/Mover.hh> // AUTO IWYU For Mover
#include <core/conformation/Residue.hh> // AUTO IWYU For Pose::Residue Residue

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::scoring;
using namespace core::chemical;
using namespace core::scoring::etable;
using namespace core::scoring::etable::count_pair;
using namespace core::scoring::lkball;
using namespace core::scoring::methods;

static basic::Tracer TR("fasol_refit");

OPT_1GRP_KEY(Real, fasol, del_ddg)
OPT_1GRP_KEY(Real, fasol, max_sasa)
OPT_1GRP_KEY(Boolean, fasol, revert)
OPT_1GRP_KEY(Boolean, fasol, classic)

class FaSolReporter : public protocols::moves::Mover {
public:
	FaSolReporter()= default;
	void apply( core::pose::Pose & pose) override {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();
		(*sf)(pose);

		// get per-atom SASA
		core::Real probe_radius=1.4;
		core::id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;
		core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius );

		EnergyGraph const & energy_graph( pose.energies().energy_graph() );
		methods::EnergyMethodOptions e_opts = sf->energy_method_options();
		Etable etable = *( ScoringManager::get_instance()->etable( e_opts ).lock() ); // copy

		core::scoring::lkball::LK_BallEnergy lkb(e_opts);

		// atomtype used for burial calcs
		// kind of hacky but use "Bsp2" as generic atom
		//   (typically not optimized, lambda=3.5, dgfree=-4, small vdw radius~2)
		core::chemical::AtomTypeSetCOP atom_set_ac( etable.atom_set() );
		int bur_type=atom_set_ac->atom_type_index("Bsp2");

		// weights
		core::Real weightS=sf->get_weight( core::scoring::fa_sol );
		core::Real weightA=sf->get_weight( core::scoring::fa_atr );
		core::Real weightR=sf->get_weight( core::scoring::fa_rep );
		core::Real weightLKbw=sf->get_weight( core::scoring::lk_ball_wtd );
		core::Real weightLKb=sf->get_weight( core::scoring::lk_ball );
		core::Real weightLKbi=sf->get_weight( core::scoring::lk_ball_iso );
		core::Real weightLKbr=sf->get_weight( core::scoring::lk_ball_bridge );
		core::Real weightLKbru=sf->get_weight( core::scoring::lk_ball_bridge_uncpl );

		for ( Size ires=1; ires<= pose.size(); ++ires ) {
			core::conformation::Residue const & rsd1( pose.residue(ires) );
			if ( !rsd1.is_protein() ) continue;
			if ( rsd_sasa[ires] != 0 ) continue;

			WaterCoords rsd1_waters, rsd2_waters;
			utility::vector1< AtomWeights > rsd1_atom_wts;
			utility::vector1< Size > rsd1_n_attached_waters, rsd2_n_attached_waters;
			utility::vector1< Size > rsd1_water_offsets, rsd2_water_offsets;

			if ( weightLKbw != 0 || weightLKb != 0 || weightLKbi != 0 || weightLKbr != 0 || weightLKbru != 0 ) {
				lkball::LKB_ResidueInfo const &lkbinfo1 = static_cast< lkball::LKB_ResidueInfo const & > (
					rsd1.data_ptr()->get( conformation::residue_datacache::LK_BALL_INFO ));
				rsd1_n_attached_waters = lkbinfo1.n_attached_waters();
				rsd1_waters = ( lkbinfo1.waters() );
				rsd1_atom_wts = ( lkbinfo1.atom_weights() );
				rsd1_water_offsets = lkbinfo1.water_offset_for_atom();
			}

			core::Real fa_sol_i = 0.0, lk_ball_i = 0.0, fa_atr_i=0, fa_rep_i=0;
			core::Real burial_i = 0.0;

			for ( Size iatm=rsd1.first_sidechain_atom(); iatm<= rsd1.natoms(); ++iatm ) {
				if ( iatm > rsd1.nheavyatoms() && iatm < rsd1.first_sidechain_hydrogen() ) continue;

				core::conformation::Atom a_i = rsd1.atom(iatm);
				a_i.type(bur_type);

				Size n_atom1_waters( 0 );
				Size atom1_offset( 0 );
				//WaterCoords atom1_waters;
				AtomWeights atom1_wts;
				numeric::xyzVector< core::Real > atom1_xyz( rsd1.xyz( iatm ) );

				if ( iatm <= rsd1.nheavyatoms() && rsd1_waters.size() > 0 ) {
					n_atom1_waters = rsd1_n_attached_waters[ iatm ];
					//atom1_waters = rsd1_waters[ iatm ];
					atom1_wts = rsd1_atom_wts[iatm];
					atom1_offset = rsd1_water_offsets[iatm];
				}

				Real const sasa_this_atom( atom_sasa[ core::id::AtomID( iatm, ires ) ] );
				if ( sasa_this_atom > option[fasol::max_sasa] ) continue;

				for ( utility::graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(ires)->const_edge_list_begin(),
						irue = energy_graph.get_node(ires)->const_edge_list_end();
						iru != irue; ++iru ) {
					auto const * edge( static_cast< EnergyEdge const *> (*iru) );
					auto jres=(int)edge->get_other_ind(ires);
					core::conformation::Residue const &rsd2( pose.residue(jres) );
					if ( !rsd2.is_protein() ) continue;

					if ( weightLKbw != 0 || weightLKb != 0 || weightLKbi != 0 || weightLKbr != 0 || weightLKbru != 0 ) {
						lkball::LKB_ResidueInfo const &lkbinfo2 = static_cast< lkball::LKB_ResidueInfo const & > (
							rsd2.data_ptr()->get( conformation::residue_datacache::LK_BALL_INFO ));
						rsd2_n_attached_waters = lkbinfo2.n_attached_waters();
						rsd2_waters = ( lkbinfo2.waters() );
						rsd2_water_offsets = lkbinfo2.water_offset_for_atom();
					}

					// count pair
					CountPairFunctionOP cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );

					for ( Size jatm=1; jatm<= rsd2.natoms(); ++jatm ) {
						numeric::xyzVector< core::Real > const & atom2_xyz( rsd2.xyz( jatm ) );

						//WaterCoords atom2_waters;
						Size n_atom2_waters = rsd2_n_attached_waters[ jatm ];
						Size atom2_offset = rsd2_water_offsets[ jatm ];
						//if ( jatm <= rsd2.nheavyatoms()  && n_atom2_waters > 0 ) {
						// atom2_waters = rsd2_waters[ jatm ];
						//}

						core::Size path_dist;
						core::Real weight = 1;
						if ( !cpfxn->count( iatm, jatm, weight, path_dist ) ) continue;

						core::Real fasol1,fasol2, ljatr, ljrep, lkjunk, dis2;

						etable.analytic_etable_evaluation( rsd1.atom(iatm), rsd2.atom(jatm), ljatr, ljrep, lkjunk, dis2);
						fa_atr_i += 0.5*(weightA*ljatr);
						fa_rep_i += 0.5*(weightR*ljrep);

						if ( iatm <= rsd1.nheavyatoms() && jatm <= rsd2.nheavyatoms() ) {
							etable.analytic_lk_energy(rsd1.atom(iatm), rsd2.atom(jatm), fasol1,fasol2 );
							fa_sol_i += weight*weightS*fasol1;

							if ( n_atom1_waters != 0 ) {
								Real const fasol1_lkball =
									fasol1 * lkb.get_lk_fractional_contribution( atom2_xyz, rsd2.atom(jatm).type(), n_atom1_waters, atom1_offset, rsd1_waters );
								lk_ball_i += weight*weightLKbw * ( atom1_wts[1] * fasol1 + atom1_wts[2] * fasol1_lkball );
								lk_ball_i += weight*weightLKb * ( fasol1_lkball );
								lk_ball_i += weight*weightLKbi * ( fasol1);

								if ( n_atom2_waters != 0 ) {
									core::Real fasol1_lkbridge = lkb.get_lkbr_fractional_contribution(
										atom1_xyz, atom2_xyz,
										n_atom1_waters, n_atom2_waters,
										atom1_offset, atom2_offset,
										rsd1_waters, rsd2_waters );
									lk_ball_i += 0.5 * weight * weightLKbr * (fasol1+fasol2) * fasol1_lkbridge;
									lk_ball_i += 0.5 * weight * weightLKbru * fasol1_lkbridge;
								}
							}

							etable.analytic_lk_energy(a_i, rsd2.atom(jatm), fasol1,fasol2 );
							burial_i += weight*fasol1 / -etable.lk_dgfree( bur_type );
							//TR << etable.lk_dgfree( bur_type ) << " " << fasol1 << std::endl;
						}
					}
				}

				// add fa_intra_* contribution
				CountPairFunctionOP cpfxn =
					CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_4 );
				for ( Size jatm=1; jatm<= rsd1.natoms(); ++jatm ) {
					core::Real weightS=sf->get_weight( core::scoring::fa_intra_sol_xover4 );
					core::Real weightA=sf->get_weight( core::scoring::fa_intra_atr_xover4 );
					core::Real weightR=sf->get_weight( core::scoring::fa_intra_rep_xover4 );
					core::Size path_dist;
					core::Real weight = 1;
					if ( cpfxn->count( iatm, jatm, weight, path_dist ) ) {
						core::Real fasol1,fasol2, ljatr, ljrep, lkjunk, dis2;

						etable.analytic_etable_evaluation( rsd1.atom(iatm), rsd1.atom(jatm), ljatr, ljrep, lkjunk, dis2);
						fa_atr_i += 0.5*(weight*weightA*ljatr);
						fa_rep_i += 0.5*(weight*weightR*ljrep);

						if ( iatm <= rsd1.nheavyatoms() && jatm <= rsd1.nheavyatoms() ) {
							etable.analytic_lk_energy(rsd1.atom(iatm), rsd1.atom(jatm), fasol1,fasol2 );
							fa_sol_i += weight*weightS*fasol1;
							etable.analytic_lk_energy(a_i, rsd1.atom(jatm), fasol1,fasol2 );
							burial_i += weight*fasol1 / -etable.lk_dgfree( bur_type );
						}
					}
				}
			}
			std::string base_name = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();

			TR << base_name << " " << pose.pdb_info()->number(ires) << " " << pose.size()
				<< " " << pose.residue(ires).name1() << " " << fa_sol_i+lk_ball_i << " " << fa_atr_i+fa_rep_i << " " << burial_i << std::endl;
		}
	}
	std::string get_name() const override {
		return "FaSolReporter";
	}
};

///////
///////
void*
my_main( void* ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( utility::pointer::make_shared< FaSolReporter >() );
	} catch (utility::excn::Exception& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return nullptr;
}



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		NEW_OPT(fasol::del_ddg, "del_ddg value", -0.2);
		NEW_OPT(fasol::max_sasa, "max_sasa", 0.0);
		NEW_OPT(fasol::revert, "revert?", false);
		NEW_OPT(fasol::classic, "use classic parameters", false);

		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
