// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/*
 * DetailedScores.cc
 *
 *  Created on: Jan 29, 2009
 *      Author: dgront
 */

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <protocols/Protocol.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <utility/vector1.hh>

static basic::Tracer tr("SaxsSampler");

void register_options() {

  OPT( in::file::native );
}

int main(int argc, char * argv[]) {

  using namespace core;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::scoring;
  using std::string;
  using utility::vector1;

  protocols::abinitio::AbrelaxApplication::register_options();
  register_options();
  devel::init(argc, argv);

  core::scoring::ScoreFunctionOP sfstd(core::scoring::ScoreFunctionFactory::create_score_function(core::scoring::STANDARD_WTS));

  pose::Pose init_pose;
  if (option[in::file::native].user()) {
    core::import_pose::pose_from_pdb(init_pose, option[in::file::native]());
  }

  (*sfstd)(init_pose);
  sfstd->show(std::cout, init_pose);

  Energies & energies(init_pose.energies());

  // the neighbor/energy links
  core::scoring::EnergyGraph & energy_graph(energies.energy_graph());

  // should be zeroed at the beginning of scoring
  EnergyMap totals; //( energies.totals() );

  std::cout << totals[fa_atr] << std::endl;

  for (Size i = 1, i_end = init_pose.total_residue(); i <= i_end; ++i) {
    conformation::Residue const & resl(init_pose.residue(i));

    for (graph::Graph::EdgeListIter iru = energy_graph.get_node(i)->upper_edge_list_begin(), irue =
        energy_graph.get_node(i)->upper_edge_list_end(); iru != irue; ++iru) {
      EnergyEdge * edge(static_cast<EnergyEdge *> (*iru));
      Size const j(edge->get_second_node_ind());

      // now i have i and j
      conformation::Residue const & resu(init_pose.residue(j));
      //      resl is i
      //      and resu is j

      // the pair energies cached in the link
      EnergyMap & emap(edge->energy_map());

      // the context-dependent guys cant be cached
      emap.zero(sfstd->cd_2b_types());
      //zero_energies( methods::cd_2b, emap );
      sfstd->eval_cd_2b(resl, resu, init_pose, emap);

      std::cout << emap[fa_atr] << std::endl;
      // accumulate energies

      EnergyMap temp;
      temp.accumulate(emap, sfstd->ci_2b_types());
      temp.accumulate(emap, sfstd->cd_2b_types());

      std::cout << temp[fa_atr] << std::endl;

    } // nbrs of i
  } // i=1,nres
}
