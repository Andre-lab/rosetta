// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/kcluster.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>
#include <map>
#include <string>
#include <utility>

// External headers
#include <boost/format.hpp>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/cmiles.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

#define ITERATIONS_PER_ASSIGNMENT 10

using core::pose::Pose;
using std::cout;
using std::endl;
using std::pair;
using utility::vector1;

class PairComparator {
 public:
  bool operator() (const pair<int, int>& a, const pair<int, int>& b) const {
    if (a.first < b.first) {
      return true;
    } else if (a.first > b.first) {
      return false;
    } else {  // a.first == b.first
      return a.second < b.second;
    }
  }
};

/// @detail Computes the distance between two structures using gdtmm
double dist(const vector1<Pose>& models, int i, int j) {
  using std::map;

  static map<pair<int, int>, double, PairComparator> cache;

  pair<int, int> key(i, j);
  pair<int, int> reverse_key(j, i);

  double d;
  if (cache.find(key) == cache.end()) {
    d = (i == j) ? 0 : 1 - core::scoring::CA_gdtmm(models[i], models[j]);
    cache[key] = d;
    cache[reverse_key] = d;
  }

  return cache[key];
}

/// @detail Selects k cluster centers using a simple rule: the next cluster
/// center to be added is as far as possible from the centers chosen so far
void partition(const int k, const vector1<Pose>& models, vector1<int>* centers) {
  const int num_models = models.size();

  centers->push_back(numeric::random::random_range(1, num_models));

  for (int c = 2; c <= k; ++c) {
    int max_index = 0;
    double max_dist = std::numeric_limits<double>::min();

    for (int i = 1; i <= num_models; ++i) {
      // distance to point i's closest cluster center
      double dist_nearest = std::numeric_limits<double>::max();

      for (vector1<int>::const_iterator j = centers->begin(); j != centers->end(); ++j) {
        double distance = dist(models, i, *j);
        if (distance < dist_nearest) {
          dist_nearest = distance;
        }
      }
      //cout << "Distance to nearest cluster center for point " << i << " is: " << dist_nearest << endl;

      if (dist_nearest > max_dist) {
        max_dist = dist_nearest;
        max_index = i;
      }
    }

    //cout << "Added cluster center " << c << ": index = " << max_index << " dist = " << max_dist << endl;
    centers->push_back(max_index);
  }
}

/// @detail Assigns each model to the closest cluster center
void assign(const vector1<Pose>& models, const vector1<int>& centers) {
  using core::io::silent::SilentFileData;
  using core::io::silent::SilentStructFactory;
  using core::io::silent::SilentStructOP;
  using std::string;

  const int num_models = models.size();
  const int num_clusters = centers.size();

  vector1<string> filenames;
  for (int i = 1; i <= num_clusters; ++i) {
    filenames.push_back(str(boost::format("c.%d.out") % i));
  }

  for (int i = 1; i <= num_models; ++i) {
    double min_dist = std::numeric_limits<double>::max();
    int min_index = 0;

    for (int j = 1; j <= num_clusters; ++j) {
      double distance = dist(models, i, centers[j]);
      if (distance < min_dist) {
        min_dist = distance;
        min_index = j;
      }
    }

    SilentStructOP silent = SilentStructFactory::get_instance()->get_silent_struct_out();
    silent->fill_struct(models[i]);

    SilentFileData sfd;
    sfd.write_silent_struct(*silent, filenames[min_index]);
  }
}

int main(int argc, char* argv[]) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::chemical::ChemicalManager;
  using core::chemical::ResidueTypeSetCAP;
  using core::import_pose::pose_stream::MetaPoseInputStream;
  devel::init(argc, argv);

  ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set("fa_standard");
  MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();

  vector1<Pose> models;

  Pose pose;
  while (input.has_another_pose()) {
    input.fill_pose(pose, *rsd_set);
    models.push_back(pose);
  }
  cout << "Read " << models.size() << " models" << endl;

  int k = option[OptionKeys::cmiles::kcluster::num_clusters]();
  vector1<int> centers;
  partition(k, models, &centers);
  assign(models, centers);
  cout << "Assigned " << models.size() << " models to " << k << " clusters" << endl;
}
