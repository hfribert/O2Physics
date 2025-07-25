// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file cascadeQa.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for vzeros
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/cascadeHistManager.h"
#include "PWGCF/FemtoUnited/Core/cascadeSelection.h"
#include "PWGCF/FemtoUnited/Core/collisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/collisionSelection.h"
#include "PWGCF/FemtoUnited/Core/modes.h"
#include "PWGCF/FemtoUnited/Core/partitions.h"
#include "PWGCF/FemtoUnited/Core/trackHistManager.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCascadesDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

struct CascadeQa {

  struct : ConfigurableGroup {
    std::string prefix = std::string("Options");

    Configurable<bool> correlatedPlots{"correlatedPlots", false, "Enable multidimensional histogramms. High memory consumption."};
  } Options;

  // setup tables
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Xis = o2::soa::Join<FUXis, FUXiMasks, FUXiExtras>;
  using Omegas = o2::soa::Join<FUOmegas, FUOmegaMasks, FUOmegaExtras>;
  using Tracks = o2::soa::Join<FUTracks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  SliceCache cache;

  // setup collisions
  colhistmanager::CollisionHistManager colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  collisionselection::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  // setup for xis
  cascadeselection::ConfXiSelection confXiSelection;
  Partition<Xis> xiPartition = MAKE_CASCADE_PARTITION(confXiSelection);
  Preslice<Xis> preColXis = aod::femtobase::stored::collisionId;

  cascadehistmanager::ConfXiBinning confXiBinning;
  cascadehistmanager::ConfXiQaBinning confXiQaBinning;
  cascadehistmanager::CascadeHistManager<cascadehistmanager::PrefixXiQa> xiHistManager;

  // setup for omegas
  cascadeselection::ConfOmegaSelection confOmegaSelection;
  Partition<Omegas> omegaPartition = MAKE_CASCADE_PARTITION(confOmegaSelection);
  Preslice<Omegas> preColOmegas = aod::femtobase::stored::collisionId;

  cascadehistmanager::ConfOmegaBinning confOmegaBinning;
  cascadehistmanager::ConfOmegaQaBinning confOmegaQaBinning;
  cascadehistmanager::CascadeHistManager<cascadehistmanager::PrefixOmegaQa> omegaHistManager;

  // setup for daughters/bachelor
  trackhistmanager::ConfCascadePosDauBinning confPosDaughterBinning;
  trackhistmanager::ConfCascadePosDauQaBinning confPosDaughterQaBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDaughterBinning;
  trackhistmanager::ConfCascadeNegDauQaBinning confNegDaughterQaBinning;
  trackhistmanager::ConfCascadeBachelorBinning confBachelorBinning;
  trackhistmanager::ConfCascadeBachelorQaBinning confBachelorQaBinning;

  trackhistmanager::TrackHistManager<trackhistmanager::PrefixCascadePosDaughterQa> posDaughterManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixCascadeNegDaughterQa> negDaughterManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixCascadeBachelorQa> bachelorManager;

  HistogramRegistry hRegistry{"FemtoCascadeQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kANALYSIS>(&hRegistry, colHistSpec);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    posDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, posDaughterHistSpec);

    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);
    negDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, negDaughterHistSpec);

    auto bachelorHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confBachelorBinning, confBachelorQaBinning);
    bachelorManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, bachelorHistSpec);

    if ((doprocessXis + doprocessOmegas) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    if (doprocessXis) {
      auto xiHistSpec = cascadehistmanager::makeCascadeQaHistSpecMap(confXiBinning, confXiQaBinning);
      xiHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, xiHistSpec);
    }

    if (doprocessOmegas) {
      auto omegaHistSpec = cascadehistmanager::makeCascadeQaHistSpecMap(confOmegaBinning, confOmegaQaBinning);
      omegaHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, omegaHistSpec);
    }
  };

  void processXis(FilteredCollision const& col, Xis const& /*xis*/, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto xiSlice = xiPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& xi : xiSlice) {
      xiHistManager.fill<modes::Mode::kANALYSIS_QA, modes::Cascade::kXi>(xi);
      auto posDaugther = xi.posDau_as<Tracks>();
      posDaughterManager.fill<modes::Mode::kANALYSIS_QA>(posDaugther);
      auto negDaugther = xi.negDau_as<Tracks>();
      negDaughterManager.fill<modes::Mode::kANALYSIS_QA>(negDaugther);
      auto bachelor = xi.bachelor_as<Tracks>();
      bachelorManager.fill<modes::Mode::kANALYSIS_QA>(bachelor);
    }
  }
  PROCESS_SWITCH(CascadeQa, processXis, "Process Xis", true);

  void processOmegas(FilteredCollision const& col, Omegas const& /*omegas*/, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto omegaSlice = omegaPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& omega : omegaSlice) {
      omegaHistManager.fill<modes::Mode::kANALYSIS_QA, modes::Cascade::kOmega>(omega);
      auto posDaugther = omega.posDau_as<Tracks>();
      posDaughterManager.fill<modes::Mode::kANALYSIS_QA>(posDaugther);
      auto negDaugther = omega.negDau_as<Tracks>();
      negDaughterManager.fill<modes::Mode::kANALYSIS_QA>(negDaugther);
      auto bachelor = omega.bachelor_as<Tracks>();
      bachelorManager.fill<modes::Mode::kANALYSIS_QA>(bachelor);
    }
  }
  PROCESS_SWITCH(CascadeQa, processOmegas, "Process Omegas", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<CascadeQa>(cfgc),
  };
  return workflow;
}
