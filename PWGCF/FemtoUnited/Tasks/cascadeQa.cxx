// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
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

#include "PWGCF/FemtoUnited/Core/CascadeHistManager.h"
#include "PWGCF/FemtoUnited/Core/CascadeSelection.h"
#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"
#include "PWGCF/FemtoUnited/Core/TrackHistManager.h"
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

  collisionselection::ConfCollisionSelection collisionSelection;
  Filter filterVtxz = femtocollisions::posZ >= collisionSelection.vtxZMin && femtocollisions::posZ <= collisionSelection.vtxZMax;
  Filter filterMult = femtocollisions::mult >= collisionSelection.multMin && femtocollisions::mult <= collisionSelection.multMax;
  Filter filterCent = femtocollisions::cent >= collisionSelection.centMin && femtocollisions::cent <= collisionSelection.centMax;
  Filter filterSpher = femtocollisions::sphericity >= collisionSelection.spherMin && femtocollisions::sphericity <= collisionSelection.spherMax;
  Filter filterMagField = femtocollisions::magField >= collisionSelection.magFieldMin && femtocollisions::magField <= collisionSelection.magFieldMax;

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  colhistmanager::ConfCollisionBinning confCollisionBinning;

  using Cascades = o2::soa::Join<FUCascades, FUCascadeMasks, FUCascadeExtras>;
  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  SliceCache cache;

  lambdaselection::ConfCascadeSelection1 confCascadeSelection;

  Partition<Cascades> CascadePartition =
    (femtobase::pt > confCascadeSelection.ptMin) &&
    (femtobase::pt < confCascadeSelection.ptMax) &&
    (femtobase::eta > confCascadeSelection.etaMin) &&
    (femtobase::eta < confCascadeSelection.etaMax) &&
    (femtobase::phi > confCascadeSelection.phiMin) &&
    (femtobase::phi < confCascadeSelection.phiMax) &&
    (femtolambdas::lambdaMass > confCascadeSelection.massMin) &&
    (femtolambdas::lambdaMass < confCascadeSelection.massMax) &&
    (femtolambdas::antiCascadeMass > confCascadeSelection.antiMassMin) &&
    (femtolambdas::antiCascadeMass < confCascadeSelection.antiMassMax) &&
    ncheckbit(femtolambdas::lambdaMask, confCascadeSelection.mask);
  Preslice<Cascades> perColReco = aod::femtobase::collisionId;

  lambdahistmanager::ConfCascadeBinning1 confCascadeBinning;
  lambdahistmanager::ConfCascadeQaBinning1 confCascadeQaBinning;

  trackhistmanager::ConfCascadePosDauBinning confPosDaughterBinning;
  trackhistmanager::ConfCascadePosDauQaBinning confPosDaughterQaBinning;
  trackhistmanager::ConfCascadeNegDauBinning confNegDaughterBinning;
  trackhistmanager::ConfCascadeNegDauQaBinning confNegDaughterQaBinning;

  HistogramRegistry hRegistry{"CascadeQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  colhistmanager::CollisionHistManager colHistManager;
  lambdahistmanager::CascadeHistManager<lambdahistmanager::PrefixCascadeQa> lambdaHistManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixCascadePosDaughterQa> posDaughterManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixCascadeNegDaughterQa> negDaughterManager;

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kANALYSIS>(&hRegistry, colHistSpec);

    auto lambdaHistSpec = lambdahistmanager::makeCascadeQaHistSpecMap(confCascadeBinning, confCascadeQaBinning);
    lambdaHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, lambdaHistSpec);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    posDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, posDaughterHistSpec);

    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);
    negDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, negDaughterHistSpec);
  };

  void process(FilteredCollision const& col, Cascades const& /*lambdas*/, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto lambdaSlice = CascadePartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      lambdaHistManager.fill<modes::Mode::kANALYSIS_QA>(lambda);
      auto posDaugther = lambda.posDauCascade_as<Tracks>();
      posDaughterManager.fill<modes::Mode::kANALYSIS_QA>(posDaugther);
      auto negDaugther = lambda.negDauCascade_as<Tracks>();
      negDaughterManager.fill<modes::Mode::kANALYSIS_QA>(negDaugther);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<CascadeQa>(cfgc),
  };
  return workflow;
}
