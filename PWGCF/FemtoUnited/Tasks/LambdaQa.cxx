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

/// \file LambdaQa.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for vzeros
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/LambdaHistManager.h"
#include "PWGCF/FemtoUnited/Core/LambdaSelection.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"
#include "PWGCF/FemtoUnited/Core/TrackHistManager.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoLambdasDerived.h"
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

struct LambdaQa {

  struct : ConfigurableGroup {
    std::string prefix = std::string("Options");
    Configurable<bool> correlatedPlots{"correlatedPlots", false, "Enable multidimensional histogramms. High memory consumption."};
  } Options;

  collisionselection::ConfCollisionSelection collisionSelection;
  Filter filterVtxz = femtocollisions::posZ >= collisionSelection.vtxZMin && femtocollisions::posZ <= collisionSelection.vtxZMax;
  Filter filterMult = femtocollisions::mult >= collisionSelection.multMin && femtocollisions::mult <= collisionSelection.multMax;
  Filter filterCent = femtocollisions::cent >= collisionSelection.centMin && femtocollisions::cent <= collisionSelection.centMax;
  Filter filterSpher = femtocollisions::sphericity >= collisionSelection.centMin && femtocollisions::sphericity <= collisionSelection.centMax;
  Filter filterMagField = femtocollisions::magField >= collisionSelection.magFieldMin && femtocollisions::magField <= collisionSelection.magFieldMax;

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  colhistmanager::ConfCollisionBinning confCollisionBinning;

  using Lambdas = o2::soa::Join<FULambdas, FULambdaMasks, FULambdaExtras>;
  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  SliceCache cache;

  lambdaselection::ConfLambdaSelection1 confLambdaSelection;

  Partition<Lambdas> LambdaPartition =
    (femtobase::pt > confLambdaSelection.ptMin) &&
    (femtobase::pt < confLambdaSelection.ptMax) &&
    (femtobase::eta > confLambdaSelection.etaMin) &&
    (femtobase::eta < confLambdaSelection.etaMax) &&
    (femtobase::phi > confLambdaSelection.phiMin) &&
    (femtobase::phi < confLambdaSelection.phiMax) &&
    (femtolambdas::lambdaMass > confLambdaSelection.massMin) &&
    (femtolambdas::lambdaMass < confLambdaSelection.massMax) &&
    (femtolambdas::antiLambdaMass > confLambdaSelection.antiMassMin) &&
    (femtolambdas::antiLambdaMass < confLambdaSelection.antiMassMax) &&
    ncheckbit(femtolambdas::lambdaMask, confLambdaSelection.mask);
  Preslice<Lambdas> perColReco = aod::femtobase::collisionId;

  lambdahistmanager::ConfLambdaBinning1 confLambdaBinning;
  lambdahistmanager::ConfLambdaQaBinning1 confLambdaQaBinning;

  trackhistmanager::ConfLambdaPosDauBinning confPosDaughterBinning;
  trackhistmanager::ConfLambdaPosDauQaBinning confPosDaughterQaBinning;
  trackhistmanager::ConfLambdaNegDauBinning confNegDaughterBinning;
  trackhistmanager::ConfLambdaNegDauQaBinning confNegDaughterQaBinning;

  HistogramRegistry hRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  colhistmanager::CollisionHistManager colHistManager;
  lambdahistmanager::LambdaHistManager<lambdahistmanager::PrefixLambdaQa> lambdaHistManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixLambdaPosDaughterQa> PosDaughterManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixLambdaNegDaughterQa> NegDaughterManager;

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kANALYSIS>(&hRegistry, colHistSpec);

    auto lambdaHistSpec = lambdahistmanager::makeLambdaQaHistSpecMap(confLambdaBinning, confLambdaQaBinning);
    lambdaHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, lambdaHistSpec);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    PosDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, posDaughterHistSpec);

    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);
    NegDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, negDaughterHistSpec);
  };

  void process(FilteredCollision const& col, Lambdas const& /*lambdas*/, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto lambdaSlice = LambdaPartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      lambdaHistManager.fill<modes::Mode::kANALYSIS_QA>(lambda);
      auto posDaugther = lambda.posDauLambda_as<Tracks>();
      PosDaughterManager.fill<modes::Mode::kANALYSIS>(posDaugther);
      auto negDaugther = lambda.negDauLambda_as<Tracks>();
      NegDaughterManager.fill<modes::Mode::kANALYSIS>(negDaugther);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<LambdaQa>(cfgc),
  };
  return workflow;
}
