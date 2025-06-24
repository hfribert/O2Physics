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
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/LambdaHistManager.h"
#include "PWGCF/FemtoUnited/Core/LambdaSelection.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoLambdasDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "CommonConstants/MathConstants.h"
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

  colhistmanager::ConfCollisionBinning collisionBinning;

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

  struct : ConfigurableGroup {
    std::string prefix = std::string("LambdaBinning");
    ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
    ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
    ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
    ConfigurableAxis mass{"mass", {{200, 1.0, 1.2}}, "Mass"};
    ConfigurableAxis dauPt{"dauPt", {{600, 0, 6}}, "daughter pt"};
    ConfigurableAxis dauEta{"dauEta", {{300, -1.5, 1.5}}, "daugher eta"};
    ConfigurableAxis dauPhi{"dauPhi", {{720, 0., 1.f * o2::constants::math::TwoPI}}, "daughter phi"};
    ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
    ConfigurableAxis dauDcaAtDecay{"dauDcaAtDecay", {{150, 0, 1.5}}, "Daughter DCA at decay vertex"};
    ConfigurableAxis decayVertex{"decayVertex", {{100, 0, 100}}, "Decay vertex"};
    ConfigurableAxis transRadius{"transRadius", {{100, 0, 100}}, "Transverse radius"};
    ConfigurableAxis kaonMass{"kaonMass", {{100, 0.45, 0.55}}, "Mass for kaon hypothesis"};
    ConfigurableAxis dauTpcCluster{"dauTpcCluster", {{153, -0.5, 152.5}}, "TPC cluster of daughters"};
    ConfigurableAxis dauP{"dauP", {{600, 0, 6}}, "Momentum binning for TPC Nsigma of daughters"};
    ConfigurableAxis dauDcaxy{"dauDcaxy", {{300, -1.5, 1.5}}, "DCAxy for daughters"};
    ConfigurableAxis dauDcaz{"dauDcaz", {{300, -1.5, 1.5}}, "Dcaz for daughters"};
    ConfigurableAxis dauDca{"dauDca", {{150, 0, 0.3}}, "Dca for daughters"};
    ConfigurableAxis dauTpcNsigma{"dauTpcNsigma", {{600, -6, 6}}, "TPC Nsigma for daughters"};
  } LambdaBinning;

  HistogramRegistry hRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  colhistmanager::CollisionHistManager colHistManager;
  lambdahistmanager::LambdaHistManager lambdaHistManager;

  void init(InitContext&)
  {
    // create a map for histogram specs
    std::map<colhistmanager::ColHist, std::vector<framework::AxisSpec>> colHistSpec = {
      {colhistmanager::kPosz, {collisionBinning.vtZ}},
      {colhistmanager::kMult, {collisionBinning.mult}},
      {colhistmanager::kCent, {collisionBinning.cent}},
      {colhistmanager::kSphericity, {collisionBinning.spher}},
      {colhistmanager::kMagField, {collisionBinning.magField}},
      {colhistmanager::kPoszVsMult, {collisionBinning.vtZ, collisionBinning.mult}},
      {colhistmanager::kPoszVsCent, {collisionBinning.vtZ, collisionBinning.cent}},
      {colhistmanager::kCentVsMult, {collisionBinning.cent, collisionBinning.mult}},
      {colhistmanager::kMultVsSphericity, {collisionBinning.mult, collisionBinning.spher}},
      {colhistmanager::kCentVsSphericity, {collisionBinning.cent, collisionBinning.spher}}};
    colHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, colHistSpec);

    std::map<lambdahistmanager::LambdaHist, std::vector<framework::AxisSpec>> lambdaHistSpec = {
      {lambdahistmanager::kPt, {LambdaBinning.pt}},
      {lambdahistmanager::kEta, {LambdaBinning.eta}},
      {lambdahistmanager::kPhi, {LambdaBinning.phi}},
      {lambdahistmanager::kMass, {LambdaBinning.mass}},
      {lambdahistmanager::kAntiMass, {LambdaBinning.mass}},
      {lambdahistmanager::kPosDauPt, {LambdaBinning.dauPt}},
      {lambdahistmanager::kPosDauEta, {LambdaBinning.dauEta}},
      {lambdahistmanager::kPosDauPhi, {LambdaBinning.dauPhi}},
      {lambdahistmanager::kNegDauPt, {LambdaBinning.dauPt}},
      {lambdahistmanager::kNegDauEta, {LambdaBinning.dauEta}},
      {lambdahistmanager::kNegDauPhi, {LambdaBinning.dauPhi}},
      {lambdahistmanager::kDecayDauDca, {LambdaBinning.dauDca}},
      {lambdahistmanager::kDecayVtxX, {LambdaBinning.decayVertex}},
      {lambdahistmanager::kDecayVtxY, {LambdaBinning.decayVertex}},
      {lambdahistmanager::kDecayVtxZ, {LambdaBinning.decayVertex}},
      {lambdahistmanager::kDecayVtx, {LambdaBinning.decayVertex}},
      {lambdahistmanager::kTransRadius, {LambdaBinning.transRadius}},
      {lambdahistmanager::kKaonMass, {LambdaBinning.kaonMass}},
      {lambdahistmanager::kPtVsEta, {LambdaBinning.pt, LambdaBinning.eta}},
      {lambdahistmanager::kPtVsPhi, {LambdaBinning.pt, LambdaBinning.phi}},
      {lambdahistmanager::kPhiVsEta, {LambdaBinning.phi, LambdaBinning.eta}},
      {lambdahistmanager::kPosDauTpcCluster, {LambdaBinning.dauTpcCluster}},
      {lambdahistmanager::kPosDauPtVsDcaxy, {LambdaBinning.dauPt, LambdaBinning.dauDcaxy}},
      {lambdahistmanager::kPosDauPtVsDcaz, {LambdaBinning.dauPt, LambdaBinning.dauDcaz}},
      {lambdahistmanager::kPosDauPtVsDca, {LambdaBinning.dauPt, LambdaBinning.dauDca}},
      {lambdahistmanager::kPosDauProtonTpcNsigma, {LambdaBinning.dauP, LambdaBinning.dauTpcNsigma}},
      {lambdahistmanager::kPosDauPionTpcNsigma, {LambdaBinning.dauP, LambdaBinning.dauTpcNsigma}},
      {lambdahistmanager::kNegDauTpcCluster, {LambdaBinning.dauTpcCluster}},
      {lambdahistmanager::kNegDauPtVsDcaxy, {LambdaBinning.dauPt, LambdaBinning.dauDcaxy}},
      {lambdahistmanager::kNegDauPtVsDcaz, {LambdaBinning.dauPt, LambdaBinning.dauDcaz}},
      {lambdahistmanager::kNegDauPtVsDca, {LambdaBinning.dauPt, LambdaBinning.dauDca}},
      {lambdahistmanager::kNegDauProtonTpcNsigma, {LambdaBinning.dauP, LambdaBinning.dauTpcNsigma}},
      {lambdahistmanager::kNegDauPionTpcNsigma, {LambdaBinning.dauP, LambdaBinning.dauTpcNsigma}},
    };

    lambdaHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, lambdaHistSpec);
  };

  void process(FilteredCollision const& col, Lambdas const& /*lambdas*/, Tracks const& tracks)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto lambdaSlice = LambdaPartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    for (auto const& lambda : lambdaSlice) {
      lambdaHistManager.fill<modes::Mode::kANALYSIS_QA>(lambda, tracks);
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
