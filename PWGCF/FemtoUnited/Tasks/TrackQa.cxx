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

/// \file TrackQa.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"
#include "PWGCF/FemtoUnited/Core/TrackHistManager.h"
#include "PWGCF/FemtoUnited/Core/TrackSelection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
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

struct TrackQa {

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

  SliceCache cache;

  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  trackselection::ConfTrackSelection1 trackSelections;

  Partition<Tracks> TrackPartition =
    (femtobase::pt > trackSelections.ptMin) &&
    (femtobase::pt < trackSelections.ptMax) &&
    (femtobase::eta > trackSelections.etaMin) &&
    (femtobase::eta < trackSelections.etaMax) &&
    (femtobase::phi > trackSelections.phiMin) &&
    (femtobase::phi < trackSelections.phiMax) &&
    ifnode(femtobase::pt * (nexp(femtobase::eta) + nexp(-1.f * femtobase::eta)) / 2.f <= trackSelections.pidThres,
           ncheckbit(femtotracks::trackMask, trackSelections.maskLowMomentum), ncheckbit(femtotracks::trackMask, trackSelections.maskHighMomentum));

  Preslice<Tracks> perColReco = aod::femtobase::collisionId;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackBinning");
    ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
    ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
    ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
    ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
    ConfigurableAxis itsCluster{"itsCluster", {{8, -0.5, 7.5}}, "ITS cluster"};
    ConfigurableAxis itsClusterIb{"itsClusterIb", {{4, -0.5, 3.5}}, "ITS cluster in inner barrel"};
    ConfigurableAxis tpcCluster{"tpcCluster", {{153, -0.5, 152.5}}, "TPC cluster"};
    ConfigurableAxis tpcClusterShared{"tpcClusterShared", {{153, -0.5, 152.5}}, "TPC cluster shared"};
    ConfigurableAxis dcaXy{"dcaXy", {{300, -0.3, 0.3}}, "DCA_xy"};
    ConfigurableAxis dcaZ{"dcaZ", {{300, -0.3, 0.3}}, "DCA_Z"};
    ConfigurableAxis dca{"dca", {{300, 0, 0.3}}, "DCA"};
  } TrackBinning;

  struct : ConfigurableGroup {
    std::string prefix = std::string("TrackPidBinning");
    ConfigurableAxis p{"p", {{300, 0, 6}}, "Momentum axis"};
    ConfigurableAxis itsSignal{"itsSignal", {{150, 0, 15}}, "its Signal"};
    ConfigurableAxis itsElectron{"itsElectron", {{300, -3, 3}}, "ITS PID for electron"};
    ConfigurableAxis itsPion{"itsPion", {{300, -3, 3}}, "ITS PID for pion"};
    ConfigurableAxis itsKaon{"itsKaon", {{300, -3, 3}}, "ITS PID for kaon"};
    ConfigurableAxis itsProton{"itsProton", {{300, -3, 3}}, "ITS PID for proton"};
    ConfigurableAxis itsDeuteron{"itsDeuteron", {{300, -3, 3}}, "ITS PID for deuteron"};
    ConfigurableAxis itsTriton{"itsTriton", {{300, -3, 3}}, "ITS PID for triton"};
    ConfigurableAxis itsHelium{"itsHelium", {{300, -3, 3}}, "ITS PID for helium"};
    ConfigurableAxis tpcSignal{"tpcSignal", {{150, 0, 150}}, "TPC Signal"};
    ConfigurableAxis tpcElectron{"tpcElectron", {{300, -3, 3}}, "TPC PID for electron"};
    ConfigurableAxis tpcPion{"tpcPion", {{300, -3, 3}}, "TPC PID for pion"};
    ConfigurableAxis tpcKaon{"tpcKaon", {{300, -3, 3}}, "TPC PID for kaon"};
    ConfigurableAxis tpcProton{"tpcProton", {{300, -3, 3}}, "TPC PID for proton"};
    ConfigurableAxis tpcDeuteron{"tpcDeuteron", {{300, -3, 3}}, "TPC PID for deuteron"};
    ConfigurableAxis tpcTriton{"tpcTriton", {{300, -3, 3}}, "TPC PID for triton"};
    ConfigurableAxis tpcHelium{"tpcHelium", {{300, -3, 3}}, "TPC PID for helium"};
    ConfigurableAxis tofBeta{"tofBeta", {{150, 0, 1.5}}, "TPC Signal"};
    ConfigurableAxis tofElectron{"tofElectron", {{300, -3, 3}}, "TOF PID for electron"};
    ConfigurableAxis tofPion{"tofPion", {{300, -3, 3}}, "TOF PID for pion"};
    ConfigurableAxis tofKaon{"tofKaon", {{300, -3, 3}}, "TOF PID for kaon"};
    ConfigurableAxis tofProton{"tofProton", {{300, -3, 3}}, "TOF PID for proton"};
    ConfigurableAxis tofDeuteron{"tofDeuteron", {{300, -3, 3}}, "TOF PID for deuteron"};
    ConfigurableAxis tofTriton{"tofTriton", {{300, -3, 3}}, "TOF PID for triton"};
    ConfigurableAxis tofHelium{"tofHelium", {{300, -3, 3}}, "TOF PID for helium"};
    ConfigurableAxis tpctofElectron{"tpctofElectron", {{300, 0, 3}}, "tpctof PID for electron"};
    ConfigurableAxis tpctofPion{"tpctofPion", {{300, 0, 3}}, "tpctof PID for pion"};
    ConfigurableAxis tpctofKaon{"tpctofKaon", {{300, 0, 3}}, "tpctof PID for kaon"};
    ConfigurableAxis tpctofProton{"tpctofProton", {{300, 0, 3}}, "tpctof PID for proton"};
    ConfigurableAxis tpctofDeuteron{"tpctofDeuteron", {{300, 0, 3}}, "tpctof PID for deuteron"};
    ConfigurableAxis tpctofTriton{"tpctofTriton", {{300, 0, 3}}, "tpctof PID for triton"};
    ConfigurableAxis tpctofHelium{"tpctofHelium", {{300, 0, 3}}, "tpctof PID for helium"};
  } TrackPidBinning;

  HistogramRegistry hRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  colhistmanager::CollisionHistManager colHistManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrackQa> trackHistManager;

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

    std::map<trackhistmanager::TrackHist, std::vector<framework::AxisSpec>> trackHistSpec = {
      {trackhistmanager::kPt, {TrackBinning.pt}},
      {trackhistmanager::kEta, {TrackBinning.eta}},
      {trackhistmanager::kPhi, {TrackBinning.phi}},
      {trackhistmanager::kSign, {TrackBinning.sign}},
      {trackhistmanager::kItsCluster, {TrackBinning.itsCluster}},
      {trackhistmanager::kItsClusterIb, {TrackBinning.itsClusterIb}},
      {trackhistmanager::kPtVsEta, {TrackBinning.pt, TrackBinning.eta}},
      {trackhistmanager::kPtVsPhi, {TrackBinning.pt, TrackBinning.phi}},
      {trackhistmanager::kPhiVsEta, {TrackBinning.phi, TrackBinning.eta}},
      {trackhistmanager::kPtVsItsCluster, {TrackBinning.pt, TrackBinning.itsCluster}},
      {trackhistmanager::kPtVsTpcCluster, {TrackBinning.pt, TrackBinning.tpcCluster}},
      {trackhistmanager::kPtVsTpcClusterShared, {TrackBinning.pt, TrackBinning.tpcClusterShared}},
      {trackhistmanager::kTpcClusterVsTpcClusterShared, {TrackBinning.tpcCluster, TrackBinning.tpcClusterShared}},
      {trackhistmanager::kTpcCluster, {TrackBinning.tpcCluster}},
      {trackhistmanager::kTpcClusterShared, {TrackBinning.tpcClusterShared}},
      {trackhistmanager::kPtVsDcaxy, {TrackBinning.pt, TrackBinning.dcaXy}},
      {trackhistmanager::kPtVsDcaz, {TrackBinning.pt, TrackBinning.dcaZ}},
      {trackhistmanager::kPtVsDca, {TrackBinning.pt, TrackBinning.dca}},
      {trackhistmanager::kItsSignal, {TrackPidBinning.p, TrackPidBinning.itsSignal}},
      {trackhistmanager::kItsElectron, {TrackPidBinning.p, TrackPidBinning.itsElectron}},
      {trackhistmanager::kItsPion, {TrackPidBinning.p, TrackPidBinning.itsPion}},
      {trackhistmanager::kItsKaon, {TrackPidBinning.p, TrackPidBinning.itsKaon}},
      {trackhistmanager::kItsProton, {TrackPidBinning.p, TrackPidBinning.itsProton}},
      {trackhistmanager::kItsDeuteron, {TrackPidBinning.p, TrackPidBinning.itsDeuteron}},
      {trackhistmanager::kItsTriton, {TrackPidBinning.p, TrackPidBinning.itsTriton}},
      {trackhistmanager::kItsHelium, {TrackPidBinning.p, TrackPidBinning.itsHelium}},
      {trackhistmanager::kTpcSignal, {TrackPidBinning.p, TrackPidBinning.tpcSignal}},
      {trackhistmanager::kTpcElectron, {TrackPidBinning.p, TrackPidBinning.tpcElectron}},
      {trackhistmanager::kTpcPion, {TrackPidBinning.p, TrackPidBinning.tpcPion}},
      {trackhistmanager::kTpcKaon, {TrackPidBinning.p, TrackPidBinning.tpcKaon}},
      {trackhistmanager::kTpcProton, {TrackPidBinning.p, TrackPidBinning.tpcProton}},
      {trackhistmanager::kTpcDeuteron, {TrackPidBinning.p, TrackPidBinning.tpcDeuteron}},
      {trackhistmanager::kTpcTriton, {TrackPidBinning.p, TrackPidBinning.tpcTriton}},
      {trackhistmanager::kTpcHelium, {TrackPidBinning.p, TrackPidBinning.tpcHelium}},
      {trackhistmanager::kTofBeta, {TrackPidBinning.p, TrackPidBinning.tofBeta}},
      {trackhistmanager::kTofElectron, {TrackPidBinning.p, TrackPidBinning.tofElectron}},
      {trackhistmanager::kTofPion, {TrackPidBinning.p, TrackPidBinning.tofPion}},
      {trackhistmanager::kTofKaon, {TrackPidBinning.p, TrackPidBinning.tofKaon}},
      {trackhistmanager::kTofProton, {TrackPidBinning.p, TrackPidBinning.tofProton}},
      {trackhistmanager::kTofDeuteron, {TrackPidBinning.p, TrackPidBinning.tofDeuteron}},
      {trackhistmanager::kTofTriton, {TrackPidBinning.p, TrackPidBinning.tofTriton}},
      {trackhistmanager::kTofHelium, {TrackPidBinning.p, TrackPidBinning.tofHelium}},
      {trackhistmanager::kTpctofElectron, {TrackPidBinning.p, TrackPidBinning.tpctofElectron}},
      {trackhistmanager::kTpctofPion, {TrackPidBinning.p, TrackPidBinning.tpctofPion}},
      {trackhistmanager::kTpctofKaon, {TrackPidBinning.p, TrackPidBinning.tpctofKaon}},
      {trackhistmanager::kTpctofProton, {TrackPidBinning.p, TrackPidBinning.tpctofProton}},
      {trackhistmanager::kTpctofDeuteron, {TrackPidBinning.p, TrackPidBinning.tpctofDeuteron}},
      {trackhistmanager::kTpctofTriton, {TrackPidBinning.p, TrackPidBinning.tpctofTriton}},
      {trackhistmanager::kTpctofHelium, {TrackPidBinning.p, TrackPidBinning.tpctofHelium}}};

    trackHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, trackHistSpec);
  };

  void process(FilteredCollision const& col, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto trackSlice = TrackPartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    for (auto const& track : trackSlice) {
      trackHistManager.fill<modes::Mode::kANALYSIS_QA>(track);
    }
  }
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TrackQa>(cfgc),
  };
  return workflow;
}
