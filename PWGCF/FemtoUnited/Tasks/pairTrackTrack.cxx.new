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

/// \file PairTrackTrack.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/ClosePairRejection.h"
#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"
#include "PWGCF/FemtoUnited/Core/PairCleaner.h"
#include "PWGCF/FemtoUnited/Core/PairHistManager.h"
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
#include <random>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

struct PairTrackTrack {

  struct : ConfigurableGroup {
    std::string prefix = std::string("Options");
    Configurable<bool> correlatedPlots{"correlatedPlots", false, "Enable multidimensional histogramms. High memory consumption."};
    Configurable<bool> sameSpecies{"sameSpecies", false, "Enable if track1 and track2 are the same particle"};
    Configurable<int> randomizePairSeed{"randomizePairSeed", -1, "Seed to randomize track1 and track2. Set to negative value to deactivate. Set to 0 to generate unique seed in time."};
  } Options;

  std::mt19937 rng;

  collisionselection::ConfCollisionSelection collisionSelection;
  Filter filterVtxz = femtocollisions::posZ >= collisionSelection.vtxZMin && femtocollisions::posZ <= collisionSelection.vtxZMax;
  Filter filterMult = femtocollisions::mult >= collisionSelection.multMin && femtocollisions::mult <= collisionSelection.multMax;
  Filter filterCent = femtocollisions::cent >= collisionSelection.centMin && femtocollisions::cent <= collisionSelection.centMax;
  Filter filterSpher = femtocollisions::sphericity >= collisionSelection.centMin && femtocollisions::sphericity <= collisionSelection.centMax;
  Filter filterMagField = femtocollisions::magField >= collisionSelection.magFieldMin && femtocollisions::magField <= collisionSelection.magFieldMax;

  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // Mixing configurables
  struct : ConfigurableGroup {
    std::string prefix = "Mixing";
    ConfigurableAxis multMixBins{"multMixBins", {VARIABLE_WIDTH, 0.0f, 4.0f, 8.0f, 12.0f, 16.0f, 20.0f, 24.0f, 28.0f, 32.0f, 36.0f, 40.0f, 44.0f, 48.0f, 52.0f, 56.0f, 60.0f, 64.0f, 68.0f, 72.0f, 76.0f, 80.0f, 84.0f, 88.0f, 92.0f, 96.0f, 100.0f, 200.0f}, "Mixing bins - multiplicity"};
    ConfigurableAxis multPercentileMixBins{"multPercentileMixBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - multiplicity percentile"};
    ConfigurableAxis vztxMixBins{"vztxMixBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    Configurable<int> depth{"depth", 5, "Number of events for mixing"};
    Configurable<int> policy{"policy", 0, "Binning policy for mixing (alywas in combination with z-vertex) -> 0: multiplicity, -> 1: multipliciy percentile, -> 2: both"};
  } Mixing;
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult> colBinningMult{{Mixing.vztxMixBins, Mixing.multMixBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> colBinningMultPercentile{{Mixing.vztxMixBins, Mixing.multPercentileMixBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> colBinningMultMultPercentile{{Mixing.vztxMixBins, Mixing.multMixBins, Mixing.multPercentileMixBins}, true};

  SliceCache cache;

  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  trackselection::ConfTrackSelection1 trackSelections1;
  Partition<Tracks> TrackPartition1 =
    ifnode(trackSelections1.sign.node() > 0, femtotracks::signedPt > 0.f, femtotracks::signedPt < 0.f) &&
    (nabs(femtotracks::signedPt) > trackSelections1.ptMin) &&
    (nabs(femtotracks::signedPt) < trackSelections1.ptMax) &&
    (femtobase::eta > trackSelections1.etaMin) &&
    (femtobase::eta < trackSelections1.etaMax) &&
    (femtobase::phi > trackSelections1.phiMin) &&
    (femtobase::phi < trackSelections1.phiMax) &&
    ifnode(nabs(femtotracks::signedPt) * (nexp(femtobase::eta) + nexp(-1.f * femtobase::eta)) / 2.f <= trackSelections1.pidThres,
           ncheckbit(femtotracks::trackMask, trackSelections1.maskLowMomentum), ncheckbit(femtotracks::trackMask, trackSelections1.maskHighMomentum));

  trackselection::ConfTrackSelection2 trackSelections2;
  Partition<Tracks> TrackPartition2 =
    ifnode(trackSelections2.sign.node() > 0, femtotracks::signedPt > 0.f, femtotracks::signedPt < 0.f) &&
    (nabs(femtotracks::signedPt) > trackSelections2.ptMin) &&
    (nabs(femtotracks::signedPt) < trackSelections2.ptMax) &&
    (femtobase::eta > trackSelections2.etaMin) &&
    (femtobase::eta < trackSelections2.etaMax) &&
    (femtobase::phi > trackSelections2.phiMin) &&
    (femtobase::phi < trackSelections2.phiMax) &&
    ifnode(nabs(femtotracks::signedPt) * (nexp(femtobase::eta) + nexp(-1.f * femtobase::eta)) / 2.f <= trackSelections2.pidThres,
           ncheckbit(femtotracks::trackMask, trackSelections2.maskLowMomentum), ncheckbit(femtotracks::trackMask, trackSelections2.maskHighMomentum));

  Preslice<Tracks> perColReco = aod::femtobase::collisionId;

  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;

  struct : ConfigurableGroup {
    std::string prefix = std::string("PairBinning");
    ConfigurableAxis kstar{"kstar", {{600, 0, 6}}, "kstar"};
    ConfigurableAxis kt{"kt", {{600, 0, 6}}, "kt"};
    ConfigurableAxis mt{"mt", {{200, 0.5, 2.5}}, "mt"};
  } PairBinning;

  struct : ConfigurableGroup {
    std::string prefix = std::string("ClosePairRejection");
    Configurable<bool> on{"on", true, "Trun on CPR"};
    Configurable<float> detaMax{"detaMax", 0.01f, "Maximium deta"};
    Configurable<float> dphistarMax{"dphistarMax", 0.01f, "Maximum dphistar"};
    ConfigurableAxis deta{"deta", {{100, -0.2, 0.2}}, "deta"};
    ConfigurableAxis dphistar{"dphistar", {{100, -0.2, 0.2}}, "dphi"};
  } ConfCpr;

  HistogramRegistry hRegistry{"TrackQA", {}, OutputObjHandlingPolicy::AnalysisObject};
  colhistmanager::CollisionHistManager colHistManager;

  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack1> trackHistManager1;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack2> trackHistManager2;

  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackSe> pairHistManagerSe;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackMe> pairHistManagerMe;

  closepairrejection::ClosePairRejection<closepairrejection::PrefixTrackTrackSe> cprSe;
  closepairrejection::ClosePairRejection<closepairrejection::PrefixTrackTrackMe> cprMe;
  paircleaner::PairCleaner pc;

  void init(InitContext&)
  {

    // setup columnpolicy for binning
    colBinningMult = {{Mixing.vztxMixBins, Mixing.multMixBins}, true};
    colBinningMultPercentile = {{Mixing.vztxMixBins, Mixing.multPercentileMixBins}, true};
    colBinningMultMultPercentile = {{Mixing.vztxMixBins, Mixing.multMixBins, Mixing.multPercentileMixBins}, true};

    // setup rng if necessary
    if (Options.randomizePairSeed.value > 0) {
      uint64_t randomSeed = 0;
      if (Options.randomizePairSeed.value == 0) {
        randomSeed = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      } else {
        randomSeed = static_cast<uint64_t>(Options.randomizePairSeed.value);
      }
      rng = std::mt19937(randomSeed);
    }

    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kANALYSIS>(&hRegistry, colHistSpec);

    std::map<trackhistmanager::TrackHist, std::vector<framework::AxisSpec>> trackHistSpec = {
      {trackhistmanager::kPt, {confTrackBinning1.pt}},
      {trackhistmanager::kEta, {confTrackBinning1.eta}},
      {trackhistmanager::kPhi, {confTrackBinning1.phi}}};

    auto trackHistSpec1 = makeTrackHistSpecMap(confTrackBinning1);
    trackHistManager1.init<modes::Mode::kANALYSIS>(&hRegistry, trackHistSpec1);
    if (!Options.sameSpecies.value) {
      auto trackHistSpec2 = makeTrackHistSpecMap(confTrackBinning2);
      trackHistManager2.init<modes::Mode::kANALYSIS>(&hRegistry, trackHistSpec2);
    }
    std::map<pairhistmanager::PairHist, std::vector<framework::AxisSpec>> pairHistSpec = {
      {pairhistmanager::kKstar, {PairBinning.kstar}},
      {pairhistmanager::kKt, {PairBinning.kt}},
      {pairhistmanager::kMt, {PairBinning.mt}},
      {pairhistmanager::kPt1VsPt2, {confTrackBinning1.pt, confTrackBinning1.pt}},
      {pairhistmanager::kPt1VsKstar, {confTrackBinning1.pt, PairBinning.kstar}},
      {pairhistmanager::kPt2VsKstar, {confTrackBinning1.pt, PairBinning.kstar}},
      {pairhistmanager::kPt1VsKt, {confTrackBinning1.pt, PairBinning.kt}},
      {pairhistmanager::kPt2VsKt, {confTrackBinning1.pt, PairBinning.kt}},
      {pairhistmanager::kPt1VsMt, {confTrackBinning1.pt, PairBinning.mt}},
      {pairhistmanager::kPt2VsMt, {confTrackBinning1.pt, PairBinning.mt}}};

    pairHistManagerSe.init<modes::Mode::kANALYSIS>(&hRegistry, pairHistSpec);
    pairHistManagerSe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);

    pairHistManagerMe.init<modes::Mode::kANALYSIS>(&hRegistry, pairHistSpec);
    pairHistManagerMe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);

    std::map<closepairrejection::CprHist, std::vector<framework::AxisSpec>> cprHistSpec = {
      {closepairrejection::kAverage, {ConfCpr.deta, ConfCpr.dphistar}},
      {closepairrejection::kRadius0, {ConfCpr.deta, ConfCpr.dphistar}},
      {closepairrejection::kRadius1, {ConfCpr.deta, ConfCpr.dphistar}}};
    cprSe.init<modes::Mode::kANALYSIS>(&hRegistry, cprHistSpec);
    cprSe.setLimits(ConfCpr.detaMax.value, ConfCpr.dphistarMax.value);
    cprMe.init<modes::Mode::kANALYSIS>(&hRegistry, cprHistSpec);
    cprMe.setLimits(ConfCpr.detaMax.value, ConfCpr.dphistarMax.value);
  };

  template <modes::Mode mode, typename P1, typename P2>
  void doSameEvent(P1 SliceTrk1, P2 SliceTrk2)
  {
    for (auto const& part : SliceTrk1) {
      trackHistManager1.fill<mode>(part);
    }

    if (!Options.sameSpecies.value) {
      for (auto const& part : SliceTrk2) {
        trackHistManager2.fill<mode>(part);
      }
    }

    /// Now build the combinations
    float rand = 0.;
    if (Options.sameSpecies.value) {
      for (auto const& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (ConfCpr.on.value) {
          cprSe.setPair(p1, p2);
          if (cprSe.isClosePair()) {
            continue;
          }
        }
        // track cleaning
        if (!pc.isCleanTrackPair(p1, p2)) {
          continue;
        }
        // randomize track1 and track2 if configured
        if (Options.randomizePairSeed.value >= 0) {
          std::uniform_real_distribution<float> dist(0, 1);
          rand = dist(rng);
        }
        if (rand <= 0.5) {
          pairHistManagerSe.setPair(p1, p2);
        } else {
          pairHistManagerSe.setPair(p2, p1);
        }
        cprSe.fill<mode>();
        pairHistManagerSe.fill<mode>();
      }
    } else {
      for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(SliceTrk1, SliceTrk2))) {
        if (ConfCpr.on.value) {
          cprSe.setPair(p1, p2);
          if (cprSe.isClosePair()) {
            continue;
          }
        }
        // pair cleaning
        if (!pc.isCleanTrackPair(p1, p2)) {
          continue;
        }
        cprSe.fill<mode>();
        pairHistManagerSe.setPair(p1, p2);
        pairHistManagerSe.fill<mode>();
      }
    }
  }

  template <modes::Mode mode, typename CT, typename TT, typename BP>
  void doMixedEvent(CT const& cols, TT const& /*parts*/, BP policy)
  {
    for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.depth.value, -1, cols, cols)) {
      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }
      auto sliceTrk1 = TrackPartition1->sliceByCached(aod::femtobase::collisionId, collision1.globalIndex(), cache);
      auto sliceTrk2 = TrackPartition1->sliceByCached(aod::femtobase::collisionId, collision2.globalIndex(), cache);
      if (sliceTrk1.size() == 0 || sliceTrk2.size() == 0) {
        continue;
      }
      for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceTrk2))) {
        if (ConfCpr.on.value) {
          cprMe.setMagField(collision1.magField());
          cprMe.setPair(p1, p2);
          if (cprMe.isClosePair()) {
            continue;
          }
        }
        // pair cleaning
        if (!pc.isCleanTrackPair(p1, p2)) {
          continue;
        }
        cprMe.fill<mode>();
        pairHistManagerMe.setPair(p1, p2);
        pairHistManagerMe.fill<mode>();
      }
    }
  }

  void processSameEvent(FilteredCollision const& col, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto trackSlice1 = TrackPartition1->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    auto trackSlice2 = TrackPartition2->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    if (trackSlice1.size() == 0 || trackSlice2.size() == 0) {
      return;
    }
    cprSe.setMagField(col.magField());
    doSameEvent<modes::Mode::kANALYSIS>(trackSlice1, trackSlice1);
  }
  PROCESS_SWITCH(PairTrackTrack, processSameEvent, "Enable processing same event", true);

  void processMixedEvent(FilteredCollisions const& cols, Tracks const& tracks)
  {
    switch (Mixing.policy.value) {
      case 0:
        doMixedEvent<modes::Mode::kANALYSIS>(cols, tracks, colBinningMult);
        break;
      case 1:
        doMixedEvent<modes::Mode::kANALYSIS>(cols, tracks, colBinningMultPercentile);
        break;
      case 2:
        doMixedEvent<modes::Mode::kANALYSIS>(cols, tracks, colBinningMultMultPercentile);
        break;
      default:
        LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
    }
  }
  PROCESS_SWITCH(PairTrackTrack, processMixedEvent, "Enable processing mixed event", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<PairTrackTrack>(cfgc),
  };
  return workflow;
}
