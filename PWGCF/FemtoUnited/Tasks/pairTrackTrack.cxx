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

/// \file pairTrackTrack.cxx
/// \brief Tasks that computes correlation between two tracks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/closePairRejection.h"
#include "PWGCF/FemtoUnited/Core/collisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/collisionSelection.h"
#include "PWGCF/FemtoUnited/Core/modes.h"
#include "PWGCF/FemtoUnited/Core/pairCleaner.h"
#include "PWGCF/FemtoUnited/Core/pairHistManager.h"
#include "PWGCF/FemtoUnited/Core/partitions.h"
#include "PWGCF/FemtoUnited/Core/trackHistManager.h"
#include "PWGCF/FemtoUnited/Core/trackSelection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

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

  // setup tables
  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  SliceCache cache;

  // setup collisions
  collisionselection::ConfCollisionSelection collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);
  colhistmanager::ConfCollisionBinning confCollisionBinning;
  colhistmanager::CollisionHistManager colHistManager;

  // setup tracks
  trackselection::ConfTrackSelection1 trackSelections1;
  trackhistmanager::ConfTrackBinning1 confTrackBinning1;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack1> trackHistManager1;
  Partition<Tracks> trackPartition1 = MAKE_TRACK_PARTITION(trackSelections1);

  trackselection::ConfTrackSelection2 trackSelections2;
  trackhistmanager::ConfTrackBinning2 confTrackBinning2;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixTrack2> trackHistManager2;
  Partition<Tracks> trackPartition2 = MAKE_TRACK_PARTITION(trackSelections1);

  Preslice<Tracks> perColReco = aod::femtobase::stored::collisionId;

  // setup pairs
  pairhistmanager::ConfPairBinning confPairBinning;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackSe> pairHistManagerSe;
  pairhistmanager::PairHistManager<pairhistmanager::PrefixTrackTrackMe> pairHistManagerMe;

  // setup mixing
  std::vector<double> defaultVtxBins{10, -10, 10};
  std::vector<double> defaultMultBins{50, 0, 200};
  std::vector<double> defaultCentBins{10, 0, 100};
  ColumnBinningPolicy<femtocollisions::PosZ, femtocollisions::Mult> mixBinningMult{{defaultVtxBins, defaultMultBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Cent> mixBinningCent{{defaultVtxBins, defaultCentBins}, true};
  ColumnBinningPolicy<aod::femtocollisions::PosZ, aod::femtocollisions::Mult, aod::femtocollisions::Cent> mixBinningMultCent{{defaultVtxBins, defaultMultBins, defaultCentBins}, true};
  pairhistmanager::ConfMixing confMixing;

  HistogramRegistry hRegistry{"FemtoTrackTrack", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::mt19937 rng;

  // setup cpr
  closepairrejection::ConfCpr confCpr;
  closepairrejection::ClosePairRejection<closepairrejection::PrefixTrackTrackSe> cprSe;
  closepairrejection::ClosePairRejection<closepairrejection::PrefixTrackTrackMe> cprMe;
  paircleaner::PairCleaner pc;

  void init(InitContext&)
  {

    // setup columnpolicy for binning
    // default values are used during instantiation, so we need to explicity update them here
    mixBinningMult = {{confMixing.vtxBins, confMixing.multBins.value}, true};
    mixBinningCent = {{confMixing.vtxBins.value, confMixing.centBins.value}, true};
    mixBinningMultCent = {{confMixing.vtxBins.value, confMixing.multBins.value, confMixing.centBins.value}, true};

    // setup histograms for tracks
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kANALYSIS>(&hRegistry, colHistSpec);

    auto trackHistSpec1 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning1);
    trackHistManager1.init<modes::Mode::kANALYSIS>(&hRegistry, trackHistSpec1);

    auto trackHistSpec2 = trackhistmanager::makeTrackHistSpecMap(confTrackBinning2);
    trackHistManager2.init<modes::Mode::kANALYSIS>(&hRegistry, trackHistSpec2);

    // setup histograms for pair
    auto pairHistSpec = pairhistmanager::makePairHistSpecMap(confPairBinning, confTrackBinning1, confTrackBinning2);
    pairHistManagerSe.init<modes::Mode::kANALYSIS>(&hRegistry, pairHistSpec);
    pairHistManagerSe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);
    pairHistManagerMe.init<modes::Mode::kANALYSIS>(&hRegistry, pairHistSpec);
    pairHistManagerMe.setMass(trackSelections1.pdgCode.value, trackSelections2.pdgCode.value);

    // setup histograms for cpr
    auto cprHistSpec = closepairrejection::makeCprHistSpecMap(confCpr);
    cprSe.init<modes::Mode::kANALYSIS>(&hRegistry, cprHistSpec);
    cprSe.setLimits(confCpr.detaMax, confCpr.dphistarMax);
    cprMe.init<modes::Mode::kANALYSIS>(&hRegistry, cprHistSpec);
    cprMe.setLimits(confCpr.detaMax, confCpr.dphistarMax);

    // setup rng if necessary
    if (Options.randomizePairSeed.value >= 0) {
      uint64_t randomSeed = 0;
      if (Options.randomizePairSeed.value == 0) {
        randomSeed = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      } else {
        randomSeed = static_cast<uint64_t>(Options.randomizePairSeed.value);
      }
      rng = std::mt19937(randomSeed);
    }
  };

  template <modes::Mode mode, typename P1, typename P2>
  void doSameEvent(P1 SliceTrk1, P2 SliceTrk2)
  {
    // fill single particle histograms
    for (auto const& part : SliceTrk1) {
      trackHistManager1.fill<mode>(part);
    }

    if (!Options.sameSpecies.value) {
      for (auto const& part : SliceTrk2) {
        trackHistManager2.fill<mode>(part);
      }
    }

    /// Now build the combinations
    std::uniform_real_distribution<float> dist(0, 1);
    float rand = 0.;
    if (Options.sameSpecies.value) {
      for (auto const& [p1, p2] : combinations(CombinationsStrictlyUpperIndexPolicy(SliceTrk1, SliceTrk1))) {

        if (confCpr.on.value) {
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
        if (confCpr.on.value) {
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

  // template <modes::Mode mode, typename CT, typename TT, typename BP>
  // void doMixedEvent(CT const& cols, TT const& /*parts*/, BP policy)
  // {
  //   for (auto const& [collision1, collision2] : soa::selfCombinations(policy, Mixing.depth.value, -1, cols, cols)) {
  //     if (collision1.globalIndex() == collision2.globalIndex()) {
  //       continue;
  //     }
  //     auto sliceTrk1 = TrackPartition1->sliceByCached(aod::femtobase::collisionId, collision1.globalIndex(), cache);
  //     auto sliceTrk2 = TrackPartition1->sliceByCached(aod::femtobase::collisionId, collision2.globalIndex(), cache);
  //     if (sliceTrk1.size() == 0 || sliceTrk2.size() == 0) {
  //       continue;
  //     }
  //     for (auto const& [p1, p2] : combinations(CombinationsFullIndexPolicy(sliceTrk1, sliceTrk2))) {
  //       if (ConfCpr.on.value) {
  //         cprMe.setMagField(collision1.magField());
  //         cprMe.setPair(p1, p2);
  //         if (cprMe.isClosePair()) {
  //           continue;
  //         }
  //       }
  //       // pair cleaning
  //       if (!pc.isCleanTrackPair(p1, p2)) {
  //         continue;
  //       }
  //       cprMe.fill<mode>();
  //       pairHistManagerMe.setPair(p1, p2);
  //       pairHistManagerMe.fill<mode>();
  //     }
  //   }
  // }

  void processSameEvent(FilteredCollision const& col, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS>(col);
    auto trackSlice1 = trackPartition1->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    auto trackSlice2 = trackPartition2->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    if (trackSlice1.size() == 0 || trackSlice2.size() == 0) {
      return;
    }
    cprSe.setMagField(col.magField());
    doSameEvent<modes::Mode::kANALYSIS>(trackSlice1, trackSlice1);
  }
  PROCESS_SWITCH(PairTrackTrack, processSameEvent, "Enable processing same event", true);

  // void processMixedEvent(FilteredCollisions const& cols, Tracks const& tracks)
  // {
  //   switch (Mixing.policy.value) {
  //     case 0:
  //       doMixedEvent<modes::Mode::kANALYSIS>(cols, tracks, colBinningMult);
  //       break;
  //     case 1:
  //       doMixedEvent<modes::Mode::kANALYSIS>(cols, tracks, colBinningMultPercentile);
  //       break;
  //     case 2:
  //       doMixedEvent<modes::Mode::kANALYSIS>(cols, tracks, colBinningMultMultPercentile);
  //       break;
  //     default:
  //       LOG(fatal) << "Invalid binning policiy specifed. Breaking...";
  //   }
  // }
  // PROCESS_SWITCH(PairTrackTrack, processMixedEvent, "Enable processing mixed event", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<PairTrackTrack>(cfgc),
  };
  return workflow;
}
