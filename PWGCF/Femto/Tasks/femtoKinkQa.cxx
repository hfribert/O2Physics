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

/// \file femtoKinkQa.cxx
/// \brief QA task for kinks
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/collisionHistManager.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/trackHistManager.h"
#include "PWGCF/Femto/Core/kinkBuilder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/InitContext.h"
#include "Framework/OutputObjHeader.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femto;

struct FemtoKinkQa {

  // setup for collisions
  collisionbuilder::ConfCollisionFilter collisionSelection;
  Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  colhistmanager::CollisionHistManager<modes::Mode::kAnalysis_Qa> colHistManager;
  colhistmanager::ConfCollisionBinning confCollisionBinning;

  // using Collisions = o2::soa::Join<FUCols, FUColPos, FUColMults, FUColCents>;
  using Collisions = FCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  // Define kink/sigma tables (joining tables for comprehensive information)
  using Sigmas = o2::soa::Join<FUSigmas, FUSigmaMasks, FUSigmaExtras>;
  using Tracks = o2::soa::Join<FTracks, FTrackDcas, FTrackExtras, FTrackPids>;
  using KinkCandidates = o2::aod::KinkCands;

  SliceCache cache;

  // setup for sigmas
  kinkbuilder::ConfSigmaSelection1 confSigmaSelection;

  // Setup partition to filter sigma particles based on selection criteria
  Partition<Sigmas> sigmaPartition = (aod::femtokinks::mask & confSigmaSelection.mask) == confSigmaSelection.mask;
  Preslice<Sigmas> perColSigmas = aod::femtobase::stored::collisionId;

  // Setup binning configurations for histogram creation
  kinkhistmanager::ConfSigmaBinning1 confSigmaBinning;
  kinkhistmanager::ConfSigmaQaBinning1 confSigmaQaBinning;
  
  // Create the histogram manager for Sigma particles
  kinkhistmanager::KinkHistManager<
    kinkhistmanager::PrefixSigmaQa,
    modes::Mode::kAnalysis_Qa,
    modes::Kink::kSigma>
    sigmaHistManager;

  // setup for daughter tracks
  trackhistmanager::ConfV0PosDauBinning confV0PosDaughterBinning;
  trackhistmanager::ConfV0PosDauQaBinning confV0PosDaughterQaBinning;

  // For tracking daughter particles
  trackhistmanager::ConfTrackBinning confChaDaughterBinning;
  trackhistmanager::ConfTrackQaBinning confChaDaughterQaBinning;

  HistogramRegistry hRegistry{"FemtoKinkQa", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    // Create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init(&hRegistry, colHistSpec);

    // Create charged daughter track histogram specs
    auto chaDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confChaDaughterBinning, confChaDaughterQaBinning);

    if ((doprocessSigma + doprocessKinkCandidates) > 1) {
      LOG(fatal) << "Only one process can be activated";
    }

    // Initialize sigma histograms
    if (doprocessSigma || doprocessKinkCandidates) {
      auto sigmaHistSpec = kinkhistmanager::makeKinkQaHistSpecMap(confSigmaBinning, confSigmaQaBinning);
      sigmaHistManager.init(&hRegistry, sigmaHistSpec);
    }
  };

  // Process function for sigma particles from femto tables
  void processSigma(FilteredCollision const& col, Sigmas const& /*sigmas*/, Tracks const& tracks)
  {
    colHistManager.fill(col);
    auto sigmaSlice = sigmaPartition->sliceByCached(femtobase::stored::collisionId, col.globalIndex(), cache);
    for (auto const& sigma : sigmaSlice) {
      sigmaHistManager.fill(sigma, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processSigma, "Process sigmas from femto tables", false);

  // Process function for kink candidates from PWGLF tables
  void processKinkCandidates(FilteredCollision const& col, KinkCandidates const& kinks, Tracks const& tracks)
  {
    colHistManager.fill(col);
    
    // Filter kinks belonging to this collision
    for (auto const& kink : kinks) {
      if (kink.collisionId() != col.globalIndex()) {
        continue;
      }
      
      // Adapter to map kinkCand properties to the format expected by histograms
      struct KinkAdapter {
        float pt() const { return kinkObj.ptMoth(); }
        float eta() const { return kinkObj.eta(); }
        float phi() const { return kinkObj.phi(); }
        float mass() const { return kinkObj.mSigmaMinus(); }
        float sign() const { return kinkObj.mothSign() < 0 ? -1.0f : 1.0f; }
        float kinkAngle() const { return kinkObj.kinkAngle(); }
        float dcaMothToPV() const { return kinkObj.dcaMothPv(); }
        float dcaDaugToPV() const { return kinkObj.dcaDaugPv(); }
        float decayVtxX() const { return kinkObj.xDecVtx(); }
        float decayVtxY() const { return kinkObj.yDecVtx(); }
        float decayVtxZ() const { return kinkObj.zDecVtx(); }
        float transRadius() const { return kinkObj.transRadius(); }
        
        KinkCandidates::iterator kinkObj;
      } adapter{kink};
      
      sigmaHistManager.fill(adapter, tracks);
    }
  }
  PROCESS_SWITCH(FemtoKinkQa, processKinkCandidates, "Process kink candidates from LF tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<FemtoKinkQa>(cfgc),
  };
  return workflow;
}