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

/// \file TwoTrackResonanceQa.cxx
/// \brief Tasks that reads the particle tables and fills QA histograms for vzeros
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/FemtoUnited/Core/CollisionHistManager.h"
#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"
#include "PWGCF/FemtoUnited/Core/TrackHistManager.h"
#include "PWGCF/FemtoUnited/Core/TwoTrackResonanceHistManager.h"
#include "PWGCF/FemtoUnited/Core/TwoTrackResonanceSelection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTwoTrackResonancesDerived.h"

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

struct TwoTrackResonanceQa {

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

  using Collisions = FUCols;
  using Collision = Collisions::iterator;

  using FilteredCollisions = o2::soa::Filtered<Collisions>;
  using FilteredCollision = FilteredCollisions::iterator;

  colhistmanager::ConfCollisionBinning confCollisionBinning;

  using Resonances = o2::soa::Join<FUResos, FUResoMasks>;
  using Tracks = o2::soa::Join<FUTracks, FUTrackMasks, FUTrackDCAs, FUTrackExtras, FUTrackPids>;

  SliceCache cache;

  twotrackresonanceselection::ConfTwoTrackResonaceSelection1 confResonanceSelection1;

  Partition<Resonances> resonancePartition =
    (femtobase::pt > confResonanceSelection1.ptMin) &&
    (femtobase::pt < confResonanceSelection1.ptMax) &&
    (femtobase::eta > confResonanceSelection1.etaMin) &&
    (femtobase::eta < confResonanceSelection1.etaMax) &&
    (femtobase::phi > confResonanceSelection1.phiMin) &&
    (femtobase::phi < confResonanceSelection1.phiMax) &&
    (femtotwotrackresonances::resonanceMass > confResonanceSelection1.massMin) &&
    (femtotwotrackresonances::resonanceMass < confResonanceSelection1.massMax);
  // ncheckbit(femtotwotrackresonances::resonanceType, confResonanceSelection1.type) &&
  // ifnode(ncheckbit(femtotwotrackresonances::resonanceType, (o2::aod::femtodatatypes::TwoTrackResonaceType)o2::analysis::femtounited::modes::TwoTrackResonace::kPosDauAboveThres),
  //        ncheckbit(femtotwotrackresonances::resonanceMask, confResonanceSelection1.maskPosDauAboveThres),
  //        ncheckbit(femtotwotrackresonances::resonanceMask, confResonanceSelection1.maskPosDauBelowThres)) &&
  //
  // ifnode(ncheckbit(femtotwotrackresonances::resonanceType, (o2::aod::femtodatatypes::TwoTrackResonaceType)o2::analysis::femtounited::modes::TwoTrackResonace::kNegDauAboveThres),
  //        ncheckbit(femtotwotrackresonances::resonanceMask, confResonanceSelection1.maskNegDauAboveThres),
  //        ncheckbit(femtotwotrackresonances::resonanceMask, confResonanceSelection1.maskNegDauBelowThres));
  // ncheckbit(femtotwotrackresonances::resonanceType, confResonanceSelection1.type) &&
  // ncheckbit(femtotwotrackresonances::resonanceMask, confResonanceSelection1.mask);

  Preslice<Resonances> perColReco = aod::femtobase::collisionId;

  twotrackresonancehistmanager::ConfTwoTrackResonanceBinning1 confResonanceBinning;

  trackhistmanager::ConfResonancePosDauBinning confPosDaughterBinning;
  trackhistmanager::ConfResonancePosDauQaBinning confPosDaughterQaBinning;
  trackhistmanager::ConfResonanceNegDauBinning confNegDaughterBinning;
  trackhistmanager::ConfResonanceNegDauQaBinning confNegDaughterQaBinning;

  HistogramRegistry hRegistry{"ResonanceQA", {}, OutputObjHandlingPolicy::AnalysisObject};

  colhistmanager::CollisionHistManager colHistManager;
  twotrackresonancehistmanager::TwoTrackResonanceHistManager<twotrackresonancehistmanager::PrefixTwoTrackResonanceQa> resonanceHistManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixResonancePosDaughterQa> PosDaughterManager;
  trackhistmanager::TrackHistManager<trackhistmanager::PrefixResonanceNegDaughterQa> NegDaughterManager;

  void init(InitContext&)
  {
    // create a map for histogram specs
    auto colHistSpec = colhistmanager::makeColHistSpecMap(confCollisionBinning);
    colHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, colHistSpec);

    auto resonanceHistSpec = twotrackresonancehistmanager::makeTwoTrackResonanceQaHistSpecMap(confResonanceBinning);
    resonanceHistManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, resonanceHistSpec);

    auto posDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confPosDaughterBinning, confPosDaughterQaBinning);
    PosDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, posDaughterHistSpec);

    auto negDaughterHistSpec = trackhistmanager::makeTrackQaHistSpecMap(confNegDaughterBinning, confNegDaughterQaBinning);
    NegDaughterManager.init<modes::Mode::kANALYSIS_QA>(&hRegistry, negDaughterHistSpec);
  };

  void process(FilteredCollision const& col, Resonances const& /*resonances*/, Tracks const& /*tracks*/)
  {
    colHistManager.fill<modes::Mode::kANALYSIS_QA>(col);
    auto resonanceSlice = resonancePartition->sliceByCached(femtobase::collisionId, col.globalIndex(), cache);
    for (auto const& resonance : resonanceSlice) {
      resonanceHistManager.fill<modes::Mode::kANALYSIS_QA>(resonance);
      auto posDaugther = resonance.posDauResonance_as<Tracks>();
      PosDaughterManager.fill<modes::Mode::kANALYSIS_QA>(posDaugther);
      auto negDaugther = resonance.negDauResonance_as<Tracks>();
      NegDaughterManager.fill<modes::Mode::kANALYSIS_QA>(negDaugther);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<TwoTrackResonanceQa>(cfgc),
  };
  return workflow;
}
