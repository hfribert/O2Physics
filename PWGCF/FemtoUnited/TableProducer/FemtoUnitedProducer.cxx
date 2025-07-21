// Copyright 2019-2024 CERN andhcopyright holders of ALICE O2.
//
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUnitedProducer.cxx
/// \brief Tasks that produces the all femto tables
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/FemtoUnited/Core/CascadeSelection.h"
#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/FemtoUtils.h"
#include "PWGCF/FemtoUnited/Core/TrackSelection.h"
// #include "PWGCF/FemtoUnited/Core/TwoTrackResonanceSelection.h"
#include "PWGCF/FemtoUnited/Core/V0Selection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCascadesDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
// #include "PWGCF/FemtoUnited/DataModel/FemtoTwoTrackResonancesDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoV0sDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "fairlogger/Logger.h"

#include <string>
#include <unordered_map>
#include <vector>

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

namespace o2::analysis::femtounited
{
namespace consumeddata
{
using Run3PpCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;
using Run3PpWithoutCentCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;

// using Run3Tracks = soa::Join<Tracks, TracksExtra, TracksDCA,
//                              pidTPCEl, pidTPCPi, pidTPCKa, pidTPCPr, pidTPCDe, pidTPCTr, pidTPCHe,
//                              pidTOFEl, pidTOFPi, pidTOFKa, pidTOFPr, pidTOFDe, pidTOFTr, pidTOFHe>;

using Run3FullPidTracks =
  soa::Join<Tracks, TracksExtra, TracksDCA,
            pidTPCFullEl, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTPCFullDe, pidTPCFullTr, pidTPCFullHe,
            pidTOFFullEl, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr, pidTOFFullDe, pidTOFFullTr, pidTOFFullHe,
            pidTOFbeta>;

using Run3PpVzeros = V0Datas;

using Run3PpCascades = CascDatas;

} // namespace consumeddata
} // namespace o2::analysis::femtounited

struct FemtoUnitedProducer {
  SliceCache cache;
  // preslicing
  Preslice<Tracks> perColTracks = track::collisionId;

  // produced objectes
  struct : ProducesGroup {
    Produces<FUCols> producedCollision;

    Produces<FUTracks> producedTracks;
    Produces<FUTrackMasks> producedTrackMasks;
    Produces<FUTrackDCAs> producedTrackDCAs;
    Produces<FUTrackExtras> producedTrackExtras;
    Produces<FUTrackPids> producedTrackPids;

    Produces<FULambdas> producedLambdas;
    Produces<FULambdaMasks> producedLambdaMasks;
    Produces<FULambdaExtras> producedLambdaExtras;

    Produces<FUKshorts> producedK0shorts;
    Produces<FUKshortMasks> producedK0shortMasks;
    Produces<FUKshortExtras> producedK0shortExtras;

    Produces<FUXis> producedXis;
    Produces<FUXiMasks> producedXiMasks;
    Produces<FUXiExtras> producedXiExtras;

    Produces<FUOmegas> producedOmegas;
    Produces<FUOmegaMasks> producedOmegaMasks;
    Produces<FUOmegaExtras> producedOmegaExtras;

    // Produces<FUPhis> producedPhis;
    // Produces<FUPhiMasks> producedPhiMasks;
    //
    // Produces<FUKstars> producedKstars;
    // Produces<FUKstarMasks> producedKstarMasks;
    //
    // Produces<FURhos> producedRhos;
    // Produces<FURhoMasks> producedRhoMasks;
  } products;

  // configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("General");
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL to ccdb"};
    Configurable<std::string> grpPath{"grpPath", "GLO/Config/GRPMagField", "Path to GRP object (Run3 -> GLO/Config/GRPMagField/Run2 -> GLO/GRP/GRP"};
    Configurable<bool> produceExtraTables{"produceExtraTables", false, "Flag to produce extra tables (for all actived tables)"};
    Configurable<bool> produceK0short{"produceK0short", false, "Flag to produce K0shorts"};
    Configurable<bool> produceLambda{"produceLambda", false, "Flag to produce Lambda"};
    Configurable<bool> produceXi{"produceXi", false, "Flag to produce Xi"};
    Configurable<bool> produceOmega{"produceOmega", false, "Flag to produce Omega"};
    Configurable<bool> producePhi{"producePhi", false, "Flag to produce Phi"};
    Configurable<bool> produceKstar{"produceKstar", false, "Flag to produce Kstar"};
    Configurable<bool> produceRho{"produceRho", false, "Flag to produce Rho"};
  } ConfOptions;

  // Event selections
  collisionselection::ConfCollisionSelection confCollisionFilter;
  Filter collisionFilter = o2::aod::collision::posZ >= confCollisionFilter.vtxZMin &&
                           o2::aod::collision::posZ <= confCollisionFilter.vtxZMax;
  collisionselection::CollisionSelection collisionSel;

  // filters for tracks
  trackselection::ConfTrackFilters confTrackFilters;
  Filter trackFilter = track::pt >= confTrackFilters.ptMin && track::pt <= confTrackFilters.ptMax &&
                       track::eta >= confTrackFilters.etaMin && track::eta <= confTrackFilters.etaMax &&
                       track::phi >= confTrackFilters.phiMin && track::phi <= confTrackFilters.phiMax;
  // track bits
  trackselection::ConfTrackBits confTrackBits;
  trackselection::TrackSelection trackSel;

  // lambda filters
  // most v0 columns are now dynamic columns, so we cannot prefilter anymore
  v0selection::ConfV0Filters confV0Filters;

  // K0short bits
  v0selection::ConfK0shortBits confK0shortBits;
  v0selection::V0Selection k0shortSel;

  // lambda bits
  v0selection::ConfLambdaBits confLambdaBits;
  v0selection::V0Selection lambdaSel;

  // cascade filters
  cascadeselection::ConfCascadeFilters confCascadeFilters;

  // xi bits
  cascadeselection::ConfXiBits confXiBits;
  cascadeselection::CascadeSelection xiSel;

  // omega bits
  cascadeselection::ConfOmegaBits confOmegaBits;
  cascadeselection::CascadeSelection omegaSel;

  // resonance filters
  // twotrackresonanceselection::ConfTwoTrackResonanceDaughterFilters confResonanceDaughterFilters;
  // twotrackresonanceselection::ConfRhoFilters confRhoFilters;
  // twotrackresonanceselection::ConfPhiFilters confPhiFilters;
  // twotrackresonanceselection::ConfKstarFilters confKstarFilters;
  // twotrackresonanceselection::ConfAntiKstarFilters confAntiKstarFilters;

  // phi bits
  // twotrackresonanceselection::ConfPhiBits confPhiBits;
  // twotrackresonanceselection::TwoTrackResonanceSelection phiSels;

  // kstar bits
  // twotrackresonanceselection::ConfKstarBits confKstarBits;
  // twotrackresonanceselection::TwoTrackResonanceSelection kstarSels;
  // twotrackresonanceselection::TwoTrackResonanceSelection kstarBarSels;

  // rho bits
  // twotrackresonanceselection::ConfRhoBits confRhoBits;
  // twotrackresonanceselection::TwoTrackResonanceSelection rhoSels;

  // Partition<Filtered<consumeddata::Run3FullPidTracks>> partitionPositiveTracks = track::signed1Pt > 0.f;
  // Partition<Filtered<consumeddata::Run3FullPidTracks>> partitionNegativeTracks = track::signed1Pt < 0.f;

  // histogramming
  HistogramRegistry hRegistry{"Producer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // data members
  int runNumber = -1;
  float magField = 0.f;
  Service<o2::ccdb::BasicCCDBManager> ccdb;            /// Accessing the CCDB
  std::unordered_map<int64_t, int64_t> indexMapTracks; // for mapping tracks to lambdas, cascades and resonances

  // functions
  void initFromCcdb(o2::aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();

    static o2::parameters::GRPMagField* grpo = nullptr;
    grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ConfOptions.grpPath.value, timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
      return;
    }
    magField = 0.1 * grpo->getNominalL3Field(); // get magnetic field in tesla
    runNumber = bc.runNumber();
  };

  void init(InitContext& /*contex*/)
  {
    // init ccdb
    ccdb->setURL(ConfOptions.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // collision selection
    collisionSel.configure(confCollisionFilter);

    // init track selection objects
    trackSel.configure(confTrackBits, confTrackFilters);

    // init v0 selection ojects
    k0shortSel.configure<modes::V0::kK0short>(confK0shortBits, confV0Filters);
    k0shortSel.configure<modes::V0::kLambda>(confLambdaBits, confV0Filters);

    // cascade selections
    xiSel.configure<modes::Cascade::kXi>(confXiBits, confCascadeFilters);
    omegaSel.configure<modes::Cascade::kOmega>(confOmegaBits, confCascadeFilters);

    // resonance selections
    // rhoSels.configure<modes::TwoTrackResonace::kRho>(confRhoBits, confRhoFilters, confResonanceDaughterFilters);
    // phiSels.configure<modes::TwoTrackResonace::kPhi>(confPhiBits, confPhiFilters, confResonanceDaughterFilters);
    // kstarSels.configure<modes::TwoTrackResonace::kKstar>(confKstarBits, confKstarFilters, confResonanceDaughterFilters);
    // kstarBarSels.configure<modes::TwoTrackResonace::kKstarBar>(confKstarBits, confKstarFilters, confResonanceDaughterFilters);
  }

  template <modes::System sys, typename T>
  void fillCollision(T const& col)
  {
    if constexpr (!modes::isFlagSet(sys, modes::System::kNoCentCal)) {
      products.producedCollision(col.posZ(),
                                 col.multNTracksPV(),
                                 col.centFT0M(),
                                 collisionSel.getSphericity(),
                                 collisionSel.getMagneticField());
    }

    if constexpr (modes::isFlagSet(sys, modes::System::kNoCentCal)) {
      products.producedCollision(col.posZ(),
                                 col.multNTracksPV(),
                                 0,
                                 collisionSel.getSphericity(),
                                 collisionSel.getMagneticField());
    }
  }

  template <modes::Mode mode, modes::Track type, typename T1>
  void fillTrack(T1 const& track)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      products.producedTracks(products.producedCollision.lastIndex(),
                              track.pt() * track.sign(),
                              track.eta(),
                              track.phi());
      if constexpr (type == modes::Track::kPrimaryTrack) {
        products.producedTrackMasks(trackSel.getBitmask());
      } else {
        products.producedTrackMasks(static_cast<femtodatatypes::TrackMaskType>(0u));
      }
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      products.producedTrackDCAs(track.dcaXY(), track.dcaZ());
      products.producedTrackExtras(track.isPVContributor(),
                                   track.itsNCls(),
                                   track.itsNClsInnerBarrel(),
                                   track.itsChi2NCl(),
                                   track.itsClusterSizes(),
                                   track.tpcSignal(),
                                   track.tpcInnerParam(),
                                   track.tpcNClsFound(),
                                   track.tpcNClsCrossedRows(),
                                   track.tpcNClsShared(),
                                   track.beta());

      if constexpr (type == modes::Track::kPrimaryTrack) {
        products.producedTrackPids(track.itsNSigmaEl(),
                                   track.itsNSigmaPi(),
                                   track.itsNSigmaKa(),
                                   track.itsNSigmaPr(),
                                   track.itsNSigmaDe(),
                                   track.itsNSigmaTr(),
                                   track.itsNSigmaHe(),
                                   track.tpcNSigmaEl(),
                                   track.tpcNSigmaPi(),
                                   track.tpcNSigmaKa(),
                                   track.tpcNSigmaPr(),
                                   track.tpcNSigmaDe(),
                                   track.tpcNSigmaTr(),
                                   track.tpcNSigmaHe(),
                                   track.tofNSigmaEl(),
                                   track.tofNSigmaPi(),
                                   track.tofNSigmaKa(),
                                   track.tofNSigmaPr(),
                                   track.tofNSigmaDe(),
                                   track.tofNSigmaTr(),
                                   track.tofNSigmaHe());
      } else {
        products.producedTrackPids(0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   track.tpcNSigmaEl(),
                                   track.tpcNSigmaPi(),
                                   track.tpcNSigmaKa(),
                                   track.tpcNSigmaPr(),
                                   track.tpcNSigmaDe(),
                                   track.tpcNSigmaTr(),
                                   track.tpcNSigmaHe(),
                                   track.tofNSigmaEl(),
                                   track.tofNSigmaPi(),
                                   track.tofNSigmaKa(),
                                   track.tofNSigmaPr(),
                                   track.tofNSigmaDe(),
                                   track.tofNSigmaTr(),
                                   track.tofNSigmaHe());
      }
    }
    indexMapTracks.emplace(track.globalIndex(), products.producedTracks.lastIndex());
  }

  template <modes::Mode mode, typename T>
  void fillTracks(T const& tracks)
  {
    for (const auto& track : tracks) {
      if (!trackSel.hasTofAboveThreshold(track)) {
        continue;
      }
      trackSel.applySelections(track);
      if (!trackSel.passesMinimalCuts()) {
        continue;
      }
      fillTrack<mode, modes::Track::kPrimaryTrack>(track);
    }
  }

  template <modes::Mode mode, typename T>
  void fillLambda(T const& v0, float sign, int posDaughterIndex, int negDaughterIndex)
  {
    float mass, massAnti;
    if (sign > 0.f) {
      mass = v0.mLambda();
      massAnti = v0.mAntiLambda();
    } else {
      mass = v0.mAntiLambda();
      massAnti = v0.mLambda();
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      products.producedLambdas(products.producedCollision.lastIndex(),
                               sign * v0.pt(),
                               v0.eta(),
                               v0.phi(),
                               mass,
                               posDaughterIndex,
                               negDaughterIndex);
      products.producedLambdaMasks(lambdaSel.getBitmask());
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      products.producedLambdaExtras(
        v0.v0cosPA(),
        v0.dcaV0daughters(),
        v0.x(),
        v0.y(),
        v0.z(),
        v0.v0radius(),
        massAnti,
        v0.mK0Short());
    }
  }

  template <modes::Mode mode, typename T>
  void fillK0short(T const& v0, int posDaughterIndex, int negDaughterIndex)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      products.producedK0shorts(products.producedCollision.lastIndex(),
                                v0.pt(),
                                v0.eta(),
                                v0.phi(),
                                v0.mK0Short(),
                                posDaughterIndex,
                                negDaughterIndex);
      products.producedLambdaMasks(k0shortSel.getBitmask());
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      products.producedK0shortExtras(
        v0.v0cosPA(),
        v0.dcaV0daughters(),
        v0.x(),
        v0.y(),
        v0.z(),
        v0.v0radius());
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fillV0s(T1 const& v0s, T2 const& fullTracks)
  {

    int64_t posDaughterIndex;
    int64_t negDaughterIndex;

    for (const auto& v0 : v0s) {
      if (v0.pt() < confV0Filters.ptMin.value || v0.pt() > confV0Filters.ptMax.value ||
          v0.eta() < confV0Filters.etaMin.value || v0.eta() > confV0Filters.etaMax.value ||
          v0.phi() < confV0Filters.phiMin.value || v0.phi() > confV0Filters.phiMax.value) {
        continue;
      }

      if (ConfOptions.produceLambda.value) {
        lambdaSel.applySelections(v0, fullTracks);
        // check if selections are passed
        if (lambdaSel.passesMinimalCuts() && (lambdaSel.checkLambdaHypothesis(v0) || lambdaSel.checkAntiLambdaHypothesis(v0))) {
          auto posDaughter = v0.template posTrack_as<T2>();
          auto negDaughter = v0.template negTrack_as<T2>();

          // get index of positive daughter
          auto resultPosDaughter = utils::getDaughterIndex(posDaughter.globalIndex(), indexMapTracks);
          if (resultPosDaughter) {
            posDaughterIndex = resultPosDaughter.value();
          } else {
            fillTrack<mode, modes::Track::kLambdaDaugher>(posDaughter);
            posDaughterIndex = products.producedTracks.lastIndex();
            indexMapTracks.emplace(posDaughter.globalIndex(), products.producedTracks.lastIndex());
          }

          // get index of negative daughter
          auto resultNegDaughter = utils::getDaughterIndex(negDaughter.globalIndex(), indexMapTracks);
          if (resultNegDaughter) {
            negDaughterIndex = resultNegDaughter.value();
          } else {
            fillTrack<mode, modes::Track::kLambdaDaugher>(negDaughter);
            negDaughterIndex = products.producedTracks.lastIndex();
            indexMapTracks.emplace(posDaughter.globalIndex(), products.producedTracks.lastIndex());
          }
          // fill lambda or antilambda
          if (lambdaSel.checkLambdaHypothesis(v0)) {
            fillLambda<mode>(v0, 1.f, posDaughterIndex, negDaughterIndex);
          }
          if (lambdaSel.checkAntiLambdaHypothesis(v0)) {
            fillLambda<mode>(v0, -1.f, posDaughterIndex, negDaughterIndex);
          }
        }
      }
      if (ConfOptions.produceK0short.value) {
        k0shortSel.applySelections(v0, fullTracks);
        if (k0shortSel.passesMinimalCuts() && lambdaSel.checkK0shortHypothesis(v0)) {
          auto posDaughter = v0.template posTrack_as<T2>();
          auto negDaughter = v0.template negTrack_as<T2>();

          // get index of positive daughter
          auto resultPosDaughter = utils::getDaughterIndex(posDaughter.globalIndex(), indexMapTracks);
          if (resultPosDaughter) {
            posDaughterIndex = resultPosDaughter.value();
          } else {
            fillTrack<mode, modes::Track::kLambdaDaugher>(posDaughter);
            posDaughterIndex = products.producedTracks.lastIndex();
            indexMapTracks.emplace(posDaughter.globalIndex(), products.producedTracks.lastIndex());
          }

          // get index of negative daughter
          auto resultNegDaughter = utils::getDaughterIndex(negDaughter.globalIndex(), indexMapTracks);
          if (resultNegDaughter) {
            negDaughterIndex = resultNegDaughter.value();
          } else {
            fillTrack<mode, modes::Track::kLambdaDaugher>(negDaughter);
            negDaughterIndex = products.producedTracks.lastIndex();
            indexMapTracks.emplace(posDaughter.globalIndex(), products.producedTracks.lastIndex());
          }
          fillK0short<mode>(v0, posDaughterIndex, negDaughterIndex);
        }
      }
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3>
  void fillCacades(T1 const& fullCascades, T2 const& fullTracks, T3 const& col)
  {
    for (const auto& cascade : fullCascades) {
      if (cascade.pt() < confCascadeFilters.ptMin.value || cascade.pt() > confCascadeFilters.ptMax.value ||
          cascade.eta() < confCascadeFilters.etaMin.value || cascade.eta() > confCascadeFilters.etaMax.value ||
          cascade.phi() < confCascadeFilters.phiMin.value || cascade.phi() > confCascadeFilters.phiMax.value ||
          cascade.mLambda() < confCascadeFilters.massLambdaMin.value || cascade.mLambda() > confCascadeFilters.massLambdaMax.value) {
        continue;
      }
      xiSel.applySelections(cascade, fullTracks, col);
      omegaSel.applySelections(cascade, fullTracks, col);

      bool xiFails = ConfOptions.produceXi.value &&
                     (!xiSel.passesMinimalCuts() || !xiSel.checkXiHypothesis(cascade));

      bool omegaFails = ConfOptions.produceOmega.value &&
                        (!omegaSel.passesMinimalCuts() || !omegaSel.checkOmegaHypothesis(cascade));

      if (xiFails && omegaFails) {
        continue;
      }

      auto bachelor = cascade.template bachelor_as<T2>();
      auto posDaughter = cascade.template posTrack_as<T2>();
      auto negDaughter = cascade.template negTrack_as<T2>();

      // daughters
      int64_t bachelorIndex;
      int64_t posDaughterIndex;
      int64_t negDaughterIndex;

      // get index of bachelor
      auto resultBachelor = utils::getDaughterIndex(bachelor.globalIndex(), indexMapTracks);
      if (resultBachelor) {
        bachelorIndex = resultBachelor.value();
      } else {
        fillTrack<mode, modes::Track::kLambdaDaugher>(bachelor);
        bachelorIndex = products.producedTracks.lastIndex();
      }

      // get index of positive daughter
      auto resultPosDaughter = utils::getDaughterIndex(posDaughter.globalIndex(), indexMapTracks);
      if (resultPosDaughter) {
        posDaughterIndex = resultPosDaughter.value();
      } else {
        fillTrack<mode, modes::Track::kLambdaDaugher>(posDaughter);
        posDaughterIndex = products.producedTracks.lastIndex();
      }

      // get index of negative daughter
      auto resultNegDaughter = utils::getDaughterIndex(negDaughter.globalIndex(), indexMapTracks);
      if (resultNegDaughter) {
        negDaughterIndex = resultNegDaughter.value();
      } else {
        fillTrack<mode, modes::Track::kLambdaDaugher>(negDaughter);
        negDaughterIndex = products.producedTracks.lastIndex();
      }
      if (!xiFails) {
        fillCascade<mode, o2::analysis::femtounited::modes::Cascade::kXi>(cascade, bachelorIndex, posDaughterIndex, negDaughterIndex, col);
      }
      if (!omegaFails) {
        fillCascade<mode, o2::analysis::femtounited::modes::Cascade::kOmega>(cascade, bachelorIndex, posDaughterIndex, negDaughterIndex, col);
      }
    }
  }

  template <modes::Mode mode, modes::Cascade C, typename T1, typename T2>
  void fillCascade(T1 const& cascade, int bachelorIndex, int posDaughterIndex, int negDaughterIndex, T2 const& col)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      if constexpr (modes::isFlagSet(C, modes::Cascade::kXi)) {
        products.producedXis(products.producedCollision.lastIndex(),
                             cascade.sign() * cascade.pt(),
                             cascade.eta(),
                             cascade.phi(),
                             cascade.mXi(),
                             bachelorIndex,
                             posDaughterIndex,
                             negDaughterIndex);
        products.producedXiMasks(xiSel.getBitmask());
      }
      if constexpr (modes::isFlagSet(C, modes::Cascade::kOmega)) {
        products.producedOmegas(products.producedCollision.lastIndex(),
                                cascade.sign() * cascade.pt(),
                                cascade.eta(),
                                cascade.phi(),
                                cascade.mOmega(),
                                bachelorIndex,
                                posDaughterIndex,
                                negDaughterIndex);
        products.producedOmegaMasks(omegaSel.getBitmask());
      }
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      if constexpr (modes::isFlagSet(C, modes::Cascade::kXi)) {
        products.producedXiExtras(
          cascade.mOmega(),
          cascade.casccosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcacascdaughters(),
          cascade.cascradius(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posY(), col.posY(), col.posZ()));
      }
      if constexpr (modes::isFlagSet(C, modes::Cascade::kOmega)) {
        products.producedOmegaExtras(
          cascade.mXi(),
          cascade.casccosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcacascdaughters(),
          cascade.cascradius(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posY(), col.posY(), col.posZ()));
      }
    }
  }

  // template <modes::Mode mode, typename T1, typename T2>
  // void fillResonance(T1 const& posDaughter, T1 const& negDaughter, T2& selectionContainer)
  // {
  //   if (!selectionContainer.hasTofAboveThreshold(posDaughter, negDaughter)) {
  //     return;
  //   }
  //
  //   selectionContainer.applySelections(posDaughter, negDaughter);
  //
  //   if (!selectionContainer.passesMinimalCuts()) {
  //     return;
  //   }
  //
  //   selectionContainer.reconstructResonance(posDaughter, negDaughter);
  //   if (!selectionContainer.checkFilters() || !selectionContainer.checkDaughterPids()) {
  //     return;
  //   }
  //
  //   int64_t posDaughterIndex;
  //   int64_t negDaughterIndex;
  //
  //   // get index of positive daughter
  //   auto resultPosDaughter = utils::getDaughterIndex(posDaughter.globalIndex(), indexMapTracks);
  //   if (resultPosDaughter) {
  //     posDaughterIndex = resultPosDaughter.value();
  //   } else {
  //     fillTrack<mode, modes::Track::kResonanceDaughter>(posDaughter);
  //     posDaughterIndex = products.producedTracks.lastIndex();
  //   }
  //
  //   // get index of negative daughter
  //   auto resultNegDaughter = utils::getDaughterIndex(negDaughter.globalIndex(), indexMapTracks);
  //   if (resultNegDaughter) {
  //     negDaughterIndex = resultNegDaughter.value();
  //   } else {
  //     fillTrack<mode, modes::Track::kResonanceDaughter>(negDaughter);
  //     negDaughterIndex = products.producedTracks.lastIndex();
  //   }
  //
  //   if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
  //     products.producedResonances(
  //       products.producedCollision.lastIndex(),
  //       selectionContainer.getPt(),
  //       selectionContainer.getEta(),
  //       selectionContainer.getPhi(),
  //       selectionContainer.getMass(),
  //       static_cast<femtodatatypes::TwoTrackResonaceType>(selectionContainer.getType()),
  //       posDaughterIndex,
  //       negDaughterIndex);
  //     products.producedResonanceMasks(selectionContainer.getBitmask());
  //   }
  // }
  //
  // template <modes::Mode mode, typename T1, typename T2>
  // void fillResonances(T1 const& col, T2 const& /*tracks*/)
  // {
  //   auto groupPositiveTracks = partitionPositiveTracks->sliceByCached(track::collisionId, col.globalIndex(), cache);
  //   auto groupNegativeTracks = partitionNegativeTracks->sliceByCached(track::collisionId, col.globalIndex(), cache);
  //   for (auto const& positiveTrack : groupPositiveTracks) {
  //     for (auto const& negativeTrack : groupNegativeTracks) {
  //       fillResonance<mode>(positiveTrack, negativeTrack, rhoSels);
  //       fillResonance<mode>(positiveTrack, negativeTrack, phiSels);
  //       fillResonance<mode>(positiveTrack, negativeTrack, kstarSels);
  //       fillResonance<mode>(positiveTrack, negativeTrack, antiKstarSels);
  //     }
  //   }
  // }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3>
  void processTracks(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks)
  {
    initFromCcdb(col.template bc_as<T2>());
    collisionSel.setMagneticField(magField);
    collisionSel.setSphericity(fullTracks);
    if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
      return;
    }
    fillCollision<system>(col);
    indexMapTracks.clear();
    fillTracks<mode>(fullTracks);
  }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processTracksV0s(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks, T4 const& fullTracksWithItsPid, T5 const& fullV0s)
  {
    initFromCcdb(col.template bc_as<T2>());
    collisionSel.setMagneticField(magField);
    collisionSel.setSphericity(fullTracks);
    if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
      return;
    }
    fillCollision<system>(col);
    indexMapTracks.clear();
    fillTracks<mode>(fullTracksWithItsPid);
    fillV0s<mode>(fullV0s, fullTracks);
  }

  // template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3>
  // void processTracksResonances(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks)
  // {
  //   initFromCcdb(col.template bc_as<T2>());
  //   collisionSel.setMagneticField(magField);
  //   collisionSel.setSphericity(fullTracks);
  //   if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
  //     return;
  //   }
  //   fillCollision<system>(col);
  //   indexMapTracks.clear();
  //   fillTracks<mode>(fullTracks);
  //   fillResonances<mode>(col, fullTracks);
  // }

  // template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  // void processTracksResonancesVzeros(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks, T4 const& fullTracksWithItsPid, T5 const& fullV0s)
  // {
  //   initFromCcdb(col.template bc_as<T2>());
  //   collisionSel.setMagneticField(magField);
  //   collisionSel.setSphericity(fullTracksWithItsPid);
  //   if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
  //     return;
  //   }
  //   fillCollision<system>(col);
  //   indexMapTracks.clear();
  //   fillTracks<mode>(fullTracksWithItsPid);
  //   fillResonances<mode>(col, fullTracks);
  //   fillV0s<mode>(fullV0s, fullTracks);
  // }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processTracksV0sCascades(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks, T4 const& fullTracksWithItsPid, T5 const& fullV0s, T6 const& fullCascades)
  {
    initFromCcdb(col.template bc_as<T2>());
    collisionSel.setMagneticField(magField);
    collisionSel.setSphericity(fullTracksWithItsPid);
    if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
      return;
    }
    fillCollision<system>(col);
    indexMapTracks.clear();
    fillTracks<mode>(fullTracksWithItsPid);
    fillV0s<mode>(fullV0s, fullTracks);
    fillCacades<mode>(fullCascades, fullTracks, col);
  }

  // proccess tracks
  void processTracksRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
                           BCsWithTimestamps const& bcs,
                           Filtered<consumeddata::Run3FullPidTracks> const& tracks)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    if (ConfOptions.produceExtraTables.value) {
      processTracks<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(col, bcs, tracksWithItsPid);
    } else {
      processTracks<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(col, bcs, tracksWithItsPid);
    }
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksRun3pp, "Process tracks", true);

  // process tracks and v0s
  void processTracksVzerosRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
                                 BCsWithTimestamps const& bcs,
                                 Filtered<consumeddata::Run3FullPidTracks> const& tracks,
                                 consumeddata::Run3PpVzeros const& v0s)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    if (ConfOptions.produceExtraTables.value) {
      processTracksV0s<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(col, bcs, tracks, tracksWithItsPid, v0s);
    } else {
      processTracksV0s<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(col, bcs, tracks, tracksWithItsPid, v0s);
    }
  };
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksVzerosRun3pp, "Process tracks and v0s", false);

  // process tracks, v0s and casacades
  void processTracksV0sCascadesRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
                                      BCsWithTimestamps const& bcs,
                                      Filtered<consumeddata::Run3FullPidTracks> const& tracks,
                                      consumeddata::Run3PpVzeros const& v0s,
                                      consumeddata::Run3PpCascades const& cascades)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    if (ConfOptions.produceExtraTables.value) {
      processTracksV0sCascades<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(col, bcs, tracks, tracksWithItsPid, v0s, cascades);
    } else {
      processTracksV0sCascades<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(col, bcs, tracks, tracksWithItsPid, v0s, cascades);
    }
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksV0sCascadesRun3pp, "Provide Tracks, V0s and Cascades for Run3", false);

  // void processTracksResonancesVzerosRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
  //                                          BCsWithTimestamps const& bcs,
  //                                          Filtered<consumeddata::Run3FullPidTracks> const& tracks,
  //                                          consumeddata::Run3PpVzeros const& v0s)
  // {
  //   // its pid information is generated dynamically, so we need to add it here
  //   auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
  //                                           pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
  //   processTracksResonancesVzeros<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(col, bcs, tracks, tracksWithItsPid, v0s);
  // }
  // PROCESS_SWITCH(FemtoUnitedProducer, processTracksResonancesVzerosRun3pp, "Provide tracks, resonances and V0s for Run3 analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoUnitedProducer>(cfgc)};
  return workflow;
}
