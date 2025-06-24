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
/// \brief Tasks that produces the track tables used for the pairing
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/FemtoUnited/Core/CollisionSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/FemtoUtils.h"
#include "PWGCF/FemtoUnited/Core/LambdaSelection.h"
#include "PWGCF/FemtoUnited/Core/ResonanceSelection.h"
#include "PWGCF/FemtoUnited/Core/TrackSelection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoLambdasDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTwoTrackResonancesDerived.h"
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

using Run3Tracks = soa::Join<Tracks, TracksExtra, TracksDCA,
                             pidTPCEl, pidTPCPi, pidTPCKa, pidTPCPr, pidTPCDe, pidTPCTr, pidTPCHe,
                             pidTOFEl, pidTOFPi, pidTOFKa, pidTOFPr, pidTOFDe, pidTOFTr, pidTOFHe>;

using Run3TracksFullPid =
  soa::Join<Tracks, TracksExtra, TracksDCA,
            pidTPCFullEl, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTPCFullDe, pidTPCFullTr, pidTPCFullHe,
            pidTOFFullEl, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr, pidTOFFullDe, pidTOFFullTr, pidTOFFullHe,
            pidTOFbeta>;

using Run3PpVzeros = V0Datas;

} // namespace consumeddata
} // namespace o2::analysis::femtounited

struct FemtoUnitedProducer {

  // preslicing
  Preslice<Tracks> perColTracks = track::collisionId;
  Preslice<V0Datas> perColV0s = v0data::collisionId;

  // produced objectes
  struct ProducesGroup {
    Produces<FUCols> producedCollision;

    Produces<FUTracks> producedTracks;
    Produces<FUTrackMasks> producedTrackMasks;
    Produces<FUTrackDCAs> producedTrackDCAs;
    Produces<FUTrackExtras> producedTrackExtras;
    Produces<FUTrackPids> producedTrackPids;

    Produces<FULambdas> producedLambdas;
    Produces<FULambdaMasks> producedLambdaMasks;
    Produces<FULambdaExtras> producedLambdaExtras;

    Produces<FUResos> producedResonances;
    Produces<FUResoMasks> producedResonanceMasks;
  } products;

  // configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("General");
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL to ccdb"};
    Configurable<std::string> grpPath{"grpPath", "GLO/Config/GRPMagField", "Path to GRP object (Run3 -> GLO/Config/GRPMagField/Run2 -> GLO/GRP/GRP"};
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
  lambdaselection::ConfLambdaFilters confLambdaFilters;

  // lambda bits
  lambdaselection::ConfLambdaBits confLambdaBits;
  lambdaselection::LambdaSelection lambdaSel;

  // resonance filters
  twotrackresonanceselection::ConfTwoTrackResonanceDaughterFilters confResonanceDaughterFilters;
  twotrackresonanceselection::ConfRhoFilters confRhoFilters;
  twotrackresonanceselection::ConfPhiFilters confPhiFilters;
  twotrackresonanceselection::ConfKstarFilters confKstarFilters;

  // resonance bits
  twotrackresonanceselection::ConfRhoBits confRhoBits;
  twotrackresonanceselection::ConfPhiBits confPhiBits;
  twotrackresonanceselection::ConfKstarBits confKstarBits;

  twotrackresonanceselection::TwoTrackResonanceSelection rhoSels;
  twotrackresonanceselection::TwoTrackResonanceSelection phiSels;
  twotrackresonanceselection::TwoTrackResonanceSelection kstarSels;
  twotrackresonanceselection::TwoTrackResonanceSelection antiKstarSels;

  // histogramming
  HistogramRegistry hRegistry{"Producer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // data members
  int runNumber = -1;
  float magField = 0.f;
  Service<o2::ccdb::BasicCCDBManager> ccdb;             /// Accessing the CCDB
  std::unordered_map<int64_t, int64_t> indexMapTracks;  // for mapping tracks to lambdas, cascades and resonances
  std::unordered_map<int64_t, int64_t> indexMapLambdas; // for mapping lambdas to cascades

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
    collisionSel.useOfflineSelection(confCollisionFilter.useOfflineSelection.value);
    collisionSel.setMagneticFieldLimits(confCollisionFilter.magFieldMin.value, confCollisionFilter.magFieldMax.value);
    collisionSel.setMultiplicityLimits(confCollisionFilter.multMin.value, confCollisionFilter.multMax.value);
    collisionSel.setCentralityLimits(confCollisionFilter.centMin.value, confCollisionFilter.centMax.value);
    collisionSel.setSphericityLimits(confCollisionFilter.spherMin.value, confCollisionFilter.spherMax.value);

    // track selections
    // require TOF above set momentum for all tracks
    trackSel.setMinimalMomentumForTof(confTrackBits.minMomentumForTof.value);
    // add selections for track quality
    trackSel.addSelection(confTrackBits.tpcClustersMin.value, trackselection::kTPCnClsMin, limits::kLowerLimit, true, true);
    trackSel.addSelection(confTrackBits.tpcCrossedRowsMin.value, trackselection::kTPCcRowsMin, limits::kLowerLimit, true, true);
    trackSel.addSelection(confTrackBits.tpcSharedClustersMax.value, trackselection::kTPCsClsMax, limits::kUpperLimit, true, true);
    trackSel.addSelection(confTrackBits.tpcSharedClusterFractionMax.value, trackselection::kTPCsClsFracMax, limits::kUpperLimit, true, true);
    trackSel.addSelection(confTrackBits.itsClustersMin.value, trackselection::kITSnClsMin, limits::kLowerLimit, true, true);
    trackSel.addSelection(confTrackBits.itsIbClustersMin.value, trackselection::kITSnClsIbMin, limits::kLowerLimit, true, true);
    trackSel.addSelection(confTrackBits.dcaxyMax.name, confTrackFilters.ptMin.value, confTrackFilters.ptMax.value, confTrackBits.dcaxyMax.value, trackselection::kDCAxyMax, limits::kAbsUpperFunctionLimit, true, true);
    trackSel.addSelection(confTrackBits.dcazMax.name, confTrackFilters.ptMin.value, confTrackFilters.ptMax.value, confTrackBits.dcazMax.value, trackselection::kDCAzMax, limits::kAbsUpperFunctionLimit, true, true);
    // add selections for its pid
    trackSel.addSelection(confTrackBits.itsElectron.value, trackselection::kItsElectron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.itsPion.value, trackselection::kItsPion, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.itsKaon.value, trackselection::kItsKaon, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.itsProton.value, trackselection::kItsProton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.itsDeuteron.value, trackselection::kItsDeuteron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.itsTriton.value, trackselection::kItsTriton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.itsHelium.value, trackselection::kItsHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tpc pid
    trackSel.addSelection(confTrackBits.tpcElectron.value, trackselection::kTpcElectron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpcPion.value, trackselection::kTpcPion, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpcKaon.value, trackselection::kTpcKaon, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpcProton.value, trackselection::kTpcProton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpcDeuteron.value, trackselection::kTpcDeuteron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpcTriton.value, trackselection::kTpcTriton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpcHelium.value, trackselection::kTpcHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tof pid
    trackSel.addSelection(confTrackBits.tofElectron.value, trackselection::kTofElectron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tofPion.value, trackselection::kTofPion, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tofKaon.value, trackselection::kTofKaon, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tofProton.value, trackselection::kTofProton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tofDeuteron.value, trackselection::kTofDeuteron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tofTriton.value, trackselection::kTofTriton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tofHelium.value, trackselection::kTofHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tpctof pid
    trackSel.addSelection(confTrackBits.tpctofElectron.value, trackselection::kTpctofElectron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpctofPion.value, trackselection::kTpctofPion, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpctofKaon.value, trackselection::kTpctofKaon, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpctofProton.value, trackselection::kTpctofProton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpctofDeuteron.value, trackselection::kTpctofDeuteron, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpctofTriton.value, trackselection::kTpctofTriton, limits::kAbsUpperLimit, false, false);
    trackSel.addSelection(confTrackBits.tpctofHelium.value, trackselection::kTpctofHelium, limits::kAbsUpperLimit, false, false);

    /// lambda selections
    lambdaSel.setK0ShortMassLimits(confLambdaBits.kaonMassRejectionLow.value, confLambdaBits.kaonMassRejectionHigh.value);
    lambdaSel.addSelection(confLambdaBits.dcaDaughMax.value, lambdaselection::kDcaDaughMax, limits::kAbsUpperLimit, true, true);
    lambdaSel.addSelection(confLambdaBits.cpaMin.value, lambdaselection::kCpaMin, limits::kLowerLimit, true, true);
    lambdaSel.addSelection(confLambdaBits.transRadMin.value, lambdaselection::kTransRadMin, limits::kLowerLimit, true, true);
    lambdaSel.addSelection(confLambdaBits.transRadMax.value, lambdaselection::kTransRadMax, limits::kUpperLimit, true, true);
    // lambda positiv daughter selections
    lambdaSel.addSelection(confLambdaBits.posDaughDcaMin.value, lambdaselection::kPosDauDcaMin, limits::kLowerLimit, true, true);
    lambdaSel.addSelection(confLambdaBits.posDaughTpcClustersMin.value, lambdaselection::kPosDauTpcClsMin, limits::kLowerLimit, true, true);
    lambdaSel.addSelection(confLambdaBits.posDaughTpcPion.value, lambdaselection::kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
    lambdaSel.addSelection(confLambdaBits.posDaughTpcProton.value, lambdaselection::kPosDaughTpcProton, limits::kAbsUpperLimit, false, false);
    // lambda negative daughter selections
    lambdaSel.addSelection(confLambdaBits.negDaughDcaMin.value, lambdaselection::kNegDauDcaMin, limits::kLowerLimit, true, true);
    lambdaSel.addSelection(confLambdaBits.negDaughTpcClustersMin.value, lambdaselection::kNegDauTpcClsMin, limits::kLowerLimit, true, true);
    lambdaSel.addSelection(confLambdaBits.negDaughTpcPion.value, lambdaselection::kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
    lambdaSel.addSelection(confLambdaBits.negDaughTpcProton.value, lambdaselection::kNegDaughTpcProton, limits::kAbsUpperLimit, false, false);

    // rho selections
    rhoSels.setResonanceType(twotrackresonanceselection::TwoTrackResonace::kRho);
    rhoSels.setMinimalMomentumForTof(confRhoBits.minMomentumForTof.value);
    rhoSels.setPtLimits(confRhoFilters.ptMin.value, confRhoFilters.ptMax.value);
    rhoSels.setEtaLimits(confRhoFilters.etaMin.value, confRhoFilters.etaMax.value);
    rhoSels.setPhiLimits(confRhoFilters.phiMin.value, confRhoFilters.phiMax.value);
    rhoSels.setMassLimits(confRhoFilters.massMin.value, confRhoFilters.massMax.value);
    // rho positive daughter selections
    rhoSels.addSelection(confRhoBits.posDaughTpcClustersMin.value, twotrackresonanceselection::kPosDauTpcClsMin, limits::kUpperLimit, true, true);
    rhoSels.addSelection(confRhoBits.posDaughDcaxyMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confRhoBits.posDaughDcaxyMax.value, twotrackresonanceselection::kPosDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    rhoSels.addSelection(confRhoBits.posDaughDcazMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confRhoBits.posDaughDcazMax.value, twotrackresonanceselection::kPosDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    rhoSels.addSelection(confRhoBits.posDaughTpcPion.value, twotrackresonanceselection::kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
    rhoSels.addSelection(confRhoBits.posDaughTofPion.value, twotrackresonanceselection::kPosDaughTofPion, limits::kAbsUpperLimit, false, false);
    rhoSels.addSelection(confRhoBits.posDaughTpctofPion.value, twotrackresonanceselection::kPosDaughTpctofPion, limits::kAbsUpperLimit, false, false);
    // rho negative daughter selections
    rhoSels.addSelection(confRhoBits.negDaughTpcClustersMin.value, twotrackresonanceselection::kNegDauTpcClsMin, limits::kUpperLimit, true, true);
    rhoSels.addSelection(confRhoBits.negDaughDcaxyMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confRhoBits.negDaughDcaxyMax.value, twotrackresonanceselection::kNegDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    rhoSels.addSelection(confRhoBits.negDaughDcazMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confRhoBits.negDaughDcazMax.value, twotrackresonanceselection::kNegDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    rhoSels.addSelection(confRhoBits.negDaughTpcPion.value, twotrackresonanceselection::kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
    rhoSels.addSelection(confRhoBits.negDaughTofPion.value, twotrackresonanceselection::kNegDaughTofPion, limits::kAbsUpperLimit, false, false);
    rhoSels.addSelection(confRhoBits.negDaughTpctofPion.value, twotrackresonanceselection::kNegDaughTpctofPion, limits::kAbsUpperLimit, false, false);

    // rho selections
    phiSels.setResonanceType(twotrackresonanceselection::TwoTrackResonace::kPhi);
    phiSels.setMinimalMomentumForTof(confPhiBits.minMomentumForTof.value);
    phiSels.setPtLimits(confPhiFilters.ptMin.value, confPhiFilters.ptMax.value);
    phiSels.setEtaLimits(confPhiFilters.etaMin.value, confPhiFilters.etaMax.value);
    phiSels.setPhiLimits(confPhiFilters.phiMin.value, confPhiFilters.phiMax.value);
    phiSels.setMassLimits(confPhiFilters.massMin.value, confPhiFilters.massMax.value);
    // phi positive daughter selections
    phiSels.addSelection(confPhiBits.posDaughTpcClustersMin.value, twotrackresonanceselection::kPosDauTpcClsMin, limits::kUpperLimit, true, true);
    phiSels.addSelection(confPhiBits.posDaughDcaxyMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confPhiBits.posDaughDcaxyMax.value, twotrackresonanceselection::kPosDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    phiSels.addSelection(confPhiBits.posDaughDcazMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confPhiBits.posDaughDcazMax.value, twotrackresonanceselection::kPosDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    phiSels.addSelection(confPhiBits.posDaughTpcPion.value, twotrackresonanceselection::kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
    phiSels.addSelection(confPhiBits.posDaughTofPion.value, twotrackresonanceselection::kPosDaughTofPion, limits::kAbsUpperLimit, false, false);
    phiSels.addSelection(confPhiBits.posDaughTpctofPion.value, twotrackresonanceselection::kPosDaughTpctofPion, limits::kAbsUpperLimit, false, false);
    // phi negative daughter selections
    phiSels.addSelection(confPhiBits.negDaughTpcClustersMin.value, twotrackresonanceselection::kNegDauTpcClsMin, limits::kUpperLimit, true, true);
    phiSels.addSelection(confPhiBits.negDaughDcaxyMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confPhiBits.negDaughDcaxyMax.value, twotrackresonanceselection::kNegDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    phiSels.addSelection(confPhiBits.negDaughDcazMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confPhiBits.negDaughDcazMax.value, twotrackresonanceselection::kNegDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    phiSels.addSelection(confPhiBits.negDaughTpcPion.value, twotrackresonanceselection::kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
    phiSels.addSelection(confPhiBits.negDaughTofPion.value, twotrackresonanceselection::kNegDaughTofPion, limits::kAbsUpperLimit, false, false);
    phiSels.addSelection(confPhiBits.negDaughTpctofPion.value, twotrackresonanceselection::kNegDaughTpctofPion, limits::kAbsUpperLimit, false, false);

    // kstar selections
    kstarSels.setResonanceType(twotrackresonanceselection::TwoTrackResonace::kKstar);
    kstarSels.setMinimalMomentumForTof(confKstarBits.minMomentumForTof.value);
    kstarSels.setPtLimits(confKstarFilters.ptMin.value, confKstarFilters.ptMax.value);
    kstarSels.setEtaLimits(confKstarFilters.etaMin.value, confKstarFilters.etaMax.value);
    kstarSels.setPhiLimits(confKstarFilters.phiMin.value, confKstarFilters.phiMax.value);
    kstarSels.setMassLimits(confKstarFilters.massMin.value, confKstarFilters.massMax.value);
    // kstar positive daughter selections
    kstarSels.addSelection(confKstarBits.posDaughTpcClustersMin.value, twotrackresonanceselection::kPosDauTpcClsMin, limits::kUpperLimit, true, true);
    kstarSels.addSelection(confKstarBits.posDaughDcaxyMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confKstarBits.posDaughDcaxyMax.value, twotrackresonanceselection::kPosDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    kstarSels.addSelection(confKstarBits.posDaughDcazMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confKstarBits.posDaughDcazMax.value, twotrackresonanceselection::kPosDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    kstarSels.addSelection(confKstarBits.posDaughTpcPion.value, twotrackresonanceselection::kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
    kstarSels.addSelection(confKstarBits.posDaughTofPion.value, twotrackresonanceselection::kPosDaughTofPion, limits::kAbsUpperLimit, false, false);
    kstarSels.addSelection(confKstarBits.posDaughTpctofPion.value, twotrackresonanceselection::kPosDaughTpctofPion, limits::kAbsUpperLimit, false, false);
    // kstar negative daughter selections
    kstarSels.addSelection(confKstarBits.negDaughTpcClustersMin.value, twotrackresonanceselection::kNegDauTpcClsMin, limits::kUpperLimit, true, true);
    kstarSels.addSelection(confKstarBits.negDaughDcaxyMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confKstarBits.negDaughDcaxyMax.value, twotrackresonanceselection::kNegDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    kstarSels.addSelection(confKstarBits.negDaughDcazMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confKstarBits.negDaughDcazMax.value, twotrackresonanceselection::kNegDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    kstarSels.addSelection(confKstarBits.negDaughTpcPion.value, twotrackresonanceselection::kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
    kstarSels.addSelection(confKstarBits.negDaughTofPion.value, twotrackresonanceselection::kNegDaughTofPion, limits::kAbsUpperLimit, false, false);
    kstarSels.addSelection(confKstarBits.negDaughTpctofPion.value, twotrackresonanceselection::kNegDaughTpctofPion, limits::kAbsUpperLimit, false, false);

    // antiKstar selections
    antiKstarSels.setResonanceType(twotrackresonanceselection::TwoTrackResonace::kAntiKstar);
    antiKstarSels.setMinimalMomentumForTof(confKstarBits.minMomentumForTof.value);
    antiKstarSels.setPtLimits(confKstarFilters.ptMin.value, confKstarFilters.ptMax.value);
    antiKstarSels.setEtaLimits(confKstarFilters.etaMin.value, confKstarFilters.etaMax.value);
    antiKstarSels.setPhiLimits(confKstarFilters.phiMin.value, confKstarFilters.phiMax.value);
    antiKstarSels.setMassLimits(confKstarFilters.massMin.value, confKstarFilters.massMax.value);
    // antiKstar positive daughter selections
    antiKstarSels.addSelection(confKstarBits.posDaughTpcClustersMin.value, twotrackresonanceselection::kPosDauTpcClsMin, limits::kUpperLimit, true, true);
    antiKstarSels.addSelection(confKstarBits.posDaughDcaxyMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confKstarBits.posDaughDcaxyMax.value, twotrackresonanceselection::kPosDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    antiKstarSels.addSelection(confKstarBits.posDaughDcazMax.name, confResonanceDaughterFilters.posDaughPtMin.value, confResonanceDaughterFilters.posDaughPtMax.value, confKstarBits.posDaughDcazMax.value, twotrackresonanceselection::kPosDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    antiKstarSels.addSelection(confKstarBits.posDaughTpcPion.value, twotrackresonanceselection::kPosDaughTpcPion, limits::kAbsUpperLimit, false, false);
    antiKstarSels.addSelection(confKstarBits.posDaughTofPion.value, twotrackresonanceselection::kPosDaughTofPion, limits::kAbsUpperLimit, false, false);
    antiKstarSels.addSelection(confKstarBits.posDaughTpctofPion.value, twotrackresonanceselection::kPosDaughTpctofPion, limits::kAbsUpperLimit, false, false);
    // antiKstar negative daughter selections
    antiKstarSels.addSelection(confKstarBits.negDaughTpcClustersMin.value, twotrackresonanceselection::kNegDauTpcClsMin, limits::kUpperLimit, true, true);
    antiKstarSels.addSelection(confKstarBits.negDaughDcaxyMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confKstarBits.negDaughDcaxyMax.value, twotrackresonanceselection::kNegDauDcaxyMax, limits::kAbsUpperFunctionLimit, true, true);
    antiKstarSels.addSelection(confKstarBits.negDaughDcazMax.name, confResonanceDaughterFilters.negDaughPtMin.value, confResonanceDaughterFilters.negDaughPtMax.value, confKstarBits.negDaughDcazMax.value, twotrackresonanceselection::kNegDauDcazMax, limits::kAbsUpperFunctionLimit, true, true);
    antiKstarSels.addSelection(confKstarBits.negDaughTpcPion.value, twotrackresonanceselection::kNegDaughTpcPion, limits::kAbsUpperLimit, false, false);
    antiKstarSels.addSelection(confKstarBits.negDaughTofPion.value, twotrackresonanceselection::kNegDaughTofPion, limits::kAbsUpperLimit, false, false);
    antiKstarSels.addSelection(confKstarBits.negDaughTpctofPion.value, twotrackresonanceselection::kNegDaughTpctofPion, limits::kAbsUpperLimit, false, false);
  }

  template <modes::System sys, typename T>
  void fillCollision(T const& col)
  {
    if constexpr (!modes::isSystemSet(sys, modes::System::kNoCentCal)) {
      products.producedCollision(col.posZ(),
                                 col.multNTracksPV(),
                                 col.centFT0M(),
                                 collisionSel.getSphericity(),
                                 collisionSel.getMagneticField());
    }

    if constexpr (modes::isSystemSet(sys, modes::System::kNoCentCal)) {
      products.producedCollision(col.posZ(),
                                 col.multNTracksPV(),
                                 0,
                                 collisionSel.getSphericity(),
                                 collisionSel.getMagneticField());
    }
  }

  template <modes::Mode mode, typename T1>
  void fillTrack(T1 const& track, bool setBitmask = true)
  {
    if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
      products.producedTracks(products.producedCollision.lastIndex(),
                              track.pt() * track.sign(),
                              track.eta(),
                              track.phi());
      if (setBitmask) {
        products.producedTrackMasks(trackSel.getBitmask());
      } else {
        // if we call with setBitmaks set to false, this not a track that passed standard selections but the daugher of some other particle type
        // for these we set the bitmask to 0, so they will never be selected by accident for other further analysis
        products.producedTrackMasks(static_cast<femtodatatypes::TrackMaskType>(0u));
      }

      if constexpr (modes::isModeSet(mode, modes::Mode::kQA)) {
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
      fillTrack<mode>(track);
    }
  }

  template <modes::Mode mode, typename T>
  void fillLambda(T const& v0, int posDaughterIndex, int negDaughterIndex)
  {
    if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
      products.producedLambdas(products.producedCollision.lastIndex(),
                               v0.pt(),
                               v0.eta(),
                               v0.phi(),
                               v0.mLambda(),
                               v0.mAntiLambda(),
                               posDaughterIndex,
                               negDaughterIndex);
      products.producedLambdaMasks(lambdaSel.getBitmask());
    }
    if constexpr (modes::isModeSet(mode, modes::Mode::kQA)) {
      products.producedLambdaExtras(
        v0.dcaV0daughters(),
        v0.x(),
        v0.y(),
        v0.z(),
        v0.v0radius(),
        v0.mK0Short());
    }
    indexMapLambdas.emplace(v0.globalIndex(), products.producedLambdas.lastIndex());
  }

  template <modes::Mode mode, typename V, typename T>
  void fillLambdas(V const& v0s, T const& tracks)
  {
    for (const auto& v0 : v0s) {
      if (v0.pt() < confLambdaFilters.ptMin.value || v0.pt() > confLambdaFilters.ptMax.value ||
          v0.eta() < confLambdaFilters.etaMin.value || v0.eta() > confLambdaFilters.etaMax.value ||
          v0.phi() < confLambdaFilters.phiMin.value || v0.phi() > confLambdaFilters.phiMax.value) {
        continue;
      }
      lambdaSel.applySelections(v0, tracks);
      if (!lambdaSel.checkK0ShortMass(v0) || !lambdaSel.passesMinimalCuts()) {
        continue;
      }
      if (!lambdaSel.checkDaughterPidsLambda() && !lambdaSel.checkDaughterPidsAntiLambda()) {
        continue;
      }
      auto posDaughter = v0.template posTrack_as<T>();
      auto negDaughter = v0.template negTrack_as<T>();
      int64_t posDaughterIndex;
      int64_t negDaughterIndex;

      // get index of positive daughter
      auto resultPosDaughter = utils::getDaughterIndex(posDaughter.globalIndex(), indexMapTracks);
      if (resultPosDaughter) {
        posDaughterIndex = resultPosDaughter.value();
      } else {
        fillTrack<mode>(posDaughter, false);
        posDaughterIndex = products.producedTracks.lastIndex();
      }

      // get index of negative daughter
      auto resultNegDaughter = utils::getDaughterIndex(negDaughter.globalIndex(), indexMapTracks);
      if (resultNegDaughter) {
        negDaughterIndex = resultNegDaughter.value();
      } else {
        fillTrack<mode>(posDaughter, false);
        negDaughterIndex = products.producedTracks.lastIndex();
      }
      fillLambda<mode>(v0, posDaughterIndex, negDaughterIndex);
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fillResonance(T1 const& posDaughter, T1 const& negDaughter, T2& selectionContainer)
  {

    int64_t posDaughterIndex;
    int64_t negDaughterIndex;

    if (!selectionContainer.passesMinimalCuts()) {
      return;
    }

    selectionContainer.reconstructResonance(posDaughter, negDaughter);
    if (!selectionContainer.checkFilters() || !selectionContainer.checkDaughterPids()) {
      return;
    }
    // get index of positive daughter
    auto resultPosDaughter = utils::getDaughterIndex(posDaughter.globalIndex(), indexMapTracks);
    if (resultPosDaughter) {
      posDaughterIndex = resultPosDaughter.value();
    } else {
      fillTrack<mode>(posDaughter);
      posDaughterIndex = products.producedTracks.lastIndex();
    }

    // get index of negative daughter
    auto resultNegDaughter = utils::getDaughterIndex(negDaughter.globalIndex(), indexMapTracks);
    if (resultNegDaughter) {
      negDaughterIndex = resultNegDaughter.value();
    } else {
      fillTrack<mode>(negDaughter);
      negDaughterIndex = products.producedTracks.lastIndex();
    }

    if constexpr (modes::isModeSet(mode, modes::Mode::kANALYSIS)) {
      products.producedResonances(
        products.producedCollision.lastIndex(),
        selectionContainer.getPt(),
        selectionContainer.getEta(),
        selectionContainer.getPhi(),
        selectionContainer.getMass(),
        static_cast<femtodatatypes::TwoTrackResonaceType>(selectionContainer.getType()),
        posDaughterIndex,
        negDaughterIndex);
      products.producedResonanceMasks(selectionContainer.getBitmask());
    }
  }

  template <modes::Mode mode, typename T>
  void fillResonances(T const& tracks)
  {
    Partition<T> positiveTracks = track::signed1Pt > 0 &&
                                  track::pt > confResonanceDaughterFilters.posDaughPtMin && track::pt < confResonanceDaughterFilters.posDaughPtMax &&
                                  track::eta > confResonanceDaughterFilters.posDaughEtaMin && track::eta < confResonanceDaughterFilters.posDaughEtaMax &&
                                  track::phi > confResonanceDaughterFilters.posDaughPhiMin && track::phi < confResonanceDaughterFilters.posDaughPhiMax;
    positiveTracks.bindTable(tracks);
    Partition<T> negativeTracks = track::signed1Pt < 0 &&
                                  track::pt > confResonanceDaughterFilters.negDaughPtMin && track::pt < confResonanceDaughterFilters.negDaughPtMax &&
                                  track::eta > confResonanceDaughterFilters.negDaughEtaMin && track::eta < confResonanceDaughterFilters.negDaughEtaMax &&
                                  track::phi > confResonanceDaughterFilters.negDaughPhiMin && track::phi < confResonanceDaughterFilters.negDaughPhiMax;
    negativeTracks.bindTable(tracks);
    for (auto const& positiveTrack : positiveTracks) {
      for (auto const& negativeTrack : negativeTracks) {
        fillResonance<mode>(positiveTrack, negativeTrack, rhoSels);
        fillResonance<mode>(positiveTrack, negativeTrack, phiSels);
        fillResonance<mode>(positiveTrack, negativeTrack, kstarSels);
        fillResonance<mode>(positiveTrack, negativeTrack, antiKstarSels);
      }
    }
  }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3>
  void processTracks(T1 const& fullCols, T2 const& /* fullBcs*/, T3 const& fullTracks)
  {
    for (const auto& col : fullCols) {
      initFromCcdb(col.template bc_as<T2>());
      collisionSel.setMagneticField(magField);
      auto tracks = fullTracks.sliceBy(perColTracks, col.globalIndex());
      collisionSel.setSphericity(tracks);
      if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
        continue;
      }
      fillCollision<system>(col);
      indexMapTracks.clear();
      fillTracks<mode>(tracks);
    }
  }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void processTracksV0s(T1 const& fullCols, T2 const& /* fullBcs*/, T3 const& fullTracks, T4 const& fullV0s)
  {
    for (const auto& col : fullCols) {
      initFromCcdb(col.template bc_as<T2>());
      collisionSel.setMagneticField(magField);
      auto tracks = fullTracks.sliceBy(perColTracks, col.globalIndex());
      collisionSel.setSphericity(tracks);
      if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
        continue;
      }
      fillCollision<system>(col);
      indexMapTracks.clear();
      fillTracks<mode>(tracks);
      indexMapLambdas.clear();
      auto v0s = fullV0s.sliceBy(perColV0s, col.globalIndex());
      fillLambdas<mode>(v0s, tracks);
    }
  }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3>
  void processTracksResonances(T1 const& fullCols, T2 const& /* fullBcs*/, T3 const& fullTracks)
  {
    for (const auto& col : fullCols) {
      initFromCcdb(col.template bc_as<T2>());
      collisionSel.setMagneticField(magField);
      auto tracks = fullTracks.sliceBy(perColTracks, col.globalIndex());
      collisionSel.setSphericity(tracks);
      if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
        continue;
      }
      fillCollision<system>(col);
      indexMapTracks.clear();
      fillTracks<mode>(tracks);
      fillResonances<mode>(tracks);
    }
  }

  // proccess functions
  // produce tracks for analysis
  void processTracksRun3pp(Filtered<consumeddata::Run3PpCollisions> const& cols,
                           BCsWithTimestamps const& bcs,
                           Filtered<consumeddata::Run3Tracks> const& tracks)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3Tracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(cols, bcs, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksRun3pp, "Provide tracks for Run3 analysis", true);

  // produce tracks for QA
  void proccessQaTracksRun3pp(Filtered<consumeddata::Run3PpCollisions> const& cols,
                              BCsWithTimestamps const& bcs,
                              Filtered<consumeddata::Run3TracksFullPid> const& tracks)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3TracksFullPid, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracks<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(cols, bcs, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, proccessQaTracksRun3pp, "Provide tracks for Run2 with QA", false);

  // produce tracks and v0s for analysis
  void processTracksVzerosRun3pp(Filtered<consumeddata::Run3PpCollisions> const& cols,
                                 BCsWithTimestamps const& bcs,
                                 Filtered<consumeddata::Run3Tracks> const& tracks,
                                 Filtered<consumeddata::Run3PpVzeros> const& v0s)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3Tracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracksV0s<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(cols, bcs, tracksWithItsPid, v0s);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksVzerosRun3pp, "Provide Tracks and V0s for Run3 with QA", false);

  // produce tracks and v0s for QA
  void processQaTracksVzerosRun3pp(Filtered<consumeddata::Run3PpCollisions> const& cols,
                                   BCsWithTimestamps const& bcs,
                                   Filtered<consumeddata::Run3TracksFullPid> const& tracks,
                                   Filtered<consumeddata::Run3PpVzeros> const& v0s)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3TracksFullPid, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracksV0s<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(cols, bcs, tracksWithItsPid, v0s);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processQaTracksVzerosRun3pp, "Provide Tracks and V0s for Run3 with QA", false);

  // produce tracks and resonances for analysis
  void processTracksResonancesRun3pp(Filtered<consumeddata::Run3PpCollisions> const& cols,
                                     BCsWithTimestamps const& bcs,
                                     Filtered<consumeddata::Run3Tracks> const& tracks)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3Tracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracksResonances<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(cols, bcs, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksResonancesRun3pp, "Provide tracks and resonances for Run3 analysis", true);

  // produce tracks and v0s for QA
  void processQaTracksResonancesRun3pp(Filtered<consumeddata::Run3PpCollisions> const& cols,
                                       BCsWithTimestamps const& bcs,
                                       Filtered<consumeddata::Run3TracksFullPid> const& tracks)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3TracksFullPid, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    processTracksResonances<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(cols, bcs, tracksWithItsPid);
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processQaTracksResonancesRun3pp, "Provide Tracks and Resonances for Run3 with QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoUnitedProducer>(cfgc)};
  return workflow;
}
