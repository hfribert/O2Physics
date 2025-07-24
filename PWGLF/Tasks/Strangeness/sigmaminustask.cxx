// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   sigmaminustask.cxx
/// \brief Example of a simple task for the analysis of the Sigma-minus
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>
/// \author Henrik Fribert <henrik.fribert@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include <cmath>
#include <numeric>
#include <array>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSel>;

struct sigmaminustask {
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaMinus{"sigmaminus", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rDebug{"debug", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaPi{"cutNSigmaPi", 4, "NSigmaTPCPion"};
  Configurable<float> cutQtAPSigma{"cutQtAPSigma", 0.15f, "Min qT_AP cut for Sigma candidates"};
  Configurable<float> cutAlphaAPSigma{"cutAlphaAPSigma", 0.0f, "Max Alpha_AP cut for Sigma candidates (for Alpha_AP < 0, set to 0.0)"};
  Configurable<float> cutMaxSigmaRadius{"cutMaxSigmaRadius", 35.0f, "Max radius for Sigma decay vertex (cm)"};
  Configurable<float> cutMaxDCAtoPVSigma{"cutMaxDCAtoPVSigma", 0.01f, "Max DCA of Sigma mother to primary vertex (cm)"};
  Configurable<float> cutMinDCAtoPVPiFromSigma{"cutMinDCAtoPVPiFromSigma", 2.0f, "Min DCA of pion from Sigma to primary vertex (cm)"};

  std::array<float, 3> forbiddenRadii = {19.6213f, 24.5597f, 34.388f};
  float forbiddenRadiusTolerance = 1.0f;

  Preslice<aod::KinkCands> mPerCol = aod::track::collisionId;

  float alphaAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    std::array<float, 3> momMissing = {momMother[0] - momKink[0], momMother[1] - momKink[1], momMother[2] - momKink[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    return (lQlP - lQlN) / (lQlP + lQlN);
  }

  float qtAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    float dp = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momKink.begin(), momKink.end(), momKink.begin(), 0.f);
    return std::sqrt(p2A - dp * dp / p2V0);
  }

  float calculateKinkAngle(const std::array<float, 3>& momMother, const std::array<float, 3>& momDaughter)
  {
    float pMoth = std::sqrt(momMother[0] * momMother[0] + momMother[1] * momMother[1] + momMother[2] * momMother[2]);
    float pDaug = std::sqrt(momDaughter[0] * momDaughter[0] + momDaughter[1] * momDaughter[1] + momDaughter[2] * momDaughter[2]);
    float dotProduct = std::inner_product(momMother.begin(), momMother.end(), momDaughter.begin(), 0.f);

    if (pMoth <= 0 || pDaug <= 0) {
      return 0.f;
    }

    float cosAngle = dotProduct / (pMoth * pDaug);
    if (cosAngle > 1.0f) cosAngle = 1.0f;
    if (cosAngle < -1.0f) cosAngle = -1.0f;

    return std::acos(cosAngle);
  }

  float calculateDecayLength(float xDecVtx, float yDecVtx, float zDecVtx)
  {
    return std::sqrt(xDecVtx * xDecVtx + yDecVtx * yDecVtx + zDecVtx * zDecVtx);
  }

  void init(InitContext const&)
  {
    const AxisSpec ptAxis{50, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec pAxis{100, 0, 10, "#it{p} (GeV/#it{c})"};
    const AxisSpec nSigmaPiAxis{100, -5, 5, "n#sigma_{#pi}"};
    const AxisSpec sigmaMassAxis{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    const AxisSpec xiMassAxis{100, 1.2, 1.6, "m_{#Xi} (GeV/#it{c}^{2})"};
    const AxisSpec pdgAxis{10001, -5000, 5000, "PDG code"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};
    const AxisSpec kinkAngleAxis{100, 0, M_PI, "Kink Angle (rad)"};
    const AxisSpec decayLengthAxis{100, 0, 100, "Decay Length (cm)"};
    const AxisSpec decayRadiusAxis{100, 0, 50, "Decay Radius (cm)"};
    const AxisSpec dEdxAxis{100, 0, 200, "TPC dE/dx (arb. units)"};
    const AxisSpec phiAxis{100, 0, 2 * M_PI, "#phi (rad)"};
    const AxisSpec etaAxis{100, -1.5, 1.5, "#eta"};
    const AxisSpec itsClusterSizeAxis{100, 0, 30, "ITS Avg. Cluster Size (Corrected)"};
    const AxisSpec decayVertexXYAxis{400, -50, 50, "Decay Vertex Position (cm)"};

    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rSigmaMinus.add("h2MassSigmaMinusPt", "h2MassSigmaMinusPt", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2SigmaMassVsXiMass", "h2SigmaMassVsXiMass", {HistType::kTH2F, {xiMassAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2NSigmaPiPt", "h2NSigmaPiPt", {HistType::kTH2F, {ptAxis, nSigmaPiAxis}});
    rSigmaMinus.add("h2MassSigmaMinusPtBeforeForbiddenRadii", "h2MassSigmaMinusPtBeforeForbiddenRadii", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});

    // Histograms for decay properties before cuts
    rDebug.add("hDecayLengthBeforeCuts", "Decay Length (Before cuts)", {HistType::kTH1F, {decayLengthAxis}});
    rDebug.add("hDecayRadiusBeforeCuts", "Decay Radius (Transverse, Before cuts)", {HistType::kTH1F, {decayRadiusAxis}});
    rDebug.add("hDecayVertexXYBeforeCuts", "Decay Vertex X-Y Position (Before cuts)", {HistType::kTH2F, {decayVertexXYAxis, decayVertexXYAxis}});

    // Debug histograms filled after cuts
    rDebug.add("hPtMother", "pT Mother (after cuts)", {HistType::kTH1F, {ptAxis}});
    rDebug.add("hPtDaughter", "pT Daughter (after cuts)", {HistType::kTH1F, {ptAxis}});
    rDebug.add("hPtDaughterNeg", "pT Daughter (Negative sign, after cuts)", {HistType::kTH1F, {ptAxis}});
    rDebug.add("hPtDaughterPos", "pT Daughter (Positive sign, after cuts)", {HistType::kTH1F, {ptAxis}});
    rDebug.add("hKinkAngle", "Kink Angle (after cuts)", {HistType::kTH1F, {kinkAngleAxis}});
    rDebug.add("hDecayLengthAfterCuts", "Decay Length (after cuts)", {HistType::kTH1F, {decayLengthAxis}});
    rDebug.add("hDecayRadiusAfterCuts", "Decay Radius (Transverse, after cuts)", {HistType::kTH1F, {decayRadiusAxis}});
    rDebug.add("hdEdxDaughter", "TPC dE/dx Daughter (after cuts)", {HistType::kTH2F, {pAxis, dEdxAxis}});
    rDebug.add("hdEdxDaughterNeg", "TPC dE/dx Daughter (Negative sign, after cuts)", {HistType::kTH2F, {pAxis, dEdxAxis}});
    rDebug.add("hdEdxDaughterPos", "TPC dE/dx Daughter (Positive sign, after cuts)", {HistType::kTH2F, {pAxis, dEdxAxis}});
    rDebug.add("hPhiMother", "Phi Mother (after cuts)", {HistType::kTH1F, {phiAxis}});
    rDebug.add("hPhiDaughter", "Phi Daughter (after cuts)", {HistType::kTH1F, {phiAxis}});
    rDebug.add("hPhiDaughterNeg", "Phi Daughter (Negative sign, after cuts)", {HistType::kTH1F, {phiAxis}});
    rDebug.add("hPhiDaughterPos", "Phi Daughter (Positive sign, after cuts)", {HistType::kTH1F, {phiAxis}});
    rDebug.add("hEtaMother", "Eta Mother (after cuts)", {HistType::kTH1F, {etaAxis}});
    rDebug.add("hEtaDaughter", "Eta Daughter (after cuts)", {HistType::kTH1F, {etaAxis}});
    rDebug.add("hEtaDaughterNeg", "Eta Daughter (Negative sign, after cuts)", {HistType::kTH1F, {etaAxis}});
    rDebug.add("hEtaDaughterPos", "Eta Daughter (Positive sign, after cuts)", {HistType::kTH1F, {etaAxis}});

    rDebug.add("hITSClusterSizeMother", "ITS Mother Avg. Cluster Size (Corrected, after cuts)", {HistType::kTH2F, {pAxis, itsClusterSizeAxis}});
    rDebug.add("hITSClusterSizeMotherNeg", "ITS Mother Avg. Cluster Size (Corrected, Neg, after cuts)", {HistType::kTH2F, {pAxis, itsClusterSizeAxis}});
    rDebug.add("hITSClusterSizeMotherPos", "ITS Mother Avg. Cluster Size (Corrected, Pos, after cuts)", {HistType::kTH2F, {pAxis, itsClusterSizeAxis}});

    rDebug.add("hITSClusterSizeDaughter", "ITS Daughter Avg. Cluster Size (Corrected, after cuts)", {HistType::kTH2F, {pAxis, itsClusterSizeAxis}});
    rDebug.add("hITSClusterSizeDaughterNeg", "ITS Daughter Avg. Cluster Size (Corrected, Neg, after cuts)", {HistType::kTH2F, {pAxis, itsClusterSizeAxis}});
    rDebug.add("hITSClusterSizeDaughterPos", "ITS Daughter Avg. Cluster Size (Corrected, Pos, after cuts)", {HistType::kTH2F, {pAxis, itsClusterSizeAxis}});

    rDebug.add("hDecayVertexXYAfterCuts", "Decay Vertex X-Y Position (after cuts)", {HistType::kTH2F, {decayVertexXYAxis, decayVertexXYAxis}});

    if (doprocessMC) {
      rSigmaMinus.add("h2MassPtMCRec", "h2MassPtMCRec", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
      rSigmaMinus.add("h2MassPtMCGen", "h2MassPtMCGen", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    }
  }

  void processData(CollisionsFull::iterator const& collision, aod::KinkCands const& KinkCands, TracksFull const& TracksFullData)
  {
    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    for (const auto& kinkCand : KinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
      auto mothTrack = kinkCand.trackMoth_as<TracksFull>();

      auto sigmaMom = std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
      auto kinkDauMom = std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};

      float currentKinkAngle = calculateKinkAngle(sigmaMom, kinkDauMom);
      float currentDecayLength = calculateDecayLength(kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx());
      float currentSigmaRadius = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx());

      // ITS Cluster Size Calculation (for Mother Track - innermost 3 layers)
      uint32_t motherClusterSizesRaw = mothTrack.itsClusterSizes();
      float sumSizesMother = 0.f;
      int layersWithClustersMother = 0;
      for (int i = 0; i < 3; ++i) {
        int layerSize = (motherClusterSizesRaw >> (i * 4)) & 0xf;
        if (layerSize > 0) {
          sumSizesMother += layerSize;
          layersWithClustersMother++;
        }
      }

      float averageClusterSizeMother = 0.f;
      if (layersWithClustersMother > 0) {
        averageClusterSizeMother = sumSizesMother / layersWithClustersMother;
      }
      
      float motherEta = mothTrack.eta();
      float motherPolarAngle = 2 * std::atan(std::exp(-std::abs(motherEta)));
      float correctedAverageClusterSizeMother = averageClusterSizeMother;
      if (motherPolarAngle > 0) {
        correctedAverageClusterSizeMother = averageClusterSizeMother * std::sin(motherPolarAngle);
      }

      // ITS Cluster Size Calculation (for Daughter Track - layers 3, 4, 5, 6)
      uint32_t daughterClusterSizesRaw = dauTrack.itsClusterSizes();
      float sumSizesDaughter = 0.f;
      int layersWithClustersDaughter = 0;
      for (int i = 3; i < 7; ++i) {
        int layerSize = (daughterClusterSizesRaw >> (i * 4)) & 0xf;
        if (layerSize > 0) {
          sumSizesDaughter += layerSize;
          layersWithClustersDaughter++;
        }
      }
      
      float averageClusterSizeDaughter = 0.f;
      if (layersWithClustersDaughter > 0) {
        averageClusterSizeDaughter = sumSizesDaughter / layersWithClustersDaughter;
      }

      float daughterEta = dauTrack.eta();
      float daughterPolarAngle = 2 * std::atan(std::exp(-std::abs(daughterEta)));
      float correctedAverageClusterSizeDaughter = averageClusterSizeDaughter;
      if (daughterPolarAngle > 0) {
        correctedAverageClusterSizeDaughter = averageClusterSizeDaughter * std::sin(daughterPolarAngle);
      }

      if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
        continue;
      }

      float currentSigmaAlphaAP = alphaAP(sigmaMom, kinkDauMom);
      float currentSigmaQtAP = qtAP(sigmaMom, kinkDauMom);
      float currentDCAtoPVSigma = std::abs(kinkCand.dcaMothPv());
      float currentDCAtoPVPiFromSigma = std::abs(kinkCand.dcaDaugPv());

      if (currentSigmaQtAP < cutQtAPSigma) {
        continue;
      }
      if (currentSigmaAlphaAP >= cutAlphaAPSigma) {
        continue;
      }
      if (currentSigmaRadius >= cutMaxSigmaRadius) {
        continue;
      }   
      if (currentDCAtoPVSigma >= cutMaxDCAtoPVSigma) {
        continue;
      }
      if (currentDCAtoPVPiFromSigma <= cutMinDCAtoPVPiFromSigma) {
        continue;
      }

      rDebug.fill(HIST("hDecayLengthBeforeCuts"), currentDecayLength);
      rDebug.fill(HIST("hDecayRadiusBeforeCuts"), currentSigmaRadius);
      rDebug.fill(HIST("hDecayVertexXYBeforeCuts"), kinkCand.xDecVtx(), kinkCand.yDecVtx());
      rSigmaMinus.fill(HIST("h2MassSigmaMinusPtBeforeForbiddenRadii"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());

      // Apply the forbidden radii cut
      bool isForbiddenRadius = false;
      for (float rForbidden : forbiddenRadii) {
        if (std::abs(currentSigmaRadius - rForbidden) < forbiddenRadiusTolerance) {
          isForbiddenRadius = true;
          break;
        }
      }
      if (isForbiddenRadius) {
        continue;
      }

      // Fill debug histograms after cuts
      rDebug.fill(HIST("hPtMother"), mothTrack.pt() * mothTrack.sign());
      rDebug.fill(HIST("hPtDaughter"), dauTrack.pt() * dauTrack.sign());
      if (dauTrack.sign() < 0) {
        rDebug.fill(HIST("hPtDaughterNeg"), dauTrack.pt());
      } else {
        rDebug.fill(HIST("hPtDaughterPos"), dauTrack.pt());
      }
      rDebug.fill(HIST("hKinkAngle"), currentKinkAngle);
      rDebug.fill(HIST("hDecayLengthAfterCuts"), currentDecayLength);
      rDebug.fill(HIST("hDecayRadiusAfterCuts"), currentSigmaRadius);
      rDebug.fill(HIST("hdEdxDaughter"), dauTrack.p(), dauTrack.tpcSignal());
      if (dauTrack.sign() < 0) {
        rDebug.fill(HIST("hdEdxDaughterNeg"), dauTrack.p(), dauTrack.tpcSignal());
      } else {
        rDebug.fill(HIST("hdEdxDaughterPos"), dauTrack.p(), dauTrack.tpcSignal());
      }
      rDebug.fill(HIST("hPhiMother"), mothTrack.phi());
      rDebug.fill(HIST("hPhiDaughter"), dauTrack.phi());
      if (dauTrack.sign() < 0) {
        rDebug.fill(HIST("hPhiDaughterNeg"), dauTrack.phi());
      } else {
        rDebug.fill(HIST("hPhiDaughterPos"), dauTrack.phi());
      }
      rDebug.fill(HIST("hEtaMother"), mothTrack.eta());
      rDebug.fill(HIST("hEtaDaughter"), dauTrack.eta());
      if (dauTrack.sign() < 0) {
        rDebug.fill(HIST("hEtaDaughterNeg"), dauTrack.eta());
      } else {
        rDebug.fill(HIST("hEtaDaughterPos"), dauTrack.eta());
      }

      // Fill ITS Cluster Size histograms only if enough hits for a meaningful average
      if (layersWithClustersMother >= 3) {
        rDebug.fill(HIST("hITSClusterSizeMother"), mothTrack.p(), correctedAverageClusterSizeMother);
        if (mothTrack.sign() < 0) {
          rDebug.fill(HIST("hITSClusterSizeMotherNeg"), mothTrack.p(), correctedAverageClusterSizeMother);
        } else {
          rDebug.fill(HIST("hITSClusterSizeMotherPos"), mothTrack.p(), correctedAverageClusterSizeMother);
        }
      }
      if (layersWithClustersDaughter >= 4) {
        rDebug.fill(HIST("hITSClusterSizeDaughter"), dauTrack.p(), correctedAverageClusterSizeDaughter);
        if (dauTrack.sign() < 0) {
          rDebug.fill(HIST("hITSClusterSizeDaughterNeg"), dauTrack.p(), correctedAverageClusterSizeDaughter);
        } else {
          rDebug.fill(HIST("hITSClusterSizeDaughterPos"), dauTrack.p(), correctedAverageClusterSizeDaughter);
        }
      }
      
      rDebug.fill(HIST("hDecayVertexXYAfterCuts"), kinkCand.xDecVtx(), kinkCand.yDecVtx());

      rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2NSigmaPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());
    }
  }
  PROCESS_SWITCH(sigmaminustask, processData, "Data processing", true);

  void processMC(CollisionsFullMC const& collisions, aod::KinkCands const& KinkCands, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, TracksFull const& TracksFullMC)
  {
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }

      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      auto kinkCandPerColl = KinkCands.sliceBy(mPerCol, collision.globalIndex());
      for (const auto& kinkCand : kinkCandPerColl) {
        auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
        auto mothTrack = kinkCand.trackMoth_as<TracksFull>();

        auto sigmaMom = std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
        auto kinkDauMom = std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};

        float currentKinkAngle = calculateKinkAngle(sigmaMom, kinkDauMom);
        float currentDecayLength = calculateDecayLength(kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx());
        float currentSigmaRadius = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx());

        // ITS Cluster Size Calculation (for Mother Track - innermost 3 layers)
        uint32_t motherClusterSizesRaw = mothTrack.itsClusterSizes();
        float sumSizesMother = 0.f;
        int layersWithClustersMother = 0;
        for (int i = 0; i < 3; ++i) {
          int layerSize = (motherClusterSizesRaw >> (i * 4)) & 0xf;
          if (layerSize > 0) {
            sumSizesMother += layerSize;
            layersWithClustersMother++;
          }
        }

        float averageClusterSizeMother = 0.f;
        if (layersWithClustersMother > 0) {
            averageClusterSizeMother = sumSizesMother / layersWithClustersMother;
        }

        float motherEta = mothTrack.eta();
        float motherPolarAngle = 2 * std::atan(std::exp(-std::abs(motherEta)));
        float correctedAverageClusterSizeMother = averageClusterSizeMother;
        if (motherPolarAngle > 0) {
          correctedAverageClusterSizeMother = averageClusterSizeMother * std::sin(motherPolarAngle);
        }

        // ITS Cluster Size Calculation (for Daughter Track - layers 3, 4, 5, 6)
        uint32_t daughterClusterSizesRaw = dauTrack.itsClusterSizes();
        float sumSizesDaughter = 0.f;
        int layersWithClustersDaughter = 0;
        for (int i = 3; i < 7; ++i) {
          int layerSize = (daughterClusterSizesRaw >> (i * 4)) & 0xf;
          if (layerSize > 0) {
            sumSizesDaughter += layerSize;
            layersWithClustersDaughter++;
          }
        }

        float averageClusterSizeDaughter = 0.f;
        if (layersWithClustersDaughter > 0) {
            averageClusterSizeDaughter = sumSizesDaughter / layersWithClustersDaughter;
        }

        float daughterEta = dauTrack.eta();
        float daughterPolarAngle = 2 * std::atan(std::exp(-std::abs(daughterEta)));
        float correctedAverageClusterSizeDaughter = averageClusterSizeDaughter;
        if (daughterPolarAngle > 0) {
          correctedAverageClusterSizeDaughter = averageClusterSizeDaughter * std::sin(daughterPolarAngle);
        }

        if (dauTrack.sign() != mothTrack.sign()) {
          continue;
        }

        if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
          continue;
        }

        float currentSigmaAlphaAP = alphaAP(sigmaMom, kinkDauMom);
        float currentSigmaQtAP = qtAP(sigmaMom, kinkDauMom);
        float currentDCAtoPVSigma = std::abs(kinkCand.dcaMothPv());
        float currentDCAtoPVPiFromSigma = std::abs(kinkCand.dcaDaugPv());

        if (currentSigmaQtAP < cutQtAPSigma) {
          continue;
        }
        if (currentSigmaAlphaAP >= cutAlphaAPSigma) {
          continue;
        }
        if (currentSigmaRadius >= cutMaxSigmaRadius) {
          continue;
        }
        if (currentDCAtoPVSigma >= cutMaxDCAtoPVSigma) {
          continue;
        }
        if (currentDCAtoPVPiFromSigma <= cutMinDCAtoPVPiFromSigma) {
          continue;
        }

        // Fill decay properties (after all cuts except forbidden radii)
        rDebug.fill(HIST("hDecayLengthBeforeCuts"), currentDecayLength);
        rDebug.fill(HIST("hDecayRadiusBeforeCuts"), currentSigmaRadius);
        rDebug.fill(HIST("hDecayVertexXYBeforeCuts"), kinkCand.xDecVtx(), kinkCand.yDecVtx());
        rSigmaMinus.fill(HIST("h2MassSigmaMinusPtBeforeForbiddenRadii"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());

        // Apply the forbidden radii cut
        bool isForbiddenRadius = false;
        for (float rForbidden : forbiddenRadii) {
          if (std::abs(currentSigmaRadius - rForbidden) < forbiddenRadiusTolerance) {
            isForbiddenRadius = true;
            break;
          }
        }
        if (isForbiddenRadius) {
          continue;
        }

        // Fill debug histograms after *all* cuts
        rDebug.fill(HIST("hPtMother"), mothTrack.pt() * mothTrack.sign());
        rDebug.fill(HIST("hPtDaughter"), dauTrack.pt() * dauTrack.sign());
        if (dauTrack.sign() < 0) {
          rDebug.fill(HIST("hPtDaughterNeg"), dauTrack.pt());
        } else {
          rDebug.fill(HIST("hPtDaughterPos"), dauTrack.pt());
        }
        rDebug.fill(HIST("hKinkAngle"), currentKinkAngle);
        rDebug.fill(HIST("hDecayLengthAfterCuts"), currentDecayLength);
        rDebug.fill(HIST("hDecayRadiusAfterCuts"), currentSigmaRadius);
        rDebug.fill(HIST("hdEdxDaughter"), dauTrack.p(), dauTrack.tpcSignal());
        if (dauTrack.sign() < 0) {
          rDebug.fill(HIST("hdEdxDaughterNeg"), dauTrack.p(), dauTrack.tpcSignal());
        } else {
          rDebug.fill(HIST("hdEdxDaughterPos"), dauTrack.p(), dauTrack.tpcSignal());
        }
        rDebug.fill(HIST("hPhiMother"), mothTrack.phi());
        rDebug.fill(HIST("hPhiDaughter"), dauTrack.phi());
        if (dauTrack.sign() < 0) {
          rDebug.fill(HIST("hPhiDaughterNeg"), dauTrack.phi());
        } else {
          rDebug.fill(HIST("hPhiDaughterPos"), dauTrack.phi());
        }
        rDebug.fill(HIST("hEtaMother"), mothTrack.eta());
        rDebug.fill(HIST("hEtaDaughter"), dauTrack.eta());
        if (dauTrack.sign() < 0) {
          rDebug.fill(HIST("hEtaDaughterNeg"), dauTrack.eta());
        } else {
          rDebug.fill(HIST("hEtaDaughterPos"), dauTrack.eta());
        }

        // Fill ITS Cluster Size histograms only if enough hits for a meaningful average
        if (layersWithClustersMother >= 3) {
            rDebug.fill(HIST("hITSClusterSizeMother"), mothTrack.p(), correctedAverageClusterSizeMother);
            if (mothTrack.sign() < 0) {
            rDebug.fill(HIST("hITSClusterSizeMotherNeg"), mothTrack.p(), correctedAverageClusterSizeMother);
            } else {
            rDebug.fill(HIST("hITSClusterSizeMotherPos"), mothTrack.p(), correctedAverageClusterSizeMother);
            }
        }
        if (layersWithClustersDaughter >= 4) {
            rDebug.fill(HIST("hITSClusterSizeDaughter"), dauTrack.p(), correctedAverageClusterSizeDaughter);
            if (dauTrack.sign() < 0) {
            rDebug.fill(HIST("hITSClusterSizeDaughterNeg"), dauTrack.p(), correctedAverageClusterSizeDaughter);
            } else {
            rDebug.fill(HIST("hITSClusterSizeDaughterPos"), dauTrack.p(), correctedAverageClusterSizeDaughter);
            }
        }

        rDebug.fill(HIST("hDecayVertexXYAfterCuts"), kinkCand.xDecVtx(), kinkCand.yDecVtx());

        rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2NSigmaPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());

        auto mcLabSigma = trackLabelsMC.rawIteratorAt(mothTrack.globalIndex());
        auto mcLabPiDau = trackLabelsMC.rawIteratorAt(dauTrack.globalIndex());
        if (mcLabSigma.has_mcParticle() && mcLabPiDau.has_mcParticle()) {
          auto mcTrackSigma = mcLabSigma.mcParticle_as<aod::McParticles>();
          auto mcTrackPiDau = mcLabPiDau.mcParticle_as<aod::McParticles>();
          if (!mcTrackPiDau.has_mothers()) {
            continue;
          }
          for (auto& piMother : mcTrackPiDau.mothers_as<aod::McParticles>()) {
            if (piMother.globalIndex() != mcTrackSigma.globalIndex()) {
              continue;
            }
            if (std::abs(mcTrackSigma.pdgCode()) != 3112 || std::abs(mcTrackPiDau.pdgCode()) != 211) {
              continue;
            }
            rSigmaMinus.fill(HIST("h2MassPtMCRec"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
          }
        }
      }
    }
    for (const auto& mcPart : particlesMC) {
      if (std::abs(mcPart.pdgCode()) != 3112 || std::abs(mcPart.y()) > 0.5) {
        continue;
      }
      if (!mcPart.has_daughters()) {
        continue;
      }
      bool hasPiDaughter = false;
      for (const auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) == 211) {
          hasPiDaughter = true;
          break;
        }
      }
      if (!hasPiDaughter) {
        continue;
      }
      float mcMass = std::sqrt(mcPart.e() * mcPart.e() - mcPart.p() * mcPart.p());
      int sigmaSign = mcPart.pdgCode() > 0 ? 1 : -1;
      rSigmaMinus.fill(HIST("h2MassPtMCGen"), sigmaSign * mcPart.pt(), mcMass);
    }
  }
  PROCESS_SWITCH(sigmaminustask, processMC, "MC processing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<sigmaminustask>(cfgc)};
}
