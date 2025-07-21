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

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"

#include <cmath>
#include <numeric>
#include <limits>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Define TracksFull with necessary PID information.
// For Sigma- -> n + pi-, we need pion PID for the daughter.
// If Sigma+ (p + pi0) or similar proton daughters were considered, pidTPCPr would also be needed.
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSel>;

// Helper functions for Armenteros-Podolanski variables
// These functions assume momMother and momKink are std::array<float, 3> for px, py, pz
float alphaAP_helper(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
{
  std::array<float, 3> momMissing = {momMother[0] - momKink[0], momMother[1] - momKink[1], momMother[2] - momKink[2]};
  float lQlP = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
  float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
  // Guard against division by zero for robustness, though rare for physical momenta
  float denominator = lQlP + lQlN;
  if (std::abs(denominator) < std::numeric_limits<float>::epsilon()) return 0.f;
  return (lQlP - lQlN) / denominator;
}

float qtAP_helper(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
{
  float dp = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
  float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f); // |P_mother|^2
  float p2A = std::inner_product(momKink.begin(), momKink.end(), momKink.begin(), 0.f);     // |P_daughter|^2
  if (p2V0 <= 0) return 0.f; // Guard against non-positive mother momentum squared
  float val = p2A - dp * dp / p2V0;
  if (val < 0) return 0.f; // Should be non-negative for physical momenta, but numerical precision can cause issues
  return std::sqrt(val);
}

struct sigmaminustask {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaMinus{"sigmaminus", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaPi{"cutNSigmaPi", 4, "NSigmaTPCPion"};
  // New configurable cuts for Armenteros-Podolanski and topological variables
  Configurable<float> cutQtAPSigma{"cutQtAPSigma", 0.15f, "Min qT_AP cut for Sigma candidates"};
  Configurable<float> cutAlphaAPSigma{"cutAlphaAPSigma", 0.0f, "Max Alpha_AP cut for Sigma candidates (for Alpha_AP < 0, set to 0.0)"};
  Configurable<float> cutMaxSigmaRadius{"cutMaxSigmaRadius", 35.0f, "Max radius for Sigma decay vertex (cm)"}; // abs(fRadiusSigma) < 35
  Configurable<float> cutMaxDCAtoPVSigma{"cutMaxDCAtoPVSigma", 0.01f, "Max DCA of Sigma mother to primary vertex (cm)"}; // abs(fDCAtoPVSigma) < 0.01
  Configurable<float> cutMinDCAtoPVPiFromSigma{"cutMinDCAtoPVPiFromSigma", 2.0f, "Min DCA of pion from Sigma to primary vertex (cm)"}; // abs(fDCAtoPVPiFromSigma) > 2

  Preslice<aod::KinkCands> mPerCol = aod::track::collisionId;

  void init(InitContext const&)
  {
    // Axes
    const AxisSpec ptAxis{50, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaPiAxis{100, -5, 5, "n#sigma_{#pi}"};
    const AxisSpec sigmaMassAxis{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    const AxisSpec xiMassAxis{100, 1.2, 1.6, "m_{#Xi} (GeV/#it{c}^{2})"};
    const AxisSpec pdgAxis{10001, -5000, 5000, "PDG code"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    // Sigma-minus reconstruction
    rSigmaMinus.add("h2MassSigmaMinusPt", "h2MassSigmaMinusPt", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2SigmaMassVsXiMass", "h2SigmaMassVsXiMass", {HistType::kTH2F, {xiMassAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2NSigmaPiPt", "h2NSigmaPiPt", {HistType::kTH2F, {ptAxis, nSigmaPiAxis}});

    if (doprocessMC) {
      // Add MC histograms if needed
      rSigmaMinus.add("h2MassPtMCRec", "h2MassPtMCRec", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
      rSigmaMinus.add("h2MassPtMCGen", "h2MassPtMCGen", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    }
  }

  void processData(CollisionsFull::iterator const& collision, aod::KinkCands const& KinkCands, TracksFull const&)
  {
    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
    for (const auto& kinkCand : KinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<TracksFull>();

      // Apply TPC nSigma cut for the pion daughter
      if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
        continue;
      }

      // Calculate Armenteros-Podolanski variables and topological properties
      auto sigmaMom = std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
      auto kinkDauMom = std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};

      float currentSigmaAlphaAP = alphaAP_helper(sigmaMom, kinkDauMom);
      float currentSigmaQtAP = qtAP_helper(sigmaMom, kinkDauMom);
      float currentSigmaRadius = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx()); // 2D decay radius
      float currentDCAtoPVSigma = std::abs(kinkCand.dcaMothPv());     // Absolute DCA of mother to PV
      float currentDCAtoPVPiFromSigma = std::abs(kinkCand.dcaDaugPv()); // Absolute DCA of daughter to PV

      // Apply the requested Armenteros-Podolanski and topological cuts
      if (currentSigmaQtAP < cutQtAPSigma) { // fQtAPSigma > 0.15
        continue;
      }
      if (currentSigmaAlphaAP >= cutAlphaAPSigma) { // fAlphaAPSigma < 0 (using cutAlphaAPSigma = 0.0)
        continue;
      }
      if (currentSigmaRadius >= cutMaxSigmaRadius) { // abs(fRadiusSigma) < 35 (using cutMaxSigmaRadius = 35.0)
        continue;
      }
      if (currentDCAtoPVSigma >= cutMaxDCAtoPVSigma) { // abs(fDCAtoPVSigma) < 0.01 (using cutMaxDCAtoPVSigma = 0.01)
        continue;
      }
      if (currentDCAtoPVPiFromSigma <= cutMinDCAtoPVPiFromSigma) { // abs(fDCAtoPVPiFromSigma) > 2 (using cutMinDCAtoPVPiFromSigma = 2.0)
        continue;
      }

      // If all cuts pass, fill histograms
      rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2NSigmaPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());
    }
  }
  PROCESS_SWITCH(sigmaminustask, processData, "Data processing", true);

  void processMC(CollisionsFullMC const& collisions, aod::KinkCands const& KinkCands, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, TracksFull const&)
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

        // Check for same sign mother and daughter, common for Sigma- decay (pi- and Sigma-)
        // This check was present in your original code.
        if (dauTrack.sign() != mothTrack.sign()) {
          // LOG(info) << "Skipping kink candidate with opposite sign daughter and mother: " << kinkCand.globalIndex();
          continue;
        }

        // Apply TPC nSigma cut for the pion daughter
        if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
          continue;
        }

        // Calculate Armenteros-Podolanski variables and topological properties
        auto sigmaMom = std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
        auto kinkDauMom = std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};

        float currentSigmaAlphaAP = alphaAP_helper(sigmaMom, kinkDauMom);
        float currentSigmaQtAP = qtAP_helper(sigmaMom, kinkDauMom);
        float currentSigmaRadius = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx());
        float currentDCAtoPVSigma = std::abs(kinkCand.dcaMothPv());
        float currentDCAtoPVPiFromSigma = std::abs(kinkCand.dcaDaugPv());

        // Apply the requested Armenteros-Podolanski and topological cuts (same as for Data)
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

        // If all cuts pass, fill histograms (reconstructed side)
        rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2NSigmaPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());

        // Do MC association (only for candidates that passed all reconstruction cuts)
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
            // Check PDG codes for true Sigma- and pion daughter
            if (std::abs(mcTrackSigma.pdgCode()) != 3112 || std::abs(mcTrackPiDau.pdgCode()) != 211) {
              continue;
            }
            rSigmaMinus.fill(HIST("h2MassPtMCRec"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
          }
        }
      }
    }
    // MC truth filling for generated particles (this loop is separate and does not use reconstructed cuts)
    for (const auto& mcPart : particlesMC) {
      if (std::abs(mcPart.pdgCode()) != 3112 || std::abs(mcPart.y()) > 0.5) {
        continue;
      }
      if (!mcPart.has_daughters()) {
        continue; // Skip if no daughters
      }
      bool hasPiDaughter = false;
      for (const auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) == 211) { // Check for pi daughter (pi- for Sigma-, pi+ for Sigma+)
          hasPiDaughter = true;
          break;
        }
      }
      if (!hasPiDaughter) {
        continue; // Skip if no pi daughter found
      }
      float mcMass = std::sqrt(mcPart.e() * mcPart.e() - mcPart.p() * mcPart.p());
      int sigmaSign = mcPart.pdgCode() > 0 ? 1 : -1; // Determine the sign (1 for Sigma-, -1 for Anti-Sigma+)
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
