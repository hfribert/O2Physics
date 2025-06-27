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
///
/// \file onTheFlyTrackerPid.cxx
///
/// \brief This task produces the PID information that can be obtained from the tracker layers (i.e. ToT and possibly cluster size).
///        So far only ToT implemented. It currently contemplates 5 particle types: electrons, muons, pions, kaons and
///        protons (deuterons, tritons and helium-3 are available if added to the event generator).
///
/// \author Berkin Ulukutlu TUM
/// \author Henrik Fribert TUM
/// \author Nicolò Jacazio Università del Piemonte Orientale
/// \since  May 22, 2025
///

#include <utility>
#include <map>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>
#include <array>
#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TString.h>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TString.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "DetectorsVertexing/HelixHelper.h"
#include "TableHelper.h"
#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"

using namespace o2;
using namespace o2::framework;

namespace TrackerConstants {
static constexpr int kMaxBarrelLayers = 11;
static constexpr std::array<float, kMaxBarrelLayers> kTrackerRadii = {0.5f, 1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f}; // Radii in cm
static constexpr int kEtaBins = 50;
static constexpr float kEtaMin = -2.5;
static constexpr float kEtaMax = 2.5;
static constexpr int kPtBins = 500;
static constexpr float kPtMin = 0.0;
static constexpr float kPtMax = 10.0;
}

class ToTLUT {
public:
  ToTLUT() = default;
  ~ToTLUT() = default;

  bool load(int pdg, const std::string& filename) {
    TFile* f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) {
      std::cerr << "ERROR: Failed to open LUT file: " << filename << std::endl;
      return false;
    }

    bool success = true;
    // Iterate through all possible layers (0 to kMaxBarrelLayers - 1)
    for (int layer = 0; layer < TrackerConstants::kMaxBarrelLayers; ++layer) {
      for (int etaBin = 0; etaBin < TrackerConstants::kEtaBins; ++etaBin) {
        float etaMin = TrackerConstants::kEtaMin + etaBin * (TrackerConstants::kEtaMax - TrackerConstants::kEtaMin) / TrackerConstants::kEtaBins;
        float etaMax = etaMin + (TrackerConstants::kEtaMax - TrackerConstants::kEtaMin) / TrackerConstants::kEtaBins;

        for (int ptBin = 0; ptBin < TrackerConstants::kPtBins; ++ptBin) {
          float ptMin = TrackerConstants::kPtMin + ptBin * (TrackerConstants::kPtMax - TrackerConstants::kPtMin) / TrackerConstants::kPtBins;
          float ptMax = ptMin + (TrackerConstants::kPtMax - TrackerConstants::kPtMin) / TrackerConstants::kPtBins;

          TString histName = Form("tot_%d_barrel%d_eta%.2f-%.2f_pt%.2f-%.2f", pdg, layer, etaMin, etaMax, ptMin, ptMax);
          TH1F* hist = dynamic_cast<TH1F*>(f->Get(histName));
          if (hist) {
            mLUT[pdg][{layer, etaBin, ptBin}] = *hist;
          }
          if (hist) delete hist;
        }
      }
    }

    f->Close();
    delete f;
    return success;
  }

  TH1F* getHistogram(int pdg, int layer, int etaBin, int ptBin) {
    auto pdgIt = mLUT.find(pdg);
    if (pdgIt == mLUT.end()) return nullptr;

    auto key = std::make_tuple(layer, etaBin, ptBin);
    auto histIt = pdgIt->second.find(key);
    return (histIt != pdgIt->second.end()) ? &histIt->second : nullptr;
  }

  // Get Most Probable Value (MPV) from the ToT distributions
  float getMPV(int pdg, int layer, int etaBin, int ptBin) {
    TH1F* hist = getHistogram(pdg, layer, etaBin, ptBin);
    if (hist && hist->GetEntries() > 0) {
      // For now returning maximum bin position, compared to fitting Landau's (could be done separately and added to the LUTs in the future)
      return hist->GetBinLowEdge(hist->GetMaximumBin());
    }
    return -1.f;
  }

  // Calculate truncated standard deviation using a cumulative fraction
  float getTruncatedStdDev(int pdg, int layer, int etaBin, int ptBin, float fraction = 0.80) {
      TH1F* hist = getHistogram(pdg, layer, etaBin, ptBin);
      if (!hist || hist->GetEntries() < 100) return -1.f;

      double total_entries_hist = hist->Integral(1, hist->GetNbinsX());
      double target_entries = total_entries_hist * fraction;
      double accumulated_entries = 0;
      int lastBin = 0;

      for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
          accumulated_entries += hist->GetBinContent(bin);
          if (accumulated_entries >= target_entries) {
              lastBin = bin;
              break;
          }
      }

      double sum_x = 0;
      double sum_x_sq = 0;
      double entries_in_range = 0;

      for (int bin = 1; bin <= lastBin; ++bin) {
          double bin_center = hist->GetBinCenter(bin);
          double bin_content = hist->GetBinContent(bin);
          sum_x += bin_center * bin_content;
          sum_x_sq += bin_center * bin_center * bin_content;
          entries_in_range += bin_content;
      }

      double mean = sum_x / entries_in_range;
      double variance = (sum_x_sq / entries_in_range) - (mean * mean);

      if (variance < 0) variance = 0;

      return static_cast<float>(std::sqrt(variance));
  }

private:
  std::map<int, std::map<std::tuple<int, int, int>, TH1F>> mLUT;
};

// On-The-Fly Tracker Pid
struct OnTheFlyTrackerPid {

private:
  float calculateNsigma(float measuredToT, float expectedToT, float resolution) {
    return (measuredToT - expectedToT) / resolution;
  }

  // Helper function adapted from OnTheFlyTofPid to calculate track length to a given radius
  /// function to calculate track length of this track up to a certain radius
  /// \param track the input track
  /// \param radius the radius of the layer you're calculating the length to
  /// \param magneticField the magnetic field to use when propagating
  float computeTrackLength(o2::track::TrackParCov track, float radius, float magneticField)
  {
    // don't make use of the track parametrization
    float length = -100;

    o2::math_utils::CircleXYf_t trcCircle;
    float sna, csa;
    track.getCircleParams(magneticField, trcCircle, sna, csa);

    // distance between circle centers (one circle is at origin -> easy)
    const float centerDistance = std::hypot(trcCircle.xC, trcCircle.yC);

    // condition of circles touching - if not satisfied returned length will be -100
    if (centerDistance < trcCircle.rC + radius && centerDistance > std::fabs(trcCircle.rC - radius)) {
      length = 0.0f;

      // base radical direction
      const float ux = trcCircle.xC / centerDistance;
      const float uy = trcCircle.yC / centerDistance;
      // calculate perpendicular vector (normalized) for +/- displacement
      const float vx = -uy;
      const float vy = +ux;
      // calculate coordinate for radical line
      const float radical = (centerDistance * centerDistance - trcCircle.rC * trcCircle.rC + radius * radius) / (2.0f * centerDistance);
      // calculate absolute displacement from center-to-center axis
      const float displace = (0.5f / centerDistance) * std::sqrt(
                                                         (-centerDistance + trcCircle.rC - radius) *
                                                         (-centerDistance - trcCircle.rC + radius) *
                                                         (-centerDistance + trcCircle.rC + radius) *
                                                         (centerDistance + trcCircle.rC + radius));

      // possible intercept points of track and TOF layer in 2D plane
      const float point1[2] = {radical * ux + displace * vx, radical * uy + displace * vy};
      const float point2[2] = {radical * ux - displace * vx, radical * uy - displace * vy};

      // decide on correct intercept point
      std::array<float, 3> mom;
      track.getPxPyPzGlo(mom);
      const float scalarProduct1 = point1[0] * mom[0] + point1[1] * mom[1];
      const float scalarProduct2 = point2[0] * mom[0] + point2[1] * mom[1];

      // get start point
      std::array<float, 3> startPoint;
      track.getXYZGlo(startPoint);

      float cosAngle = -1000, modulus = -1000;

      if (scalarProduct1 > scalarProduct2) {
        modulus = std::hypot(point1[0] - trcCircle.xC, point1[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point1[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point1[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
      } else {
        modulus = std::hypot(point2[0] - trcCircle.xC, point2[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point2[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point2[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
      }
      cosAngle /= modulus;
      length = trcCircle.rC * std::acos(cosAngle);
      length *= std::sqrt(1.0f + track.getTgl() * track.getTgl());
    }
    return length;
  }

public:
  Produces<aod::UpgradeTrkPidSignals> tableUpgradeTrkPidSignals;
  Produces<aod::UpgradeTrkPids> tableUpgradeTrkPids;

  Service<o2::framework::O2DatabasePDG> pdg;
  ToTLUT mToTLUT;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> lutTotEl{"lutTotEl", "lut_tot_11.root", "ToT LUT for electrons"};
  Configurable<std::string> lutTotMu{"lutTotMu", "lut_tot_13.root", "ToT LUT for muons"};
  Configurable<std::string> lutTotPi{"lutTotPi", "lut_tot_211.root", "ToT LUT for pions"};
  Configurable<std::string> lutTotKa{"lutTotKa", "lut_tot_321.root", "ToT LUT for kaons"};
  Configurable<std::string> lutTotPr{"lutTotPr", "lut_tot_2212.root", "ToT LUT for protons"};
  Configurable<std::string> lutTotDe{"lutTotDe", "lut_tot_1000010020.root", "ToT LUT for deuteron"};
  Configurable<std::string> lutTotTr{"lutTotTr", "lut_tot_1000010030.root", "ToT LUT for triton"};
  Configurable<std::string> lutTotHe{"lutTotHe", "lut_tot_1000020030.root", "ToT LUT for helium-3"};
  Configurable<float> truncationFraction{"truncationFraction", 0.80f, "Fraction of lower entries to consider for truncated standard deviation (e.g., 0.80 for lower 80%)"};
  Configurable<float> dBz{"dBz", 20, "magnetic field (kilogauss) for track propagation"};

  std::vector<double> mLogBins;
  static constexpr int kNumLogBins = 200;

  std::map<int, std::string> mPdgToName = {
      {11, "e"}, {13, "mu"}, {211, "pi"}, {321, "ka"}, {2212, "pr"},
      {1000010020, "de"}, {1000010030, "tr"}, {1000020030, "he"}
  };

  std::vector<int> mPIDSortedPdgCodes = {
      11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030
  };

  void init(o2::framework::InitContext&) {
    bool loaded = true;

    loaded &= mToTLUT.load(11, lutTotEl.value);
    loaded &= mToTLUT.load(13, lutTotMu.value);
    loaded &= mToTLUT.load(211, lutTotPi.value);
    loaded &= mToTLUT.load(321, lutTotKa.value);
    loaded &= mToTLUT.load(2212, lutTotPr.value);
    loaded &= mToTLUT.load(1000010020, lutTotDe.value);
    loaded &= mToTLUT.load(1000010030, lutTotTr.value);
    loaded &= mToTLUT.load(1000020030, lutTotHe.value);

    if (!loaded) {
        std::cerr << "WARNING: Failed to load one or more ToT LUTs. PID results might be incomplete." << std::endl;
    }

    // Logarithmic momentum bins
    mLogBins.clear();
    double pMin = 0.05;
    double pMax = 10;
    double logMin = std::log10(pMin);
    double logMax = std::log10(pMax);
    double dLog = (logMax - logMin) / kNumLogBins;
    for (int i = 0; i <= kNumLogBins; ++i) {
      mLogBins.push_back(std::pow(10, logMin + i * dLog));
    }

    // To be done: loop through particles/histograms instead of creating every hist by hand
    histos.add("hToTvsP", "ToT vs #it{p}; #it{p} (GeV/#it{c}); ToT (#mus/10#mum)",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(400, 0., 200.)});
    histos.add("hToTvsPt", "ToT vs #it{p}_{T}; #it{p}_{T} (GeV/#it{c}); ToT (#mus/10#mum)",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(400, 0., 200.)});
    histos.add("hHitLayers", "Number of hits on each detector layer;Layer;Counts",
               kTH1F, {AxisSpec(TrackerConstants::kMaxBarrelLayers, -0.5, TrackerConstants::kMaxBarrelLayers - 0.5)});
    histos.add("hHitMultiplicity", "Hit multiplicity along the track; # Hits per track;Counts",
               kTH1F, {AxisSpec(TrackerConstants::kMaxBarrelLayers, -0.5, TrackerConstants::kMaxBarrelLayers - 0.5)});

    histos.add("h2dBarrelNsigmaTrueElecVsElecHypothesis", "Nsigma (True e vs Hyp e); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueMuonVsMuonHypothesis", "Nsigma (True #mu vs Hyp #mu); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTruePionVsPionHypothesis", "Nsigma (True #pi vs Hyp #pi); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueKaonVsKaonHypothesis", "Nsigma (True K vs Hyp K); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueProtVsProtHypothesis", "Nsigma (True p vs Hyp p); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    // histos.add("h2dBarrelNsigmaTrueDeutVsDeutHypothesis", "Nsigma (True d vs Hyp d); #it{p} (GeV/#it{c}); N#sigma",
    //            kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    // histos.add("h2dBarrelNsigmaTrueTritVsTritHypothesis", "Nsigma (True t vs Hyp t); #it{p} (GeV/#it{c}); N#sigma",
    //            kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    // histos.add("h2dBarrelNsigmaTrueHe3VsHe3Hypothesis", "Nsigma (True He vs Hyp He); #it{p} (GeV/#it{c}); N#sigma",
    //            kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});

    histos.add("h2dBarrelNsigmaTrueElecVsMuonHypothesis", "Nsigma (True e vs Hyp #mu); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueMuonVsElecHypothesis", "Nsigma (True #mu vs Hyp e); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueElecVsPionHypothesis", "Nsigma (True e vs Hyp #pi); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTruePionVsElecHypothesis", "Nsigma (True #pi vs Hyp e); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueElecVsKaonHypothesis", "Nsigma (True e vs Hyp K); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueKaonVsElecHypothesis", "Nsigma (True K vs Hyp e); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueElecVsProtHypothesis", "Nsigma (True e vs Hyp p); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueProtVsElecHypothesis", "Nsigma (True p vs Hyp e); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});

    histos.add("h2dBarrelNsigmaTrueMuonVsPionHypothesis", "Nsigma (True #mu vs Hyp #pi); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTruePionVsMuonHypothesis", "Nsigma (True #pi vs Hyp #mu); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueMuonVsKaonHypothesis", "Nsigma (True #mu vs Hyp K); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueKaonVsMuonHypothesis", "Nsigma (True K vs Hyp #mu); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueMuonVsProtHypothesis", "Nsigma (True #mu vs Hyp p); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueProtVsMuonHypothesis", "Nsigma (True p vs Hyp #mu); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});

    histos.add("h2dBarrelNsigmaTruePionVsKaonHypothesis", "Nsigma (True #pi vs Hyp K); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueKaonVsPionHypothesis", "Nsigma (True K vs Hyp #pi); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTruePionVsProtHypothesis", "Nsigma (True #pi vs Hyp p); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueProtVsPionHypothesis", "Nsigma (True p vs Hyp #pi); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});

    histos.add("h2dBarrelNsigmaTrueKaonVsProtHypothesis", "Nsigma (True K vs Hyp p); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
    histos.add("h2dBarrelNsigmaTrueProtVsKaonHypothesis", "Nsigma (True p vs Hyp K); #it{p} (GeV/#it{c}); N#sigma",
               kTH2F, {AxisSpec(mLogBins), AxisSpec(200, -10., 10.)});
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
               soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks,
               aod::McParticles const&,
               aod::McCollisions const&) {

    std::array<float, 6> mcPvCov = {0.};
    o2::dataformats::VertexBase mcPvVtx({0.0f, 0.0f, 0.0f}, mcPvCov);

    std::array<float, TrackerConstants::kMaxBarrelLayers> timeOverThresholdBarrel;

    std::array<int, 8> hypothesisPdgCodes = {
        11,         // Electron
        13,         // Muon
        211,        // Pion
        321,        // Kaon
        2212,       // Proton
        1000010020, // Deuteron
        1000010030, // Triton
        1000020030  // Helium-3
    };

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }

      const auto& mcParticle = track.mcParticle();
      const auto& pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        continue;
      }

      const float pt = mcParticle.pt();
      const float p = mcParticle.p();
      const float eta = mcParticle.eta();
      const int truePdgCode = std::abs(mcParticle.pdgCode());

      // Get momentum bin index for histograms
      auto it_p = std::upper_bound(mLogBins.begin(), mLogBins.end(), p);
      int pLogBin = std::distance(mLogBins.begin(), it_p) - 1;
      pLogBin = std::max(0, std::min(pLogBin, kNumLogBins - 1));

      // Compute bin indices in (pt, eta) for accessing the appropriate ToT LUT histogram
      const float clampedPt = std::max(TrackerConstants::kPtMin, std::min(pt, TrackerConstants::kPtMax - 1e-6f));
      const float ptFraction = (clampedPt - TrackerConstants::kPtMin)/(TrackerConstants::kPtMax - TrackerConstants::kPtMin);
      const int binnedPt = std::min(static_cast<int>(ptFraction * TrackerConstants::kPtBins), TrackerConstants::kPtBins-1);

      const float clampedEta = std::max(TrackerConstants::kEtaMin, std::min(eta, TrackerConstants::kEtaMax - 1e-6f));
      const float etaFraction = (clampedEta - TrackerConstants::kEtaMin)/(TrackerConstants::kEtaMax - TrackerConstants::kEtaMin);
      const int binnedEta = std::min(static_cast<int>(etaFraction * TrackerConstants::kEtaBins), TrackerConstants::kEtaBins-1);

      // Hit determination
      uint16_t hitMap = 0;
      int nHitLayers = 0;
      o2::track::TrackParCov o2track = o2::upgrade::convertMCParticleToO2Track(mcParticle, pdg);
      // Propagate the MC track to DCA for proper starting point for length calculation
      float xPv = -100.f;
      static constexpr float kTrkXThreshold = -99.f; // Threshold to consider a good propagation of the track
      if (o2track.propagateToDCA(mcPvVtx, dBz)) {
        xPv = o2track.getX();
      }

      if (xPv > kTrkXThreshold) { // Only proceed if propagation to vertex was successful
        for (int layer = 0; layer < TrackerConstants::kMaxBarrelLayers; ++layer) {
          float layerRadius = TrackerConstants::kTrackerRadii[layer];
          float trackLength = computeTrackLength(o2track, layerRadius, dBz);
          
          if (trackLength > 0) { // If a valid length is returned, the layer is considered "hit"
            hitMap |= (1 << layer); // Set the bit for this layer
            histos.fill(HIST("hHitLayers"), layer);
            nHitLayers++;
          }
        }
      }

      histos.fill(HIST("hHitMultiplicity"), nHitLayers);

      std::vector<float> validToTs;

      // Simulate measured ToT for the track by sampling from the true particle's LUTs
      for (int layer = 3; layer < TrackerConstants::kMaxBarrelLayers; ++layer) { // don't consider vertex layers for ToT
        timeOverThresholdBarrel[layer] = -1;
        if ((hitMap >> layer) & 0x1) { // Check if this layer was hit geometrically
          TH1F* totHist = mToTLUT.getHistogram(truePdgCode, layer, binnedEta, binnedPt);
          if (totHist && totHist->GetEntries() > 100) {
            float val = totHist->GetRandom();
            timeOverThresholdBarrel[layer] = val;
            validToTs.push_back(val);
          }
        }
      }

      float truncatedMeanToT = -1.0f;
      const size_t nValid = validToTs.size();
      size_t nUse = 0;

      // Number of points to be used for truncated (geometric) mean calculation
      if (nValid == 4) nUse = 2;
      else if (nValid == 5 || nValid == 6) nUse = 3;
      else if (nValid >= 7) nUse = 4; // Since only tracks with > 6 hits, this line would be enough for now

      if (nUse > 0 && nValid >= nUse) {
        std::sort(validToTs.begin(), validToTs.end());
        float product = 1.0f;
        for (size_t i = 0; i < nUse; ++i) {
          product *= validToTs[i];
        }
        truncatedMeanToT = std::pow(product, 1.0f / static_cast<float>(nUse)); // geometric mean

        histos.fill(HIST("hToTvsPt"), pt, truncatedMeanToT);
        histos.fill(HIST("hToTvsP"), p, truncatedMeanToT);
      }

      std::array<float, 8> nSigmaValues;
      nSigmaValues.fill(999.f);

      if (truncatedMeanToT > 0) {
        for (size_t i = 0; i < hypothesisPdgCodes.size(); ++i) {
          int hypPdgCode = hypothesisPdgCodes[i];

          float sumExpectedToTMPV = 0.f;
          float sumExpectedToTResolution = 0.f;
          int numValidHistograms = 0;

          // Iterate through layers based on the geometrically determined 'hitMap'
          for (int layer = 0; layer < TrackerConstants::kMaxBarrelLayers; ++layer) {
              if ((hitMap >> layer) & 0x1) { // Check if this layer was hit geometrically
                  TH1F* hypToTHist = mToTLUT.getHistogram(hypPdgCode, layer, binnedEta, binnedPt);
                  if (hypToTHist && hypToTHist->GetEntries() > 100) {
                      float mpv = mToTLUT.getMPV(hypPdgCode, layer, binnedEta, binnedPt);
                      float resolution = mToTLUT.getTruncatedStdDev(hypPdgCode, layer, binnedEta, binnedPt, truncationFraction.value);

                      if (mpv > 0 && resolution > 0) {
                          sumExpectedToTMPV += mpv;
                          sumExpectedToTResolution += resolution;
                          numValidHistograms++;
                      }
                  }
              }
          }

          float expectedToT = (numValidHistograms > 0) ? sumExpectedToTMPV / numValidHistograms : -1.f;
          float averagedResolution = (numValidHistograms > 0) ? sumExpectedToTResolution / numValidHistograms : -1.f;

          if (expectedToT > 0 && averagedResolution > 0) {
              nSigmaValues[i] = calculateNsigma(truncatedMeanToT, expectedToT, averagedResolution);

              // Same particle
              if (truePdgCode == 11 && hypPdgCode == 11) histos.fill(HIST("h2dBarrelNsigmaTrueElecVsElecHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 13 && hypPdgCode == 13) histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsMuonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 211 && hypPdgCode == 211) histos.fill(HIST("h2dBarrelNsigmaTruePionVsPionHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 321 && hypPdgCode == 321) histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsKaonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 2212 && hypPdgCode == 2212) histos.fill(HIST("h2dBarrelNsigmaTrueProtVsProtHypothesis"), p, nSigmaValues[i]);
              // else if (truePdgCode == 1000010020 && hypPdgCode == 1000010020) histos.fill(HIST("h2dBarrelNsigmaTrueDeutVsDeutHypothesis"), p, nSigmaValues[i]);
              // else if (truePdgCode == 1000010030 && hypPdgCode == 1000010030) histos.fill(HIST("h2dBarrelNsigmaTrueTritVsTritHypothesis"), p, nSigmaValues[i]);
              // else if (truePdgCode == 1000020030 && hypPdgCode == 1000020030) histos.fill(HIST("h2dBarrelNsigmaTrueHe3VsHe3Hypothesis"), p, nSigmaValues[i]);

              // Electron combinations
              else if (truePdgCode == 11 && hypPdgCode == 13) histos.fill(HIST("h2dBarrelNsigmaTrueElecVsMuonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 13 && hypPdgCode == 11) histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsElecHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 11 && hypPdgCode == 211) histos.fill(HIST("h2dBarrelNsigmaTrueElecVsPionHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 211 && hypPdgCode == 11) histos.fill(HIST("h2dBarrelNsigmaTruePionVsElecHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 11 && hypPdgCode == 321) histos.fill(HIST("h2dBarrelNsigmaTrueElecVsKaonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 321 && hypPdgCode == 11) histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsElecHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 11 && hypPdgCode == 2212) histos.fill(HIST("h2dBarrelNsigmaTrueElecVsProtHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 2212 && hypPdgCode == 11) histos.fill(HIST("h2dBarrelNsigmaTrueProtVsElecHypothesis"), p, nSigmaValues[i]);

              // Muon combinations
              else if (truePdgCode == 13 && hypPdgCode == 211) histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsPionHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 211 && hypPdgCode == 13) histos.fill(HIST("h2dBarrelNsigmaTruePionVsMuonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 13 && hypPdgCode == 321) histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsKaonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 321 && hypPdgCode == 13) histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsMuonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 13 && hypPdgCode == 2212) histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsProtHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 2212 && hypPdgCode == 13) histos.fill(HIST("h2dBarrelNsigmaTrueProtVsMuonHypothesis"), p, nSigmaValues[i]);

              // Pion combinations
              else if (truePdgCode == 211 && hypPdgCode == 321) histos.fill(HIST("h2dBarrelNsigmaTruePionVsKaonHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 321 && hypPdgCode == 211) histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsPionHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 211 && hypPdgCode == 2212) histos.fill(HIST("h2dBarrelNsigmaTruePionVsProtHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 2212 && hypPdgCode == 211) histos.fill(HIST("h2dBarrelNsigmaTrueProtVsPionHypothesis"), p, nSigmaValues[i]);

              // Kaon-Proton
              else if (truePdgCode == 321 && hypPdgCode == 2212) histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsProtHypothesis"), p, nSigmaValues[i]);
              else if (truePdgCode == 2212 && hypPdgCode == 321) histos.fill(HIST("h2dBarrelNsigmaTrueProtVsKaonHypothesis"), p, nSigmaValues[i]);

          } else {
              nSigmaValues[i] = 999.f; // Dummy holder
          }
        }
      }

      tableUpgradeTrkPidSignals(truncatedMeanToT);
      tableUpgradeTrkPids(nSigmaValues[0], nSigmaValues[1], nSigmaValues[2], nSigmaValues[3],
                         nSigmaValues[4], nSigmaValues[5], nSigmaValues[6], nSigmaValues[7]);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTrackerPid>(cfgc)};
}