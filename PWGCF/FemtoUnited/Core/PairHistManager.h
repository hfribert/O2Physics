//
// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PairHistManager.h
/// \brief histogram manager class for pair tasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_PAIRHISTMANAGER_H_
#define PWGCF_FEMTOUNITED_CORE_PAIRHISTMANAGER_H_

#include "PWGCF/FemtoUnited/Core/FemtoUtils.h"
#include "PWGCF/FemtoUnited/Core/HistManager.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "Framework/HistogramRegistry.h"

#include "Math/Boost.h"
#include "Math/Vector4D.h"
#include "TMath.h"

#include <array>
#include <map>
#include <string>
#include <vector>

namespace o2::analysis::femtounited
{
namespace pairhistmanager
{
// enum for track histograms
enum PairHist {
  // kinemtics
  kKstar,
  kKt,
  kMt,
  // 2d qa
  kPt1VsPt2,
  kPt1VsKstar,
  kPt2VsKstar,
  kPt1VsKt,
  kPt2VsKt,
  kPt1VsMt,
  kPt2VsMt,
  kPairHistogramLast
};

constexpr char PrefixTrackTrackSe[] = "TrackTrack/SE/";
constexpr char PrefixTrackTrackMe[] = "TrackTrack/ME/";
constexpr char PrefixTrackVzeroSe[] = "TrackVzero/ME/";
constexpr char PrefixTrackVzeroMe[] = "TrackVzero/ME/";

constexpr std::string_view AnalysisDir = "Analysis/";
constexpr std::string_view QaDir = "QA/";

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<PairHist>, kPairHistogramLast> HistTable = {
  {{kKstar, o2::framework::kTH1F, "hKstar", "k*; k* (GeV/#it{c}); Entries"},
   {kKt, o2::framework::kTH1F, "hKt", "transverse momentum; k_{T} (GeV/#it{c}); Entries"},
   {kMt, o2::framework::kTH1F, "hMt", "transverse mass; m_{T} (GeV/#it{c}^{2}); Entries"},
   {kPt1VsPt2, o2::framework::kTH2F, "hPt1VsPt2", "track1 p_{T} vs track2 p_{T}; track1 p_T (GeV/#it{c}); track2 p_{T} (GeV/#it{c})"},
   {kPt1VsKstar, o2::framework::kTH2F, "hPt1VsKstar", "track1 p_{T} vs k*; track1 p_{T} (GeV/#it{c}); k* (GeV/#it{c})"},
   {kPt2VsKstar, o2::framework::kTH2F, "hPt2VsKstar", "track2 p_{T} vs k*; track2 p_{T} (GeV/#it{c}); k* (GeV/#it{c})"},
   {kPt1VsKt, o2::framework::kTH2F, "hPt1VsKt", "track1 p_{T} vs k_{T}; track1 p_{T} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
   {kPt2VsKt, o2::framework::kTH2F, "hPt2VsKt", "track2 p_{T} vs k_{T}; track2 p_{T} (GeV/#it{c}); k_{T} (GeV/#it{c})"},
   {kPt1VsMt, o2::framework::kTH2F, "hPt1VsMt", "track1 p_{T} vs m_{T}; track1 p_{T} (GeV/#it{c}); m_{T} (GeV/#it{c})"},
   {kPt2VsMt, o2::framework::kTH2F, "hPt2VsMt", "track1 p_{T} vs m_{T}; track2 p_{T} (GeV/#it{c}); m_{T} (GeV/#it{c})"}}};

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix>
class PairHistManager
{
 public:
  /// Destructor
  virtual ~PairHistManager() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry, std::map<PairHist, std::vector<o2::framework::AxisSpec>> Specs)
  {
    mHistogramRegistry = registry;

    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + GetHistNamev2(kKstar, HistTable), GetHistDesc(kKstar, HistTable), GetHistType(kKstar, HistTable), {Specs[kKstar]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kMt, HistTable), GetHistDesc(kMt, HistTable), GetHistType(kMt, HistTable), {Specs[kMt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt1VsPt2, HistTable), GetHistDesc(kPt1VsPt2, HistTable), GetHistType(kPt1VsPt2, HistTable), {Specs[kPt1VsPt2]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt1VsKstar, HistTable), GetHistDesc(kPt1VsKstar, HistTable), GetHistType(kPt1VsKstar, HistTable), {Specs[kPt1VsKstar]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt2VsKstar, HistTable), GetHistDesc(kPt2VsKstar, HistTable), GetHistType(kPt2VsKstar, HistTable), {Specs[kPt2VsKstar]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt1VsKt, HistTable), GetHistDesc(kPt1VsKt, HistTable), GetHistType(kPt1VsKt, HistTable), {Specs[kPt1VsKt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt2VsKt, HistTable), GetHistDesc(kPt2VsKt, HistTable), GetHistType(kPt2VsKt, HistTable), {Specs[kPt2VsKt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt1VsMt, HistTable), GetHistDesc(kPt1VsMt, HistTable), GetHistType(kPt1VsMt, HistTable), {Specs[kPt1VsMt]});
      mHistogramRegistry->add(analysisDir + GetHistNamev2(kPt2VsMt, HistTable), GetHistDesc(kPt2VsMt, HistTable), GetHistType(kPt2VsMt, HistTable), {Specs[kPt2VsMt]});
    }

    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    // std::string qaDir = std::string(prefix) + std::string(QaDir);
    // }
  }

  void setMass(int PdgParticle1, int PdgParticle2)
  {
    mMass1 = o2::analysis::femtounited::utils::getMass(PdgParticle1);
    mMass1 = o2::analysis::femtounited::utils::getMass(PdgParticle2);
  }

  template <typename T1, typename T2>
  void setPair(T1 particle1, T2 particle2)
  {
    mTrack1 = ROOT::Math::PtEtaPhiMVector(particle1.pt(), particle1.eta(), particle1.phi(), mMass1);
    mTrack2 = ROOT::Math::PtEtaPhiMVector(particle2.pt(), particle2.eta(), particle2.phi(), mMass2);
    auto partSum = mTrack1 + mTrack1;

    // set kT
    mKt = partSum.Pt();

    // set mT
    float averageMass = (mMass1 + mMass2) / 2.;
    mMt = std::hypot(mKt, averageMass);

    // boost to pair rest frame to get kstar
    ROOT::Math::PxPyPzEVector track1Cms(mTrack1);
    ROOT::Math::PxPyPzEVector track2Cms(mTrack2);
    ROOT::Math::Boost boostPrf = ROOT::Math::Boost(partSum.BoostToCM());
    track1Cms = boostPrf(mTrack1);
    track2Cms = boostPrf(mTrack2);
    auto trackRel = track1Cms - track2Cms;

    // set kstar
    mKstar = 0.5 * trackRel.P();
  }

  template <modes::Mode mode>
  void fill()
  {
    if constexpr (isFlagSet(mode, modes::Mode::kANALYSIS)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kKstar, HistTable)), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kMt, HistTable)), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt1VsPt2, HistTable)), mTrack1.Pt(), mTrack2.Pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt1VsKstar, HistTable)), mTrack1.Pt(), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt1VsMt, HistTable)), mTrack1.Pt(), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt1VsKt, HistTable)), mTrack1.Pt(), mKt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt2VsKstar, HistTable)), mTrack2.Pt(), mKstar);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt2VsMt, HistTable)), mTrack2.Pt(), mMt);
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(GetHistName(kPt2VsKt, HistTable)), mTrack2.Pt(), mKt);
    }

    // if constexpr (isFlagSet(mode, modes::Mode::kQA)) {
    //  mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(GetHistName(kPtVsDcaz, HistTable)), track.pt(), track.dcaZ());
    // }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry;
  float mMass1 = 0.f;
  float mMass2 = 0.f;
  ROOT::Math::PtEtaPhiMVector mTrack1{};
  ROOT::Math::PtEtaPhiMVector mTrack2{};
  float mKstar = 0.f;
  float mKt = 0.f;
  float mMt = 0.f;
};
}; // namespace pairhistmanager
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_PAIRHISTMANAGER_H_
